function regOverTime_latShift_dF(xlFile, sheet, beh, varargin)

p = inputParser;
% default parameters if none given
p.addParameter('bernFlag', 1);
p.addParameter('shiftForward', 1);
p.addParameter('modelName', 'sixParam_absPePeAN_bi')
p.addParameter('vars', [{'O'}]);
p.addParameter('modelVars', [{'peBar'} {'pePe'}]);
p.addParameter('coefWindow', []);
p.parse(varargin{:});

[root, sep] = currComputer();

[numbers, sessionCellList, ~] = xlsread(xlFile, sheet);
revForFlag = numbers(:,1);
intanFlag = numbers(:,2);
cellList = sessionCellList(2:end, 1);
sessionList = sessionCellList(2:end, 2);

allVarNames = [p.Results.vars p.Results.modelVars];

if size(numbers,2) > 2
    trialList = numbers(:,3:4);
else
    trialList = nan(length(cellList),2);
end 

%set time window for spike analyses
tB = -1500;
tF = 5000;
tW = 500;
tS = 100;
timeToPlot = tB:tS:tF - tW;

smoothKern = normpdf(0:5000, 0, 250);
smoothKern = smoothKern/sum(smoothKern);

mdlStruct = [];
sessionName = [];
lastwarn('');

for currCell = 1:length(sessionList)
    
    fprintf('Analyzing cell %d of %d \n', currCell, length(sessionList));
    
    if strcmp(sessionName, sessionList{currCell}) == 0
        sessionName = sessionList{currCell};
        [animalName, date] = strtok(sessionName, 'd'); 
        animalName = animalName(2:end);

        %load spike data
        if isstrprop(sessionName(end), 'alpha')
            sortedFolderLocation = [root animalName sep sessionName(1:end-1) sep 'sorted' sep 'session ' sessionName(end) sep];
        else
            sortedFolderLocation = [root animalName sep sessionName sep 'sorted' sep 'session' sep];
        end
        sortedFolder = dir(sortedFolderLocation);

        if revForFlag(currCell) == 1
            if any(~cellfun(@isempty,strfind({sortedFolder.name},'_intan.mat'))) == 1
                sessionDataInd = ~cellfun(@isempty,strfind({sortedFolder.name},'_intan.mat'));
                load([sortedFolderLocation sortedFolder(sessionDataInd).name])
            else
                [sessionData] = generateSessionData_intan_operantMatching(sessionName);
            end
        else
            if any(~cellfun(@isempty,strfind({sortedFolder.name},'_nL.mat'))) == 1
                sessionDataInd = ~cellfun(@isempty,strfind({sortedFolder.name},'_nL.mat'));
                load([sortedFolderLocation sortedFolder(sessionDataInd).name])
            elseif any(~cellfun(@isempty,strfind({sortedFolder.name},'_intan.mat'))) == 1
                sessionDataInd = ~cellfun(@isempty,strfind({sortedFolder.name},'_intan.mat'));
                load([sortedFolderLocation sortedFolder(sessionDataInd).name])
            else
                if intanFlag(currCell) == 1
                    [sessionData] = generateSessionData_intan_operantMatching(sessionName);
                else
                    [sessionData] = generateSessionData_nL_operantMatching(sessionName);
                end
            end
        end
        [s] = behDataForStruct_opMD(sessionName, 'revForFlag', revForFlag(currCell));
        O = s.allRewardsBinary;
        C = s.allChoices;
        X = O.*C;
        pO = [nan O(1:end-1)];
        
        spikeFields = fields(sessionData);
        cellInds = find(~cellfun(@isempty,strfind(spikeFields, 'TT')));
        sessionTime = (sessionData(end).CSon + tF) - (sessionData(1).CSon + tB);
        csTimes = [sessionData.CSon] -  (sessionData(1).CSon + tB);
        csTimes = [csTimes(s.responseInds) csTimes(end) + tF];
        tEndTimes = [0 csTimes(1:end-1)+1500];
        
        %get model terms
        if p.Results.bernFlag
            modelPath = [root animalName sep animalName 'sorted' sep 'stan' sep 'bernoulli' sep p.Results.modelName sep animalName...
            beh '_' p.Results.modelName '.mat'];
        else
            modelPath = [root animalName sep animalName 'sorted' sep 'stan' sep p.Results.modelName sep animalName...
            beh '_' p.Results.modelName '.mat'];
        end
        t = generateStanModelTerms_opMD(p.Results.modelName, modelPath, sessionName, 1, revForFlag(currCell));
        sumQ = [zscore(t.Q(:,1) + t.Q(:,2))]';
        diffQ = [zscore(t.Q(:,1) - t.Q(:,2))]';
        confQ = [zscore(abs(diffQ))]';
        pe = [zscore(abs(t.pe))]';
     
        %extract relevant variables for use as regressors
        vars = [];
        for termInd = 1:length(p.Results.vars)
            vars(termInd,:) = eval([p.Results.vars{termInd} ';']);
        end
        modelVars = [];
        for termInd = 1:length(p.Results.modelVars)
            modelVars(termInd,:) = zscore(eval(['t.' p.Results.modelVars{termInd} ';']));
        end
        
        %add variables that are shifted by one trial
        if isempty(vars)
            vars_s = [];
        elseif p.Results.shiftForward
            vars_s = [nan(size(vars,1), 1) vars(:,1:end-1)];
        else
            vars_s = [vars(:,2:end) nan(size(vars,1), 1)];
        end
        if isempty(modelVars)
            modelVars_s = [];
        elseif p.Results.shiftForward
            modelVars_s = [nan(size(modelVars,1),1) modelVars(:,1:end-1)];
        else
            modelVars_s = [modelVars(:,2:end) nan(size(modelVars,1), 1)];
        end
        
        allVars = [vars; modelVars];
        allVars_s = [vars_s; modelVars_s];
    end

    %% smooth and z-score spikes
    cellInd = find(~cellfun(@isempty,strfind(spikeFields, cellList{currCell})));
    cellNum = find(cellInds == cellInd);
    spikeTimes = sessionData(cellNum).allSpikes - (sessionData(1).CSon + tB);
    spikeTimes = spikeTimes(spikeTimes > 0 & spikeTimes <= sessionTime);
    sessionSpikes = zeros(1, sessionTime);
    sessionSpikes(spikeTimes) = 1;
    smoothSpikes = conv(sessionSpikes, smoothKern);
    smoothSpikes = smoothSpikes(1:(end-length(smoothKern)+1));
    zSpikes = zscore(smoothSpikes);

    for regTime = timeToPlot
        spikeTmp = zeros(1,length(s.responseInds));
        for trialInd = 1:length(s.responseInds)
            winInds = [(csTimes(trialInd) + regTime + 1) : (csTimes(trialInd) + regTime + tW)];
            winInds = winInds(winInds > tEndTimes(trialInd) & winInds <  csTimes(trialInd + 1));
            if ~isempty(winInds)
                spikeTmp(trialInd) = nanmean(zSpikes(winInds));
            end
        end
        
        glm = fitlm([allVars' allVars_s'], spikeTmp);
        
        regName = ['reg_' num2str(regTime)];
        regName = strrep(regName, '-', 'min');
        mdlStruct(currCell).(regName).coef = glm.Coefficients.Estimate(2:end);
        mdlStruct(currCell).(regName).pVal = glm.Coefficients.pValue(2:end);
        mdlStruct(currCell).(regName).tStat = glm.Coefficients.tStat(2:end);
        mdlStruct(currCell).(regName).rSqr = glm.Rsquared(1).Ordinary(1);
        mdlStruct(currCell).(regName).warning = lastwarn;
        lastwarn('');
    end
end

for termInd = 1:length(allVarNames)
    eval([char(allVarNames(termInd)) '_sig = [];']);
    eval([char(allVarNames(termInd)) '_coef = [];']);
    eval([char(allVarNames(termInd)) '_s_sig = [];']);
    eval([char(allVarNames(termInd)) '_s_coef = [];']);
end

rSqr = [];
warn_inds = logical([]);
numTerms = length(allVarNames);
for n = 1:length(sessionList)
    for t = timeToPlot
        regName = ['reg_' num2str(t)];
        regName = strrep(regName, '-', 'min');

        t_ind = find(timeToPlot == t);
        for termInd = 1:numTerms
            eval([char(allVarNames(termInd)) '_sig(t_ind, n) = mdlStruct(n).(regName).pVal(termInd);']);
            eval([char(allVarNames(termInd)) '_coef(t_ind, n) = mdlStruct(n).(regName).coef(termInd);']);
            eval([char(allVarNames(termInd)) '_s_sig(t_ind, n) = mdlStruct(n).(regName).pVal(termInd+numTerms);']);
            eval([char(allVarNames(termInd)) '_s_coef(t_ind, n) = mdlStruct(n).(regName).coef(termInd+numTerms);']);
        end
        
        rSqr(t_ind, n) = mdlStruct(n).(regName).rSqr;
        warn_inds(t_ind, n) = ~isempty(mdlStruct(n).(regName).warning);
    end
end

for termInd = 1:numTerms
    eval([char(allVarNames(termInd)) '_sig(warn_inds) = NaN;']);
    eval([char(allVarNames(termInd)) '_coef(warn_inds) = NaN;']);
    eval([char(allVarNames(termInd)) '_s_sig(warn_inds) = NaN;']);
    eval([char(allVarNames(termInd)) '_s_coef(warn_inds) = NaN;']);
end

rSqr(warn_inds) = NaN;

%plot signicance across populatioon
figure; hold on;
colors = hsv(numTerms*2);
pCut = 0.05;
lgTxt = [];
for termInd = 1:numTerms
    eval(['plot(timeToPlot + tW, nanmean(' char(allVarNames(termInd))...
        '_sig < pCut, 2), ''linewidth'', 2, ''color'', colors(termInd,:))']);
    eval(['plot(timeToPlot + tW, nanmean(' char(allVarNames(termInd))...
        '_s_sig < pCut, 2), ''linewidth'', 2, ''color'', colors(termInd+numTerms,:))']);
    lgTxt = [lgTxt allVarNames(termInd), strcat(allVarNames(termInd), ' shifted')];
end
xlim([timeToPlot(1)+tW timeToPlot(end)+tW])
legend(lgTxt)
xlabel('Time from cue (ms)')
ylabel('Fraction of significant neurons (P < 0.05)')
ylim([0 1])
set(gca, 'tickdir', 'out')

%plot coefficient values
figure; hold on;
colors = hsv(size(rSqr,2));
for currCell = 1:size(rSqr,2)
    plot(timeToPlot + tW, rSqr(:,currCell), 'color', colors(currCell,:))
end
xlabel('Time from cue (ms)')
ylabel('R^2')
set(gca, 'tickdir', 'out')

titleTxt = strrep([sheet ' - ' p.Results.modelName], '_', ' ');
suptitle(titleTxt)
set(gcf, 'renderer', 'painters')


%plot coefficient distributions
if ~isempty(p.Results.coefWindow)
    winInd = find(timeToPlot == 0);
else
    winInd = 1:length(timeToPlot);
end
numRows = ceil(numTerms*2/4);
colors = hsv(numTerms*2);

figure; hold on;
for currVar = 1:numTerms
    subplot(numRows,4,currVar); hold on; 
    eval([allVarNames{currVar} '_coef = ' allVarNames{currVar} '_coef(winInd,:);']);
    eval(['sigInds = find(' allVarNames{currVar} '_sig(winInd,:) < pCut);']);
    eval(['histogram(' allVarNames{currVar} '_coef, 20, ''FaceColor'', ''k'', ''FaceAlpha'', 0.3)']);
    eval(['histogram(' allVarNames{currVar} '_coef(sigInds), 20, ''FaceColor'', colors(currVar,:))']);
    title(allVarNames{currVar})
    set(gca, 'tickdir', 'out')
    
    subplot(numRows,4,currVar+numTerms); hold on; 
    eval([allVarNames{currVar} '_s_coef = ' allVarNames{currVar} '_s_coef(winInd,:);']);
    eval(['sigInds = find(' allVarNames{currVar} '_sig(winInd,:) < pCut);']);
    eval(['histogram(' allVarNames{currVar} '_s_coef, 20, ''FaceColor'', ''k'', ''FaceAlpha'', 0.3)']);
    eval(['histogram(' allVarNames{currVar} '_s_coef(sigInds), 20, ''FaceColor'', colors(currVar+numTerms,:))']);
    title([allVarNames{currVar} ' shift'])
    set(gca, 'tickdir', 'out')
end


set(gcf, 'renderer', 'painters')


