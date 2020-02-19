function regOverTime_lat_dF(xlFile, sheet, beh, varargin)

p = inputParser;
% default parameters if none given
p.addParameter('bernFlag', 1);
p.addParameter('modelName', 'sixParam_absPePeAN_bi')
p.addParameter('allVars', [{'O'}, {'C'}, {'X'}, {'pO'}, {'diffQ'}, {'sumQ'}, {'confQ'}, {'pe'}]);
p.addParameter('modelVars', [{'peBar'} {'pePe'}]);
p.addParameter('coefWindow', []);
p.parse(varargin{:});

[root, sep] = currComputer();

[numbers, sessionCellList, ~] = xlsread(xlFile, sheet);
revForFlag = numbers(:,1);
intanFlag = numbers(:,2);
cellList = sessionCellList(2:end, 1);
sessionList = sessionCellList(2:end, 2);

if size(numbers,2) > 2
    trialList = numbers(:,3:4);
else
    trialList = nan(length(cellList),2);
end 

%set time window for spike analyses
tB = -3000;
tF = 5000;
tW = 2000;
tS = 400;
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
     
        allVars = [];
        for termInd = 1:length(p.Results.allVars)
            allVars(termInd,:) = eval([p.Results.allVars{termInd} ';']);
        end
        modelVars = [];
        for termInd = 1:length(p.Results.modelVars)
            modelVars(termInd,:) = zscore(eval(['t.' p.Results.modelVars{termInd} ';']));
        end
        modelVars = [modelVars(2:end) NaN];
    end

    %% Sort all spikes into a raster-able matrix
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
        
        glm = fitlm([allVars' modelVars'], spikeTmp);
        
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

for termInd = 1:length(p.Results.allVars)
    eval([char(p.Results.allVars(termInd)) '_sig = [];']);
    eval([char(p.Results.allVars(termInd)) '_coef = [];']);
end
for termInd = 1:length(p.Results.modelVars)
    eval([char(p.Results.modelVars(termInd)) '_sig = [];']);
    eval([char(p.Results.modelVars(termInd)) '_coef = [];']);
end
rSqr = [];
warn_inds = logical([]);
numTerms = length(p.Results.allVars);
for n = 1:length(sessionList)
    for t = timeToPlot
        regName = ['reg_' num2str(t)];
        regName = strrep(regName, '-', 'min');
        t_ind = find(timeToPlot == t);
        for termInd = 1:length(p.Results.allVars)
            eval([char(p.Results.allVars(termInd)) '_sig(t_ind, n) = mdlStruct(n).(regName).pVal(termInd);']);
            eval([char(p.Results.allVars(termInd)) '_coef(t_ind, n) = mdlStruct(n).(regName).coef(termInd);']);
        end
        for termInd = 1:length(p.Results.modelVars)
            eval([char(p.Results.modelVars(termInd)) '_sig(t_ind, n) = mdlStruct(n).(regName).pVal(numTerms+termInd);']);
            eval([char(p.Results.modelVars(termInd)) '_coef(t_ind, n) = mdlStruct(n).(regName).coef(numTerms+termInd);']);
        end
        
        rSqr(t_ind, n) = mdlStruct(n).(regName).rSqr;
        warn_inds(t_ind, n) = ~isempty(mdlStruct(n).(regName).warning);
    end
end

for termInd = 1:length(p.Results.allVars)
    eval([char(p.Results.allVars(termInd)) '_sig(warn_inds) = NaN;']);
    eval([char(p.Results.allVars(termInd)) '_coef(warn_inds) = NaN;']);
end
for termInd = 1:length(p.Results.modelVars)
    eval([char(p.Results.modelVars(termInd)) '_sig(warn_inds) = NaN;']);
    eval([char(p.Results.modelVars(termInd)) '_coef(warn_inds) = NaN;']);
end

rSqr(warn_inds) = NaN;

%plot signicance across populatioon
figure; 
subplot(1,2,1); hold on;
totalNumTerms = length(p.Results.allVars)+length(p.Results.modelVars); 
colors = hsv(totalNumTerms);
pCut = 0.05;
for termInd = 1:length(p.Results.allVars)
    eval(['plot(timeToPlot + tW, nanmean(' char(p.Results.allVars(termInd))...
        '_sig < pCut, 2), ''linewidth'', 2, ''color'', colors(termInd,:))']);
end
for termInd = 1:length(p.Results.modelVars)
    eval(['plot(timeToPlot + tW, nanmean(' char(p.Results.modelVars(termInd))...
        '_sig < pCut, 2), ''linewidth'', 2, ''color'', colors(numTerms+termInd,:))']);
end
xlim([timeToPlot(1)+tW timeToPlot(end)+tW])
lgTxt = [p.Results.allVars, p.Results.modelVars];
legend(lgTxt)
xlabel('Time from cue (ms)')
ylabel('Fraction of significant neurons (P < 0.05)')
ylim([0 1])
set(gca, 'tickdir', 'out')

subplot(1,2,2); hold on;
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
numRows = ceil(totalNumTerms/4);
colors = hsv(totalNumTerms);

figure; hold on;
for currVar = 1:length(p.Results.allVars)
    subplot(numRows,4,currVar); hold on; 
    eval([p.Results.allVars{currVar} '_coef = ' p.Results.allVars{currVar} '_coef(winInd,:);']);
    eval(['sigInds = find(' p.Results.allVars{currVar} '_sig(winInd,:) < pCut);']);
    eval(['histogram(' p.Results.allVars{currVar} '_coef, 20, ''FaceColor'', ''k'', ''FaceAlpha'', 0.3)']);
    eval(['histogram(' p.Results.allVars{currVar} '_coef(sigInds), 20, ''FaceColor'', colors(currVar,:))']);
    title(p.Results.allVars{currVar})
    set(gca, 'tickdir', 'out')
end
for currVar = 1:length(p.Results.modelVars)
    subplot(numRows,4,numTerms+currVar); hold on; 
    eval([p.Results.modelVars{currVar} '_coef = ' p.Results.modelVars{currVar} '_coef(winInd,:);']);
    eval(['sigInds = find(' p.Results.modelVars{currVar} '_sig(winInd,:) < pCut);']);
    eval(['histogram(' p.Results.modelVars{currVar} '_coef, 20, ''FaceColor'', ''k'', ''FaceAlpha'', 0.3)']);
    eval(['histogram(' p.Results.modelVars{currVar} '_coef(sigInds), 20, ''FaceColor'', colors(numTerms+currVar,:))']);
    title(p.Results.modelVars{currVar})
    set(gca, 'tickdir', 'out')
end

set(gcf, 'renderer', 'painters')


