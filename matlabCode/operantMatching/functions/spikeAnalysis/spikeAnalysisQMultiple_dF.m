function spikeAnalysisQMultiple_dF(xlFile, sheet, beh, varargin)


p = inputParser;
% default parameters if none given
p.addParameter('biasFlag',0);
p.addParameter('bernFlag', 1);
p.addParameter('modelName', 'sixParam_absPePeAN_bi')
p.addParameter('regNames', [{'sumQ'} {'diffQ'} {'confQ'} {'peBar'} {'pePe'}]);
p.addParameter('modelVars', [{'peBar'} {'pePe'}]);
p.parse(varargin{:});


% path info
[root,sep] = currComputer();

% load data from excel sheet
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
tb = 1.5;
tf = 5;
time = -1000*tb:1000*tf;
trialBeg = tb*1000;
CSoff = tb*1000 + 500;
smoothWin = 250;
respWin = 2000;

rBarFRbins = [];

sessionName = [];
rho = nan(1,length(cellList));
pVal = nan(1,length(cellList));

for currCell = 1:length(cellList)
    
    fprintf('Analyzing cell %d of %d \n', currCell, length(sessionList));
    
    if strcmp(sessionName, sessionList{currCell}) == 0
        sessionName = sessionList{currCell};
        [animalName, date] = strtok(sessionName, 'd'); 
        animalName = animalName(2:end);

        [animalName] = strtok(sessionName, 'd');
        animalName = animalName(2:end);

        if isstrprop(sessionName(end), 'alpha')
            sortedFolderLocation = [root animalName sep sessionName(1:end-1) sep 'sorted' sep 'session ' sessionName(end) sep];
        else
            sortedFolderLocation = [root animalName sep sessionName sep 'sorted' sep 'session' sep];
        end
        sortedFolder = dir(sortedFolderLocation);


        if intanFlag(currCell)
            if any(~cellfun(@isempty,strfind({sortedFolder.name},'_intan.mat'))) == 1
                sessionDataInd = ~cellfun(@isempty,strfind({sortedFolder.name},'_intan.mat'));
                load([sortedFolderLocation sortedFolder(sessionDataInd).name])
            else
                [sessionData] = generateSessionData_intan_operantMatching(sessionName);
            end
        else
%             if any(~cellfun(@isempty,strfind({sortedFolder.name},'_nL.mat'))) == 1
%                 sessionDataInd = ~cellfun(@isempty,strfind({sortedFolder.name},'_nL.mat'));
%                 load([sortedFolderLocation sortedFolder(sessionDataInd).name])
%             else
                [sessionData] = generateSessionData_nL_operantMatching(sessionName);
%             end
        end
        [s] = behAnalysisNoPlot_opMD(sessionName, 'revForFlag', revForFlag(currCell));
        prevRwd = [NaN s.allRewardsBinary(1:end-1)]';
        
        oTime = nan(1, length(s.responseInds));
        for currT = 1:length(s.responseInds)
            if isnan(sessionData(s.responseInds(currT)).rewardR)
                if length(sessionData(s.responseInds(currT)).licksL) > 1 & ...
                        sessionData(s.responseInds(currT)).licksL(2) - sessionData(s.responseInds(currT)).CSon <  respWin;
                    oTime(currT) = sessionData(s.responseInds(currT)).licksL(2) - sessionData(s.responseInds(currT)).CSon;
                else
                    oTime(currT) = sessionData(s.responseInds(currT)).licksL(1) - sessionData(s.responseInds(currT)).CSon;
                end
            else
                if length(sessionData(s.responseInds(currT)).licksR) > 1 & ...
                        sessionData(s.responseInds(currT)).licksR(2) - sessionData(s.responseInds(currT)).CSon <  respWin;
                    oTime(currT) = sessionData(s.responseInds(currT)).licksR(2) - sessionData(s.responseInds(currT)).CSon;
                else
                    oTime(currT) = sessionData(s.responseInds(currT)).licksR(1) - sessionData(s.responseInds(currT)).CSon;
                end
            end
        end
        eInds = [];
        if any(oTime > 2000)
            eInds = find(oTime > 2000)
            oTime(eInds) = [];
            s.responseInds(eInds) = [];
        end
                
            
        if p.Results.bernFlag
            modelPath = [root animalName sep animalName 'sorted' sep 'stan' sep 'bernoulli' sep p.Results.modelName sep animalName...
            beh '_' p.Results.modelName '.mat'];
        else
            modelPath = [root animalName sep animalName 'sorted' sep 'stan' sep p.Results.modelName sep animalName...
            beh '_' p.Results.modelName '.mat'];
        end
        t = generateStanModelTerms_opMD(p.Results.modelName, modelPath, sessionName, 1, revForFlag(currCell));
        sumQ = t.Q(:,1) + t.Q(:,2);
        diffQ = t.Q(:,1) - t.Q(:,2);
        confQ = [0; abs(diffQ)];
        pe = t.pe(1:end-1);

        for termInd = 1:length(p.Results.modelVars)
            eval([p.Results.modelVars{termInd} '= zscore(t.' p.Results.modelVars{termInd} ');']);
        end
        modelVars = [];
        for regInd = 1:length(p.Results.regNames)
            modelVars = [modelVars eval(p.Results.regNames{regInd})];
        end
        if ~isempty('eInds')
            modelVars(eInds, :) = [];
        end
    end


    %% Sort all spikes into a raster-able matrix

    spikeFields = fields(sessionData);
    cellInd = find(~cellfun(@isempty,strfind(spikeFields,cellList(currCell))));

    allTrial_spike = {};
    for k = 1:length(sessionData)
        if k == 1
            prevTrial_spike = [];
        else
            prevTrial_spikeInd = sessionData(k-1).(spikeFields{cellInd}) > sessionData(k-1).trialEnd-tb*1000;
            prevTrial_spike = sessionData(k-1).(spikeFields{cellInd})(prevTrial_spikeInd) - sessionData(k).CSon;
        end

        currTrial_spikeInd = sessionData(k).(spikeFields{cellInd}) < sessionData(k).CSon+tf*1000;
        currTrial_spike = sessionData(k).(spikeFields{cellInd})(currTrial_spikeInd) - sessionData(k).CSon;

        allTrial_spike{k} = [prevTrial_spike currTrial_spike];
    end

    % sometimes no licks/spikes are considered 1x0 and sometimes they are []
    % plotSpikeRaster does not place nicely with [] so this converts all empty indices to 1x0
    allTrial_spike(cellfun(@isempty,allTrial_spike)) = {zeros(1,0)}; 


    %% set time window and smoothing parameters, run analysis for all cells


    for i = 1:length(s.behSessionData)
        trialDurDiff(i) = (s.behSessionData(i).trialEnd - s.behSessionData(i).CSon)- tf*1000;
    end


    allTrial_spikeMatx = NaN(length(sessionData),length(time));
    for j = 1:length(allTrial_spike)
        tempSpike = allTrial_spike{j};
        tempSpike = tempSpike + tb*1000; % add this to pad time for SDF
        allTrial_spikeMatx(j,tempSpike) = 1;
        if trialDurDiff(j) < 0
            allTrial_spikeMatx(j, isnan(allTrial_spikeMatx(j, 1:end+trialDurDiff(j)))) = 0;  %converts within trial duration NaNs to 0's
        else
            allTrial_spikeMatx(j, isnan(allTrial_spikeMatx(j,:))) = 0;
        end
        if sum(allTrial_spikeMatx(j,:)) == 0     %if there is no spike data for this trial, don't count it
            allTrial_spikeMatx(j,:) = NaN;
        end
    end


%% find features of spike rate on each trial   
    preCSspikeCount = nan(1, length(s.responseInds));
    postCSspikeCount = nan(1, length(s.responseInds));
    postOspikeCount = nan(1, length(s.responseInds));
    maxFRcs = nan(1, length(s.responseInds));
    minFRcs = nan(1, length(s.responseInds));
 
    for j = 1:length(s.responseInds)
        if ~isempty(allTrial_spikeMatx(i,j))
            preCSspikeCount(j) = nansum(allTrial_spikeMatx(s.responseInds(j), 1:trialBeg));              %find total spikes before CS on
            postCSspikeCount(j) = nansum(allTrial_spikeMatx(s.responseInds(j), trialBeg:oTime(j)+trialBeg));
            postOspikeCount(j) = nansum(allTrial_spikeMatx(s.responseInds(j), oTime(j)+trialBeg:end));

            spikeTemp = fastsmooth(allTrial_spikeMatx(s.responseInds(j),:)*1000, smoothWin, 3);         %smooth raw spikes to find features of spike rate
            maxFRcs(j) = max(spikeTemp(trialBeg:CSoff));
            minFRcs(j) = min(spikeTemp(trialBeg:CSoff));
            if ~isnan(maxFRcs(j))
                maxFRtime(j) = find(spikeTemp == max(spikeTemp(trialBeg:CSoff)), 1);
            else
                maxFRtime(j) = NaN;
            end
        else
            preCSspikeCount(j) = 0;
            postCSspikeCount(j) = 0;
            maxFRcs(j) = 0;
            minFRcs(j) = 0;
            maxFRtime(j) = 0;
        end
    end
    
    preCSspikeCount = zscore(preCSspikeCount);
    postCSspikeCount = zscore(postCSspikeCount);
    maxFRcs = zscore(maxFRcs);
    minFRcs = zscore(minFRcs);

%% create models and put them into the structure
    mdl_postCS = fitlm([modelVars], postCSspikeCount);
    mdl_preCS = fitlm([modelVars], postOspikeCount);
    
    postCS_p(currCell,:) = mdl_postCS.Coefficients.pValue(2:end);
    postCS_rSqr(currCell,:) = [mdl_postCS.Rsquared.Ordinary mdl_postCS.Rsquared.Adjusted];
    preCS_p(currCell,:) = mdl_preCS.Coefficients.pValue(2:end);
    preCS_rSqr(currCell,:) = [mdl_preCS.Rsquared.Ordinary mdl_preCS.Rsquared.Adjusted];
    
    
    [rho(currCell) pVal(currCell)] = corr(modelVars, postCSspikeCount');
    if pVal(currCell) < 0.05
        figure; plot(postCSspikeCount, '-k', 'linewidth', 1.3)
        hold on; plot(modelVars, '-b', 'linewidth', 1.3)
    end
end

regNames = p.Results.regNames;
numRegs = length(regNames);

blue = [0 1 1];
purp = [0.7 0 1];
colors = [linspace(blue(1),purp(1),numRegs)', linspace(blue(2),purp(2),numRegs)', linspace(blue(3),purp(3),numRegs)'];

sigRates = [];
eitherSig = [];
figure;
for i = 1:numRegs
    subplot(2,3,i); hold on; title(regNames{i});
    
    tmpInds = find(postCS_p(:,i) >= 0.05 & preCS_p(:,i) >= 0.05);
    scatter(postCS_p(tmpInds,i), preCS_p(tmpInds,i), [], 'k', 'filled')
    
    tmpInds = find(postCS_p(:,i) < 0.05 | preCS_p(:,i) < 0.05);
    scatter(postCS_p(tmpInds,i), preCS_p(tmpInds,i), [], colors(i,:), 'filled')
    
    set(gca, 'TickDir', 'out')
    xlabel('CS firing rate P Values')
    ylabel('Pre-CS firing rate P Values')
    
    tmpInds_post = find(postCS_p(:,i) < 0.05);
    sigRatePost = length(tmpInds_post) / length(postCS_p);
    tmpInds_pre = find(preCS_p(:,i) < 0.05);
    sigRatePre = length(tmpInds_pre) / length(preCS_p);
    sigRates = [sigRates; sigRatePost sigRatePre];
    eitherSig = [eitherSig; length(unique([tmpInds_post; tmpInds_pre]))];
end

subplot(2,3,4)
colormap(gray)
bar(sigRates)
set(gca, 'Box', 'off', 'TickDir', 'out')
xticks([1:numRegs]); xticklabels(regNames)
ylabel('Significance Rate')
legend('CS firing rate', 'Pre-CS firing rate')

subplot(2,3,5)
scatter(postCS_rSqr(:,1), preCS_rSqr(:,1), 'm', 'filled')
xlabel('CS model R^2 Value')
ylabel('Pre-CS model R^2 Value')


sigCount_pre = preCS_p < 0.05;
sigCount_post = postCS_p < 0.05;
sigCount = sigCount_pre + sigCount_post;
sigCount(sigCount==2) = 1;
sigRates = sum(sigCount);
figure;
colormap(gray)
bar(sigRates/length(cellList))
xticks([1:numRegs]); xticklabels(regNames)
set(gca, 'Box', 'off', 'TickDir', 'out')
ylabel('Fraction of cells with significant regression coefficient')
ylim([0 1])


% tbl_maxCS = table(sessionList, cellList, postCS_rSqr(:,1), postCS_rSqr(:,2), postCS_t(:,1), postCS_t(:,2), postCS_t(:,3),...
%     postCS_t(:,4), postCS_t(:,5), postCS_t(:,6), 'VariableNames', {'Session' 'Cell' 'rSqr_ord' 'rSqr_adj'...
%     'Qsum' 'Qdiff' 'Qconf' 'Qchoice' 'rBar' 'beta'});
% 
% tbl_preCS = table(sessionList, cellList, preCS_rSqr(:,1), preCS_rSqr(:,2), preCS_t(:,1), preCS_t(:,2), preCS_t(:,3),...
%     preCS_t(:,4), preCS_t(:,5), preCS_t(:,6), 'VariableNames', {'Session' 'Cell' 'rSqr_ord' 'rSqr_adj'...
%     'Qsum' 'Qdiff' 'Qconf' 'Qchoice' 'rBar' 'beta'});

    