function [mdlStruct] = linRegSpikes_dF(xlFile, sheet, varargin)

p = inputParser;
% default parameters if none given
p.addParameter('trialFlag', 1);
p.addParameter('figFlag', 1)
p.addParameter('plotFlag', 1)
p.addParameter('intanFlag', 0)
p.parse(varargin{:});

[root, sep] = currComputer();

[revForFlagList, sessionCellList, ~] = xlsread(xlFile, sheet);
cellList = sessionCellList(2:end, 1);
sessionList = sessionCellList(2:end, 2);
revForFlagList = sessionCellList(:,1);

timeMax = 181000;
binSize = 30000;
timeBinEdges = [1000:binSize:timeMax];  %no trials shorter than 1s between outcome and CS on
tMax = length(timeBinEdges) - 1;

for i = 1:length(sessionList)
    
    fprintf('Analyzing cell %d of %d \n', i, length(sessionList))
    
    %get spike information
    spikeStruct = [];
    [spikeStruct, behSessionData] = spikeProps_opMD(sessionList{i}, cellList{i});
    if p.Results.intanFlag
        [s] = behAnalysisNoPlot_opMD(sessionList{i}, 'revForFlag', revForFlagList{i});
        behSessionData = s.behSessionData;
    end
        
    %load behavioral data
    [animalName, date] = strtok(sessionList{i}, 'd'); 
    animalName = animalName(2:end);
    date = date(1:9);
    sessionFolder = ['m' animalName date];


    %create arrays for choices and rewards
    responseInds = find(~isnan([behSessionData.rewardTime])); % find CS+ trials with a response in the lick window
    allReward_R = [behSessionData(responseInds).rewardR]; 
    allReward_L = [behSessionData(responseInds).rewardL]; 
    allChoices = NaN(1,length(behSessionData(responseInds)));
    allChoices(~isnan(allReward_R)) = 1;
    allChoices(~isnan(allReward_L)) = -1;
    allChoice_R = double(allChoices == 1);
    allChoice_L = double(allChoices == -1);

    allReward_R(isnan(allReward_R)) = 0;
    allReward_L(isnan(allReward_L)) = 0;
    allRewards = zeros(1,length(allChoices));
    allRewards(logical(allReward_R)) = 1;
    allRewards(logical(allReward_L)) = 1;
    
    allNoRewards = allChoices;
    allNoRewards(logical(allReward_R)) = 0;
    allNoRewards(logical(allReward_L)) = 0;
    
    %create outcome matrices for trialwise analysis
    rwdTmpMatxTrial = [];
    noRwdTmpMatxTrial = [];
    for j = 1:tMax
        rwdTmpMatxTrial(j,:) = [NaN(1,j) allRewards(1:end-j)];
        noRwdTmpMatxTrial(j,:) = [NaN(1,j) allNoRewards(1:end-j)];
    end

    %create binned outcome matrices for timewise analysis
    rwdTmpMatx = NaN(tMax, length(responseInds));     %initialize matrices for number of response trials x number of time bins
    noRwdTmpMatx = NaN(tMax, length(responseInds));
    for j = 2:length(responseInds)          
        k = 1;
        %find time between "current" choice and previous rewards, up to timeMax in the past 
        timeTmpL = []; timeTmpR = []; nTimeTmpL = []; nTimeTmpR = [];
        while j-k > 0 & behSessionData(responseInds(j)).rewardTime - behSessionData(responseInds(j-k)).rewardTime < timeMax
            if behSessionData(responseInds(j-k)).rewardL == 1
                timeTmpL = [timeTmpL (behSessionData(responseInds(j)).rewardTime - behSessionData(responseInds(j-k)).rewardTime)];
            end
            if behSessionData(responseInds(j-k)).rewardR == 1
                timeTmpR = [timeTmpR (behSessionData(responseInds(j)).rewardTime - behSessionData(responseInds(j-k)).rewardTime)];
            end
            if behSessionData(responseInds(j-k)).rewardL == 0
                nTimeTmpL = [nTimeTmpL (behSessionData(responseInds(j)).rewardTime - behSessionData(responseInds(j-k)).rewardTime)];
            end
            if behSessionData(responseInds(j-k)).rewardR == 0
                nTtimeTmpR = [nTimeTmpR (behSessionData(responseInds(j)).rewardTime - behSessionData(responseInds(j-k)).rewardTime)];
            end
            k = k + 1;
        end
        %bin outcome times and use to fill matrices
        if ~isempty(timeTmpL)
            binnedRwds = discretize(timeTmpL,timeBinEdges);
            for k = 1:tMax
                if ~isempty(binnedRwds == k)
                    rwdTmpMatx(k,j) = sum(binnedRwds == k);
                else
                    rwdTmpMatx(k,j) = 0;
                end
            end
        end
        if ~isempty(timeTmpR)
            binnedRwds = discretize(timeTmpR,timeBinEdges);
            for k = 1:tMax
                if ~isempty(binnedRwds == k) & isnan(rwdTmpMatx(k,j))
                    rwdTmpMatx(k,j) = sum(binnedRwds == k);
                elseif ~isempty(binnedRwds == k) & ~isnan(rwdTmpMatx(k,j))
                    rwdTmpMatx(k,j) = rwdTmpMatx(k,j) + sum(binnedRwds == k);
                else
                    rwdTmpMatx(k,j) = 0;
                end
            end
        end
        if isempty(timeTmpL) & isempty(timeTmpR)
            rwdTmpMatx(:,j) = 0;
        end
        if ~isempty(nTimeTmpL)
            binnedNoRwds = discretize(nTimeTmpL,timeBinEdges);
            for k = 1:tMax
                if ~isempty(binnedNoRwds == k)
                    noRwdTmpMatx(k,j) = sum(binnedNoRwds == k);
                else
                    noRwdTmpMatx(k,j) = 0;
                end
            end
        end
        if ~isempty(nTimeTmpR)
            binnedNoRwds = discretize(nTimeTmpR,timeBinEdges);
            for k = 1:tMax
                if ~isempty(binnedNoRwds == k) & isnan(noRwdTmpMatx(k,j))
                    noRwdTmpMatx(k,j) = sum(binnedNoRwds == k);
                elseif ~isempty(binnedNoRwds == k) & ~isnan(noRwdTmpMatx(k,j))
                    noRwdTmpMatx(k,j) = noRwdTmpMatx(k,j) + sum(binnedNoRwds == k);
                else
                    noRwdTmpMatx(k,j) = 0;
                end
            end
        end
        if isempty(nTimeTmpL) & isempty(nTimeTmpR)
            noRwdTmpMatx(:,j) = 0;
        end
    end
    
    % linear regression models
    time_preCS = fitlm([rwdTmpMatx' noRwdTmpMatx'], spikeStruct.preCScount);
    time_postCS = fitlm([rwdTmpMatx' noRwdTmpMatx'], spikeStruct.postCScount);
    time_postCSrate = fitlm([rwdTmpMatx' noRwdTmpMatx'], spikeStruct.maxCSrate);
    trial_preCS = fitlm([rwdTmpMatxTrial' noRwdTmpMatxTrial'], spikeStruct.preCScount);
    trial_postCS = fitlm([rwdTmpMatxTrial' noRwdTmpMatxTrial'], spikeStruct.postCScount);
    trial_postCSrate = fitlm([rwdTmpMatxTrial' noRwdTmpMatxTrial'], spikeStruct.maxCSrate);
    
    tmp = [sessionList{i} '_' cellList{i}];
    mdlStruct.(tmp).time_preCS = time_preCS;
    mdlStruct.(tmp).time_postCS = time_postCS;
    mdlStruct.(tmp).time_postCSrate = time_postCSrate;
    mdlStruct.(tmp).trial_preCS = trial_preCS;
    mdlStruct.(tmp).trial_postCS = trial_postCS;
    mdlStruct.(tmp).trial_postCSrate = trial_postCSrate;
    
    
    
    
    if p.Results.plotFlag
        % plot beta coefficients from models
        figure; 
        subplot(2,3,1); title('pre CS'); hold on;
        relevInds = 2:tMax+1;
        coefVals = time_preCS.Coefficients.Estimate(relevInds);
        CIbands = coefCI(time_preCS);
        errorL = abs(coefVals - CIbands(relevInds,1));
        errorU = abs(coefVals - CIbands(relevInds,2));
        errorbar(((1:tMax)*binSize/1000),coefVals,errorL,errorU,'Color', [0.7 0 1],'linewidth',2)
        relevInds = tMax+2:tMax*2+1;
        coefVals = time_preCS.Coefficients.Estimate(relevInds);
        CIbands = coefCI(time_preCS);
        errorL = abs(coefVals - CIbands(relevInds,1));
        errorU = abs(coefVals - CIbands(relevInds,2));
        errorbar(((1:tMax)*binSize/1000),coefVals,errorL,errorU,'-b', 'linewidth',2)
        xlabel('Outcome n seconds back')
        ylabel('\beta Coefficient')
        xlim([0 tMax*binSize/1000 + binSize/1000])

        subplot(2,3,2); title('post CS'); hold on;
        relevInds = 2:tMax+1;
        coefVals = time_postCS.Coefficients.Estimate(relevInds);
        CIbands = coefCI(time_postCS);
        errorL = abs(coefVals - CIbands(relevInds,1));
        errorU = abs(coefVals - CIbands(relevInds,2));
        errorbar(((1:tMax)*binSize/1000),coefVals,errorL,errorU,'Color', [0.7 0 1],'linewidth',2)
        relevInds = tMax+2:tMax*2+1;
        coefVals = time_postCS.Coefficients.Estimate(relevInds);
        CIbands = coefCI(time_postCS);
        errorL = abs(coefVals - CIbands(relevInds,1));
        errorU = abs(coefVals - CIbands(relevInds,2));
        errorbar(((1:tMax)*binSize/1000),coefVals,errorL,errorU,'-b', 'linewidth',2)
        xlabel('Outcome n seconds back')
        ylabel('\beta Coefficient')
        xlim([0 tMax*binSize/1000 + binSize/1000])

        subplot(2,3,3); title('post CS rate'); hold on;
        relevInds = 2:tMax+1;
        coefVals = time_postCSrate.Coefficients.Estimate(relevInds);
        CIbands = coefCI(time_postCSrate);
        errorL = abs(coefVals - CIbands(relevInds,1));
        errorU = abs(coefVals - CIbands(relevInds,2));
        errorbar(((1:tMax)*binSize/1000),coefVals,errorL,errorU,'Color', [0.7 0 1],'linewidth',2)
        relevInds = tMax+2:tMax*2+1;
        coefVals = time_postCSrate.Coefficients.Estimate(relevInds);
        CIbands = coefCI(time_postCSrate);
        errorL = abs(coefVals - CIbands(relevInds,1));
        errorU = abs(coefVals - CIbands(relevInds,2));
        errorbar(((1:tMax)*binSize/1000),coefVals,errorL,errorU,'-b', 'linewidth',2)
        xlabel('Outcome n seconds back')
        ylabel('\beta Coefficient')
        xlim([0 tMax*binSize/1000 + binSize/1000])

        subplot(2,3,4); title('pre CS'); hold on;
        relevInds = 2:tMax+1;
        coefVals = trial_preCS.Coefficients.Estimate(relevInds);
        CIbands = coefCI(trial_preCS);
        errorL = abs(coefVals - CIbands(relevInds,1));
        errorU = abs(coefVals - CIbands(relevInds,2));
        errorbar([1:tMax],coefVals,errorL,errorU,'Color', [0.7 0 1],'linewidth',2)
        relevInds = tMax+2:tMax*2+1;
        coefVals = trial_preCS.Coefficients.Estimate(relevInds);
        CIbands = coefCI(trial_preCS);
        errorL = abs(coefVals - CIbands(relevInds,1));
        errorU = abs(coefVals - CIbands(relevInds,2));
        errorbar([1:tMax],coefVals,errorL,errorU,'-b', 'linewidth',2)
        xlabel('Outcome n trials back')
        ylabel('\beta Coefficient')
        xlim([0 tMax+1]);

        subplot(2,3,5); title('post CS'); hold on;
        relevInds = 2:tMax+1;
        coefVals = trial_postCS.Coefficients.Estimate(relevInds);
        CIbands = coefCI(trial_postCS);
        errorL = abs(coefVals - CIbands(relevInds,1));
        errorU = abs(coefVals - CIbands(relevInds,2));
        errorbar([1:tMax],coefVals,errorL,errorU,'Color', [0.7 0 1],'linewidth',2)
        relevInds = tMax+2:tMax*2+1;
        coefVals = trial_postCS.Coefficients.Estimate(relevInds);
        CIbands = coefCI(trial_postCS);
        errorL = abs(coefVals - CIbands(relevInds,1));
        errorU = abs(coefVals - CIbands(relevInds,2));
        errorbar([1:tMax],coefVals,errorL,errorU,'-b', 'linewidth',2)
        xlabel('Outcome n trials back')
        ylabel('\beta Coefficient')
        xlim([0 tMax+1]);

        subplot(2,3,6); title('post CS rate'); hold on;
        relevInds = 2:tMax+1;
        coefVals = trial_postCSrate.Coefficients.Estimate(relevInds);
        CIbands = coefCI(trial_postCSrate);
        errorL = abs(coefVals - CIbands(relevInds,1));
        errorU = abs(coefVals - CIbands(relevInds,2));
        errorbar([1:tMax],coefVals,errorL,errorU,'Color', [0.7 0 1],'linewidth',2)
        relevInds = tMax+2:tMax*2+1;
        coefVals = trial_postCSrate.Coefficients.Estimate(relevInds);
        CIbands = coefCI(trial_postCSrate);
        errorL = abs(coefVals - CIbands(relevInds,1));
        errorU = abs(coefVals - CIbands(relevInds,2));
        errorbar([1:tMax],coefVals,errorL,errorU,'-b', 'linewidth',2)
        xlabel('Outcome n trials back')
        ylabel('\beta Coefficient')
        xlim([0 tMax+1]);

        legend('Reward', 'No Reward')    
        set(0, 'DefaulttextInterpreter', 'none')
        suptitle([sessionList{i} ' ' cellList{i}])
        set(gcf, 'Position', [-1596 335 1271 532])
        set(0, 'DefaulttextInterpreter', 'tex')
    end
    
end

