function [glm_rwdLick, stayLickLat, switchLickLat, binSize, timeMax, combinedITIlicks] = combineLinRegLickLatTime_opMD(xlFile, animal, category, varargin)

p = inputParser;
% default parameters if none given
p.addParameter('revForFlag',0)
p.addParameter('plotFlag', 0)
p.addParameter('binSize', 10000);
p.addParameter('numBins', 9);
p.parse(varargin{:});

[root, sep] = currComputer();

[weights, dayList, ~] = xlsread(xlFile, animal);
[~,col] = find(~cellfun(@isempty,strfind(dayList, category)) == 1);
dayList = dayList(2:end,col);
endInd = find(cellfun(@isempty,dayList),1);
if ~isempty(endInd)
    dayList = dayList(1:endInd-1,:);
end

binSize = p.Results.binSize;
timeMax = binSize * p.Results.numBins + 1000;
timeBinEdges = [1000:binSize:timeMax];  %no trials shorter than 1s between outcome and CS on
tMax = length(timeBinEdges) - 1;
rwdMatx = [];
noRwdMatx = [];
combinedLickLat = [];
stayLickLat = [];
switchLickLat = [];
combinedITIlicks = [];
combinedTimeInSesh = [];
combinedChangeChoice = [];
latAndTime = [];


for i = 1: length(dayList)
    sessionName = dayList{i};
    [animalName, date] = strtok(sessionName, 'd'); 
    animalName = animalName(2:end);
    date = date(1:9);
    sessionFolder = ['m' animalName date];

    if isstrprop(sessionName(end), 'alpha')
        sessionDataPath = [root animalName sep sessionFolder sep 'sorted' sep 'session' sep sessionName(end) sep sessionName '_sessionData_behav.mat'];
    else
        sessionDataPath = [root animalName sep sessionFolder sep 'sorted' sep 'session' sep sessionName '_sessionData_behav.mat'];
    end

    if exist(sessionDataPath,'file')
        load(sessionDataPath)
        if p.Results.revForFlag
            behSessionData = sessionData;
        end
    else
        [behSessionData, ~] = generateSessionData_operantMatchingDecoupled(sessionName);
    end
    
    %%generate reward matrix for tMax trials
    responseInds = find(~isnan([behSessionData.rewardTime])); % find CS+ trials with a response in the lick window
    allReward_R = [behSessionData(responseInds).rewardR]; 
    allReward_L = [behSessionData(responseInds).rewardL]; 
    allChoices = NaN(1,length(behSessionData(responseInds)));
    allChoices(~isnan(allReward_R)) = 1;
    allChoices(~isnan(allReward_L)) = -1;
    allReward_R(isnan(allReward_R)) = 0;
    allReward_L(isnan(allReward_L)) = 0;
    allRewards = zeros(1,length(allChoices));
    allRewards(logical(allReward_R)) = 1;
    allRewards(logical(allReward_L)) = 1;
    timeInSesh = ([behSessionData(responseInds).CSon] - behSessionData(1).CSon) / 1000;
    changeChoice = [0 abs(diff(allChoices)) > 0];

    %create binned outcome matrices
    rwdTmpMatx = zeros(tMax, length(responseInds));     %initialize matrices for number of response trials x number of time bins
    noRwdTmpMatx = zeros(tMax, length(responseInds));
    for j = 2:length(responseInds)          
        k = 1;
        %find time between "current" choice and previous rewards, up to timeMax in the past 
        timeTmp = [];
        timeTmpNoRwd = [];
        while j-k > 0 & behSessionData(responseInds(j)).rewardTime - behSessionData(responseInds(j-k)).rewardTime < timeMax
            if behSessionData(responseInds(j-k)).rewardL == 1 || behSessionData(responseInds(j-k)).rewardR == 1
                timeTmp = [timeTmp (behSessionData(responseInds(j)).rewardTime - behSessionData(responseInds(j-k)).rewardTime)];
            end
            if behSessionData(responseInds(j-k)).rewardL == 0 || behSessionData(responseInds(j-k)).rewardR == 0
                timeTmpNoRwd = [timeTmpNoRwd (behSessionData(responseInds(j)).rewardTime - behSessionData(responseInds(j-k)).rewardTime)];
            end
            k = k + 1;
        end
        %bin outcome times and use to fill matrices
        if ~isempty(timeTmp)
            binnedRwds = discretize(timeTmp,timeBinEdges);
            for k = 1:tMax
                if ~isempty(find(binnedRwds == k))
                    rwdTmpMatx(k,j) = sum(binnedRwds == k);
                end
            end
        end
        if ~isempty(timeTmpNoRwd)
            binnedNoRwds = discretize(timeTmpNoRwd,timeBinEdges);
            for k = 1:tMax
                if ~isempty(find(binnedNoRwds == k))
                    noRwdTmpMatx(k,j) = sum(binnedNoRwds == k);
                end
            end
        end
    end
    
    %fill in NaNs at beginning of session
    j = 2;
    while behSessionData(responseInds(j)).rewardTime - behSessionData(responseInds(1)).rewardTime < timeMax
        tmpDiff = behSessionData(responseInds(j)).rewardTime - behSessionData(responseInds(1)).rewardTime;
        binnedDiff = discretize(tmpDiff, timeBinEdges);
        rwdTmpMatx(binnedDiff:tMax,j) = NaN;
        noRwdTmpMatx(binnedDiff:tMax,j) = NaN;
        j = j+1;
    end
    %concatenate temp matrix with combined matrix
    rwdTmpMatx(:,1) = NaN;
    rwdMatx = [rwdMatx NaN(length(timeBinEdges)-1, 100) rwdTmpMatx];
    noRwdTmpMatx(:,1) = NaN;
    noRwdMatx = [noRwdMatx NaN(length(timeBinEdges)-1, 100) noRwdTmpMatx];
    combinedChangeChoice = [combinedChangeChoice NaN(1, 100) changeChoice];
    combinedTimeInSesh = [combinedTimeInSesh NaN(1, 100) timeInSesh];
    
    
    %% determine lick latency distributions for each spout
    lickLat = [behSessionData(responseInds).rewardTime] - [behSessionData(responseInds).CSon];
    indsR = find(allChoices == 1);
    indsL = find(allChoices == -1);
    lickLat_R = zscore(lickLat(indsR));
    lickLat_L = zscore(lickLat(indsL));
    lickLatZ = NaN(1, length(allChoices));
    lickLatZ(indsR) = lickLat_R;
    lickLatZ(indsL) = lickLat_L;
    
    combinedLickLat = [combinedLickLat NaN(1,101) lickLatZ(2:end)];
    
    %% determine lick latency for stay v switch trials
    changeChoice = [false abs(diff(allChoices)) > 0];
    stayLickLat = [stayLickLat lickLatZ(~changeChoice)]; 
    switchLickLat = [switchLickLat lickLatZ(changeChoice)];
    
    %% combine lick latency and time in session
    tmp = [lickLatZ; timeInSesh];
    latAndTime = [latAndTime tmp];
    
    %% determine iti lick rates on each trial
    itiLicks = zeros(1,length(behSessionData));
    for currT = 1:length(behSessionData)
        itiLicksTmp = [behSessionData(currT).licksL(behSessionData(currT).licksL > behSessionData(currT).CSon + 2000)...
                        behSessionData(currT).licksR(behSessionData(currT).licksR > behSessionData(currT).CSon + 2000)];
        if ~isempty(itiLicksTmp)
            itiLicks(currT) = length(itiLicksTmp) / (behSessionData(currT).trialEnd - behSessionData(currT).CSon - 2000);
        end
    end
    combinedITIlicks = [combinedITIlicks itiLicks];
end

%linear regression model
glm_rwdLick = fitlm([rwdMatx' noRwdMatx'], combinedLickLat);
glm_rwdLickAll = fitlm([rwdMatx' noRwdMatx' combinedTimeInSesh' combinedChangeChoice'], combinedLickLat);

if p.Results.plotFlag
    figure; suptitle([animal ' ' category])
    colors = cool(5);
    subplot(1,3,1); hold on
    relevInds = 2:tMax+1;
    coefVals = glm_rwdLickAll.Coefficients.Estimate(relevInds);
    CIbands = coefCI(glm_rwdLickAll);
    errorL = abs(coefVals - CIbands(relevInds,1));
    errorU = abs(coefVals - CIbands(relevInds,2));
    errorbar(((1:tMax)*binSize/1000),coefVals,errorL,errorU,'Color', colors(1,:),'linewidth',2)
    
    relevInds = tMax+2:2*tMax+1;
    coefVals = glm_rwdLickAll.Coefficients.Estimate(relevInds);
    CIbands = coefCI(glm_rwdLickAll);
    errorL = abs(coefVals - CIbands(relevInds,1));
    errorU = abs(coefVals - CIbands(relevInds,2));
    errorbar(((1:tMax)*binSize/1000),coefVals,errorL,errorU,'Color', colors(4,:),'linewidth',2)
    
    legend('reward', 'no reward')
    xlabel('reward n seconds back')
    ylabel('\beta Coefficient')
    xlim([0 (tMax*binSize/1000 + 5)])
    set(gca, 'tickdir', 'out')
        
    subplot(1,3,2); hold on
    histogram(stayLickLat,20,'FaceColor', colors(2,:), 'Normalization', 'Probability')
    histogram(switchLickLat,20,'FaceColor', colors(5,:), 'Normalization', 'Probability')
    ylabel('probability')
    xlabel('lick latency (ms)')
    legend('stay', 'switch')
    set(gca, 'tickdir', 'out')
    
    subplot(1,3,3); hold on
    numBins = 6;
    sortInds = discretize(latAndTime(2,:), numBins);
    for currBin = 1:numBins
        tmp = latAndTime(1, sortInds==currBin);
        meanLat(currBin) = mean(tmp);
        semLat(currBin) = std(tmp) / sqrt(length(tmp));
    end
    errorbar([1:numBins], meanLat, semLat, 'Color', colors(3,:),'linewidth',2);
    xlim([0 numBins+1]);
    xticks([])
    xlabel('time in session')
    ylabel('z-scored lick latency')
    set(gca, 'tickdir', 'out')
    set(gcf, 'renderer', 'painters', 'position', [-1704 409 1552 420])
    
end