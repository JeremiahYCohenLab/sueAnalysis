function [glm_rwdLick, stayLickLat, switchLickLat, tMax, combinedITIlicks] = combineLinRegLickLat_opMD(xlFile, animal, category, varargin)

p = inputParser;
% default parameters if none given
p.addParameter('revForFlag',0)
p.addParameter('plotFlag', 1)
p.parse(varargin{:});

[root, sep] = currComputer();

[~, dayList, ~] = xlsread([root xlFile], animal);
[~,col] = find(~cellfun(@isempty,strfind(dayList, category)) == 1);
dayList = dayList(2:end,col);
endInd = find(cellfun(@isempty,dayList),1);
if ~isempty(endInd)
    dayList = dayList(1:endInd-1,:);
end
combinedRewardsMatx = [];
combinedLickLat = [];
stayLickLat = []; 
switchLickLat = [];
tMax = 12;
combinedITIlicks = [];


for i = 1: length(dayList)
    sessionName = dayList{i};
    [animalName, date] = strtok(sessionName, 'd'); 
    animalName = animalName(2:end);
    date = date(1:9);
    sessionFolder = ['m' animalName date];

    if isstrprop(sessionName(end), 'alpha')
        sessionDataPath = [root animalName sep sessionFolder sep 'sorted' sep 'session' sep sessionName(end) sep sessionName sep '_sessionData_behav.mat'];
    else
        sessionDataPath = [root animalName sep sessionFolder sep 'sorted' sep 'session' sep sessionName '_sessionData_behav.mat'];
    end

    if exist(sessionDataPath,'file')
        load(sessionDataPath)
        if p.Results.revForFlag
            behSessionData = sessionData;
        end
    else
        [behSessionData, ~] = generateSessionData_operantMatchingDecoupledRwdDelay(sessionName);
    end
    
    %%generate reward matrix for tMax trials
    responseInds = find(~isnan([behSessionData.respondTime])); % find CS+ trials with a response in the lick window
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

    
    rwdMatxTmp = [];
    for j = 1:tMax
        rwdMatxTmp(j,:) = [NaN(1,j) allRewards(1:end-j)];
    end

    combinedRewardsMatx = [combinedRewardsMatx NaN(tMax,100) rwdMatxTmp(:,1:end-1)];
    
    
    %% determine and plot lick latency distributions for each spout
    lickLat = [behSessionData(responseInds).respondTime] - [behSessionData(responseInds).CSon];
    indsR = find(allChoices == 1);
    indsL = find(allChoices == -1);
    lickLat_R = zscore(lickLat(indsR));
    lickLat_L = zscore(lickLat(indsL));
    lickLat = NaN(1, length(allChoices));
    lickLat(indsR) = lickLat_R;
    lickLat(indsL) = lickLat_L;
    
    combinedLickLat = [combinedLickLat NaN(1,100) lickLat(2:end)];
    
    %% determine lick latency for stay v switch trials
    changeChoice = [false abs(diff(allChoices)) > 0];
    stayLickLat = [stayLickLat lickLat(~changeChoice)]; 
    switchLickLat = [switchLickLat lickLat(changeChoice)];
    
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
glm_rwdLick = fitlm([combinedRewardsMatx]', combinedLickLat);

if p.Results.plotFlag
    figure; hold on
    relevInds = 2:tMax+1;
    coefVals = glm_rwdLick.Coefficients.Estimate(relevInds);
    CIbands = coefCI(glm_rwdLick);
    errorL = abs(coefVals - CIbands(relevInds,1));
    errorU = abs(coefVals - CIbands(relevInds,2));
    errorbar((1:tMax)+0.2,coefVals,errorL,errorU,'Color', [0.7 0 1],'linewidth',2)

    xlabel('Reward n Trials Back')
    ylabel('\beta Coefficient')
    xlim([0.5 tMax+0.5])

    suptitle([animal ' ' category])
end