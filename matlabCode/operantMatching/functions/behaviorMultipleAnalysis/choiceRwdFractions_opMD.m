function [choiceRwdFractions] = choiceRwdFractions_opMD(xlFile, sheet, category)

[root, sep] = currComputer();
[~, dayList, ~] = xlsread(xlFile, sheet);
[~,col] = find(~cellfun(@isempty,strfind(dayList, category)) == 1);
dayList = dayList(2:end,col);

endInd = find(cellfun(@isempty,dayList),1);
if ~isempty(endInd)
    dayList = dayList(1:endInd-1,:);
end
choiceFractions = [];
rwdFractions = [];
choiceRwdFractions = [];

for i = 1: length(dayList)
    sessionName = dayList{i};
    [animalName, date] = strtok(sessionName, 'd'); 
    animalName = animalName(2:end);
    date = date(1:9);
    sessionFolder = ['m' animalName date];

    if isstrprop(sessionName(end), 'alpha')
        sessionDataPath = [root animalName sep sessionFolder sep 'sorted' sep 'session ' sessionName(end) sep sessionName '_sessionData_behav.mat'];
    else
        sessionDataPath = [root animalName sep sessionFolder sep 'sorted' sep 'session' sep sessionName '_sessionData_behav.mat'];
    end

    if exist(sessionDataPath,'file')
        load(sessionDataPath)
    else
        [behSessionData, blockSwitch, blockSwitchR, blockSwitchR] = generateSessionData_operantMatchingDecoupled(sessionName);
    end

    %% Break session down into CS+ trials where animal responded

    responseInds = find(~isnan([behSessionData.rewardTime])); % find CS+ trials with a response in the lick window
    omitInds = isnan([behSessionData.rewardTime]); 
    
    blockSwitchR = [blockSwitchR length(behSessionData)];
    tempBlockSwitchR = blockSwitchR;
    for j = 2:length(blockSwitchR)
        subVal = sum(omitInds(tempBlockSwitchR(j-1):tempBlockSwitchR(j)));
        blockSwitchR(j:end) = blockSwitchR(j:end) - subVal;
    end

    allReward_R = [behSessionData(responseInds).rewardR]; 
    allReward_L = [behSessionData(responseInds).rewardL]; 
    rewProb_R = [behSessionData(responseInds).rewardProbR]; 
    rewProb_L = [behSessionData(responseInds).rewardProbL]; 
    allChoices = NaN(1,length(behSessionData(responseInds)));
    allChoices(~isnan(allReward_R)) = 1;
    allChoices(~isnan(allReward_L)) = -1;

    allReward_R(isnan(allReward_R)) = 0;
    allReward_L(isnan(allReward_L)) = 0;
    allChoice_R = double(allChoices == 1);
    allChoice_L = double(allChoices == -1);

    allRewards = zeros(1,length(allChoices));
    allRewards(logical(allReward_R)) = 1;
    allRewards(logical(allReward_L)) = 1;
    
    
    for j = 2:length(blockSwitchR)
        choiceFractionTemp = sum(allChoice_R((blockSwitchR(j-1)+1):blockSwitchR(j))) / (blockSwitchR(j)-blockSwitchR(j-1));
        choiceFractions = [choiceFractions; choiceFractionTemp]; 
        rwdFractionTemp = sum(allReward_R((blockSwitchR(j-1)+1):blockSwitchR(j))) / sum(allRewards((blockSwitchR(j-1)+1):blockSwitchR(j)));
        rwdFractions = [rwdFractions; rwdFractionTemp];
    end
end
rmvInds = [find(isnan(choiceFractions))' find(isnan(rwdFractions))'];
choiceFractions(rmvInds) = [];
rwdFractions(rmvInds) = [];
choiceRwdFractions = [rwdFractions choiceFractions];