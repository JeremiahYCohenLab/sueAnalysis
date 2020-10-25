function [combinedChoices, combinedRewards,combinedSwitches] = combineBeh_opMD(xlFile, animal, category)

[root, sep] = currComputer();

[weights, dayList, ~] = xlsread(xlFile, animal);
[~,col] = find(~cellfun(@isempty,strfind(dayList, category)) == 1);
dayList = dayList(2:end,col);
endInd = find(cellfun(@isempty,dayList),1);
if ~isempty(endInd)
    dayList = dayList(1:endInd-1,:);
end
combinedChoices = []; 
combinedRewards = [];
combinedSwitches = [];

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
        [behSessionData, blockSwitch, blockSwitchL, blockSwitchR] = generateSessionData_operantMatchingDecoupled(sessionName);
    end
    
    
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
    
    combinedChoices = [combinedChoices allChoices];
    combinedRewards = [combinedRewards allRewards];
    
    if i == 1
        combinedSwitches = blockSwitch;
        combinedProbsR = allProbsR;
        combinedProbsL = allProbsL;
    else
        blockSwitchTemp = blockSwitch + combinedSwitches(end);
        combinedSwitches = [combinedSwitches blockSwitchTemp];
        combinedProbsR = [combinedProbsR allProbsR];
        combinedProbsL = [combinedProbsR allProbsL];
    end
    
end





