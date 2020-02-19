function [wS_rwdHx, lS_rwdHx, wS, lS] = wslsRwdCount_opMD(xlFile, sheet, category, varargin)

p = inputParser;
% default parameters if none given
p.addParameter('tMax', 10);
p.addParameter('revForFlag', 0);
p.addParameter('figFlag', 0);
p.parse(varargin{:});

[root, sep] = currComputer();
tMax = p.Results.tMax;

[~, dayList, ~] = xlsread(xlFile, sheet);
[~,col] = find(~cellfun(@isempty,strfind(dayList, category)) == 1);
dayList = dayList(2:end,col);
endInd = find(cellfun(@isempty,dayList),1);
if ~isempty(endInd)
    dayList = dayList(1:endInd-1,:);
end


combinedRwdHx = [];
combinedRwds = [];
combinedChangeChoice = [];

for seshInd = 1: length(dayList)
    sessionName = dayList{seshInd};
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
        if p.Results.revForFlag
            behSessionData = sessionData;
        end
    elseif p.Resulst.revForFlag                                    %otherwise generate the struct
        [behSessionData, ~] = generateSessionData_operantMatching(sessionName);
    else
        [behSessionData, ~, ~, ~] = generateSessionData_operantMatchingDecoupled(sessionName);
    end

    responseInds = find(~isnan([behSessionData.rewardTime])); % find CS+ trials with a response in the lick window
    allReward_R = [behSessionData(responseInds).rewardR]; 
    allReward_L = [behSessionData(responseInds).rewardL]; 
    allChoices = NaN(1,length(behSessionData(responseInds)));
    allChoices(~isnan(allReward_R)) = 1;
    allChoices(~isnan(allReward_L)) = -1;

    allReward_R(isnan(allReward_R)) = 0;
    allReward_L(isnan(allReward_L)) = 0;
    allChoice_R = double(allChoices == 1);
    allChoice_L = double(allChoices == -1);

    allRewards = zeros(1,length(allChoices));
    allRewards(logical(allReward_R)) = 1;
    allRewards(logical(allReward_L)) = -1;
    allRewardsBin = allRewards;
    allRewardsBin(allRewards == -1) = 1;

    rwdMatx = [];
    for j=1:tMax
        rwdMatx(j,:) = [nan(1,j) allRewardsBin(1:end-j)];
    end
    rwdHx = nansum(rwdMatx, 1);
    rwdHx(1) = 0;
        
    combinedRwdHx = [combinedRwdHx rwdHx(1:end-2)];
    combinedRwds = [combinedRwds allRewardsBin(2:end-1)];
    changeChoice = [abs(diff(allChoices)) > 0];
    combinedChangeChoice = [combinedChangeChoice changeChoice(2:end)];

end

wS = 1 - (sum(combinedChangeChoice(combinedRwds==1))/sum(combinedRwds==1));
lS = sum(combinedChangeChoice(combinedRwds==0))/sum(combinedRwds==0);

for i = 1:tMax+1
    tmpInds = logical(combinedRwdHx == i-1);
    wS_rwdHx(i) = 1 - (sum(combinedChangeChoice(combinedRwds==1 & tmpInds))/sum(combinedRwds==1 & tmpInds));
    lS_rwdHx(i) = sum(combinedChangeChoice(combinedRwds==0 & tmpInds))/sum(combinedRwds==0 & tmpInds);
end 


if p.Results.figFlag
    figure;
    subplot(1,2,1); hold on;
    scatter([0:tMax], wS_rwdHx, 'filled')
    xlim([-0.5 tMax+0.5])
    xlabel('number of rewards in last 10 trials')
    ylabel('probability')
    title('win-stay')
    set(gca, 'tickdir', 'out')
    subplot(1,2,2); hold on;
    scatter([0:tMax], lS_rwdHx, 'filled')
    xlim([-0.5 tMax+0.5])
    xlabel('number of rewards in last 10 trials')
    ylabel('probability')
    title('lose-shift')
    set(gca, 'tickdir', 'out')
    set(gcf, 'renderer', 'painters', 'position', [-1919 41 1920 963])
end

end