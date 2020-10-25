function combine = combineHist_opMD(xlFile, animal, category, varargin)
%task and model parameters
p = inputParser;
% default parameters if none given
p.addParameter('numBins', 10)
p.addParameter('plotFlag', 1)
p.addParameter('maxTrials', 200)
p.addParameter('revForFlag', 0)

p.parse(varargin{:})

[root, sep] = currComputer();
workbookFile = [root xlFile];
[~, dayList, ~] = xlsread(workbookFile, animal);
col = contains(dayList(1,:),category);
dayList = dayList(2:end,col);
endInd = find(cellfun(@isempty,dayList),1);
if ~isempty(endInd)
    dayList = dayList(1:endInd-1,:);
end
 
%%
combine = struct;

%%

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
        load(sessionDataPath);
    else
        [behSessionData, ~, ~, ~] = generateSessionData_operantMatchingDecoupledRwdDelay(sessionName);
    end
    behSessionData = behSessionData(1:min(p.Results.maxTrials, length(behSessionData)));
 %% 
    responseInds = find(~isnan([behSessionData.rewardTime])); % find CS+ trials with a response in the lick window
    omitInds = isnan([behSessionData.rewardTime]); 
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
    
    allNoRewards = allChoices;
    allNoRewards(logical(allReward_R)) = 0;
    allNoRewards(logical(allReward_L)) = 0;
    %% p(stay|rwd)
    rwd_stay = find(allChoices(2:end).* allChoices(1:end-1) > 0 & allRewards(1:end-1) ~= 0)+1;
    rwd = find(allRewards(1:end-1) ~= 0);
    combine(i).pstRwd = length(rwd_stay)/length(rwd); 
    %% p(swtich|nrwd)
    nrwd_switch = find(allChoices(2:end).* allChoices(1:end-1) < 0 & allRewards(1:end-1) == 0)+1;
    nrwd = find(allRewards(1:end-1) == 0);
    combine(i).pswNrwd = length(nrwd_switch)/length(nrwd);    
    %% p(stay|switch,rwd)
    switch_rwd_stay = find(allChoices(3:end).*allChoices(2:end-1) > 0 & allChoices(1:end-2).*allChoices(2:end-1) < 0 & allRewards(2:end-1) ~= 0);
    switch_rwd = find(allChoices(1:end-2).*allChoices(2:end-1) < 0 & allRewards(2:end-1) ~= 0);
    if ~isempty(switch_rwd)
        combine(i).pstaySwRwd = length(switch_rwd_stay)/length(switch_rwd); 
    else
        combine(i).pstaySwRwd = NaN;
    end
%     %% hmm
%     % explore portion
%     combine(i).rePo = length(find(states == 1))/length(states);
%     % transtion to explore
%     combine(i).itToRe = 0.5 * (trans_fit(2,1) + trans_fit(3,1));
%     % transition to exploit
%     combine(i).reToIt = 1 - trans_fit(1,1);
%     % explore length
%     erNo = length(find(states(2:end) == 1 & states(1:end - 1) ~= 1));
%     if states(1) == 1
%         erNo = erNo + 1;
%     end   
%     combine(i).reLen = length(find(states == 1))/erNo;
end
