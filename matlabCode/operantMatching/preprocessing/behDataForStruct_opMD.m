function [s] = behDataForStruct_opMD(sessionName, varargin)

p = inputParser;
% default parameters if none given
p.addParameter('revForFlag',0)
p.parse(varargin{:});


[root, sep] = currComputer();

[animalName, date] = strtok(sessionName, 'd'); 
animalName = animalName(2:end);
date = date(1:9);
sessionFolder = ['m' animalName date];

if isstrprop(sessionName(end), 'alpha')
    sessionDataPath = [root animalName sep sessionFolder sep 'sorted' sep 'session ' sessionName(end) sep sessionName '_sessionData_behav.mat'];
else
    sessionDataPath = [root animalName sep sessionFolder sep 'sorted' sep 'session' sep sessionName '_sessionData_behav.mat'];
end

if p.Results.revForFlag
    if exist(sessionDataPath,'file')
        load(sessionDataPath)
        if ~exist('behSessionData')
            behSessionData = sessionData;
        end
    else
        [behSessionData, blockSwitch, blockProbs] = generateSessionData_behav_operantMatching(sessionName);
    end
else
    if exist(sessionDataPath,'file')
        load(sessionDataPath)
    else
        [behSessionData, blockSwitch, blockSwitchL, blockSwitchR] = generateSessionData_operantMatchingDecoupled(sessionName);
    end
end

%% Break session down into CS+ trials where animal responded

responseInds = find(~isnan([behSessionData.rewardTime])); % find CS+ trials with a response in the lick window
omitInds = isnan([behSessionData.rewardTime]); 

origBlockSwitch = blockSwitch;
tempBlockSwitch = blockSwitch;
for i = 2:length(blockSwitch)
    subVal = sum(omitInds(tempBlockSwitch(i-1):tempBlockSwitch(i)));
    blockSwitch(i:end) = blockSwitch(i:end) - subVal;
end

allReward_R = [behSessionData(responseInds).rewardR]; 
allReward_L = [behSessionData(responseInds).rewardL]; 
allChoices = NaN(1,length(behSessionData(responseInds)));
allChoices(~isnan(allReward_R)) = 1;
allChoices(~isnan(allReward_L)) = -1;
changeChoice = [false abs(diff(allChoices)) > 0];

allReward_R(isnan(allReward_R)) = 0;
allReward_L(isnan(allReward_L)) = 0;
allChoice_R = double(allChoices == 1);
allChoice_L = double(allChoices == -1);

allRewards = zeros(1,length(allChoices));
allRewards(logical(allReward_R)) = 1;
allRewards(logical(allReward_L)) = -1;
allRewardsBinary = allRewards;                      %make all rewards have the same value
allRewardsBinary(find(allRewards==-1)) = 1;
rewardsList =  allRewards(find(allRewards~=0));

if blockSwitch(end) == length(allChoices)
    blockSwitch = blockSwitch(1:end-1);
end


%% determine lick latency
lickLat = [];       lickRate = [];
lickLat_L = [];     lickRate_L = [];
lickLat_R = [];     lickRate_R = [];
for i = 1:length(behSessionData)
    if ~isempty(behSessionData(i).rewardTime)
        lickLat = [lickLat behSessionData(i).rewardTime - behSessionData(i).CSon];
        if ~isnan(behSessionData(i).rewardL)
            lickLat_L = [lickLat_L behSessionData(i).rewardTime - behSessionData(i).CSon];
            if behSessionData(i).rewardL == 1
                if length(behSessionData(i).licksL) > 1
                    lickRateTemp = 1000/(min(diff(behSessionData(i).licksL)));
                    lickRate = [lickRate lickRateTemp];
                    lickRate_L = [lickRate_L lickRateTemp];
                else
                   lickRate = [lickRate 0];
                   lickRate_L = [lickRate_L 0]; 
                end
            end
        elseif ~isnan(behSessionData(i).rewardR)
            lickLat_R = [lickLat_R behSessionData(i).rewardTime - behSessionData(i).CSon];      
            if behSessionData(i).rewardR == 1
                if length(behSessionData(i).licksR) > 1     
                    lickRateTemp = 1000/(min(diff(behSessionData(i).licksR)));
                    lickRate = [lickRate lickRateTemp];
                    lickRate_R = [lickRate_R lickRateTemp];  
                else                                                                    %make single licks zeros for easier indexing
                    lickRate = [lickRate 0];
                    lickRate_R = [lickRate_R 0];
                end
            end
        end
    end
end

%% Z-scored lick latency (gets rid of preemptive licks)

lickLatResp = lickLat(responseInds);                    %remove NaNs from lickLat array
lickLatResp = lickLatResp(2:end);                       %shift for comparison to rwd history
lickLatInds = find(lickLatResp > 250);                  %find indices of non-preemptive licks (limit to normal distribution)

if ~isnan(behSessionData(responseInds(1)).rewardR)      %remove first response for shift to compare to rwd hist
    responseLat_R = lickLat_R(2:end);
    responseLat_L = lickLat_L;
else
    responseLat_R = lickLat_R;
    responseLat_L = lickLat_L(2:end);
end

responseLat_R = responseLat_R(responseLat_R > 250);        %remove lick latencies outside of normal distribution
responseLat_L = responseLat_L(responseLat_L > 250);
responseLat_R  = zscore(responseLat_R);                   %get z scores for lick latencies based on spout side average
responseLat_L  = zscore(responseLat_L);
choicesLick = allChoices(2:end);                        %make shifted choice array without preemptive licks
choicesLick = choicesLick(lickLatInds);

L = 1;
R = 1;
for j = 1:length(choicesLick)                     %put z scored lick latencies back in trial order
    if choicesLick(j) == 1
        responseLat(j) = responseLat_R(R);
        R = R + 1;
    else
        responseLat(j) = responseLat_L(L);
        L = L + 1;
    end
end

%% Z scored lick rate (eliminates trials with impossible lick rates)

rewardsLick = allRewards(allRewards == 1 | allRewards == -1);
responseRateInds = find(lickRate < 15);

corrLickRate = lickRate(lickRate < 15);
corrLickRate_R = lickRate_R(lickRate_R < 15);
corrLickRate_L = lickRate_L(lickRate_L < 15);
corrLickRate_R = zscore(corrLickRate_R);
corrLickRate_L = zscore(corrLickRate_L);
rewardsLick = rewardsLick(responseRateInds);

L = 1;
R = 1;
for j = 1:length(rewardsLick)                     %put z scored lick rates back in trial order
    if rewardsLick(j) == 1
        corrLickRate(j) = corrLickRate_R(R);
        R = R + 1;
    elseif rewardsLick(j) == -1
        corrLickRate(j) = corrLickRate_L(L);
        L = L + 1;
    end
end


%% create rwds array in time

%find all rwd times in terms of seconds
choiceTimes = ceil(([behSessionData(responseInds).rewardTime] - behSessionData(1).CSon)/1000); %baseline to start time and convert to s from ms
if choiceTimes(1) == 0
    choiceTimes(1) = 1;
end
rwdTimes = choiceTimes(logical(allRewardsBinary));

sessionTime = ceil((behSessionData(end).CSon + 3000 - behSessionData(1).CSon)/1000);     % find total session time and pad time for reward on last trial
sessionRwds = zeros(1,sessionTime);
sessionRwds(rwdTimes) = 1;



%% make output struct 
s = struct;

s.allChoices = allChoices;
s.allRewards = allRewards;
s.allRewardsBinary = allRewardsBinary;
s.behSessionData = behSessionData;
s.blockSwitch = blockSwitch;
s.choiceTimes = choiceTimes;
s.corrLickRate = corrLickRate;
s.lickLat = lickLat;
s.lickLatInds = lickLatInds;
s.responseInds = responseInds;
s.responseLat = responseLat;
s.responseRateInds = responseRateInds; 
s.sessionRwds = sessionRwds;

if p.Results.revForFlag
    s.blockProbs = blockProbs;
else
    if exist('blockSwitchL')
        s.blockSwitchL = blockSwitchL;
        s.blockSwitchR = blockSwitchR;
    end
end


