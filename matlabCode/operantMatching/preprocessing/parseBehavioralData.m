function outputStruct = parseBehavioralData(sessionData,maxTrial)

minLickLat = 250;
responseInds = find(~isnan([sessionData.rewardTime]));
omitInds = isnan([sessionData.rewardTime]);

allReward_R = [sessionData(responseInds).rewardR]; 
allReward_L = [sessionData(responseInds).rewardL]; 
allChoices = NaN(1,length(sessionData(responseInds)));
allChoices(~isnan(allReward_R)) = 1;
allChoices(~isnan(allReward_L)) = -1;

allReward_R(isnan(allReward_R)) = 0;
allReward_L(isnan(allReward_L)) = 0;
allChoice_R = double(allChoices == 1);
allChoice_L = double(allChoices == -1);

allRewards = zeros(1,length(allChoices));
allRewards(logical(allReward_R)) = 1;
allRewards(logical(allReward_L)) = -1;

timeBtwn = [sessionData(responseInds(2:end)).CSon] - [sessionData(responseInds(1:end-1)).rewardTime];
timeBtwn(timeBtwn < 0 ) = 0;
timeBtwn = [0 timeBtwn(~isnan(timeBtwn))/1000];

% calculate lickLat
if isempty(find(contains(fields(sessionData), 'respondTime') == 1, 1))
   lickLat = [sessionData(responseInds).rewardTime] - [sessionData(responseInds).CSon];
else
   lickLat = [sessionData(responseInds).respondTime] - [sessionData(responseInds).CSon];
end

%% conver to outcome
outputStruct.responseInds = responseInds;
outputStruct.omitInds = omitInds;
outputStruct.allRewards = allRewards;
outputStruct.allReward_R = allReward_R;
outputStruct.allReward_L = allReward_L;
outputStruct.allChoices = allChoices;
outputStruct.allChoice_R = allChoice_R;
outputStruct.allChoice_L = allChoice_L;
outputStruct.timeBtwn = timeBtwn;
outputStruct.lickLat = lickLat;
%% zscore
lickLat = lickLat(1:min(length(lickLat), maxTrial));
realID = lickLat>minLickLat;
allChoices = allChoices(1:min(length(lickLat), maxTrial));
lickLatZ = NaN(1, length(allChoices));
lickLatLog = NaN(1,length(lickLat));
lickLatLogZ = NaN(1,length(lickLat));
indsR = allChoices == 1 & realID;
indsL = allChoices == -1 & realID;
lickLat_R = zscore(lickLat(indsR));
lickLat_L = zscore(lickLat(indsL));
lickLat_Rlog = log(lickLat(indsR));
lickLat_Llog = log(lickLat(indsL));
lickLatZ(indsR) = lickLat_R;
lickLatZ(indsL) = lickLat_L;
lickLatLog(indsR) = lickLat_Rlog;
lickLatLog(indsL) = lickLat_Llog;
lickLatLogZ(indsR) = zscore(lickLat_Rlog);
lickLatLogZ(indsL) = zscore(lickLat_Llog);
%% mean shift
lickLat_R = lickLat(indsR);
lickLat_L = lickLat(indsL);
R = mean(lickLat_R);
L = mean(lickLat_L);
lickLatM = NaN(1, length(allChoices));
lickLatM(indsR) = lickLat_R + 0.5*(L-R);
lickLatM(indsL) = lickLat_L + 0.5*(R-L);
%% mean log shift
lickLat_R = log(lickLat(indsR));
lickLat_L = log(lickLat(indsL));
R = mean(lickLat_R);
L = mean(lickLat_L);
lickLatMlog = NaN(1, length(allChoices));
lickLatMlog(indsR) = lickLat_R + 0.5*(L-R);
lickLatMlog(indsL) = lickLat_L + 0.5*(R-L);
%% Convert to outputStruct
outputStruct.lickLatZ = lickLatZ;
outputStruct.lickLatLog = lickLatLog;
outputStruct.lickLatLogZ = lickLatLogZ;
outputStruct.lickLatM = lickLatM;
outputStruct.lickLatMlog = lickLatMlog;
outputStruct.lickReal = realID;