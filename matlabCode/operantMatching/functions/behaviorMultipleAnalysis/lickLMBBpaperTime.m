function [rwdMatx, lickLatMatx, combinedLickLat,sessionID] = lickLMBBpaperTime(xlFile, animal, category, varargin)

p = inputParser;
% default parameters if none given)
p.addParameter('plotFlag', 0)
p.addParameter('binSize', 20000);
p.addParameter('numBins', 10);
p.addParameter('maxTrials', 200);
p.parse(varargin{:});

[root, ~] = currComputer();
maxTrial = p.Results.maxTrials;
colors = cool(5);

[~, dayList, ~] = xlsread([root xlFile], animal);
col = contains(dayList(1,:),category);
dayList = dayList(2:end,col);
endInd = find(cellfun(@isempty,dayList),1);
if ~isempty(endInd)
    dayList = dayList(1:endInd-1,:);
end

binSize = p.Results.binSize;
timeMax = binSize * p.Results.numBins + 1000;
timeBinEdges = 1000:binSize:timeMax;  %no trials shorter than 1s between outcome and CS on
tMax = length(timeBinEdges) - 1;
rwdMatx = [];
nrwdMatx = [];
lickLatMatx = [];
combinedLickLat = [];
combinedTimeInSesh = [];
combinedChangeChoice = [];
sessionLen = zeros(size(dayList));

for i = 1: length(dayList)
 %  disp(dayList{i});
    behSessionData = loadBehavioralData([dayList{i} '.asc']);
    behavStruct = parseBehavioralData(behSessionData,maxTrial);
    
    %%generate lickLat, time in session
    responseInds = behavStruct.responseInds;
    responseInds = responseInds(1:min(maxTrial,length(responseInds)));
    behSessionData = behSessionData(1:responseInds(end));
    timeInSesh = ([behSessionData(responseInds).CSon] - behSessionData(1).CSon)/1000;
   % timeInSesh = ([behSessionData(responseInds).CSon] - behSessionData(1).CSon) / (behSessionData(responseInds(end)).CSon - behSessionData(responseInds(1)).CSon);
    changeChoice = [0 abs(diff(behavStruct.allChoices(1:min(maxTrial,length(responseInds))))) > 0];
    outcomes = abs(behavStruct.allRewards(1:min(maxTrial,length(behavStruct.allRewards))));
    %% create binned outcome matrices
    rwdTmpMatx = zeros(tMax, length(responseInds));
    nrwdTmpMatx = zeros(tMax, length(responseInds)); %initialize matrices for number of response trials x number of time bins
    lickLatTempMatx = nan(tMax, length(responseInds));
    for j = 2:length(responseInds)  
        k = 1;
        %find time between "current" choice and previous rewards, up to timeMax in the past 
        timeTmp = [];
        while j-k > 0 && behSessionData(responseInds(j)).respondTime - behSessionData(responseInds(j-k)).rewardTime < timeMax
            timeTmp = [timeTmp (behSessionData(responseInds(j)).respondTime - behSessionData(responseInds(j-k)).rewardTime)];
            k = k + 1;
        end
        %bin outcome times and use to fill matrices
        if ~isempty(timeTmp)
            binnedRwds = discretize(timeTmp,timeBinEdges);
            respIdx = (j-1):-1:(j-k+1);
            for k = 1:tMax
                if ~isempty(find(binnedRwds == k, 1))
                    rwdTmpMatx(k,j) = sum(outcomes(respIdx(binnedRwds == k)));
                    nrwdTmpMatx(k,j) = sum(binnedRwds == k) - rwdTmpMatx(k,j);
                end
            end
        end
        %find licks in previous time bins
        k = 1;
        timeTmp = [];
        while j-k > 0 && behSessionData(responseInds(j)).respondTime - behSessionData(responseInds(j-k)).respondTime < timeMax
            timeTmp = [timeTmp (behSessionData(responseInds(j)).respondTime - behSessionData(responseInds(j-k)).rewardTime)];
            k = k + 1;
        end
        %bin lick times and use to fill matrices
        if ~isempty(timeTmp)
            binnedRwds = discretize(timeTmp,timeBinEdges);
            respIdx = (j-1):-1:(j-k+1);
            for k = 1:tMax
                if ~isempty(find(binnedRwds == k, 1))
                    lickLatTempMatx(k,j) = nanmean(behavStruct.lickLatMlog(respIdx(binnedRwds == k)));
                end
            end
        end
    end
    
    %fill in NaNs at beginning of session for rwdMatrix
    j = 2;
    while behSessionData(responseInds(j)).respondTime - behSessionData(responseInds(1)).respondTime < timeMax
        tmpDiff = behSessionData(responseInds(j)).respondTime - behSessionData(responseInds(1)).respondTime;
        binnedDiff = discretize(tmpDiff, timeBinEdges);
        rwdTmpMatx(binnedDiff+1:tMax,j) = NaN;
        nrwdTmpMatx(binnedDiff+1:tMax,j) = NaN;
        j = j+1;
    end
    rwdTmpMatx(:,1) = NaN;
    %concatenate temp matrix with combined matrix 
    rwdMatx = [rwdMatx rwdTmpMatx];
    nrwdMatx = [nrwdMatx nrwdTmpMatx];
    lickLatMatx = [lickLatMatx lickLatTempMatx];
    combinedTimeInSesh = [combinedTimeInSesh timeInSesh];
    combinedLickLat = [combinedLickLat behavStruct.lickLatMlog];
    combinedChangeChoice = [combinedChangeChoice changeChoice];
    sessionLen(i)=length(timeInSesh);
end

sessionID = zeros(1,length(combinedLickLat));
for i = 1:length(dayList)
    if i > 1
        sessionID(1,sum(sessionLen(1:i-1))+1:sum(sessionLen(1:i-1))+sessionLen(i))=i;
    else
        sessionID(1,1:sessionLen(i))=i;
    end

end

%linear regression model
glm_LickRwd = fitlm([rwdMatx' lickLatMatx'], combinedLickLat);
glm_Rwd = fitlm([rwdMatx' combinedTimeInSesh'], combinedLickLat);
%% plot everything
c = 1;
b = 2;
a = 1;
figure2; hold on;
subplot(1,3,1); hold on;
CIbands = coefCI(glm_LickRwd);
relevInds = 2:tMax+1;
coefVals = glm_LickRwd.Coefficients.Estimate(relevInds);
errorL = abs(coefVals - CIbands(relevInds,1));
errorU = abs(coefVals - CIbands(relevInds,2));
errorbar(((1:tMax)*binSize/1000),coefVals,errorL,errorU,'Color', colors(c,:),'linewidth',2)
line([0 (tMax*binSize/1000 + 5)],[0 0], 'Color', [0.5 0.5 0.5], 'LineStyle','--');
xlim([0 (tMax*binSize/1000 + 5)])
xlabel('Rewards back in Time')
ylabel('\beta Coefficient')
title('rewards')

subplot(1,3,2); hold on;
relevInds = tMax + 2:2*tMax+1;
coefVals = glm_LickRwd.Coefficients.Estimate(relevInds);
errorL = abs(coefVals - CIbands(relevInds,1));
errorU = abs(coefVals - CIbands(relevInds,2));
errorbar(((1:tMax)*binSize/1000),coefVals,errorL,errorU,'Color', colors(c,:),'linewidth',2)
line([0 (tMax*binSize/1000 + 5)],[0 0], 'Color', [0.5 0.5 0.5], 'LineStyle','--');
xlim([0 (tMax*binSize/1000 + 5)])
xlabel('LickLats back in Time')
ylabel('\beta Coefficient')
title('lickLats')
% 
% subplot(2,2,3); hold on;
% CIbands = coefCI(glm_preLickRwd);
% relevInds = 2:tMax;
% coefVals = glm_preLickRwd.Coefficients.Estimate(relevInds);
% errorL = abs(coefVals - CIbands(relevInds,1));
% errorU = abs(coefVals - CIbands(relevInds,2));
% errorbar(((1:tMax-1)*binSize/1000),coefVals,errorL,errorU,'Color', colors(5,:),'linewidth',2)
% line([0 (tMax*binSize/1000 + 5)],[0 0], 'Color', [0.5 0.5 0.5], 'LineStyle','--');
% xlim([0 (tMax*binSize/1000 + 5)])
% xlabel('Rewards back in Time')
% ylabel('\beta Coefficient')
% title('rewards')

subplot(1,3,3); hold on;
int = glm_LickRwd.Coefficients.Estimate(1);
bar(a, int, 'FaceColor', colors(b,:), 'EdgeColor', colors(c,:))
errorL = abs(int - CIbands(1,1));
errorU = abs(int - CIbands(1,2));
errorbar(a, int,errorL,errorU,'Color', colors(c,:),'linewidth',2)
xlim([0 3])
xlabel('intercept')
ylabel('\beta Coefficient')
title('intercept')
suptitle([animal '  ' category])