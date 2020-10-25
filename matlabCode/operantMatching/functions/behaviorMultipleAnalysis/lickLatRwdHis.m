function [combineLickHx, combineLickHxTrial,combineLickLat] = lickLatRwdHis(xlFile, animal, category, varargin)
%task and model parameters
p = inputParser;
% default parameters if none given
p.addParameter('revForFlag', 0)
p.addParameter('binSize', 5000)
p.addParameter('numBins', 10)
p.addParameter('plotFlag', 1);
p.addParameter('maxTrials', 350);
p.parse(varargin{:});

%determine root for file location
[root, sep] = currComputer();
% compute time  kernel
[t,fitresult] = combineLogRegTime_opMD(xlFile, animal, category, p.Results, 'plotFlag', 0);
kernel = fitresult.a*exp(-(1/(fitresult.b))*(0:t.timeMax/1000));
kernel = kernel./sum(kernel);
[root, sep] = currComputer();

% compute trial kernel
 tMax = p.Results.numBins;
 [~, ~, ~, fitresult] = combineLogReg_opMD(xlFile, animal, category, 'maxTrials', p.Results.maxTrials, 'plotFlag', 0);
 kernelTrial = fitresult.a*exp(-(1/(fitresult.b))*(0:tMax));
 kernelTrial = kernelTrial./sum(kernelTrial);

%import behavior session titles for desired category
[~, dayList, ~] = xlsread([root xlFile], animal);
[~,col] = find(~cellfun(@isempty,strfind(dayList, category)) == 1);
dayList = dayList(2:end,col);
endInd = find(cellfun(@isempty,dayList),1);
if ~isempty(endInd)
    dayList = dayList(1:endInd-1,:);
end

combineLickLat = [];
combineLickHx = [];
combineLickHxTrial = [];
combineQuitHx = [];
combineQuitHxT = [];
combineQuitTrial = [];

for i = 1: length(dayList)              
    sessionName = dayList{i};                       %extract relevant info from session title
    [animalName, date] = strtok(sessionName, 'd'); 
    animalName = animalName(2:end);
    date = date(1:9);
    sessionFolder = ['m' animalName date];

    if isstrprop(sessionName(end), 'alpha')         %define appropriate data path
        sessionDataPath = [root animalName sep sessionFolder sep 'sorted' sep 'session ' sessionName(end) sep sessionName '_sessionData_behav.mat'];
    else
        sessionDataPath = [root animalName sep sessionFolder sep 'sorted' sep 'session' sep sessionName '_sessionData_behav.mat'];
    end

    if exist(sessionDataPath,'file')        %load preprocessed struct if there is one
        load(sessionDataPath)
        if p.Results.revForFlag
            behSessionData = sessionData;
        end
    elseif p.Results.revForFlag                                    %otherwise generate the struct
        [behSessionData, ~] = generateSessionData_operantMatching(sessionName);
    else
        [behSessionData, ~, ~, ~] = generateSessionData_operantMatchingDecoupled(sessionName);
    end
    behSessionData = behSessionData(1:min(length(behSessionData),p.Results.maxTrials));
   
    %create reward env
    
    responseInds = find(~isnan([behSessionData.rewardTime])); % find CS+ trials with a response in the lick window
    lengthSec = round((behSessionData(length(behSessionData)).CSon + 5000 - behSessionData(1).CSon)/1000) + 1; 
    allRewardsBinary = zeros(1,lengthSec);
    
    %kernel reward env
    for j = 1: length(responseInds)
        if ~isnan(behSessionData(responseInds(j)).rewardL) || ~isnan(behSessionData(responseInds(j)).rewardR)
            secInd = fix((behSessionData(responseInds(j)).rewardTime - behSessionData(1).CSon)/1000) + 1;
            allRewardsBinary(secInd) = 1;
        end
    end
    
    rwdHx = conv(allRewardsBinary,kernel);              %convolve with exponential decay to give weighted moving average
    rwdHx = rwdHx(1:end-(length(kernel)-1));
    lickHxsecInd = fix(([behSessionData(responseInds).CSon] - behSessionData(1).CSon)/1000) + 1;
    lickHx = rwdHx(lickHxsecInd);
    
    % rwdHx by trial
    allReward_R = [behSessionData(responseInds).rewardR]; 
    allReward_L = [behSessionData(responseInds).rewardL];  
    allChoices = NaN(1,length(behSessionData(responseInds)));
    allChoices(~isnan(allReward_R)) = 1;
    allChoices(~isnan(allReward_L)) = -1;
    
    allReward_R(isnan(allReward_R)) = 0;
    allReward_L(isnan(allReward_L)) = 0;

    allRewards = zeros(1,length(responseInds));
    allRewards(logical(allReward_R)) = 1;
    allRewards(logical(allReward_L)) = -1;
    
    allRewardsBinary = allRewards;                       %makewe all rewards have the same value
    allRewardsBinary((allRewards == -1)) = 1;
    
    % rwd env by trial
    rwdHxTrial = conv(allRewardsBinary,kernelTrial);              %convolve with exponential decay to give weighted moving average
    rwdHxTrial = rwdHxTrial(1:end-(length(kernelTrial)-1));
    rwdHxTrial = [0, rwdHxTrial(2:end)];
    %to account for convolution padding

    % calculate lickLat
    lickLat = [behSessionData(responseInds).respondTime] - [behSessionData(responseInds).CSon];
    lickLatlog = zeros(1,length(lickLat));
    indsR = find(allChoices == 1);
    indsL = find(allChoices == -1);
    lickLat_R = zscore(lickLat(indsR));
    lickLat_L = zscore(lickLat(indsL));
    lickLat_Rlog = log(lickLat(indsR));
    lickLat_Llog = log(lickLat(indsL));
    lickLatZ = NaN(1, length(allChoices));
    lickLatZ(indsR) = lickLat_R;
    lickLatZ(indsL) = lickLat_L;
    lickLatlog(indsR) = lickLat_Rlog;
    lickLatlog(indsL) = lickLat_Llog;
    %find quit trials
    

    A = cellfun(@(x) strcmp(x,'CSplus'),{behSessionData.trialType});
    B = cellfun(@isnan, {behSessionData.rewardTime});
    quitTrial = find(A>0 & B>0);
    quitTime = fix(([behSessionData(quitTrial).CSon] - behSessionData(1).CSon)/1000) + 1;
    quitHx = rwdHx(quitTime);
    quitHxT = NaN(size(quitTrial));
    for j = 1:length(quitTrial)
        m = 1;
        while j - m >0 && (~strcmp(behSessionData(quitTrial(j)-m).trialType, 'CSplus') || isnan(behSessionData(quitTrial(j)-m).rewardTime))
            m = m+1;
        end
        if quitTrial(j) - m > 0
            if ~isnan(find(responseInds == (quitTrial(j) - m)))
                quitHxT(j) = find(responseInds == (quitTrial(j) - m));
            end
        end
    end
    
    quitHxT = rwdHxTrial(quitHxT(~isnan(quitHxT)));

    combineLickLat = [combineLickLat lickLat];
    combineLickHx = [combineLickHx lickHx];
    combineLickHxTrial = [combineLickHxTrial rwdHxTrial];
    combineQuitHx = [combineQuitHx quitHx];
    combineQuitHxT = [combineQuitHxT quitHxT];
    combineQuitTrial = [combineQuitTrial quitTrial];
end

% group lickLats by rwd env
[group,~] = discretize(combineLickHx,p.Results.numBins);
[groupT,~] = discretize(combineLickHxTrial,p.Results.numBins);

len = zeros(1,p.Results.numBins);
lenT = zeros(1,p.Results.numBins);

for i = 1:tMax
    len(i) = length(find(group==i));
    lenT(i) = length(find(groupT==i));
end

small = find(len<5);
smallT = find(lenT<5);
if ~isempty(small)&&length(small)==1
    if small < 0.5*tMax
        for j = small:1
            group(group==j) = j+1;
        end
    else 
        for j = small:tMax
            group(group==j) = j-1;
        end
    end
end

if ~isempty(smallT)&&length(smallT)==1
    if smallT < 0.5*tMax
        for j = smallT:1
            groupT(groupT==j) = j+1;
        end
    else 
        for j = smallT:tMax
            groupT(groupT==j) = j-1;
        end
    end
end
lickLatM = zeros(1,max(group));
lickLatSEM = zeros(1,max(group));
binMean = zeros(1,max(group));

lickLatTM = zeros(1,max(groupT));
lickLatTSEM = zeros(1,max(groupT));
binMeanT = zeros(1,max(groupT));
for i = 1:max(group)
    lickLatM(i) = mean(combineLickLat(group == i));
    lickLatSEM(i) = std(combineLickLat(group == i))/sqrt(len(i));
    binMean(i) = mean(combineLickHx(group == i));
end
for i = 1:max(groupT)
    lickLatTM(i) = mean(combineLickLat(groupT == i));
    lickLatTSEM(i) = std(combineLickLat(groupT == i))/sqrt(lenT(i));
    binMeanT(i) = mean(combineLickHxTrial(groupT == i));
end

%% plot everything

figure2('Position', [1 1 800 400]); 
suptitle([animal '  ' category])
subplot(2,2,1); hold on;
scatter(combineLickHx, combineLickLat, 3, 'MarkerFaceColor',[.7 .7 .7], 'MarkerEdgeColor', 'none');
[R,P,RL,RU] = corrcoef(combineLickHx, combineLickLat);

errorbar(binMean, lickLatM, lickLatSEM, 'Color', [0.9 0.3 0.3],'linewidth',1.5)
%ylim([200 1500]) 
ylabel('lickLat')
xlabel('rewardHx by time')

subplot(2,2,2); hold on;
scatter(combineLickHxTrial, combineLickLat, 3, 'MarkerFaceColor',[.7 .7 .7], 'MarkerEdgeColor', 'none');
%errorbar(binMeanT, lickLatTM, lickLatTSEM, 'Color', [0.3 0.3 0.9],'linewidth',1.5)
errorbar(binMeanT, lickLatTM, lickLatTSEM, 'Color', [0.9 0.3 0.3],'linewidth',1.5)
%ylim([200 1500])
ylabel('lickLat')
xlabel('rewardHx by trial')

subplot(2,3,4); hold on;
cdfplot(combineQuitHx)
xlabel('rewardEnvinTime')
ylabel('pQuit')
subplot(2,3,5); hold on;
cdfplot(combineQuitHxT)
xlabel('rewardEnvinTrial')
ylabel('pQuit')
subplot(2,3,6); hold on;
cdfplot(combineQuitTrial)
xlabel('Trial')
ylabel('pQuit')























