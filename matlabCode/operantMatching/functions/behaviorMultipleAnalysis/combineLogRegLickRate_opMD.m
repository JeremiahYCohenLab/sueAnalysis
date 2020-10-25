function [glm_rwdNoRwd, tMax] = combineLogRegLickRate_opMD(xlFile, animal, category, varargin)

%task and model parameters
p = inputParser;
% default parameters if none given
p.addParameter('revForFlag', 0)
p.addParameter('numBins', 10)
p.addParameter('plotFlag', 1);
p.addParameter('maxTrials', 350);
p.parse(varargin{:});


[root, sep] = currComputer();

[~ , dayList, ~] = xlsread([root xlFile], animal);
[~,col] = find(~cellfun(@isempty,strfind(dayList, category)) == 1);
dayList = dayList(2:end,col);
endInd = find(cellfun(@isempty,dayList),1);
if ~isempty(endInd)
    dayList = dayList(1:endInd-1,:);
end
rwdRateMatx = [];
combinedChoicesMatx = []; 
combinedRewardsMatx = [];
combinedNoRewardsMatx = [];
combineAntiLicksMatx = [];
combinedTimesMatx = [];
combinedAllChoice_R = [];
seshLength = [];
tMax = p.Results.numBins;
 

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
        if p.Results.revForFlag
            behSessionData = sessionData;
        end
    elseif p.Results.revForFlag                                    %otherwise generate the struct
        [behSessionData, ~] = generateSessionData_operantMatching(sessionName);
    else
        [behSessionData, ~, ~, ~] = generateSessionData_operantMatchingDecoupled(sessionName);
    end
    behSessionData = behSessionData(1:min(length(behSessionData), p.Results.maxTrials));
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
    
    allAntiLicks = NaN(1,length(responseInds));
    for j = 1:length(responseInds)
        if ~isnan(behSessionData(responseInds(j)).rewardL)
            allAntiLicks(j) = length(find(behSessionData(responseInds(j)).licksL < behSessionData(responseInds(j)).rewardTime));
        else
            if ~isnan(behSessionData(responseInds(j)).rewardR)
            allAntiLicks(j) = length(find(behSessionData(responseInds(j)).licksR < behSessionData(responseInds(j)).rewardTime));
            end
        end
    end
    
    outcomeTimes = [behSessionData(responseInds).rewardTime] - behSessionData(responseInds(1)).rewardTime;
    outcomeTimes = [diff(outcomeTimes) NaN];
    
    rwdMatxTmp = [];
    choiceMatxTmp = [];
    noRwdMatxTmp = [];
    antiLicksMatxTmp = [];
    
    for j = 1:tMax
        rwdMatxTmp(j,:) = [NaN(1,j) allRewards(1:end-j)];
        choiceMatxTmp(j,:) = [NaN(1,j) allChoices(1:end-j)];
        noRwdMatxTmp(j,:) = [NaN(1,j) allNoRewards(1:end-j)];
        antiLicksMatxTmp(j,:) = [NaN(1,j) allAntiLicks(1:end-j)];
    end

    timeTmp = NaN(tMax,length(allRewards)); 
    for j = 1:tMax
        for k = 1:length(outcomeTimes)-j
            timeTmp(j,k+j) = sum(outcomeTimes(k:k+j-1));
        end
    end
    
%    allRewards(allRewards == -1) = 1;
    rwdsTmp = NaN(tMax,length(allRewards)); 
    for j = 1:tMax
        for k = 1:length(outcomeTimes)-j
            rwdsTmp(j,k+j) = sum(allRewards(k:k+j-1));
        end
    end
    
%   antiLicksMatxTmp = [];
    
    combinedRewardsMatx = [combinedRewardsMatx NaN(tMax,100) rwdMatxTmp];
    combinedNoRewardsMatx = [combinedNoRewardsMatx NaN(tMax,100) noRwdMatxTmp];
    combinedChoicesMatx = [combinedChoicesMatx NaN(tMax,100) choiceMatxTmp];
    combinedTimesMatx = [combinedTimesMatx NaN(tMax, 100) timeTmp];
    combineAntiLicksMatx = [combineAntiLicksMatx NaN(tMax,100) antiLicksMatxTmp];
    combinedAllChoice_R = [combinedAllChoice_R NaN(1,100) allChoice_R];
    seshLength = [seshLength length(rwdMatxTmp)];
end

seshInd = zeros(length(combinedRewardsMatx), length(dayList));
tmpInd = 101;
for i = 1:length(dayList)-1
    seshInd(tmpInd:tmpInd+seshLength(i)-1,i) = 1;
    tmpInd = tmpInd + seshLength(i) + 100;
end



%logistic regression models
%glm_rwd = fitglm([combinedRewardsMatx]', combinedAllChoice_R,'distribution','binomial','link','logit'); rsq{1} = num2str(round(glm_rwd.Rsquared.Adjusted*100)/100);
%glm_rwdANDchoice = fitglm([combinedRewardsMatx' combinedChoicesMatx'], combinedAllChoice_R, 'distribution','binomial','link','logit'); rsq{2} = num2str(round(glm_rwdANDchoice.Rsquared.Adjusted*100)/100);
%glm_time = fitglm([combinedTimesMatx]', combinedAllChoice_R,'distribution','binomial','link','logit'); rsq{4} = num2str(round(glm_time.Rsquared.Adjusted*100)/100);
%glm_rwdANDtime = fitglm([combinedRewardsMatx' combinedTimesMatx'], combinedAllChoice_R,'distribution','binomial','link','logit'); rsq{5} = num2str(round(glm_rwdANDtime.Rsquared.Adjusted*100)/100);
%glm_rwdRate = fitglm([rwdRateMatx]', combinedAllChoice_R,'distribution','binomial','link','logit'); rsq{6} = num2str(round(glm_rwd.Rsquared.Adjusted*100)/100);
glm_rwdNoRwd = fitglm([combinedRewardsMatx' combinedNoRewardsMatx'], combinedAllChoice_R,'distribution','binomial','link','logit'); rsq{7} = num2str(round(glm_rwdNoRwd.Rsquared.Adjusted*100)/100);
glm_seshInd = fitglm([combinedRewardsMatx' combinedNoRewardsMatx' seshInd], combinedAllChoice_R,'distribution','binomial','link','logit');
%glm_all = fitglm([combinedRewardsMatx' combinedNoRewardsMatx' combinedChoicesMatx'], combinedAllChoice_R, 'distribution','binomial','link','logit');


if p.Results.plotFlag
    figure; hold on;
    relevInds = 2:tMax+1;
    coefVals = glm_rwdNoRwd.Coefficients.Estimate(relevInds);
    CIbands = coefCI(glm_rwdNoRwd);
    errorL = abs(coefVals - CIbands(relevInds,1));
    errorU = abs(coefVals - CIbands(relevInds,2));
    errorbar((1:tMax)+0.2,coefVals,errorL,errorU,'Color', [0.7 0 1],'linewidth',2)

    relevInds = tMax+2:length(glm_rwdNoRwd.Coefficients.Estimate);
    coefVals = glm_rwdNoRwd.Coefficients.Estimate(relevInds);
    CIbands = coefCI(glm_rwdNoRwd);
    errorL = abs(coefVals - CIbands(relevInds,1));
    errorU = abs(coefVals - CIbands(relevInds,2));
    errorbar((1:tMax)+0.2,coefVals,errorL,errorU,'b','linewidth',2)
    
    line([0.5 tMax+0.5], [0 0], 'Color','red','LineStyle','--')

    xlabel('Reward n Trials Back')
    ylabel('\beta Coefficient')
    xlim([0.5 tMax+0.5])
    legend('rwd', 'no rwd')
    title([animal ' ' category])
    set(gca, 'tickdir', 'out')
    set(gcf, 'renderer', 'painters')
end
