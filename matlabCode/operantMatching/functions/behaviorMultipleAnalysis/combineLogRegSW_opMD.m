function [glm_rwdNoRwd ,combinedRewardsMatx, combinedNoRewardsMatx, fitresult] = combineLogRegSW_opMD(xlFile, animal, category, varargin)

%task and model parameters
p = inputParser;
% default parameters if none given
p.addParameter('revForFlag', 0)
p.addParameter('numBins', 10)
p.addParameter('plotFlag', 1);
p.addParameter('maxTrials', 600);
p.addParameter('modelName', '5params');
p.addParameter('numSamps', 200);
p.parse(varargin{:});


[root, sep] = currComputer();

[~ , dayList, ~] = xlsread([root xlFile], animal);
[~,col] = find(~cellfun(@isempty,strfind(dayList, category)) == 1);
dayList = dayList(2:end,col);
endInd = find(cellfun(@isempty,dayList),1);
if ~isempty(endInd)
    dayList = dayList(1:endInd-1,:);
end

combinedRewardsMatx = [];
combinedNoRewardsMatx = [];
combinedChoicesMatx = [];
combinedTimesMatx = [];
combinedAllChoice_R = [];
combinedPrefAlt = [];

combinedSwMatx = [];
combinedSw = [];
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
        [behSessionData, ~, ~, ~] = generateSessionData_operantMatchingDecoupledRwdDelay(sessionName);
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
    
    outcomeTimes = [behSessionData(responseInds).rewardTime] - behSessionData(responseInds(1)).rewardTime;
    outcomeTimes = [diff(outcomeTimes) NaN];
    
    sw = [0 double(allChoices(1:end-1)~=allChoices(2:end))];
    
    rwdMatxTmp = [];
    choiceMatxTmp = [];
    noRwdMatxTmp = [];
    swMatxTmp = [];
    for j = 1:tMax
        rwdMatxTmp(j,:) = [NaN(1,j) allRewards(1:end-j)];
        choiceMatxTmp(j,:) = [NaN(1,j) allChoices(1:end-j)];
        noRwdMatxTmp(j,:) = [NaN(1,j) allNoRewards(1:end-j)];
        swMatxTmp(j,:) = [NaN(1,j) sw(1:end-j)];
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
   
    [params, modelName, ~, noSession] = getStanModelParams_sampsOnly(animalName, category, p.Results.modelName, p.Results.numSamps, 'sessionParamsFlag', 1, 'sessionName', dayList{i});
    if noSession
        fprintf([dayList{i} ' no good behavior \n'])
        continue
    end
    t = inferModelVar(dayList{i}, params, modelName);
    paramNames = getParamNames_dF(p.Results.modelName, 1);
    betaMean = mean(t.params(:,cellfun(@(x) strcmp(x, 'beta'), paramNames)));
    biasMean = mean(t.params(:,cellfun(@(x) strcmp(x, 'bias'), paramNames)));
    pref = betaMean*(t.Q(:,2)-t.Q(:,1)) + biasMean;
    
    prefAltTemp = [pref(2:end); NaN]; % preference of right choice on next trial
    
    prefAltTemp(allChoices > 0) = -prefAltTemp(allChoices > 0); % flip is current trial is right
    
    
    combinedRewardsMatx = [combinedRewardsMatx NaN(tMax,100) rwdMatxTmp];
    combinedNoRewardsMatx = [combinedNoRewardsMatx NaN(tMax,100) noRwdMatxTmp];
    combinedChoicesMatx = [combinedChoicesMatx NaN(tMax,100) choiceMatxTmp];
    combinedTimesMatx = [combinedTimesMatx NaN(tMax, 100) timeTmp];
    combinedAllChoice_R = [combinedAllChoice_R NaN(1,100) allChoice_R];
    combinedSwMatx = [combinedSwMatx NaN(tMax,100) swMatxTmp];
    combinedSw = [combinedSw NaN(1,100) sw];
    combinedPrefAlt = [combinedPrefAlt; NaN(100,1); prefAltTemp];
     
    seshLength = [seshLength length(rwdMatxTmp)];
end

seshInd = zeros(length(combinedRewardsMatx), length(dayList));
tmpInd = 101;
% for i = 1:length(dayList)-1
%     seshInd(tmpInd:tmpInd+seshLength(i)-1,i) = 1;
%     tmpInd = tmpInd + seshLength(i) + 100;
% end



%logistic regression models
%glm_rwd = fitglm([combinedRewardsMatx]', combinedAllChoice_R,'distribution','binomial','link','logit'); rsq{1} = num2str(round(glm_rwd.Rsquared.Adjusted*100)/100);
%glm_rwdANDchoice = fitglm([combinedRewardsMatx' combinedChoicesMatx'], combinedAllChoice_R, 'distribution','binomial','link','logit'); rsq{2} = num2str(round(glm_rwdANDchoice.Rsquared.Adjusted*100)/100);
%glm_time = fitglm([combinedTimesMatx]', combinedAllChoice_R,'distribution','binomial','link','logit'); rsq{4} = num2str(round(glm_time.Rsquared.Adjusted*100)/100);
%glm_rwdANDtime = fitglm([combinedRewardsMatx' combinedTimesMatx'], combinedAllChoice_R,'distribution','binomial','link','logit'); rsq{5} = num2str(round(glm_rwdANDtime.Rsquared.Adjusted*100)/100);
%glm_rwdRate = fitglm([rwdRateMatx]', combinedAllChoice_R,'distribution','binomial','link','logit'); rsq{6} = num2str(round(glm_rwd.Rsquared.Adjusted*100)/100);
glm_rwdNoRwd = fitglm([combinedRewardsMatx' combinedNoRewardsMatx'], combinedAllChoice_R,'distribution','binomial','link','logit'); rsq{7} = num2str(round(glm_rwdNoRwd.Rsquared.Adjusted*100)/100);
% glm_auto = fitglm([combinedRewardsMatx'+combinedNoRewardsMatx'], combinedAllChoice_R,'distribution','binomial','link','logit');
% glm_seshInd = fitglm([combinedRewardsMatx' combinedNoRewardsMatx' seshInd], combinedAllChoice_R,'distribution','binomial','link','logit');
%glm_all = fitglm([combinedRewardsMatx' combinedNoRewardsMatx' combinedChoicesMatx'], combinedAllChoice_R, 'distribution','binomial','link','logit');
[fitresult, gof] = singleExpFit(glm_rwdNoRwd.Coefficients.Estimate(2:tMax+1), (1:tMax)');
glm_sw = fitglm([combinedSwMatx', combinedPrefAlt], combinedSw, 'distribution','binomial','link','logit');
if p.Results.plotFlag
    figure; hold on;
    relevInds = 2:tMax+1;
    coefVals = glm_sw.Coefficients.Estimate(relevInds);
    CIbands = coefCI(glm_sw);
    errorL = abs(coefVals - CIbands(relevInds,1));
    errorU = abs(coefVals - CIbands(relevInds,2));
    errorbar((1:tMax),coefVals,errorL,errorU,'Color', [0.7 0 1],'linewidth',2)
    line([0.5 tMax+0.5], [0 0], 'Color',[0.5 0.5 0.5],'LineStyle','--')
    xlabel('Switch n Trials Back')
    ylabel('\beta Coefficient')
    xlim([0.5 tMax+0.5])
%     legend('rwd', 'no rwd', 'fit')
    title([animal ' ' category])
    set(gca, 'tickdir', 'out')
    set(gcf, 'renderer', 'painters')
end
