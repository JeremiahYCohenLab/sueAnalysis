function [mdl, s] = simulationLogRegSession(xlFile, sheet, category, varargin)


%task and model parameters
p = inputParser;
% default parameters if none given
p.addParameter('modelName', 'fiveParam_rBeta_scale');
p.addParameter('params', []);
p.addParameter('taskType', 'coupled');
p.addParameter('rwdProbs', [90 50 10]);
p.addParameter('maxTrials', 300);
p.addParameter('randomSeed', 27);
p.addParameter('figFlag', 0);
p.addParameter('revForFlag', 0);
p.parse(varargin{:});

[root, sep] = currComputer();
tMax = 12;
rwdMatx = [];
noRwdMatx = [];
combinedAllChoice_R = [];
allRewardsBin = [];
changeChoice = [];
rSeed = p.Results.randomSeed;

[~, dayList, ~] = xlsread(xlFile, sheet);
[~,col] = find(~cellfun(@isempty,strfind(dayList, category)) == 1);
dayList = dayList(2:end,col);
delInds = any(cellfun(@(x) isempty(x) || (ischar(x) && all(x==' ')),dayList),2); 
dayList(delInds,:) = [];
runs = length(dayList);

if regexp(p.Results.taskType, 'switch')
    [actual, tMax] = combineLogReg_opMD(xlFile, sheet, 'switch', p.Results.revForFlag);
else
    [actual, tMax] = combineLogReg_opMD(xlFile, sheet, category, p.Results.revForFlag);
end


for runInd = 1:runs
    fprintf('Run: %i of %i \n', runInd, runs);
    rSeed = rSeed + 1;
    
    if isempty(p.Results.params)
        modelPath = [root sheet sep sheet 'sorted' sep 'stan' sep p.Results.modelName sep sheet...
                    category '_' p.Results.modelName '.mat'];
        t = generateStanModelTerms_opMD(p.Results.modelName, modelPath, dayList{runInd}, 1, p.Results.revForFlag);
        params = t.params;
    else
        params = p.Results.params;
    end

    switch p.Results.modelName
        case 'fourParam'
            [allRewards, allChoices, blockProbs, blockSwitch] = qLearningModel_simNoPlot('params', params,...
                'maxTrials', p.Results.maxTrials, 'randomSeed', rSeed, 'rwdProbs', p.Results.rwdProbs, 'taskType', p.Results.taskType);
        case 'fiveParam_rBeta_scale'
            [~, allRewards, allChoices, blockProbs, blockSwitch] = qLearningModel_rBeta_simNoPlot('params', params,...
                'maxTrials', p.Results.maxTrials, 'randomSeed', rSeed, 'rwdProbs', p.Results.rwdProbs, 'taskType', p.Results.taskType);
        case 'fiveParam_rBeta_kappa'
            [~, allRewards, allChoices, blockProbs, blockSwitch] = qLearningModel_rBetaKappa_simNoPlot('params', params,...
                'maxTrials', p.Results.maxTrials, 'randomSeed', rSeed, 'rwdProbs', p.Results.rwdProbs, 'taskType', p.Results.taskType);
    end
    
    allRewards(allChoices == -1) = -1;
    allNoRewards = allChoices;
    allNoRewards(find(allRewards == 1)) = 0;
    allNoRewards(find(allRewards == -1)) = 0;
    
    allChoice_R = allChoices;
    allChoice_R(find(allChoice_R == -1)) = 0;
    
    rwdMatxTmp = [];
    noRwdMatxTmp = [];
    for tInd = 1:tMax
        rwdMatxTmp(tInd,:) = [NaN(1,tInd) allRewards(1:end-tInd)];
        noRwdMatxTmp(tInd,:) = [NaN(1,tInd) allNoRewards(1:end-tInd)];
    end
    
   
    rwdMatx = [rwdMatx NaN(tMax,100) rwdMatxTmp];
    noRwdMatx = [noRwdMatx NaN(tMax,100) noRwdMatxTmp];
    combinedAllChoice_R = [combinedAllChoice_R NaN(1,100) allChoice_R];
    
    allRewardsBin = allRewards;
    allRewardsBin(allRewardsBin == -1) = 1;
    allRewardsBin = allRewardsBin(1:end-1);
    changeChoice = [false abs(diff(allChoices)) > 0];
    changeChoice = changeChoice(2:end);
    s.probSwitchNoRwd(runInd) = sum(changeChoice(allRewardsBin(1:end-1)==0))/sum(allRewardsBin(1:end-1)==0);
    s.probStayRwd(runInd) = 1 - (sum(changeChoice(allRewardsBin(1:end-1)==1))/sum(allRewardsBin(1:end-1)==1));
end


s.ws(1,1)= mean(s.probStayRwd);
s.ws(1,2) = std(s.probStayRwd)/sqrt(length((s.probStayRwd)));
s.ls(1,1)= mean(s.probSwitchNoRwd);
s.ls(1,2) = std(s.probSwitchNoRwd)/sqrt(length((s.probSwitchNoRwd)));

%run the glm on the simulated data
mdl = fitglm([rwdMatx' noRwdMatx'], combinedAllChoice_R,'distribution','binomial','link','logit');

%plot the simulated beta coefficients against the actual
if p.Results.figFlag
    figure;
    subplot(1,2,1); hold on;
    relevInds = 2:tMax+1;
    coefVals = mdl.Coefficients.Estimate(relevInds);            %model coeffs
    CIbands = coefCI(mdl);
    errorL = abs(coefVals - CIbands(relevInds,1));
    errorU = abs(coefVals - CIbands(relevInds,2));
    errorbar((1:tMax)+0.2,coefVals,errorL,errorU,'Color',[0.7 0 1],'linewidth',2)
    
    coefVals = actual.Coefficients.Estimate(relevInds);         %actual coeffs
    CIbands = coefCI(actual);
    errorL = abs(coefVals - CIbands(relevInds,1));
    errorU = abs(coefVals - CIbands(relevInds,2));
    errorbar((1:tMax)+0.2,coefVals,errorL,errorU,'Color',[0.7 0.5 1],'linewidth',2)
    legend('simulated', 'actual')
    xlabel('Reward n Trials Back')
    ylabel('\beta Coefficient')
    xlim([0.5 tMax+0.5])

    subplot(1,2,2); hold on;
    relevInds = tMax+2:length(mdl.Coefficients.Estimate);
    coefVals = mdl.Coefficients.Estimate(relevInds);
    CIbands = coefCI(mdl);
    errorL = abs(coefVals - CIbands(relevInds,1));
    errorU = abs(coefVals - CIbands(relevInds,2));
    errorbar((1:tMax)+0.2,coefVals,errorL,errorU,'b','linewidth',2)

    coefVals = actual.Coefficients.Estimate(relevInds);
    CIbands = coefCI(actual);
    errorL = abs(coefVals - CIbands(relevInds,1));
    errorU = abs(coefVals - CIbands(relevInds,2));
    errorbar((1:tMax)+0.2,coefVals,errorL,errorU,'Color',[0.5 0.5 1],'linewidth',2)
    legend('simulated', 'actual')
    xlabel('No Reward n Trials Back')
    ylabel('\beta Coefficient')
    xlim([0.5 tMax+0.5])
    titleTmp = strrep([sheet ': ' category ' ' p.Results.modelName], '_', ' ');
    suptitle(titleTmp)
end
