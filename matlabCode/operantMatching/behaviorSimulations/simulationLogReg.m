function [mdl, s] = simulationLogReg(xlFile, sheet, category, varargin)


%task and model parameters
p = inputParser;
% default parameters if none given
p.addParameter('modelName', 'fiveParam_rBeta_scale');
p.addParameter('params', []);
p.addParameter('taskType', 'coupled');
p.addParameter('rwdProbs', [90 50 10]);
p.addParameter('runs', 5);
p.addParameter('maxTrials', 300);
p.addParameter('randomSeed', 27);
p.addParameter('figFlag', 0);
p.addParameter('revForFlag', 0);
p.addParameter('compareFlag', 1);
p.parse(varargin{:});

if isempty(p.Results.params)
    [root, sep] = currComputer();
    modelPath = [root sheet sep sheet 'sorted' sep 'stan' sep p.Results.modelName sep sheet...
                category '_' p.Results.modelName '.mat'];
    t = generateStanModelTerms_opMD(p.Results.modelName, modelPath, [], 0, p.Results.revForFlag);
    params = t.params;
else
    params = p.Results.params;
end

tMax = 12;
rwdMatx = [];
noRwdMatx = [];
combinedAllChoice_R = [];
allRewardsBin = [];
changeChoice = [];
rSeed = p.Results.randomSeed;

if p.Results.compareFlag
    if regexp(p.Results.taskType, 'switch')
        [actual, tMax] = combineLogReg_opMD(xlFile, sheet, 'switch', p.Results.revForFlag);
    else
        [actual, tMax] = combineLogReg_opMD(xlFile, sheet, category, p.Results.revForFlag);
    end
end

for runInd = 1:p.Results.runs
    if rem(runInd,20) == 0
        fprintf('Running simulation %d of %d \n', runInd, p.Results.runs);
    end
    rSeed = rSeed + 1;
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
        case 'sixParam_absPePeAN_bi'
            [~, allRewards, allChoices, blockProbs, blockSwitch] = qLearningModel_absPePeAN_bi_simNoPlot('params', params,...
                'maxTrials', p.Results.maxTrials, 'randomSeed', rSeed, 'rwdProbs', p.Results.rwdProbs);
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


mdl = fitglm([rwdMatx' noRwdMatx'], combinedAllChoice_R,'distribution','binomial','link','logit');
if p.Results.figFlag
    figure;
    subplot(1,2,1); hold on;
    relevInds = 2:tMax+1;
    coefVals = mdl.Coefficients.Estimate(relevInds);
    CIbands = coefCI(mdl);
    errorL = abs(coefVals - CIbands(relevInds,1));
    errorU = abs(coefVals - CIbands(relevInds,2));
    errorbar((1:tMax)+0.2,coefVals,errorL,errorU,'Color',[0.7 0 1],'linewidth',2)
    if p.Results.compareFlag
        coefVals = actual.Coefficients.Estimate(relevInds);
        CIbands = coefCI(actual);
        errorL = abs(coefVals - CIbands(relevInds,1));
        errorU = abs(coefVals - CIbands(relevInds,2));
        errorbar((1:tMax)+0.2,coefVals,errorL,errorU,'Color',[0.7 0.5 1],'linewidth',2)
        legend('simulated', 'actual')
    end
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
    if p.Results.compareFlag
        coefVals = actual.Coefficients.Estimate(relevInds);
        CIbands = coefCI(actual);
        errorL = abs(coefVals - CIbands(relevInds,1));
        errorU = abs(coefVals - CIbands(relevInds,2));
        errorbar((1:tMax)+0.2,coefVals,errorL,errorU,'Color',[0.5 0.5 1],'linewidth',2)
        legend('simulated', 'actual')
    end
    xlabel('No Reward n Trials Back')
    ylabel('\beta Coefficient')
    xlim([0.5 tMax+0.5])
%    title(sprintf(['v = %i'], p.Results.params(5)));
end
