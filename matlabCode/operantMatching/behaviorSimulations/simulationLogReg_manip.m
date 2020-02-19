function [mdl_pre, s] = simulationLogReg_manip(xlFile, sheet, pre, post, varargin)


%task and model parameters
p = inputParser;
% default parameters if none given
p.addParameter('modelName', 'sixParam_absPePeAN_bi');
p.addParameter('params', []);
p.addParameter('taskType', 'coupled');
p.addParameter('rwdProbs', [90 50 10]);
p.addParameter('runs', 1000);
p.addParameter('maxTrials', 300);
p.addParameter('randomSeed', 27);
p.addParameter('figFlag', 1);
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

rwdMatx_pre = [];
noRwdMatx_pre = [];
combinedAllChoice_R_pre = [];

rwdMatx_post = [];
noRwdMatx_post = [];
combinedAllChoice_R_post = [];

if p.Results.compareFlag
    [actual_pre, tMax] = combineLogReg_opMD(xlFile, sheet, pre, p.Results.revForFlag);
    [actual_post, ~] = combineLogReg_opMD(xlFile, sheet, post, p.Results.revForFlag);
end

rSeed = p.Results.randomSeed;
for runInd = 1:p.Results.runs
    if rem(runInd,20) == 0
        fprintf('Running simulation %d of %d \n', runInd, p.Results.runs);
    end
    rSeed = rSeed + 1;
    switch p.Results.modelName
        case 'fourParam'
            [allRewards, allChoices, ~, ~] = qLearningModel_simNoPlot('params', params(1,:),...
                'maxTrials', p.Results.maxTrials, 'randomSeed', rSeed, 'rwdProbs', p.Results.rwdProbs, 'taskType', p.Results.taskType);
        case 'fiveParam_rBeta_scale'
            [~, allRewards, allChoices, ~, ~] = qLearningModel_rBeta_simNoPlot('params', params(1,:),...
                'maxTrials', p.Results.maxTrials, 'randomSeed', rSeed, 'rwdProbs', p.Results.rwdProbs, 'taskType', p.Results.taskType);
        case 'fiveParam_rBeta_kappa'
            [~, allRewards, allChoices, ~, ~] = qLearningModel_rBetaKappa_simNoPlot('params', params(1,:),...
                'maxTrials', p.Results.maxTrials, 'randomSeed', rSeed, 'rwdProbs', p.Results.rwdProbs, 'taskType', p.Results.taskType);
        case 'sixParam_absPePeAN_bi'
            [~, allRewards, allChoices, ~, ~] = qLearningModel_absPePeAN_bi_simNoPlot('params', params(1,:),...
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
    
   
    rwdMatx_pre = [rwdMatx_pre NaN(tMax,100) rwdMatxTmp];
    noRwdMatx_pre = [noRwdMatx_pre NaN(tMax,100) noRwdMatxTmp];
    combinedAllChoice_R_pre = [combinedAllChoice_R_pre NaN(1,100) allChoice_R];
    
    allRewardsBin = allRewards;
    allRewardsBin(allRewardsBin == -1) = 1;
    allRewardsBin = allRewardsBin(1:end-1);
    changeChoice = [false abs(diff(allChoices)) > 0];
    changeChoice = changeChoice(2:end);
    s.probSwitchNoRwd_pre(runInd) = sum(changeChoice(allRewardsBin(1:end-1)==0))/sum(allRewardsBin(1:end-1)==0);
    s.probStayRwd_pre(runInd) = 1 - (sum(changeChoice(allRewardsBin(1:end-1)==1))/sum(allRewardsBin(1:end-1)==1));
end


s.ws_pre(1,1)= mean(s.probStayRwd_pre);
s.ws_pre(1,2) = std(s.probStayRwd_pre)/sqrt(length((s.probStayRwd_pre)));
s.ls_pre(1,1)= mean(s.probSwitchNoRwd_pre);
s.ls_pre(1,2) = std(s.probSwitchNoRwd_pre)/sqrt(length((s.probSwitchNoRwd_pre)));

mdl_pre = fitglm([rwdMatx_pre' noRwdMatx_pre'], combinedAllChoice_R_pre,'distribution','binomial','link','logit');


rSeed = p.Results.randomSeed;
for runInd = 1:p.Results.runs
    if rem(runInd,20) == 0
        fprintf('Running simulation %d of %d \n', runInd, p.Results.runs);
    end
    rSeed = rSeed + 1;
    switch p.Results.modelName
        case 'fourParam'
            [allRewards, allChoices, ~, ~] = qLearningModel_simNoPlot('params', params(2,:),...
                'maxTrials', p.Results.maxTrials, 'randomSeed', rSeed, 'rwdProbs', p.Results.rwdProbs, 'taskType', p.Results.taskType);
        case 'fiveParam_rBeta_scale'
            [~, allRewards, allChoices, ~, ~] = qLearningModel_rBeta_simNoPlot('params', params(2,:),...
                'maxTrials', p.Results.maxTrials, 'randomSeed', rSeed, 'rwdProbs', p.Results.rwdProbs, 'taskType', p.Results.taskType);
        case 'fiveParam_rBeta_kappa'
            [~, allRewards, allChoices, ~, ~] = qLearningModel_rBetaKappa_simNoPlot('params', params(2,:),...
                'maxTrials', p.Results.maxTrials, 'randomSeed', rSeed, 'rwdProbs', p.Results.rwdProbs, 'taskType', p.Results.taskType);
        case 'sixParam_absPePeAN_bi'
            [~, allRewards, allChoices, ~, ~] = qLearningModel_absPePeAN_bi_simNoPlot('params', params(2,:),...
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
    
   
    rwdMatx_post = [rwdMatx_post NaN(tMax,100) rwdMatxTmp];
    noRwdMatx_post = [noRwdMatx_post NaN(tMax,100) noRwdMatxTmp];
    combinedAllChoice_R_post = [combinedAllChoice_R_post NaN(1,100) allChoice_R];
    
    allRewardsBin = allRewards;
    allRewardsBin(allRewardsBin == -1) = 1;
    allRewardsBin = allRewardsBin(1:end-1);
    changeChoice = [false abs(diff(allChoices)) > 0];
    changeChoice = changeChoice(2:end);
    s.probSwitchNoRwd_post(runInd) = sum(changeChoice(allRewardsBin(1:end-1)==0))/sum(allRewardsBin(1:end-1)==0);
    s.probStayRwd_post(runInd) = 1 - (sum(changeChoice(allRewardsBin(1:end-1)==1))/sum(allRewardsBin(1:end-1)==1));
end


s.ws_post(1,1)= mean(s.probStayRwd_post);
s.ws_post(1,2) = std(s.probStayRwd_post)/sqrt(length((s.probStayRwd_post)));
s.ls_post(1,1)= mean(s.probSwitchNoRwd_post);
s.ls_post(1,2) = std(s.probSwitchNoRwd_post)/sqrt(length((s.probSwitchNoRwd_post)));

mdl_post = fitglm([rwdMatx_post' noRwdMatx_post'], combinedAllChoice_R_post,'distribution','binomial','link','logit');


colors = cool(8);

if p.Results.figFlag
    figure;
    if p.Results.compareFlag
        subplot(2,2,1); hold on;
        relevInds = 2:tMax+1;
        coefVals = actual_pre.Coefficients.Estimate(relevInds);
        CIbands = coefCI(actual_pre);
        errorL = abs(coefVals - CIbands(relevInds,1));
        errorU = abs(coefVals - CIbands(relevInds,2));
        errorbar((1:tMax)+0.2,coefVals,errorL,errorU,'Color',colors(1,:),'linewidth',2)
        
        coefVals = actual_post.Coefficients.Estimate(relevInds);
        CIbands = coefCI(actual_post);
        errorL = abs(coefVals - CIbands(relevInds,1));
        errorU = abs(coefVals - CIbands(relevInds,2));
        errorbar((1:tMax)+0.2,coefVals,errorL,errorU,'Color',colors(2,:),'linewidth',2)
        legend('pre', 'post')
        title('actual behavior')
        xlabel('Reward n Trials Back')
        ylabel('\beta Coefficient')
        xlim([0.5 tMax+0.5])
        set(gca, 'tickdir', 'out')
        
        subplot(2,2,2); hold on;
        relevInds = tMax+2:length(actual_pre.Coefficients.Estimate);
        coefVals = actual_pre.Coefficients.Estimate(relevInds);
        CIbands = coefCI(actual_pre);
        errorL = abs(coefVals - CIbands(relevInds,1));
        errorU = abs(coefVals - CIbands(relevInds,2));
        errorbar((1:tMax)+0.2,coefVals,errorL,errorU,'color',colors(5,:),'linewidth',2)
        
        coefVals = actual_post.Coefficients.Estimate(relevInds);
        CIbands = coefCI(actual_post);
        errorL = abs(coefVals - CIbands(relevInds,1));
        errorU = abs(coefVals - CIbands(relevInds,2));
        errorbar((1:tMax)+0.2,coefVals,errorL,errorU,'color',colors(6,:),'linewidth',2)
        legend('pre', 'post')
        title('actual behavior')
        xlabel('No Reward n Trials Back')
        ylabel('\beta Coefficient')
        xlim([0.5 tMax+0.5])
        set(gca, 'tickdir', 'out')
    end
    
    if p.Results.compareFlag
        subplot(2,2,3); hold on;
    else
        subplot(1,2,1); hold on;
    end
        
    relevInds = 2:tMax+1;
    coefVals = mdl_pre.Coefficients.Estimate(relevInds);
    CIbands = coefCI(mdl_pre);
    errorL = abs(coefVals - CIbands(relevInds,1));
    errorU = abs(coefVals - CIbands(relevInds,2));
    errorbar((1:tMax)+0.2,coefVals,errorL,errorU,'color',colors(3,:),'linewidth',2)

    coefVals = mdl_post.Coefficients.Estimate(relevInds);
    CIbands = coefCI(mdl_post);
    errorL = abs(coefVals - CIbands(relevInds,1));
    errorU = abs(coefVals - CIbands(relevInds,2));
    errorbar((1:tMax)+0.2,coefVals,errorL,errorU,'color',colors(4,:),'linewidth',2)
    legend('pre', 'post')
    title('simulated behavior')
    xlabel('Reward n Trials Back')
    ylabel('\beta Coefficient')
    xlim([0.5 tMax+0.5])
    set(gca, 'tickdir', 'out')

    if p.Results.compareFlag
        subplot(2,2,4); hold on;
    else
        subplot(1,2,2); hold on;
    end
    
    relevInds = tMax+2:length(mdl_pre.Coefficients.Estimate);
    coefVals = mdl_pre.Coefficients.Estimate(relevInds);
    CIbands = coefCI(mdl_pre);
    errorL = abs(coefVals - CIbands(relevInds,1));
    errorU = abs(coefVals - CIbands(relevInds,2));
    errorbar((1:tMax)+0.2,coefVals,errorL,errorU,'color',colors(7,:),'linewidth',2)
    
    coefVals = mdl_post.Coefficients.Estimate(relevInds);
    CIbands = coefCI(mdl_post);
    errorL = abs(coefVals - CIbands(relevInds,1));
    errorU = abs(coefVals - CIbands(relevInds,2));
    errorbar((1:tMax)+0.2,coefVals,errorL,errorU,'color',colors(8,:),'linewidth',2)
    legend('pre', 'post')
    title('simulated behavior')
    xlabel('No Reward n Trials Back')
    ylabel('\beta Coefficient')
    xlim([0.5 tMax+0.5])
    set(gca, 'tickdir', 'out')  
    
    t = suptitle([sheet ' ' p.Results.modelName]);
    t.Interpreter = 'none';
end


end
