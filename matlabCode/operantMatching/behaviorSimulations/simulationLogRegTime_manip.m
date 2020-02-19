function [mdl_pre, s] = simulationLogRegTime_manip(xlFile, sheet, pre, post, varargin)


%task and model parameters
p = inputParser;
% default parameters if none given
p.addParameter('revForFlag', 0);
p.addParameter('binSize', 6000);
p.addParameter('numBins', 10);
p.addParameter('modelName', 'sixParam_absPePeAN_bi');
p.addParameter('params', []);
p.addParameter('taskType', 'coupled');
p.addParameter('rwdProbs', [90 50 10]);
p.addParameter('runs', 1000);
p.addParameter('maxTrials', 300);
p.addParameter('randomSeed', 9872);
p.addParameter('figFlag', 1);
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

rwdMatx_pre = [];
noRwdMatx_pre = [];
combinedAllChoice_R_pre = [];

rwdMatx_post = [];
noRwdMatx_post = [];
combinedAllChoice_R_post = [];

if p.Results.compareFlag
    [actual_pre, t] = combineLogRegTime_opMD(xlFile, sheet, pre, 'revForFlag', p.Results.revForFlag,...
       'binSize', p.Results.binSize, 'numBins', p.Results.numBins);
    [actual_post, ~] = combineLogRegTime_opMD(xlFile, sheet, post, 'revForFlag', p.Results.revForFlag,...
       'binSize', p.Results.binSize, 'numBins', p.Results.numBins);
end

rSeed = p.Results.randomSeed;
for runInd = 1:p.Results.runs
    if rem(runInd,20) == 0
        fprintf('Running simulation %d of %d \n', runInd, p.Results.runs);
    end
    rSeed = rSeed + 1;
    switch p.Results.modelName
        case 'fourParam'
            [allRewards, allChoices, ~, ~,ITI] = qLearningModel_simNoPlot('params', params(1,:),...
                'maxTrials', p.Results.maxTrials, 'randomSeed', rSeed, 'rwdProbs', p.Results.rwdProbs, 'taskType', p.Results.taskType);
        case 'fiveParam_rBeta_scale'
            [~, allRewards, allChoices, ~, ~,ITI] = qLearningModel_rBeta_simNoPlot('params', params(1,:),...
                'maxTrials', p.Results.maxTrials, 'randomSeed', rSeed, 'rwdProbs', p.Results.rwdProbs, 'taskType', p.Results.taskType);
        case 'fiveParam_rBeta_kappa'
            [~, allRewards, allChoices, ~, ~,ITI] = qLearningModel_rBetaKappa_simNoPlot('params', params(1,:),...
                'maxTrials', p.Results.maxTrials, 'randomSeed', rSeed, 'rwdProbs', p.Results.rwdProbs, 'taskType', p.Results.taskType);
        case 'sixParam_absPePeAN_bi'
            [~, allRewards, allChoices, ~, ~,ITI] = qLearningModel_absPePeAN_bi_simNoPlot('params', params(1,:),...
                'maxTrials', p.Results.maxTrials, 'randomSeed', rSeed, 'rwdProbs', p.Results.rwdProbs);
    end
    
    ITI = ITI + 4; %add time for trial, fixed ITI, and no lick window
    rTime = zeros(1, length(ITI));
    for currI = 2:length(ITI)
        rTime(currI) = ITI(currI-1) + rTime(currI-1);
    end
    rTime = rTime*1000;
        
    allNoRewards = allChoices;
    allNoRewards(find(allRewards == 1)) = 0;
    allNoRewards(find(allRewards == -1)) = 0;
    
    allChoice_R = allChoices;
    allChoice_R(find(allChoice_R == -1)) = 0;
    
    %create binned outcome matrices
    rwdTmpMatx = NaN(t.tMax, length(allRewards));     %initialize matrices for number of response trials x number of time bins
    noRwdTmpMatx = NaN(t.tMax, length(allRewards));
    for j = 2:length(allRewards)          
        k = 1;
        %find time between "current" choice and previous rewards, up to timeMax in the past 
        timeTmpL = []; timeTmpR = []; nTimeTmpL = []; nTimeTmpR = [];
        while j-k > 0 & rTime(j) - rTime(j-k) < t.timeMax
            if allRewards(j-k) == -1
                timeTmpL = [timeTmpL rTime(j) - rTime(j-k)];
            end
            if allRewards(j-k) == 1
                timeTmpR = [timeTmpR rTime(j) - rTime(j-k)];
            end
            if allNoRewards(j-k) == -1
                nTimeTmpL = [nTimeTmpL rTime(j) - rTime(j-k)];
            end
            if allNoRewards(j-k) == 1
                nTtimeTmpR = [nTimeTmpR rTime(j) - rTime(j-k)];
            end
            k = k + 1;
        end
        %bin outcome times and use to fill matrices
        if ~isempty(timeTmpL)
            binnedRwds = discretize(timeTmpL, t.timeBinEdges);
            for k = 1:t.tMax
                if ~isempty(binnedRwds == k)
                    rwdTmpMatx(k,j) = -1*sum(binnedRwds == k);
                else
                    rwdTmpMatx(k,j) = 0;
                end
            end
        end
        if ~isempty(timeTmpR)
            binnedRwds = discretize(timeTmpR, t.timeBinEdges);
            for k = 1:t.tMax
                if ~isempty(binnedRwds == k) & isnan(rwdTmpMatx(k,j))
                    rwdTmpMatx(k,j) = sum(binnedRwds == k);
                elseif ~isempty(binnedRwds == k) & ~isnan(rwdTmpMatx(k,j))
                    rwdTmpMatx(k,j) = rwdTmpMatx(k,j) + sum(binnedRwds == k);
                else
                    rwdTmpMatx(k,j) = 0;
                end
            end
        end
        if isempty(timeTmpL) & isempty(timeTmpR)
            rwdTmpMatx(:,j) = 0;
        end
        if ~isempty(nTimeTmpL)
            binnedNoRwds = discretize(nTimeTmpL, t.timeBinEdges);
            for k = 1:t.tMax
                if ~isempty(binnedNoRwds == k)
                    noRwdTmpMatx(k,j) = -1*sum(binnedNoRwds == k);
                else
                    noRwdTmpMatx(k,j) = 0;
                end
            end
        end
        if ~isempty(nTimeTmpR)
            binnedNoRwds = discretize(nTimeTmpR, t.timeBinEdges);
            for k = 1:t.tMax
                if ~isempty(binnedNoRwds == k) & isnan(noRwdTmpMatx(k,j))
                    noRwdTmpMatx(k,j) = sum(binnedNoRwds == k);
                elseif ~isempty(binnedNoRwds == k) & ~isnan(noRwdTmpMatx(k,j))
                    noRwdTmpMatx(k,j) = noRwdTmpMatx(k,j) + sum(binnedNoRwds == k);
                else
                    noRwdTmpMatx(k,j) = 0;
                end
            end
        end
        if isempty(nTimeTmpL) & isempty(nTimeTmpR)
            noRwdTmpMatx(:,j) = 0;
        end
    end
    
    %concatenate temp matrix with combined matrix
    rwdTmpMatx(:,1) = NaN;
    rwdMatx_pre = [rwdMatx_pre NaN(length(t.timeBinEdges)-1, 100) rwdTmpMatx];
    noRwdMatx_pre = [noRwdMatx_pre NaN(length(t.timeBinEdges)-1, 100) noRwdTmpMatx];
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
            [allRewards, allChoices, ~, ~,ITI] = qLearningModel_simNoPlot('params', params(2,:),...
                'maxTrials', p.Results.maxTrials, 'randomSeed', rSeed, 'rwdProbs', p.Results.rwdProbs, 'taskType', p.Results.taskType);
        case 'fiveParam_rBeta_scale'
            [~, allRewards, allChoices, ~, ~,ITI] = qLearningModel_rBeta_simNoPlot('params', params(2,:),...
                'maxTrials', p.Results.maxTrials, 'randomSeed', rSeed, 'rwdProbs', p.Results.rwdProbs, 'taskType', p.Results.taskType);
        case 'fiveParam_rBeta_kappa'
            [~, allRewards, allChoices, ~, ~,ITI] = qLearningModel_rBetaKappa_simNoPlot('params', params(2,:),...
                'maxTrials', p.Results.maxTrials, 'randomSeed', rSeed, 'rwdProbs', p.Results.rwdProbs, 'taskType', p.Results.taskType);
        case 'sixParam_absPePeAN_bi'
            [~, allRewards, allChoices, ~, ~,ITI] = qLearningModel_absPePeAN_bi_simNoPlot('params', params(2,:),...
                'maxTrials', p.Results.maxTrials, 'randomSeed', rSeed, 'rwdProbs', p.Results.rwdProbs);
    end
    
    ITI = ITI + 4; %add time for trial, fixed ITI, and no lick window
    rTime = zeros(1, length(ITI));
    for currI = 2:length(ITI)
        rTime(currI) = ITI(currI-1) + rTime(currI-1);
    end
    rTime = rTime*1000;
        
    allNoRewards = allChoices;
    allNoRewards(find(allRewards == 1)) = 0;
    allNoRewards(find(allRewards == -1)) = 0;
    
    allChoice_R = allChoices;
    allChoice_R(find(allChoice_R == -1)) = 0;
    
    %create binned outcome matrices
    rwdTmpMatx = NaN(t.tMax, length(allRewards));     %initialize matrices for number of response trials x number of time bins
    noRwdTmpMatx = NaN(t.tMax, length(allRewards));
    for j = 2:length(allRewards)          
        k = 1;
        %find time between "current" choice and postvious rewards, up to timeMax in the past 
        timeTmpL = []; timeTmpR = []; nTimeTmpL = []; nTimeTmpR = [];
        while j-k > 0 & rTime(j) - rTime(j-k) < t.timeMax
            if allRewards(j-k) == -1
                timeTmpL = [timeTmpL rTime(j) - rTime(j-k)];
            end
            if allRewards(j-k) == 1
                timeTmpR = [timeTmpR rTime(j) - rTime(j-k)];
            end
            if allNoRewards(j-k) == -1
                nTimeTmpL = [nTimeTmpL rTime(j) - rTime(j-k)];
            end
            if allNoRewards(j-k) == 1
                nTtimeTmpR = [nTimeTmpR rTime(j) - rTime(j-k)];
            end
            k = k + 1;
        end
        %bin outcome times and use to fill matrices
        if ~isempty(timeTmpL)
            binnedRwds = discretize(timeTmpL, t.timeBinEdges);
            for k = 1:t.tMax
                if ~isempty(binnedRwds == k)
                    rwdTmpMatx(k,j) = -1*sum(binnedRwds == k);
                else
                    rwdTmpMatx(k,j) = 0;
                end
            end
        end
        if ~isempty(timeTmpR)
            binnedRwds = discretize(timeTmpR, t.timeBinEdges);
            for k = 1:t.tMax
                if ~isempty(binnedRwds == k) & isnan(rwdTmpMatx(k,j))
                    rwdTmpMatx(k,j) = sum(binnedRwds == k);
                elseif ~isempty(binnedRwds == k) & ~isnan(rwdTmpMatx(k,j))
                    rwdTmpMatx(k,j) = rwdTmpMatx(k,j) + sum(binnedRwds == k);
                else
                    rwdTmpMatx(k,j) = 0;
                end
            end
        end
        if isempty(timeTmpL) & isempty(timeTmpR)
            rwdTmpMatx(:,j) = 0;
        end
        if ~isempty(nTimeTmpL)
            binnedNoRwds = discretize(nTimeTmpL, t.timeBinEdges);
            for k = 1:t.tMax
                if ~isempty(binnedNoRwds == k)
                    noRwdTmpMatx(k,j) = -1*sum(binnedNoRwds == k);
                else
                    noRwdTmpMatx(k,j) = 0;
                end
            end
        end
        if ~isempty(nTimeTmpR)
            binnedNoRwds = discretize(nTimeTmpR, t.timeBinEdges);
            for k = 1:t.tMax
                if ~isempty(binnedNoRwds == k) & isnan(noRwdTmpMatx(k,j))
                    noRwdTmpMatx(k,j) = sum(binnedNoRwds == k);
                elseif ~isempty(binnedNoRwds == k) & ~isnan(noRwdTmpMatx(k,j))
                    noRwdTmpMatx(k,j) = noRwdTmpMatx(k,j) + sum(binnedNoRwds == k);
                else
                    noRwdTmpMatx(k,j) = 0;
                end
            end
        end
        if isempty(nTimeTmpL) & isempty(nTimeTmpR)
            noRwdTmpMatx(:,j) = 0;
        end
    end
    
    %concatenate temp matrix with combined matrix
    rwdTmpMatx(:,1) = NaN;
    rwdMatx_post = [rwdMatx_post NaN(length(t.timeBinEdges)-1, 100) rwdTmpMatx];
    noRwdMatx_post = [noRwdMatx_post NaN(length(t.timeBinEdges)-1, 100) noRwdTmpMatx];
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
        relevInds = 2:t.tMax+1;
        coefVals = actual_pre.Coefficients.Estimate(relevInds);
        CIbands = coefCI(actual_pre);
        errorL = abs(coefVals - CIbands(relevInds,1));
        errorU = abs(coefVals - CIbands(relevInds,2));
        errorbar(((1:t.tMax)*t.binSize/1000),coefVals,errorL,errorU,'Color',colors(1,:),'linewidth',2)
        
        coefVals = actual_post.Coefficients.Estimate(relevInds);
        CIbands = coefCI(actual_post);
        errorL = abs(coefVals - CIbands(relevInds,1));
        errorU = abs(coefVals - CIbands(relevInds,2));
        errorbar(((1:t.tMax)*t.binSize/1000),coefVals,errorL,errorU,'Color',colors(2,:),'linewidth',2)
        legend('pre', 'post')
        title('actual behavior')
        xlabel('Reward n seconds back')
        ylabel('\beta Coefficient')
        xlim([0 t.tMax*t.binSize/1000 + t.binSize/1000])
        set(gca, 'tickdir', 'out')
        
        subplot(2,2,2); hold on;
        relevInds =  t.tMax+2:t.tMax*2+1;
        coefVals = actual_pre.Coefficients.Estimate(relevInds);
        CIbands = coefCI(actual_pre);
        errorL = abs(coefVals - CIbands(relevInds,1));
        errorU = abs(coefVals - CIbands(relevInds,2));
        errorbar(((1:t.tMax)*t.binSize/1000),coefVals,errorL,errorU,'color',colors(5,:),'linewidth',2)
        
        coefVals = actual_post.Coefficients.Estimate(relevInds);
        CIbands = coefCI(actual_post);
        errorL = abs(coefVals - CIbands(relevInds,1));
        errorU = abs(coefVals - CIbands(relevInds,2));
        errorbar(((1:t.tMax)*t.binSize/1000),coefVals,errorL,errorU,'color',colors(6,:),'linewidth',2)
        legend('pre', 'post')
        title('actual behavior')
        xlabel('No reward n seconds back')
        ylabel('\beta Coefficient')
        xlim([0 t.tMax*t.binSize/1000 + t.binSize/1000])
        set(gca, 'tickdir', 'out')
    end
    
    if p.Results.compareFlag
        subplot(2,2,3); hold on;
    else
        subplot(1,2,1); hold on;
    end
        
    relevInds = 2:t.tMax+1;
    coefVals = mdl_pre.Coefficients.Estimate(relevInds);
    CIbands = coefCI(mdl_pre);
    errorL = abs(coefVals - CIbands(relevInds,1));
    errorU = abs(coefVals - CIbands(relevInds,2));
    errorbar(((1:t.tMax)*t.binSize/1000),coefVals,errorL,errorU,'color',colors(3,:),'linewidth',2)

    coefVals = mdl_post.Coefficients.Estimate(relevInds);
    CIbands = coefCI(mdl_post);
    errorL = abs(coefVals - CIbands(relevInds,1));
    errorU = abs(coefVals - CIbands(relevInds,2));
    errorbar(((1:t.tMax)*t.binSize/1000),coefVals,errorL,errorU,'color',colors(4,:),'linewidth',2)
    legend('pre', 'post')
    title('simulated behavior')
    xlabel('Reward n seconds back')
    ylabel('\beta Coefficient')
    xlim([0 t.tMax*t.binSize/1000 + t.binSize/1000])
    set(gca, 'tickdir', 'out')

    if p.Results.compareFlag
        subplot(2,2,4); hold on;
    else
        subplot(1,2,2); hold on;
    end
    
    relevInds = t.tMax+2:t.tMax*2+1;
    coefVals = mdl_pre.Coefficients.Estimate(relevInds);
    CIbands = coefCI(mdl_pre);
    errorL = abs(coefVals - CIbands(relevInds,1));
    errorU = abs(coefVals - CIbands(relevInds,2));
    errorbar(((1:t.tMax)*t.binSize/1000),coefVals,errorL,errorU,'color',colors(7,:),'linewidth',2)
    
    coefVals = mdl_post.Coefficients.Estimate(relevInds);
    CIbands = coefCI(mdl_post);
    errorL = abs(coefVals - CIbands(relevInds,1));
    errorU = abs(coefVals - CIbands(relevInds,2));
    errorbar(((1:t.tMax)*t.binSize/1000),coefVals,errorL,errorU,'color',colors(8,:),'linewidth',2)
    legend('pre', 'post')
    title('simulated behavior')
    xlabel('No reward n seconds back')
    ylabel('\beta Coefficient')
    xlim([0 t.tMax*t.binSize/1000 + t.binSize/1000])
    set(gca, 'tickdir', 'out')  
    
    t = suptitle([sheet ' ' p.Results.modelName]);
    t.Interpreter = 'none';
end


end
