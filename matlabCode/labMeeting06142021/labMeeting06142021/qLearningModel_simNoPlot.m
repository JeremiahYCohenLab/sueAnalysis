function [allRewards, allChoices, blockProbs, blockSwitch] = qLearningModel_simNoPlot(varargin)
%
%
% Simulate dynamic foraging task with temporally-forgetting Q learning model
% 
%
%task and model parameters
a = inputParser;
% default parameters if none given
a.addParameter('taskType', 'decoupled');
a.addParameter('maxTrials', 1000);
a.addParameter('blockLength', [20 35]);
a.addParameter('rwdProbs', [90 50 10]);
a.addParameter('ITIparam', 0.3);
a.addParameter('params', [0.0596149,0.305917,0.642195,3.31916,0.1]);
a.addParameter('tForgetFlag', false);
a.addParameter('randomSeed', 1);
a.addParameter('biasFlag',1);
a.addParameter('plotFlag',1);
a.parse(varargin{:});

alphaNPE = a.Results.params(1);
alphaPPE = a.Results.params(2);
if a.Results.tForgetFlag == true
    tForget = a.Results.params(3);
else
    alphaForget = a.Results.params(3);
end
beta = a.Results.params(4);
if a.Results.biasFlag
    bias = a.Results.params(5);
else
    bias = 0;
end
%initialize task class
switch a.Results.taskType
    case 'coupled'
        p = RestlessBandit('RandomSeed',a.Results.randomSeed,'BlockLength', a.Results.blockLength,'MaxTrials', a.Results.maxTrials,...
            'RewardProbabilities', a.Results.rwdProbs);
    case 'switch'
        p = RestlessBanditSwitch('RandomSeed', a.Results.randomSeed,'BlockLength', a.Results.blockLength,'MaxTrials', a.Results.maxTrials);
    case 'decoupled'
        p = RestlessBanditDecoupled('RandomSeed',a.Results.randomSeed,'BlockLength', a.Results.blockLength,'MaxTrials', a.Results.maxTrials,...
            'RewardProbabilities', a.Results.rwdProbs);
end
    
% [left, right]; these are Q values going INTO that trial, before making a decision
Q = [0 0; NaN(a.Results.maxTrials-1, 2)]; % initialize Q values as 0
%rBar values, initialized by specified input

allChoices = ones(1, a.Results.maxTrials);
allRewards = zeros(1, a.Results.maxTrials);

for currT = 1:p.MaxTrials - 1    
    % Select action
    pLeft = 1/(1 + exp(-beta*diff(Q(currT, :) - bias)));
    if binornd(1, pLeft) == 0 % left choice selected probabilistically
        p.inputChoice([1 0]);
        allChoices(currT) = -1;
        allRewards(currT) = p.AllRewards(currT, 1) * -1;
        rpe = p.AllRewards(currT, 1) - Q(currT, 1);
        if rpe >= 0
            Q(currT + 1, 1) = Q(currT, 1) + alphaPPE*rpe;
        else
            Q(currT + 1, 1) = Q(currT, 1) + alphaNPE*rpe;
        end
        Q(currT + 1, 2) = Q(currT, 2)*alphaForget;
    else
        p.inputChoice([0 1]);
        allChoices(currT) = 1;
        allRewards(currT) = p.AllRewards(currT, 2);
        rpe = p.AllRewards(currT, 2) - Q(currT, 2);
        if rpe >= 0
            Q(currT + 1, 2) = Q(currT, 2) + alphaPPE*rpe;
        else
            Q(currT + 1, 2) = Q(currT, 2) + alphaNPE*rpe;
        end
        Q(currT + 1, 1) = Q(currT, 1)*alphaForget;
    end
end



switch a.Results.taskType
    case 'decoupled'
        blockSwitch = sort(unique([p.BlockSwitchL p.BlockSwitchR]));
        blockSwitch = blockSwitch(blockSwitch < p.MaxTrials-1) + 1;
        blockProbs = p.BlockProbs;
    case 'coupled'
        blockSwitch = p.BlockSwitch(p.BlockSwitch < p.MaxTrials-1) + 1;
        for i =1:length(blockSwitch)
            if rem(i,2) == 1
                blockProbs(i,:) = a.Results.rwdProbs;
            else
                blockProbs(i,:) = fliplr(a.Results.rwdProbs);
            end
        end
    case 'switch'
        blockSwitch = p.BlockSwitch(p.BlockSwitch < p.MaxTrials-1) + 1;
        blockProbs = p.BlockProbs;
end
    

if a.Results.plotFlag
    figure;
    subplot(2, 4, 1:3); hold on; ylabel('<--- L          R --->');
    plot([1:a.Results.maxTrials; 1:a.Results.maxTrials], [zeros(1,a.Results.maxTrials); (allChoices + allRewards)./2], 'k')
    subplot(2, 4, 5:7); hold on; title('Q values'); xlabel('Trials'); ylabel('Q values');
    plot(1:a.Results.maxTrials, Q(:,1), 'c', 'LineWidth', 2);
    plot(1:a.Results.maxTrials, Q(:,2), 'm', 'LineWidth', 2);
    legend('left', 'right');
    
    subplot(2,4, [4 8]); hold on;
    tMax = 10;
    allChoice_R = allChoices;
    allChoice_R(allChoice_R == -1) = 0;
    rwdMatx = [];
    for i = 1:tMax
        rwdMatx(i,:) = [NaN(1,i) allRewards(1:end-i)];
    end

    allNoRewards = allChoices;
    allNoRewards(allRewards~=0) = 0;
    noRwdMatx = [];
    for i = 1:tMax
        noRwdMatx(i,:) = [NaN(1,i) allNoRewards(1:end-i)];
    end

    glm_rwdANDnoRwd = fitglm([rwdMatx; noRwdMatx]', allChoice_R, 'distribution','binomial','link','logit'); rsq = num2str(round(glm_rwdANDnoRwd.Rsquared.Adjusted*100)/100);

    relevInds = 2:tMax+1;
    coefVals = glm_rwdANDnoRwd.Coefficients.Estimate(relevInds);
    CIbands = coefCI(glm_rwdANDnoRwd);
    errorL = abs(coefVals - CIbands(relevInds,1));
    errorU = abs(coefVals - CIbands(relevInds,2));
    errorbar((1:tMax)+0.2,coefVals,errorL,errorU,'Color','c','linewidth',2)

    relevInds = tMax+2:length(glm_rwdANDnoRwd.Coefficients.Estimate);
    coefVals = glm_rwdANDnoRwd.Coefficients.Estimate(relevInds);
    CIbands = coefCI(glm_rwdANDnoRwd);
    errorL = abs(coefVals - CIbands(relevInds,1));
    errorU = abs(coefVals - CIbands(relevInds,2));
    errorbar((1:tMax)+0.2,coefVals,errorL,errorU,'Color','m','linewidth',2)
    xlabel('Outcome n Trials Back')
    ylabel('\beta Coefficient')
    xlim([0.5 tMax+0.5])
    plot([0 tMax],[0 0],'k--')
    legend('rwd', [sprintf('\n%s\n%s%s',['no rwd'], ['R^2' rsq ' | '], ['Int: ' num2str(round(100*glm_rwdANDnoRwd.Coefficients.Estimate(1))/100)])], ...
           'location','northeast')
    
    
end

end