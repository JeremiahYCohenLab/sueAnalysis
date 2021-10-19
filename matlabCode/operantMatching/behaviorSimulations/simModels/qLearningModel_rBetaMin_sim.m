function [rBar, allRewards, allChoices, blockProbs, blockSwitch] = qLearningModel_rBetaMin_sim(varargin)
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
a.addParameter('rwdProb', [90 50 10]);
a.addParameter('ITIparam', 0.3);
a.addParameter('bestParams', [0.0596149,0.305917,0.642195,3.31916,0.1]);
a.addParameter('tForgetFlag', false);
a.addParameter('randomSeed', 1);
a.parse(varargin{:});

alphaNPE = a.Results.bestParams(1);
alphaPPE = a.Results.bestParams(2);
if a.Results.tForgetFlag == true;
    tForget = a.Results.bestParams(3);
else
    aF = a.Results.bestParams(3);
end
v = a.Results.bestParams(4);
betaScale = a.Results.bestParams(5);
betaMin = a.Results.bestParams(6);

%initialize task class
switch a.Results.taskType
    case 'coupled'
        p = RestlessBandit('RandomSeed',a.Results.randomSeed,'BlockLength', a.Results.blockLength,'MaxTrials', a.Results.maxTrials,...
            'RewardProbabilities', a.Results.rwdProb);
    case 'switch'
        p = RestlessBanditSwitch('RandomSeed', a.Results.randomSeed,'BlockLength', a.Results.blockLength,'MaxTrials', a.Results.maxTrials);
    case 'decoupled'
        p = RestlessBanditDecoupled('RandomSeed',a.Results.randomSeed,'BlockLength', a.Results.blockLength,'MaxTrials', a.Results.maxTrials,...
            'RewardProbabilities', a.Results.rwdProb);
end
    
% [left, right]; these are Q values going INTO that trial, before making a decision
Q = [0 0; NaN(a.Results.maxTrials-1, 2)]; % initialize Q values as 0
%rBar values, initialized by specified input
rBar = [0; NaN(a.Results.maxTrials-1, 1)];
beta =[0; NaN(a.Results.maxTrials-1, 1)];

% plot parameters
figure;
rawData_plot = subplot(2, 1, 1); hold on;  ylabel('<--- L          R --->');
qValue_plot = subplot(2, 1, 2); hold on; title('Q values');  xlabel('Trials'); ylabel('Q values');

allChoices = ones(1, a.Results.maxTrials);
allRewards = zeros(1, a.Results.maxTrials);

for currT = 1:p.MaxTrials - 1
        
    % Select action
    pLeft = 1/(1 + exp(-beta(currT)*diff(Q(currT, :))));
    subplot(qValue_plot)
    title(sprintf('Probability of right choice %d%%', round(pLeft*100)))
    if binornd(1, pLeft) == 0 % left choice selected probabilistically
        p = p.inputChoice([1 0]);
        rpe = p.AllRewards(currT, 1) - Q(currT, 1);
        if rpe >= 0
            Q(currT + 1, 1) = Q(currT, 1) + alphaPPE*rpe;
        else
            Q(currT + 1, 1) = Q(currT, 1) + alphaNPE*rpe;
        end
        Q(currT + 1, 2) = Q(currT, 2) * aF;
        rBar(currT + 1) = v*p.AllRewards(currT, 1) + (1-v)*rBar(currT);
    else
        p = p.inputChoice([0 1]);
        rpe = p.AllRewards(currT, 2) - Q(currT, 2);
        if rpe >= 0
            Q(currT + 1, 2) = Q(currT, 2) + alphaPPE*rpe;
        else
            Q(currT + 1, 2) = Q(currT, 2) + alphaNPE*rpe;
        end
        Q(currT + 1, 1) = Q(currT, 1) * aF;
        rBar(currT + 1) = v*p.AllRewards(currT, 2) + (1-v)*rBar(currT);
    end
    beta(currT + 1) = rBar(currT + 1) * betaScale + betaMin;
    
    ITI = exprnd(1/a.Results.ITIparam);
    if a.Results.tForgetFlag == true
        Q(currT + 1, :) = Q(currT + 1, :)*exp(-tForget*ITI);
    else
        Q(currT + 1, :) = Q(currT + 1, :)*aF;
    end

    subplot(rawData_plot); xlim([0 currT+1]); ylim([-1 1])
    title(sprintf('ITI of %2.1f seconds', round(ITI*10)/10))
    if p.AllChoices(currT, 1) == 1 % left choice
        allChoices(currT) = -1;
        if p.AllRewards(currT, 1) == 1 % reward
            plot([currT, currT], [0 -1], 'k')
            allRewards(currT) = -1;
        else
            plot([currT, currT], [0 -0.5], 'k')
        end
    elseif p.AllChoices(currT, 2) == 1 % right choice
        if p.AllRewards(currT, 2) == 1 % reward
            plot([currT, currT], [0 1], 'k')
            allRewards(currT) = 1;
        else
            plot([currT, currT], [0 0.5], 'k')
        end
    end    
    
    % plot block switches with reward probabilities
    if p.BlockSwitch_Flag == true
        plot([currT currT], [-1 1], '--c');
        if strcmp(a.Results.taskType, 'decoupled')
            if rem(length([p.BlockSwitchL p.BlockSwitchR]), 2) == 0
                labelOffset = 1.12;
            else
                labelOffset = 1.04;
            end    
        else
            labelOffset = 1.04;
        end
        label = [num2str(p.RewardProbabilities(1,1)) '/' num2str(p.RewardProbabilities(1,2))];
        text(currT,labelOffset,label);
        set(text,'FontSize',3);
    end
    
    subplot(qValue_plot); xlim([0 currT + 1]); ylim([-0.1 1.1]);
    plot(Q(1:currT + 1, 1),'c')
    plot(Q(1:currT + 1, 2),'m')
    legend('Left','Right')
end
subplot(qValue_plot); ylim([0 max(max(Q))]);
suptitle('Q learning opponency model simulated behavior');


switch a.Results.taskType
    case 'decoupled'
        blockSwitch = sort(unique([p.BlockSwitchL p.BlockSwitchR]));
        blockSwitch = blockSwitch(blockSwitch < p.MaxTrials-1) + 1;
        blockProbs = p.BlockProbs;
    case 'coupled'
        blockSwitch = p.BlockSwitch(p.BlockSwitch < p.MaxTrials-1) + 1;
        for i =1:length(blockSwitch)
            if rem(i,2) == 1
                blockProbs(i,:) = a.Results.rwdProb;
            else
                blockProbs(i,:) = fliplr(a.Results.rwdProb);
            end
        end
    case 'switch'
        blockSwitch = p.BlockSwitch(p.BlockSwitch < p.MaxTrials-1) + 1;
        blockProbs = p.BlockProbs;
end
    
% tMax = 10;
% allChoice_R = allChoices;
% allChoice_R(allChoice_R == -1) = 0;
% rwdMatx = [];
% for i = 1:tMax
%     rwdMatx(i,:) = [NaN(1,i) allRewards(1:end-i)];
% end
% 
% allNoRewards = allChoices;
% allNoRewards(allRewards~=0) = 0;
% noRwdMatx = [];
% for i = 1:tMax
%     noRwdMatx(i,:) = [NaN(1,i) allNoRewards(1:end-i)];
% end
% 
% glm_rwdANDnoRwd = fitglm([rwdMatx; noRwdMatx]', allChoice_R, 'distribution','binomial','link','logit'); rsq = num2str(round(glm_rwdANDnoRwd.Rsquared.Adjusted*100)/100);
% 
% figure; hold on
% relevInds = 2:tMax+1;
% coefVals = glm_rwdANDnoRwd.Coefficients.Estimate(relevInds);
% CIbands = coefCI(glm_rwdANDnoRwd);
% errorL = abs(coefVals - CIbands(relevInds,1));
% errorU = abs(coefVals - CIbands(relevInds,2));
% errorbar((1:tMax)+0.2,coefVals,errorL,errorU,'Color','c','linewidth',2)
% 
% relevInds = tMax+2:length(glm_rwdANDnoRwd.Coefficients.Estimate);
% coefVals = glm_rwdANDnoRwd.Coefficients.Estimate(relevInds);
% CIbands = coefCI(glm_rwdANDnoRwd);
% errorL = abs(coefVals - CIbands(relevInds,1));
% errorU = abs(coefVals - CIbands(relevInds,2));
% errorbar((1:tMax)+0.2,coefVals,errorL,errorU,'Color','m','linewidth',2)
% xlabel('Outcome n Trials Back')
% ylabel('\beta Coefficient')
% legend('rwd', [sprintf('\n%s\n%s%s',['no rwd'], ['R^2' rsq ' | '], ['Int: ' num2str(round(100*glm_rwdANDnoRwd.Coefficients.Estimate(1))/100)])], ...
%        'location','northeast')
% xlim([0.5 tMax+0.5])
% plot([0 tMax],[0 0],'k--')

end