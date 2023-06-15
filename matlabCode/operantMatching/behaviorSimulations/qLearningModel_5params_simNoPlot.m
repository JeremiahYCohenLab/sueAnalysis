function [t, allRewards, allChoices, blockProbs, blockSwitch] = qLearningModel_5params_simNoPlot(params, maxTrial, randomSeed, varargin)
%
%
% Simulate dynamic foraging task with temporally-forgetting Q learning model
% 
%
%task and model parameters
a = inputParser;
% default parameters if none given
a.addParameter('taskType', 'decoupled');
a.addParameter('blockLength', [20 35]);
a.addParameter('rwdProbs', [90 50 10]);
a.addParameter('ITIparam', 0.3);
a.parse(varargin{:});

alphaNPE = params(1);
alphaPPE = params(2);
alphaForget = params(3);
beta = params(4);
bias = params(5);

%initialize task class
switch a.Results.taskType
    case 'coupled'
        p = RestlessBandit('RandomSeed',randomSeed,'BlockLength', a.Results.blockLength,'MaxTrials', maxTrial,...
            'RewardProbabilities', a.Results.rwdProbs);
    case 'switch'
        p = RestlessBanditSwitch('RandomSeed', randomSeed,'BlockLength', a.Results.blockLength,'MaxTrials', maxTrial);
    case 'decoupled'
        p = RestlessBanditDecoupled('RandomSeed', randomSeed,'BlockLength', a.Results.blockLength,'MaxTrials', maxTrial,...
            'RewardProbabilities', a.Results.rwdProbs);
end
    
% [left, right]; these are Q values going INTO that trial, before making a decision
Q = [0 0; NaN(maxTrial-1, 2)]; % initialize Q values as 0
%rBar values, initialized by specified input

allChoices = ones(1, maxTrial);
allRewards = zeros(1, maxTrial);
probChosen = zeros(1, maxTrial);
for currT = 1:maxTrial
    % Select action
    pRight = logistic(beta*diff(Q(currT, :)) + bias);
    if binornd(1, pRight) == 0 % left choice selected probabilistically
        probChosen(currT) = 1 - pRight;
        p.inputChoice([1 0]);
        allChoices(currT) = -1;
        allRewards(currT) = p.AllRewards(currT, 1) * -1;
        rpe = p.AllRewards(currT, 1) - Q(currT, 1);
        pe(currT) = rpe;
        if rpe >= 0
            Q(currT + 1, 1) = Q(currT, 1) + alphaPPE*rpe;
        else
            Q(currT + 1, 1) = Q(currT, 1) + alphaNPE*rpe;
        end
        Q(currT + 1, 2) = Q(currT, 2)*alphaForget;
    else
        probChosen(currT) = pRight;
        p.inputChoice([0 1]);
        allChoices(currT) = 1;
        allRewards(currT) = p.AllRewards(currT, 2);
        rpe = p.AllRewards(currT, 2) - Q(currT, 2);
        pe(currT) = rpe;
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
        blockSwitch = blockSwitch(blockSwitch < maxTrial) + 1;
        blockSwitch(1) = 1;
        blockProbs = p.BlockProbs;
    case 'coupled'
        blockSwitch = p.BlockSwitch(p.BlockSwitch < maxTrial-1) + 1;
        for i =1:length(blockSwitch)
            if rem(i,2) == 1
                blockProbs(i,:) = a.Results.rwdProbs;
            else
                blockProbs(i,:) = fliplr(a.Results.rwdProbs);
            end
        end
    case 'switch'
        blockSwitch = p.BlockSwitch(p.BlockSwitch < maxTrial-1) + 1;
        blockProbs = p.BlockProbs;
end

t.Q = Q;
t.probChoice = probChosen;
t.pe = pe;

end