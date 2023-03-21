function [t, allRewards, allChoices, blockProbs, blockSwitch, laser] = qLearningModel_5params_k_bias_LaserDisengage_simNoPlot(params, maxTrial, randomSeed, varargin)
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
a.addParameter('laser', 0.3);
a.parse(varargin{:});

alpha = params(1);
alphaForget = params(2);
beta = params(3);
kappa = params(4);
scale = params(5);
bias = params(end);


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
laser = zeros(1, maxTrial);
pe = zeros(maxTrial,1);
pRight = zeros(maxTrial, 1);

prevChoice = 0;
for currT = 1:maxTrial
    % Select action
    if currT > 1
        if laser(currT-1) == 1
            pRight(currT) = logistic(scale*beta*(Q(currT, 2) - Q(currT, 1)) + bias + kappa*prevChoice);
        else
            pRight(currT) = logistic(beta*(Q(currT, 2) - Q(currT, 1)) + bias + kappa*prevChoice);
        end
    else
        pRight(currT) = logistic(beta*(Q(currT, 2) - Q(currT, 1)) + bias + kappa*prevChoice);
    end
            
    if binornd(1, a.Results.laser) == 1
        laser(currT) = 1;
    end    
    if binornd(1, pRight(currT)) == 0 % left choice selected probabilistically
        p.inputChoice([1 0]);
        allChoices(currT) = -1;
        allRewards(currT) = p.AllRewards(currT, 1) * -1;
        rpe = p.AllRewards(currT, 1) - Q(currT, 1);
        pe(currT) = rpe;
        Q(currT + 1, 1) = Q(currT, 1) + alpha*rpe;
        Q(currT + 1, 2) = Q(currT, 2)*alphaForget;
        prevChoice = -1;
    else
        p.inputChoice([0 1]);
        allChoices(currT) = 1;
        allRewards(currT) = p.AllRewards(currT, 2);
        rpe = p.AllRewards(currT, 2) - Q(currT, 2);
        pe(currT) = rpe;
        Q(currT + 1, 2) = Q(currT, 2) + alpha*rpe;
        Q(currT + 1, 1) = Q(currT, 1)*alphaForget;
        prevChoice = 1;
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
t.pe = pe;
t.pRight = pRight;

end