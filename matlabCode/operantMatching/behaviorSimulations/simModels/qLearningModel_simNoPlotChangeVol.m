function [t, allRewards, allChoices, blockProbs, blockSwitch, vol] = qLearningModel_simNoPlotChangeVol(varargin)
%
%
% Simulate dynamic foraging task with temporally-forgetting Q learning model
% 
%
%task and model parameters
a = inputParser;
% default parameters if none given
a.addParameter('taskType', 'coupled');
a.addParameter('maxTrials', 2000);
a.addParameter('blockLengthL', [40 50]);
a.addParameter('blockLengthH', [10 20]);
a.addParameter('rwdProbs', [70 10]);
a.addParameter('ITIparam', 0.3);
a.addParameter('params', [0.0596149,0.705917,0.342195,3.31916,0.1]);
a.addParameter('initState', 1);
a.addParameter('randomSeed', 1);
a.parse(varargin{:});

alphaNPE = a.Results.params(1);
alphaPPE = a.Results.params(2);
alphaForget = a.Results.params(3);
beta = a.Results.params(4);
if length(a.Results.params) == 5
    bias = a.Results.params(5);
else
    bias = 0;
end

if a.Results.initState==1
    blockLength = a.Results.blockLengthH;
else
    blockLength = a.Results.blockLengthL;
end

%initialize task class
switch a.Results.taskType
    case 'coupled'
        p = RestlessBandit('RandomSeed',a.Results.randomSeed,'BlockLength', blockLength,'MaxTrials', a.Results.maxTrials,...
            'RewardProbabilities', a.Results.rwdProbs);
    case 'switch'
        p = RestlessBanditSwitch('RandomSeed', a.Results.randomSeed,'BlockLength', blockLength,'MaxTrials', a.Results.maxTrials);
    case 'decoupled'
        p = RestlessBanditDecoupled('RandomSeed',a.Results.randomSeed,'BlockLength', blockLength,'MaxTrials', a.Results.maxTrials,...
            'RewardProbabilities', a.Results.rwdProbs);
end
    
% [left, right]; these are Q values going INTO that trial, before making a decision
Q = [0 0; NaN(a.Results.maxTrials-1, 2)]; % initialize Q values as 0
%rBar values, initialized by specified input

allChoices = ones(1, a.Results.maxTrials);
allRewards = zeros(1, a.Results.maxTrials);
rwdProbs = zeros(a.Results.maxTrials,2);
vol = zeros(1, a.Results.maxTrials);
currVol = a.Results.initState;
switched = false;

for currT = 1:p.MaxTrials
    % Select action
    pRight = logistic(beta*diff(Q(currT, :)) + bias);
    if binornd(1, pRight) == 0 % left choice selected probabilistically
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
    
    rwdProbs(currT,:) = p.RewardProbabilities;
    if p.NewBlock_Flag && currT >= 0.5*a.Results.maxTrials && ~switched
        if a.Results.initState==0
            p.BlockLength = a.Results.blockLengthH;
        else
            p.BlockLength = a.Results.blockLengthL;
        end
        currVol = 1-currVol;
        switched = true;
    end
        
    vol(currT) = currVol;
end



switch a.Results.taskType
    case 'decoupled'
        blockSwitch = sort(unique([p.BlockSwitchL p.BlockSwitchR]));
        blockSwitch = blockSwitch(blockSwitch < p.MaxTrials) + 1;
        blockSwitch(1) = 1;
        blockProbs = p.BlockProbs;
    case 'coupled'
        blockSwitch = p.BlockSwitch(p.BlockSwitch < p.MaxTrials-1) + 1;
        blockSwitch(1) = 1;
        blockProbs = rwdProbs;
    case 'switch'
        blockSwitch = p.BlockSwitch(p.BlockSwitch < p.MaxTrials-1) + 1;
        blockProbs = p.BlockProbs;
end

t.Q = Q;

end