function [t, allRewards, allChoices, blockProbs, blockSwitch, ITI] = qLearningModel_absPePe_scale_bias_simNoPlot(varargin)
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
a.addParameter('randomSeed', 1);
a.parse(varargin{:});

aN = a.Results.params(1);
aP = a.Results.params(2);
aF = a.Results.params(3);
aPE = a.Results.params(4);
beta = a.Results.params(5);
if length(a.Results.params) == 6
    bias = a.Results.params(6);
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

peBar = [0; NaN(a.Results.maxTrials-1, 1)];
pePe = NaN(a.Results.maxTrials, 1);
rpe = NaN(a.Results.maxTrials, 1);

allChoices = ones(1, a.Results.maxTrials);
allRewards = zeros(1, a.Results.maxTrials);
ITI = zeros(1, a.Results.maxTrials);

for currT = 1:p.MaxTrials - 1
    ITI(currT) = p.ITI;
        
    % Select action
    pRight = logistic(beta*diff(Q(currT, :)) + bias);
    if binornd(1, pRight) == 0 % left choice selected probabilistically
        p.inputChoice([1 0]);
        allChoices(currT) = -1;
        allRewards(currT) = p.AllRewards(currT, 1) * -1;
        rpe(currT) = p.AllRewards(currT, 1) - Q(currT, 1);
        pePe(currT) = abs(rpe(currT)) - peBar(currT);
        peBar(currT + 1) = peBar(currT) + aPE * pePe(currT);
        if rpe(currT) >= 0
            Q(currT + 1, 1) = Q(currT, 1) + aP*rpe(currT);
        else
            Q(currT + 1, 1) = Q(currT, 1) + aN*rpe(currT) * (1 - peBar(currT + 1));
        end
        Q(currT + 1, 2) = Q(currT, 2)*aF;
    else
        p.inputChoice([0 1]);
        allRewards(currT) = p.AllRewards(currT, 2);
        rpe(currT) = p.AllRewards(currT, 2) - Q(currT, 2);
        pePe(currT) = abs(rpe(currT)) - peBar(currT);
        peBar(currT + 1) = peBar(currT) + aPE * pePe(currT);
        if rpe(currT) >= 0
            Q(currT + 1, 2) = Q(currT, 2) + aP*rpe(currT);
        else
            Q(currT + 1, 2) = Q(currT, 2) + aN*rpe(currT) * (1 - peBar(currT + 1));
        end
        Q(currT + 1, 1) = Q(currT, 1)*aF;
    end
end


switch a.Results.taskType
    case 'decoupled'
        blockSwitch = sort(unique([p.BlockSwitchL p.BlockSwitchR]));
        blockSwitch = blockSwitch(blockSwitch < p.MaxTrials-1) + 1;
        blockSwitch(1) = 1;
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

t.Q = Q;
t.pe = rpe;
t.peBar = peBar;
t.pePe = pePe;

end