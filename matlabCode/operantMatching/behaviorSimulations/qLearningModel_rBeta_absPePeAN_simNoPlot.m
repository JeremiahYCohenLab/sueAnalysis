function [peBar, allRewards, allChoices, blockProbs, blockSwitch] = qLearningModel_rBeta_absPePeAN_simNoPlot(varargin)
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
a.parse(varargin{:});

aNscale = a.Results.params(1);
aNmin = a.Results.params(2);
aP = a.Results.params(3);
aF = a.Results.params(4);
aPE = a.Results.params(5);
v = a.Results.params(6);
betaScale = a.Results.params(7);
betaMin = a.Results.params(8);

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
peBar = [0; NaN(a.Results.maxTrials-1, 1)];
rBar = [0; NaN(a.Results.maxTrials-1, 1)];
aN = aNmin;
beta = [betaMin; NaN(a.Results.maxTrials-1, 1)];

allChoices = ones(1, a.Results.maxTrials);
allRewards = zeros(1, a.Results.maxTrials);

for currT = 1:p.MaxTrials - 1
        
    % Select action
    pLeft = 1/(1 + exp(-beta(currT)*diff(Q(currT, :))));
    if binornd(1, pLeft) == 0 % left choice selected probabilistically
        p.inputChoice([1 0]);
        allChoices(currT) = -1;
        allRewards(currT) = p.AllRewards(currT, 1) * -1;
        rpe = p.AllRewards(currT, 1) - Q(currT, 1);
        if rpe >= 0
            Q(currT + 1, 1) = Q(currT, 1) + aP*rpe;
        else
            Q(currT + 1, 1) = Q(currT, 1) + aN*rpe;
        end
        Q(currT + 1, 2) = Q(currT, 2)*aF;
        rBar(currT + 1) = v*p.AllRewards(currT, 1) + (1-v)*rBar(currT);
    else
        p.inputChoice([0 1]);
        allChoices(currT) = 1;
        allRewards(currT) = p.AllRewards(currT, 2);
        rpe = p.AllRewards(currT, 2) - Q(currT, 2);
        peBar(currT + 1) = v * abs(rpe) + (1-v) * peBar(currT);
        aN = peBar(currT + 1) * aNscale + aNmin;
        if rpe >= 0
            Q(currT + 1, 2) = Q(currT, 2) + aP*rpe;
        else
            Q(currT + 1, 2) = Q(currT, 2) + aN*rpe;
        end
        Q(currT + 1, 1) = Q(currT, 1)*aF;
        rBar(currT + 1) = v*p.AllRewards(currT, 2) + (1-v)*rBar(currT);
    end
    pePe = abs(abs(rpe) - peBar(currT));
    aN = (1 - pePe) * aNscale + aNmin;
    peBar(currT + 1) = peBar(currT) + aPE * (abs(rpe) - peBar(currT));
    beta(currT + 1) = rBar(currT + 1) * betaScale + betaMin;
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
    


end