function [t, allRewards, allChoices, blockProbs, blockSwitch, laser] = qLearningModel_5paramsLaserNegRPE_simNoPlot(params, maxTrial, randomSeed, varargin)
%
% Simulate dynamic foraging task with temporally-forgetting Q learning model
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

alphaNPE = params(1);
alphaPPE = params(2);
alphaForget = params(3);
beta = params(4);
rpeShift = params(5);
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

allChoices = zeros(maxTrial,1);
allRewards = zeros(maxTrial,1);
laser = zeros(maxTrial,1);
pe = zeros(maxTrial,1);
peChanged = zeros(maxTrial,1);

for currT = 1:maxTrial
    % Select action
    pRight(currT) = logistic(beta*(Q(currT, 2) - Q(currT, 1)) + bias);
    if binornd(1, a.Results.laser) == 1
        laser(currT) = 1;
    end
    if binornd(1, pRight(currT)) == 0 % left choice selected probabilistically
        p.inputChoice([1 0]);
        allChoices(currT) = -1;
        allRewards(currT) = p.AllRewards(currT, 1) * -1;
        rpe = p.AllRewards(currT, 1) - Q(currT, 1);
        pe(currT) = rpe;
        peChanged(currT) = rpe;
        if rpe < 0
            if laser(currT) == 1
                peChanged(currT) = rpe + rpeShift;
                rpe = rpe + rpeShift;
                if rpe > 0
                    rpe = 0;
                end
            end
            Q(currT + 1, 1) = Q(currT, 1) + alphaNPE*rpe;
        else
            Q(currT + 1, 1) = Q(currT, 1) + alphaPPE*rpe;
        end       
        Q(currT + 1, 2) = Q(currT, 2)*alphaForget;
    else
        p.inputChoice([0 1]);
        allChoices(currT) = 1;
        allRewards(currT) = p.AllRewards(currT, 2);
        rpe = p.AllRewards(currT, 2) - Q(currT, 2);
        pe(currT) = rpe;
        peChanged(currT) = rpe;
        if rpe < 0
           if laser(currT) == 1
                peChanged(currT) = rpe + rpeShift;
                rpe = rpe + rpeShift;
                if rpe > 0
                    rpe = 0;
                end
           end
           Q(currT + 1, 2) = Q(currT, 2) + alphaNPE*rpe;
        else
           Q(currT + 1, 2) = Q(currT, 2) + alphaPPE*rpe;
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
t.pe = pe;
t.peChanged = peChanged;
t.pRight = pRight;

end