function [allRewards, allChoices, blockProbs, blockSwitch, correctArray] = optimizedVKF_simNoPlot(varargin)
%
%
% Simulate dynamic foraging task with temporally-forgetting Q learning model
% 
%
%task and model parameters
a = inputParser;
% default parameters if none given
a.addParameter('taskType', 'decoupled')
a.addParameter('maxTrials', 1000)
a.addParameter('blockLength', [20 35])
a.addParameter('rwdProbs', [90 50 10])
a.addParameter('ITIparam', 0.3)
a.addParameter('tForgetFlag', false)
a.addParameter('randomSeed', 1)
a.addParameter('plotFlag', 0)
a.parse(varargin{:});

lambda = 0.549245864306393;  % volatility update rate
v0 = 0.422311804752216;      % initial volatility
omega = 0.975697922267773;   % observation noise
beta = 1.341230225800368;    % inverse-temp parameter for softmax decision function

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
correctArray = zeros(1, a.Results.maxTrials);


w = zeros(a.Results.maxTrials+1, 2);
v = [v0 v0; zeros(a.Results.maxTrials, 2)];
m = zeros(a.Results.maxTrials+1, 2);

for t = 1:p.MaxTrials
    % Select action
    mDiff = m(t, 2) -  m(t, 1);
    pRight = logistic(beta*mDiff);
    choiceProb(t,:) = [1-pRight pRight];

    if binornd(1, pRight) == 0 % left choice selected probabilistically
        p.inputChoice([1 0]);
        o = p.AllRewards(t, 1);
        allChoices(t) = -1;
        allRewards(t) = p.AllRewards(t, 1) * -1;
        
        k = (w(t, 1) + v(t, 1)) / (w(t, 1) + v(t, 1) + omega);
        alpha = sqrt(w(t, 1) + v(t, 1));
        m(t+1, 1) = m(t, 1) + alpha * (o - logistic(m(t, 1)));
        w(t+1, 1) = (1 - k) * (w(t, 1) + v(t, 1));
        wcov = (1 - k) * w(t, 1);
        v(t+1, 1) = v(t, 1) + lambda * ( (m(t+1, 1) - m(t, 1))^2 + w(t, 1) + w(t+1, 1) - 2 * wcov - v(t, 1) );

        m(t+1, 2) = m(t, 2); 
        w(t+1, 2) = w(t, 2); 
        v(t+1, 2) = v(t, 2); 

        if p.RewardProbabilities(1) >= p.RewardProbabilities(2)
            correctArray(t) = 1;
        end

    else
        p.inputChoice([0 1]);
        o = p.AllRewards(t, 2);
        allChoices(t) = 1;
        allRewards(t) = p.AllRewards(t, 2);

        k = (w(t, 2) + v(t, 2)) / (w(t, 2) + v(t, 2) + omega);
        alpha = sqrt(w(t, 2) + v(t, 2));
        m(t+1, 2) = m(t, 2) + alpha * (o - logistic(m(t, 2)));
        w(t+1, 2) = (1 - k) * (w(t, 2) + v(t, 2));
        wcov = (1 - k)* w(t, 2);
        v(t+1, 2) = v(t, 2) + lambda * ( (m(t+1, 2) - m(t, 2))^2 + w(t, 2) + w(t+1, 2) - 2 * wcov - v(t, 2) );

        m(t+1, 1) = m(t, 1); 
        w(t+1, 1) = w(t, 1); 
        v(t+1, 1) = v(t, 1);

        if p.RewardProbabilities(2) >= p.RewardProbabilities(1)
            correctArray(t) = 1;
        end
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
    figure; hold on
    
    plot(v(:,1), '-c', 'linewidth', 1.5)
    plot(v(:,2), '-m', 'linewidth', 1.5)
    legend('volatility L', 'volatility R')
    yU = max(max(v));
    yL = min(min(v));
    
    for i = 1:length(blockSwitch)
        bs_loc = blockSwitch(i);
        plot([bs_loc bs_loc],[yL yU+0.5],'--k','linewidth',1)
        if rem(i,2) == 0
            labelOffset = yU + 0.42;
        else
            labelOffset = yU + 0.34;
        end
        b = num2str(blockProbs(i,1));
        c = '/';
        d = num2str(blockProbs(i,2));
        label = strcat(b,c,d);
        text(bs_loc,labelOffset,label);
        set(text,'FontSize',3);
    end

    text(0,1.12,'L/R');

    rMag = 0.3;
    nrMag = rMag/2;

    % trial plot
    j = 1;
    for i = 1:length(allChoices)
        if allChoices(i) == 1
            if allRewards(i) == 1 % R side rewarded
                plot([i i],[yU yU+rMag],'k')
            else
                plot([i i],[yU yU+nrMag],'color', [0.4 0.4 0.4]) % R side not rewarded
            end
        else
            if allRewards(i) == -1 % L side rewarded
                plot([i i],[-1*rMag + yL yL],'k')
            else
                plot([i i],[-1*nrMag + yL yL], 'color', [0.4 0.4 0.4])
            end
        end
    end

    xlim([0 a.Results.maxTrials])
    ylim([-1*rMag + yL yU+0.5])
    set(gca,'tickdir', 'out')
    set(gcf, 'renderer', 'painters', 'position', [-1697 322 1461 531])
end

end