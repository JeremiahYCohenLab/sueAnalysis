function [s, outcome, choice, pRight, blockProbs, blockSwitch] = vkfSim_w0(varargin)
%task and model parameters
a = inputParser;
% default parameters if none given
a.addParameter('maxTrials', 1000);
a.addParameter('blockLength', [20 35]);
a.addParameter('rwdProbs', [90 50 10]);
a.addParameter('ITIparam', 0.3);
a.addParameter('params', [0.6,0.9,1.5,0.8]);
a.addParameter('randomSeed', 1);
a.addParameter('plotFlag', 1);
a.parse(varargin{:});




lambda = a.Results.params(1);  % volatility update rate
w0 = a.Results.params(2);      % initial volatility
omega = a.Results.params(3);   % observation noise
beta = a.Results.params(4);    % inverse-temp parameter for softmax decision function

p = RestlessBanditDecoupled('RandomSeed',a.Results.randomSeed,'BlockLength', a.Results.blockLength,'MaxTrials', a.Results.maxTrials,...
    'RewardProbabilities', a.Results.rwdProbs);
numT = a.Results.maxTrials;
v = zeros(numT+1, 2);
w = [w0 w0; zeros(numT, 2)];
m = zeros(numT+1, 2);

choice = zeros(numT,1);
outcome = zeros(numT,1);
pRight = zeros(numT,1);
for t = 1:numT
    % Select action
    mDiff = m(t, 2) -  m(t, 1);
    pRight(t) = logistic(beta*mDiff);

    if binornd(1, pRight(t)) == 1 % right choice selected probabilistically
        p.inputChoice([0 1]);
        outcome(t) = p.AllRewards(t, 2);
        choice(t) = 1;
        k = (w(t, 2) + v(t, 2)) / (w(t, 2) + v(t, 2) + omega);
        alpha = sqrt(w(t, 2) + v(t, 2));
        m(t+1, 2) = m(t, 2) + alpha * (outcome(t) - logistic(m(t, 2)));
        w(t+1, 2) = (1 - k) * (w(t, 2) + v(t, 2));
        wcov = (1 - k) * w(t, 2);
        v(t+1, 2) = v(t, 2) + lambda * ( (m(t+1, 2) - m(t, 2))^2 + w(t, 2) + w(t+1, 2) - 2 * wcov - v(t, 2) );

        m(t+1, 1) = m(t, 1); 
        w(t+1, 1) = w(t, 1); 
        v(t+1, 1) = v(t, 1); 

    else
        p.inputChoice([1 0]);
        outcome(t) = p.AllRewards(t, 1);
        choice(t) = 0;
        k = (w(t, 1) + v(t, 1)) / (w(t, 1) + v(t, 1) + omega);
        alpha = sqrt(w(t, 1) + v(t, 1));
        m(t+1, 1) = m(t, 1) + alpha * (outcome(t) - logistic(m(t, 1)));
        w(t+1, 1) = (1 - k) * (w(t, 1) + v(t, 1));
        wcov = (1 - k)* w(t, 1);
        v(t+1, 1) = v(t, 1) + lambda * ( (m(t+1, 1) - m(t, 1))^2 + w(t, 1) + w(t+1, 1) - 2 * wcov - v(t, 1) );

        m(t+1, 2) = m(t, 2); 
        w(t+1, 2) = w(t, 2); 
        v(t+1, 2) = v(t, 2);

    end
    
    s.m = m(1:numT,:);
    s.v = v(1:numT,:);
    s.w = w(1:numT,:);

end

blockSwitch = sort(unique([p.BlockSwitchL p.BlockSwitchR]));
blockSwitch = blockSwitch(blockSwitch < p.MaxTrials) + 1;
blockSwitch(1) = 1;
blockProbs = p.BlockProbs;

if a.Results.plotFlag
    figure; hold on
    smoothKern = ones(1,5);
    smoothKernSize = length(smoothKern);
    smoothKern = smoothKern/length(smoothKern);

    allChoices = ones(1, numT);
    allChoices(p.AllChoices(:,1) == 1) = -1;
    allRewards = zeros(1, numT);
    allRewards(p.AllRewards(:,1) == 1) = -1;
    allRewards(p.AllRewards(:,2) == 1) = 1;

    %smooth choices and outcomes
    smoothChoice = conv(allChoices, smoothKern);
    smoothRewards = conv(allRewards, smoothKern);

    %plot smoothed curves
    plot(smoothChoice,'k','linewidth',2);
    plot(smoothRewards,'-','Color',[30 144 255]./255,'linewidth',2)
    xlabel('Trials')
    ylabel('<-- Left       Right -->')
    xlim([1 numT])

    blockSwitch = unique(sort([p.BlockSwitchL p.BlockSwitchR]));
    blockSwitch = blockSwitch(blockSwitch < numT);

    for i = 1:length(blockSwitch)
        bs_loc = blockSwitch(i);
        plot([bs_loc bs_loc],[-1 1.5],'--k','linewidth',1)
        if rem(i,2) == 0
            labelOffset = 1.42;
        else
            labelOffset = 1.34;
        end
        a = num2str(p.BlockProbs(i,1));
        b = '/';
        c = num2str(p.BlockProbs(i,2));
        label = strcat(a,b,c);
        text(bs_loc,labelOffset,label);
        set(text,'FontSize',3);
    end

    text(0,1.12,'L/R');

    rMag = 0.3;
    nrMag = rMag/2;

    % trial plot
    for i = 1:numT
        if allChoices(i) == 1
            if allRewards(i) == 1 % R side rewarded
                plot([i i],[1 1+rMag],'k')
            else
                plot([i i],[1 1+nrMag],'color', [0.4 0.4 0.4]) % R side not rewarded
            end
        else
            if allRewards(i) == -1 % L side rewarded
                plot([i i],[-1*rMag - 1 -1],'k')
            else
                plot([i i],[-1*nrMag - 1 -1], 'color', [0.4 0.4 0.4])
            end
        end
    end

    set(gca,'tickdir', 'out')
    set(gcf, 'renderer', 'painters', 'position', [-1919          41        1920         963])
end




