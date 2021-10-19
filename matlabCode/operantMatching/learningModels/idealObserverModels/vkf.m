function [s,choice,outcome] = vkfSim(varargin)
%task and model parameters
a = inputParser;
% default parameters if none given
a.addParameter('maxTrials', 1000);
a.addParameter('blockLength', [20 35]);
a.addParameter('rwdProbs', [90 50 10]);
a.addParameter('ITIparam', 0.3);
a.addParameter('params', [0.0596149,0.305917,0.642195,3.31916,0.1]);
a.addParameter('randomSeed', 1);
a.parse(varargin{:});




lambda = a.params(1);  % volatility update rate
v0 = a.params(2);      % initial volatility
omega = a.params(3);   % observation noise
beta = a.params(4);    % inverse-temp parameter for softmax decision function

p = RestlessBanditDecoupled('RandomSeed',a.Results.randomSeed,'BlockLength', a.Results.blockLength,'MaxTrials', a.Results.maxTrials,...
    'RewardProbabilities', a.Results.rwdProbs);
numT = a.Results.maxTrials;
w = zeros(numT+1, 2);
v = [v0 v0; zeros(numT, 2)];
m = zeros(numT+1, 2);

choice = zeros(numT,1);
outcome = zeros(numT,1);

for t = 1:numT
    % Select action
    mDiff = m(t, 2) -  m(t, 1);
    pRight = logistic(beta*mDiff);

    if binornd(1, pRight) == 1 % right choice selected probabilistically
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


