function [LH, probChosen, m, pe, v, w, k, alpha, probChoice] = vkf_LL_kappa(startValues, choice, outcome, plotFlag)

lambda = startValues(1);  % volatility update rate
v0 = startValues(2);      % initial volatility
omega = startValues(3);   % observation noise
beta = startValues(4);    % inverse-temp parameter for softmax decision function
kappa = startValues(5);   % choice autoCorrelation

numT = length(choice);
w = omega*ones(numT, 2);
v = [v0 v0; zeros(numT-1, 2)];
m = zeros(numT, 2);
pe = zeros(numT, 1);
k = zeros(numT, 1);
alpha = zeros(numT, 1);
prevChoice = zeros(numT+1,1);
for t = 1:(length(choice)-1)
    if choice(t) == 1 % right choice selected probabilistically
        prevChoice(t+1) = 1;
        k(t) = (w(t, 2) + v(t, 2)) / (w(t, 2) + v(t, 2) + omega);
        alpha(t) = sqrt(w(t, 2) + v(t, 2));
        pe(t) = outcome(t) - logistic(m(t, 2));
        m(t+1, 2) = m(t, 2) + alpha(t) * (outcome(t) - logistic(m(t, 2)));
        w(t+1, 2) = (1 - k(t)) * (w(t, 2) + v(t, 2));
        wcov = (1 - k(t)) * w(t, 2);
        v(t+1, 2) = v(t, 2) + lambda * ( (m(t+1, 2) - m(t, 2))^2 + w(t, 2) + w(t+1, 2) - 2 * wcov - v(t, 2) );

        m(t+1, 1) = m(t, 1); 
        w(t+1, 1) = w(t, 1); 
        v(t+1, 1) = v(t, 1); 

    else
        prevChoice(t+1) = -1;
        k(t) = (w(t, 1) + v(t, 1)) / (w(t, 1) + v(t, 1) + omega);
        alpha(t) = sqrt(w(t, 1) + v(t, 1));
        pe(t) = outcome(t) - logistic(m(t, 1));
        m(t+1, 1) = m(t, 1) + alpha(t) * (outcome(t) - logistic(m(t, 1)));
        w(t+1, 1) = (1 - k(t)) * (w(t, 1) + v(t, 1));
        wcov = (1 - k(t))* w(t, 1);
        v(t+1, 1) = v(t, 1) + lambda * ( (m(t+1, 1) - m(t, 1))^2 + w(t, 1) + w(t+1, 1) - 2 * wcov - v(t, 1) );

        m(t+1, 2) = m(t, 2); 
        w(t+1, 2) = w(t, 2); 
        v(t+1, 2) = v(t, 2);

    end

end
% last trial
if choice(end) == 1
    pe(end) = outcome(end) - logistic(m(end, 2));
else
    pe(end) = outcome(end) - logistic(m(end, 1));
end
% p(right)
probChoice = logistic(beta.*(m(:,2) - m(:,1)) + kappa*prevChoice(1:end-1));

% To calculate likelihood:
if size(choice,1)==1
    LH = likelihood(choice',probChoice);
else
    LH = likelihood(choice,probChoice);
end

% probChosenChoice

probChosen = probChoice;
probChosen(choice == 0) = 1 - probChoice(choice==0);

if plotFlag
    figure; hold on
    smoothKern = ones(1,5);
    smoothKernSize = length(smoothKern);
    smoothKern = smoothKern/length(smoothKern);

    allChoices = ones(1, numT);
    allChoices(choice == 0) = -1;
    allRewards = zeros(1, numT);
    allRewards(outcome == 1 & choice == 1) = 1;
    allRewards(outcome == 1 & choice == 0) = -1;

    %smooth choices and outcomes
    smoothChoice = conv(allChoices, smoothKern);
    smoothRewards = conv(allRewards, smoothKern);

    %plot smoothed curves
    plot(smoothChoice,'k','linewidth',2);
    plot(smoothRewards,'-','Color',[30 144 255]./255,'linewidth',2)
    xlabel('Trials')
    ylabel('<-- Left       Right -->')
    xlim([1 numT])

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
