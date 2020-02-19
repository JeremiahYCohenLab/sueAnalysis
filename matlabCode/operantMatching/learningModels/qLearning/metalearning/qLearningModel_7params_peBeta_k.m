function [LH, probChoice, Q, pe, beta, peBar] = qLearningModel_6params_peBeta(startValues, choice, outcome)

aN = startValues(1);
aP = startValues(2);
aF = startValues(3);
v = startValues(4);
betaMin = startValues(5);
betaScale = startValues(6);
k = startValues(7);

trials = length(choice);
Q = zeros(trials,2);
beta = zeros(trials,1);
peBar = zeros(trials,1);

% Call learning rule
for t = 1 : (trials-1)
    if choice(t, 1) == 1 % right choice
        Q(t+1, 2) = aF*Q(t, 2);
        pe(t) = outcome(t, 1) - Q(t, 1);
        if pe(t) < 0
            Q(t+1, 1) = Q(t, 1) + aN * pe(t);
        else
            Q(t+1, 1) = Q(t, 1) + aP * pe(t);
        end
    else % left choice
        Q(t+1, 1) = aF*Q(t, 1);
        pe(t) = outcome(t, 2) - Q(t, 2);
        if pe(t) < 0
            Q(t+1, 2) = Q(t, 2) + aN * pe(t);
        else
            Q(t+1, 2) = Q(t, 2) + aP * pe(t);
        end
    end
    peBar(t+1) = v * abs(pe(t)) + (1-v) * peBar(t);
    beta(t+1) = betaMin + betaScale * peBar(t);
end

if choice(t, 1) == 1
    pe(trials) = outcome(end, 1) - Q(end, 1);
else
    pe(trials) = outcome(end, 2) - Q(end, 2);
end


% Call softmax  rule
prevChoice_rl = [0; choice(1:end-1,1) - choice(1:end-1,2)];
prevChoice_lr = -prevChoice_rl;
probChoice = logistic([beta.*(Q(:, 1)-Q(:, 2)) + k*prevChoice_lr, ...
                       beta.*(Q(:, 2)-Q(:, 1)) + k*prevChoice_rl]);

% To calculate likelihood:
LH = likelihood(choice,probChoice);
end