function [LH, probChoice, Q, pe, cQ] = qLearningModel_4params(startValues, choice, outcome)

alphaNPE = startValues(1);
alphaPPE = startValues(2);
alphaForget = startValues(3);
beta = startValues(4);

trials = length(choice);
Q = zeros(trials+1,2);
pe = zeros(trials,1);
cQ = zeros(trials+1,1);

% Call learning rule
for t = 1 : (trials)
    if choice(t) == 1 % right choice
        cQ(t) = Q(t,2);
        Q(t+1, 1) = alphaForget*Q(t, 1);
        pe(t) = outcome(t) - Q(t, 2);
        if pe(t) < 0
            Q(t+1, 2) = Q(t, 2) + alphaNPE * pe(t);
        else
            Q(t+1, 2) = Q(t, 2) + alphaPPE * pe(t);
        end
    else % left choice
        cQ(t) = Q(t,1);
        Q(t+1, 2) = alphaForget*Q(t, 2);
        pe(t) = outcome(t) - Q(t, 1);
        if pe(t) < 0
            Q(t+1, 1) = Q(t, 1) + alphaNPE * pe(t);
        else
            Q(t+1, 1) = Q(t, 1) + alphaPPE * pe(t);
        end
    end
end

Q = Q(1:end-1,:);
% Call softmax  rule

probChoice = logistic(beta*(Q(:, 2)-Q(:, 1)));

% To calculate likelihood:
LH = likelihood(choice',probChoice);
end