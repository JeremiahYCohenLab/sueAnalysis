function [LH, probChosen, Q, pe] = qLearningModel_5params_k(startValues, choice, outcome)

alpha = startValues(1);
aF = startValues(2);
beta = startValues(3);
k = startValues(4);
bias = startValues(5);


trials = length(choice);
Q = zeros(trials,2);
kChoice = zeros(trials,1);
pe = zeros(trials,1);
% Call learning rule
for t = 1 : (trials-1)
    if choice(t) == 1 % right choice
        Q(t+1, 1) = aF * Q(t, 1);
        pe(t) = outcome(t) - Q(t, 2);
        Q(t+1, 2) = Q(t, 2) + alpha * pe(t);
        kChoice(t+1) = k;
    else % left choice
        Q(t+1, 2) = aF * Q(t, 2);
        pe(t) = outcome(t) - Q(t, 1);
        Q(t+1, 1) = Q(t, 1) + alpha * pe(t);
        kChoice(t+1) = -k;
    end
end

if choice(t) == 1
    pe(trials) = outcome(end) - Q(end, 2);
else
    pe(trials) = outcome(end) - Q(end, 1);
end


% Call softmax  rule

probChoice = logistic(beta.*[Q(:, 2)-Q(:, 1)] + kChoice + bias);

% To calculate likelihood:
LH = likelihood(choice,probChoice);

% probChosenChoice

probChosen = probChoice;
probChosen(choice == 0) = 1 - probChoice(choice==0);
end