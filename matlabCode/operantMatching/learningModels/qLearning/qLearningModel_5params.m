function [LH, probChosen, Q, pe] = qLearningModel_5params(startValues, choice, outcome)

alphaNPE = startValues(1);
alphaPPE = startValues(2);
alphaForget = startValues(3);
beta = startValues(4);
bias = startValues(5);

trials = length(choice);
Q = zeros(trials,2);
pe = zeros(trials,1);
% Call learning rule
for t = 1 : (trials-1)
    if choice(t) == 1 % right choice
        Q(t+1, 1) = alphaForget*Q(t, 1);
        pe(t) = outcome(t) - Q(t, 2);
        if pe(t) < 0
            Q(t+1, 2) = Q(t, 2) + alphaNPE * pe(t);
        else
            Q(t+1, 2) = Q(t, 2) + alphaPPE * pe(t);
        end
    else % left choice
        Q(t+1, 2) = alphaForget*Q(t, 2);
        pe(t) = outcome(t) - Q(t, 1);
        if pe(t) < 0
            Q(t+1, 1) = Q(t, 1) + alphaNPE * pe(t);
        else
            Q(t+1, 1) = Q(t, 1) + alphaPPE * pe(t);
        end
    end
end

if choice(t, 1) == 1
    pe(trials) = outcome(end) - Q(end, 2);
else
    pe(trials) = outcome(end) - Q(end, 1);
end

% Call softmax  rule

probChoice = logistic(beta*(Q(:, 2)-Q(:, 1)) + bias);

% To calculate likelihood:
LH = likelihood(choice,probChoice);

% probChosenChoice

probChosen = probChoice;
probChosen(choice == 0) = 1 - probChoice(choice==0);  
end