function [LH, probChosen, Q, pe, peOri] = qLearningModel_5params_k_bias_LaserDisengageScale(startValues, choice, outcome, laser)

alpha = startValues(1);
aF = startValues(2);
beta = startValues(3);
k = startValues(4);
diff = startValues(5);
bias = startValues(end);


trials = length(choice);
Q = zeros(trials,2);
kChoice = zeros(trials,1);
pe = zeros(trials,1);
peOri = zeros(trials,1);

% Call learning rule
for t = 1 : (trials-1)
    if choice(t) == 1 % right choice
        Q(t+1, 1) = aF * Q(t, 1);
        pe(t) = outcome(t) - Q(t, 2);
        peOri(t) = pe(t);
        Q(t+1, 2) = Q(t, 2) + alpha * pe(t);
        kChoice(t+1) = k;
    else % left choice
        Q(t+1, 2) = aF * Q(t, 2);
        pe(t) = outcome(t) - Q(t, 1);
        peOri(t) = pe(t);
        Q(t+1, 1) = Q(t, 1) + alpha * pe(t);
        kChoice(t+1) = -k;
    end
end



% Call softmax  rule
betaCurr = beta * ones(1, length(choice));
betaCurr(intersect(find(laser==1)+1, 1:length(choice))) = diff*betaCurr(intersect(find(laser==1)+1, 1:length(choice)));
probChoice = logistic(betaCurr'.*(Q(:, 2)-Q(:, 1)) + kChoice + bias);

% To calculate likelihood:
LH = likelihood(choice,probChoice);

% probChosenChoice

probChosen = probChoice;
probChosen(choice == 0) = 1 - probChoice(choice==0);
end