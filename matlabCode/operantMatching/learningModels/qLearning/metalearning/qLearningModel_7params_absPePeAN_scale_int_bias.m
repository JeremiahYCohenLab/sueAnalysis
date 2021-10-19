function [LH, probChosen, Q, pe, pePe, aN, peBar] = qLearningModel_7params_absPePeAN_scale_int_bias(startValues, choice, outcome)

aNmin = startValues(1);
aP = startValues(2);
aF = startValues(3);
aPE = startValues(4);
v = startValues(5);
beta = startValues(6);
if length(startValues) == 7
    bias = startValues(7);
else
    bias = 0;
end

trials = length(choice);
Q = zeros(trials+1,2);
aN = [aNmin; zeros(trials,1)];
pe = zeros(trials,1);
pePe = zeros(trials,1);
peBar = zeros(trials+1,1);

% Call learning rule
for t = 1 : trials
    if choice(t) == 1 % right choice
        Q(t+1, 1) = aF*Q(t, 1);
        pe(t) = outcome(t) - Q(t, 2);
        pePe(t) = abs(pe(t)) - peBar(t);
        if pe(t) < 0
            aN(t+1) = v * (pePe(t) + aNmin) + (1 - v) * aN(t);
            if aN(t+1) < 0
                aN(t+1) = 0;
            end
            if aN(t+1) > 1
                aN(t+1) = 1;
            end
            Q(t+1, 2) = Q(t, 2) + aN(t+1) * pe(t) * (1 - peBar(t));
        else
            Q(t+1, 2) = Q(t, 2) + aP * pe(t) * (1 - peBar(t));
            aN(t+1) = aN(t);
        end
    else % left choice
        Q(t+1, 2) = aF*Q(t, 2);
        pe(t) = outcome(t) - Q(t, 1);
        pePe(t) = abs(pe(t)) - peBar(t);
        if pe(t) < 0
            aN(t+1) = v * (pePe(t) + aNmin) + (1 - v) * aN(t);
            if aN(t+1) < 0
                aN(t+1) = 0;
            end
            if aN(t+1) > 1
                aN(t+1) = 1;
            end
            Q(t+1, 1) = Q(t, 1) + aN(t+1) * pe(t) * (1 - peBar(t));
        else
            Q(t+1, 1) = Q(t, 1) + aP * pe(t) * (1 - peBar(t));
            aN(t+1) = aN(t);
        end
    end
    peBar(t+1) = peBar(t) + aPE * pePe(t);
end

aN = aN(2:end);
Q = Q(1:end-1,:);
peBar = peBar(1:end-1);
% Call softmax  rule, p(right)
probChoice = logistic(beta.*(Q(:, 2)-Q(:, 1)) + bias);

% To calculate likelihood:
LH = likelihood(choice,probChoice);

% probChosenChoice

probChosen = probChoice;
probChosen(choice == 0) = 1 - probChoice(choice==0);

end