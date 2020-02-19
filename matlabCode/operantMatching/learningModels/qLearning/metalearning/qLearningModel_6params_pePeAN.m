function [LH, probChoice, Q, pe, pePe, aN, peBar] = qLearningModel_6params_pePeAN(startValues, choice, outcome)

aNscale = startValues(1);
aNmin = startValues(2);
aP = startValues(3);
aF = startValues(4);
aPE = startValues(5);
beta = startValues(6);


trials = length(choice);
Q = zeros(trials,2);
aN = [aNmin; zeros(trials-1,1)];
peBar = zeros(trials,1);
pePe = zeros(trials,1);

% Call learning rule
for t = 1 : (trials-1)
    if choice(t, 1) == 1 % right choice
        Q(t+1, 2) = aF*Q(t, 2);
        pe(t) = outcome(t, 1) - Q(t, 1);
        if pe(t) < 0
            Q(t+1, 1) = Q(t, 1) + aN(t) * pe(t);
        else
            Q(t+1, 1) = Q(t, 1) + aP * pe(t);
        end
    else % left choice
        Q(t+1, 1) = aF*Q(t, 1);
        pe(t) = outcome(t, 2) - Q(t, 2);
        if pe(t) < 0
            Q(t+1, 2) = Q(t, 2) + aN(t) * pe(t);
        else
            Q(t+1, 2) = Q(t, 2) + aP * pe(t);
        end
    end
    pePe(t) = abs(pe(t)) - peBar(t);
    peBar(t+1) = peBar(t) + aPE * (pePe(t));
    aN(t+1) = (1 - peBar(t+1)) * aNscale + aNmin;
end

if choice(t, 1) == 1
    pe(trials) = outcome(end, 1) - Q(end, 1);
else
    pe(trials) = outcome(end, 2) - Q(end, 2);
end
pePe(trials) = pe(trials) - peBar(trials-1);


% Call softmax  rule

probChoice = logistic([beta.*(Q(:, 1)-Q(:, 2)), ...
                       beta.*(Q(:, 2)-Q(:, 1))]);

% To calculate likelihood:
LH = likelihood(choice,probChoice);
end