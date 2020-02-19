function [LH, probChoice, Q, pe, pePe, peBar] = qLearningModel_7params_absPePeLR(startValues, choice, outcome)

aNscale = startValues(1);
aNmin = startValues(2);
aPscale = startValues(3);
aPmin = startValues(4);
aF = startValues(5);
aPE = startValues(6);
beta = startValues(7);

trials = length(choice);
Q = zeros(trials,2);
peBar = zeros(trials,1);
pePe = zeros(trials,1);

% Call learning rule
for t = 1 : (trials-1)
    if choice(t, 1) == 1 % right choice
        Q(t+1, 2) = aF*Q(t, 2);
        pe(t) = outcome(t, 1) - Q(t, 1);
        pePe(t) = abs(pe(t)) - peBar(t);
        if pe(t) < 0
            aN = pePe(t) * aNscale + aNmin;
            Q(t+1, 1) = Q(t, 1) + aN * pe(t);
        else
            aP = pePe(t) * aPscale + aPmin;
            Q(t+1, 1) = Q(t, 1) + aP * pe(t);
        end
    else % left choice
        Q(t+1, 1) = aF*Q(t, 1);
        pe(t) = outcome(t, 2) - Q(t, 2);
        pePe(t) = abs(pe(t)) - peBar(t);
        if pe(t) < 0
            aN = pePe(t) * aNscale + aNmin;
            Q(t+1, 2) = Q(t, 2) + aN * pe(t);
        else
            aP = pePe(t) * aPscale + aPmin;
            Q(t+1, 2) = Q(t, 2) + aP * pe(t);
        end
    end
    peBar(t+1) = peBar(t) + aPE * pePe(currT);
end

if choice(t, 1) == 1
    pe(trials) = outcome(end, 1) - Q(end, 1);
else
    pe(trials) = outcome(end, 2) - Q(end, 2);
end
pePe(trials) = abs(peBar(trials) - abs(pe(trials)));

% Call softmax  rule

probChoice = logistic([beta.*(Q(:, 1)-Q(:, 2)), ...
                       beta.*(Q(:, 2)-Q(:, 1))]);

% To calculate likelihood:
LH = likelihood(choice,probChoice);
end