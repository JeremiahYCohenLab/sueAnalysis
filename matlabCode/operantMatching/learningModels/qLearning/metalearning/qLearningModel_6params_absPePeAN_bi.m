function [LH, probChoice, Q, pe, pePe, aN, peBar] = qLearningModel_6params_absPePeAN_bi(startValues, choice, outcome)

aNscale = startValues(1);
aNmin = startValues(2);
aP = startValues(3);
aF = startValues(4);
aPE = startValues(5);
beta = startValues(6);

trials = length(choice);
Q = zeros(trials,2);
aN = zeros(trials,1);
pe = zeros(trials,1);
pePe = zeros(trials,1);
peBar = zeros(trials,1);

% Call learning rule
for t = 1 : (trials-1)
    if choice(t, 1) == 1 % right choice
        Q(t+1, 2) = aF*Q(t, 2);
        pe(t) = outcome(t, 1) - Q(t, 1);
        pePe(t) = abs(pe(t)) - peBar(t);
        aN(t) = pePe(t) * aNscale + aNmin;
        if aN(t) < 0
            aN(t) = 0;
        end
        if pe(t) < 0
            Q(t+1, 1) = Q(t, 1) + aN(t) * pe(t);
        else
            Q(t+1, 1) = Q(t, 1) + aP * pe(t);
        end
    else % left choice
        Q(t+1, 1) = aF*Q(t, 1);
        pe(t) = outcome(t, 2) - Q(t, 2);
        pePe(t) = abs(pe(t)) - peBar(t);
        aN(t) = pePe(t) * aNscale + aNmin;
        if aN(t) < 0
            aN(t) = 0;
        end
        if pe(t) < 0
            Q(t+1, 2) = Q(t, 2) + aN(t) * pe(t);
        else
            Q(t+1, 2) = Q(t, 2) + aP * pe(t);
        end
    end
    peBar(t+1) = peBar(t) + aPE * pePe(t);
end

if choice(t, 1) == 1
    pe(end) = outcome(end, 1) - Q(end, 1);
else
    pe(end) = outcome(end, 2) - Q(end, 2);
end
pePe(end) =  abs(pe(end)) - peBar(end);
aN(end) = pePe(end) * aNscale + aNmin;
if aN(end) < 0
    aN(end) = 0;
end


% Call softmax  rule
probChoice = logistic([beta.*(Q(:, 1)-Q(:, 2)), ...
                       beta.*(Q(:, 2)-Q(:, 1))]);

% To calculate likelihood:
LH = likelihood(choice,probChoice);
end