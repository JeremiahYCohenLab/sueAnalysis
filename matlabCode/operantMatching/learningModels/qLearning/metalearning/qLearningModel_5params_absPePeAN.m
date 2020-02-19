function [LH, probChoice, Q, pe, pePe, aN, peBar] = qLearningModel_5params_absPePeAN(startValues, choice, outcome)

aNscale = startValues(1);
aP = startValues(2);
aF = startValues(3);
aPE = startValues(4);
beta = startValues(5);


trials = length(choice);
Q = zeros(trials,2);
aN = zeros(trials,1);
peBar = zeros(trials,1);
pePe = zeros(trials,1);

% Call learning rule
for t = 1 : (trials-1)
    if choice(t, 1) == 1 % right choice
        Q(t+1, 2) = aF*Q(t, 2);
        pe(t) = outcome(t, 1) - Q(t, 1);
        pePe(t) = abs(peBar(t) - abs(pe(t)));
        aN(t+1) = pePe(t) * aNscale;
        if pe(t) < 0
            Q(t+1, 1) = Q(t, 1) + aN(t+1) * pe(t);
        else
            Q(t+1, 1) = Q(t, 1) + aP * pe(t);
        end
    else % left choice
        Q(t+1, 1) = aF*Q(t, 1);
        pe(t) = outcome(t, 2) - Q(t, 2);
        pePe(t) = abs(peBar(t) - abs(pe(t)));
        aN(t+1) = pePe(t) * aNscale;
        if pe(t) < 0
            Q(t+1, 2) = Q(t, 2) + aN(t+1) * pe(t);
        else
            Q(t+1, 2) = Q(t, 2) + aP * pe(t);
        end
    end
    peBar(t+1) = peBar(t) + aPE * (abs(pe(t)) - peBar(t));
end

if choice(t, 1) == 1
    pe(trials) = outcome(end, 1) - Q(end, 1);
else
    pe(trials) = outcome(end, 2) - Q(end, 2);
end
pePe(trials) = abs(peBar(trials) - abs(pe(trials)));
aN = aN(2:end);
aN(trials) = pePe(trials) * aNscale;

% Call softmax  rule

probChoice = logistic([beta.*(Q(:, 1)-Q(:, 2)), ...
                       beta.*(Q(:, 2)-Q(:, 1))]);

% To calculate likelihood:
LH = likelihood(choice,probChoice);
end