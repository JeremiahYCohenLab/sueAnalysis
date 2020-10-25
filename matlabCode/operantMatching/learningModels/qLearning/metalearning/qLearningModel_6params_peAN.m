function [LH, probChoice, Q, pe, pePe, aN, peBar] = qLearningModel_6params_peAN(startValues, choice, outcome)

aNscale = startValues(1);
aNmin = startValues(2);
aP = startValues(3);
aF = startValues(4);
aPE = startValues(5);
beta = startValues(6);

trials = length(choice);
Q = zeros(trials+1,2);
aN = zeros(trials,1);
pe = zeros(trials,1);
pePe = zeros(trials,1);
peBar = zeros(trials+1,1);

% Call learning rule
for t = 1 : trials
    if choice(t, 1) == 1 % right choice
        Q(t+1, 2) = aF*Q(t, 2);
        pe(t) = outcome(t, 1) - Q(t, 1);
        pePe(t) = abs(pe(t)) - peBar(t);
        peBar(t+1) = peBar(t) + aPE * pePe(t);
        aN(t) = peBar(t+1) * aNscale + aNmin;
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
        peBar(t+1) = peBar(t) + aPE * pePe(t);
        aN(t) = peBar(t+1) * aNscale + aNmin;
        if aN(t) < 0
            aN(t) = 0;
        end
        if pe(t) < 0
            Q(t+1, 2) = Q(t, 2) + aN(t) * pe(t);
        else
            Q(t+1, 2) = Q(t, 2) + aP * pe(t);
        end
    end
end


% Call softmax  rule
probChoice = logistic([beta.*(Q(:, 1)-Q(:, 2)), ...
                       beta.*(Q(:, 2)-Q(:, 1))]);

% To calculate likelihood:
LH = likelihood(choice,probChoice(1:end-1,:));
end