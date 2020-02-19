function [LH, probChoice, Q, pe, pePe, aN, peBar, peBar_L, peBar_R] = qLearningModel_6params_absPePeAN_biSep(startValues, choice, outcome)

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
peBar_L = zeros(trials,1);
peBar_R = zeros(trials,1);

% Call learning rule
for t = 1 : (trials-1)
    if choice(t, 1) == 1 % right choice
        Q(t+1, 2) = aF*Q(t, 2);
        pe(t) = outcome(t, 1) - Q(t, 1);
        pePe(t) = abs(pe(t)) - peBar_R(t);
        aN(t) = pePe(t) * aNscale + aNmin;
        if pe(t) < 0
            Q(t+1, 1) = Q(t, 1) + aN(t) * pe(t);
        else
            Q(t+1, 1) = Q(t, 1) + aP * pe(t);
        end
        peBar_R(t+1) = peBar_R(t) + aPE * pePe(t);
        peBar_L(t+1) = peBar_L(t);
        peBar(t) = peBar_R(t);
    else % left choice
        Q(t+1, 1) = aF*Q(t, 1);
        pe(t) = outcome(t, 2) - Q(t, 2);
        pePe(t) = abs(pe(t)) - peBar_L(t);
        aN(t) = pePe(t) * aNscale + aNmin;
        if pe(t) < 0
            Q(t+1, 2) = Q(t, 2) + aN(t) * pe(t);
        else
            Q(t+1, 2) = Q(t, 2) + aP * pe(t);
        end
        peBar_L(t+1) = peBar_L(t) + aPE * pePe(t);
        peBar_R(t+1) = peBar_R(t);
        peBar(t) = peBar_L(t);
    end
end

if choice(t, 1) == 1
    pe(end) = outcome(end, 1) - Q(end, 1);
else
    pe(end) = outcome(end, 2) - Q(end, 2);
end
pePe(end) =  abs(pe(end)) - peBar(end);
aN(end) = pePe(end) * aNscale + aNmin;


% Call softmax  rule
probChoice = logistic([beta.*(Q(:, 1)-Q(:, 2)), ...
                       beta.*(Q(:, 2)-Q(:, 1))]);

% To calculate likelihood:
LH = likelihood(choice,probChoice);
end