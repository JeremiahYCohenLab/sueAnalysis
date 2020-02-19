function [LH, probChoice, Q, pe, R, aN, aP] = qLearningModel_7params_peLR(startValues, choice, outcome)

aNscale = startValues(1);
aNmin = startValues(2);
aPscale = startValues(3);
aPmin = startValues(4);
aF = startValues(5);
v = startValues(6);
beta = startValues(7);


trials = length(choice);
Q = zeros(trials,2);
aN = [aNmin; zeros(trials-1,1)];
aP = [aPmin; zeros(trials-1,1)];
R = zeros(trials,1);

% Call learning rule
for t = 1 : (trials-1)
    if choice(t, 1) == 1 % right choice
        Q(t+1, 2) = aF*Q(t, 2);
        pe(t) = outcome(t, 1) - Q(t, 1);
        if pe(t) < 0
            Q(t+1, 1) = Q(t, 1) + aN(t) * pe(t);
        else
            Q(t+1, 1) = Q(t, 1) + aP(t) * pe(t);
        end
    else % left choice
        Q(t+1, 1) = aF*Q(t, 1);
        pe(t) = outcome(t, 2) - Q(t, 2);
        if pe(t) < 0
            Q(t+1, 2) = Q(t, 2) + aN(t) * pe(t);
        else
            Q(t+1, 2) = Q(t, 2) + aP(t) * pe(t);
        end
    end
    R(t+1) = v * abs(pe(t)) + (1-v) * R(t);
    aN(t+1) = R(t+1) * aNscale + aNmin;
    aP(t+1) = R(t+1) * aPscale + aPmin;
end

if choice(t, 1) == 1
    pe(trials) = outcome(end, 1) - Q(end, 1);
else
    pe(trials) = outcome(end, 2) - Q(end, 2);
end


% Call softmax  rule

probChoice = logistic([beta.*(Q(:, 1)-Q(:, 2)), ...
                       beta.*(Q(:, 2)-Q(:, 1))]);

% To calculate likelihood:
LH = likelihood(choice,probChoice);
end