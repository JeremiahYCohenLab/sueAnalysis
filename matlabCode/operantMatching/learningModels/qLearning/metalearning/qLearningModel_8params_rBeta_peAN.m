function [LH, probChoice, Q, pe, beta, R, peBar] = qLearningModel_8params_rBeta_peAN(startValues, choice, outcome)

aNscale = startValues(1);
aNmin = startValues(2);
aP = startValues(3);
aF = startValues(4);
aPE = startValues(5);
v = startValues(6);
betaScale = startValues(7);
betaMin = startValues(8);


trials = length(choice);
Q = zeros(trials,2);
beta = [betaMin; zeros(trials-1,1)];
R = zeros(trials,1);
aN = [aNmin; zeros(trials-1,1)];
peBar = zeros(trials,1);


% Call learning rule
for t = 1 : (trials-1)
    if choice(t, 1) == 1 % right choice
        Q(t+1, 2) = aF*Q(t, 2);
        pe(t) = outcome(t, 1) - Q(t, 1);
        peBar(t+1) = aPE * abs(pe(t)) + (1-aPE) * peBar(t);
        aN(t+1)= peBar(t+1) * aNscale + aNmin;
        if pe(t) < 0
            Q(t+1, 1) = Q(t, 1) + aN(t+1) * pe(t);
        else
            Q(t+1, 1) = Q(t, 1) + aP * pe(t);
        end
        R(t+1) = v * outcome(t,1) + (1-v)*R(t);
    else % left choice
        Q(t+1, 1) = aF*Q(t, 1);
        pe(t) = outcome(t, 2) - Q(t, 2);
        peBar(t+1) = aPE * abs(pe(t)) + (1-aPE) * peBar(t);
        aN(t+1)= peBar(t+1) * aNscale + aNmin;
        if pe(t) < 0
            Q(t+1, 2) = Q(t, 2) + aN(t+1) * pe(t);
        else
            Q(t+1, 2) = Q(t, 2) + aP * pe(t);
        end
        R(t+1) = v * outcome(t,2) + (1-v)*R(t);
    end
    beta(t+1) = R(t+1) * betaScale + betaMin;
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