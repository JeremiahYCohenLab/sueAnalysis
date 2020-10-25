function [LH, probChoice, Q, pe, pePe, aN, peBar] = qLearningModel_8params_absPePeAN_bi_bias_k(startValues, choice, outcome)

aNscale = startValues(1);
aNmin = startValues(2);
aP = startValues(3);
aF = startValues(4);
aPE = startValues(5);
beta = startValues(6);
k = startValues(7);
if length(startValues) == 8
    bias = startValues(8);
else
    bias = 0;
end

trials = length(choice);
Q = zeros(trials+1,2);
aN = zeros(trials,1);
pePe = zeros(trials,1);
peBar = zeros(trials+1,1);

% Call learning rule
for t = 1 : trials
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


% Call softmax  rule
prevChoice_rl = [0; choice(:,1) - choice(:,2)];
prevChoice_lr = -prevChoice_rl;
probChoice = logistic([beta.*(Q(:, 1)-Q(:, 2)) + k*prevChoice_lr - bias, ...
                       beta.*(Q(:, 2)-Q(:, 1)) + k*prevChoice_rl + bias]);

% To calculate likelihood:
LH = likelihood(choice,probChoice(1:end-1,:));
end