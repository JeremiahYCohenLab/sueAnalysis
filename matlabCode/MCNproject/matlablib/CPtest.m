% simulate environment and a RL-like agent to do task
[rewardPs] = taskSimulate(0.9, 0.2, 400, [60, 60],1);
beh = taskPerform(rewardPs,1,'sigma',0.07);
%% beginning of hacky bayesian online change point detection
rewardPs = beh.rewardsP;
aPrior = [1, 1]; % prior of beta distribution, a flat line
bPrior = [1, 1];
T = size(rewardPs,1);
prob = zeros(size(rewardPs,2),2); % prob each side need to shift to update from pior 
phase_estimates = zeros(size(rewardPs,2),1);
p_estimates = zeros(size(rewardPs));
a = ones(size(rewardPs));
b = ones(size(rewardPs));
for t = 1:T
    if beh.choices(t) > 0
        id = 2;
    else
        id = 1;
    end
    
    if t == 1
        [a(t,id), b(t,id)] = posterior_updateBi(beh.rewards(t), aPrior(id), bPrior(id));
    else
        marg_prior = marg_likelihoodBi(beh.rewards(t), aPrior(id), bPrior(id));
        marg_posterior = marg_likelihoodBi(beh.rewards(t), a(t-1,id), b(t-1,id));

        prob(t, id) = marg_prior / (marg_prior + marg_posterior);

        % compute ratios
        % if marg_prior >= marg_posterior:
        if prob(t, id) >= 0.60 % i tried different values here, in case of 0.6, it worked a bit in the beginning of the session
            fprintf(['hi' '\n'])
            [a(t,id), b(t,id)] = posterior_updateBi(beh.rewards(t), aPrior(id), bPrior(id));
            phase_estimates(t) = phase_estimates(t - 1) + 1;
        else
            [a(t,id), b(t,id)] = posterior_updateBi(beh.rewards(t), a(t-1,id), b(t-1,id));
            phase_estimates(t) = phase_estimates(t - 1);
        end
    end
    p_estimates(t,id) = a(t,id)/(a(t,id)+b(t,id));
    if t>1
        p_estimates(t,3-id) = p_estimates(t-1,3-id);
        a(t,3-id) = a(t-1,3-id);
        b(t,3-id) = b(t-1,3-id);
    else
        p_estimates(t,3-id) = a(t,3-id)/(a(t,3-id)+b(t,3-id));
    end
end
%%
figure; hold on;
subplot(2,1,1); hold on;
plot(1:T, rewardPs(:,1));
plot(1:T, p_estimates(:,1),'linestyle','--','color', 'm');
plot(1:T, p_estimates(:,2),'linestyle','--','color', 'c');

subplot(2,1,2); hold on;
scatter(1:T, [0; diff(phase_estimates)], 10, 'filled', 'k');
%%
