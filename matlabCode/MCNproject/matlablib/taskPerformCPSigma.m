function behSim = taskPerformCPSigma(rewardsP, plotFlag,varargin)
p = inputParser;
p.addParameter('swThresh',0.65);
%p.addParameter('confThresh', 0.99);
p.parse(varargin{:});

choices = [];
rewards = [];
t = 1;
wrongTrials = 0;
extension = 0;
aPrior = [1, 1]; % prior of beta distribution, a flat line
bPrior = [1, 1];
prob = []; % prob each side need to shift to update from pior 
a = ones(size(rewardsP));
b = ones(size(rewardsP));
p_estimates = [];
changePs = [];
while t<=size(rewardsP,1)
    
    choices(t) = thompson_sampling(a(t,1),b(t,1),a(t,2),b(t,2));
    rewards(t) = binornd(1, rewardsP(t,choices(t)+1));

    %% update for meta-learning
    id = choices(t)+1;
    if t == 1
        [a(t+1,id), b(t+1,id)] = posterior_updateBi(rewards(t), aPrior(id), bPrior(id));
    else
        marg_prior = marg_likelihoodBi(rewards(t), aPrior(id), bPrior(id));
        marg_posterior = marg_likelihoodBi(rewards(t), a(t,id), b(t,id));

        prob(t, id) = marg_prior / (marg_prior + marg_posterior);

        % compute ratios
        % if marg_prior >= marg_posterior:
        if prob(t, id) >= p.Results.swThresh % 
            [a(t+1,id), b(t+1,id)] = posterior_updateBi(rewards(t), aPrior(id), bPrior(id));
            changePs(t) = 1;
        else
            [a(t+1,id), b(t+1,id)] = posterior_updateBi(rewards(t), a(t,id), b(t,id));
            changePs(t) = 0;
        end
    end
    p_estimates(t+1,id) = a(t+1,id)/(a(t+1,id)+b(t+1,id));
    
    % update unchosen side
    
    a(t+1,3-id) = a(t,3-id);
    b(t+1,3-id) = b(t,3-id);
    p_estimates(t+1,3-id) = a(t+1,3-id)/(a(t+1,3-id)+b(t+1,3-id));

    
    %%

    if sign(choices(t)-0.5) * sign(rewardsP(t,2) - rewardsP(t,1))<0 && rewards(t)==0
        wrongTrials = wrongTrials + 1;
    else
        wrongTrials = 0;
    end
    
    if wrongTrials >= 4
        extension = extension + 1;
        repeat = repmat(rewardsP(t,:),1,1);
        rewardsP = [rewardsP(1:t,:); repeat; rewardsP(t+1:end,:)];
    end
    
    t = t+1;
end

behSim.choices = choices;
behSim.rewards = rewards;
behSim.rewardsP = rewardsP;
behSim.changePs = changePs;
behSim.p_estimates = p_estimates;

if plotFlag
    figure; 
    subplot(2,1,1); hold on;
    plot(p_estimates(:,1),'color','r','linewidth',1.5);
    plot(p_estimates(:,2),'color','g','linewidth',1.5);
    legend({'red', 'green'})
    title('p estimates')
    subplot(2,1,2); hold on;
    tickLen = 1;
    for t = 1:length(choices)
        line([t t], [0 sign(choices(t)-0.5)*(0.5*tickLen + 0.5*rewards(t)*tickLen)]);
    end
    scatter(1:length(choices), changePs, 10, [1 0.3 0.3]);
    plot(rewardsP(:,2), 'color', [0.5 0.5 0.5], 'linewidth', 1.5)
    title([num2str(extension)]);
    
end
end