function behSim = taskPerformDynSigma(rewardsP, plotFlag,varargin)
p = inputParser;
p.addParameter('sigma', 0.10);
p.addParameter('aP', 0.06);
p.addParameter('aN', 0.06);
p.addParameter('aF', 0.90);
p.addParameter('aConf', 0.25);
p.addParameter('aFConf', 0.9);
p.addParameter('aSup',1);
%p.addParameter('confThresh', 0.99);
p.parse(varargin{:});

c = [0 0];
choices = [];
pChoices = [];
rewards = [];
t = 1;
wrongTrials = 0;
extension = 0;
choiceConf = [];
sigmaCurr = [];
confBar = 0;
confSup = 0;
a = 0.15;
while t<=size(rewardsP,1)
    if t > 1
        sigmaCurr(t) = FIcurve(confSup(t-1), 'a', a, 'b', a*p.Results.aSup, 'd', 25) + p.Results.sigma;
    else
        sigmaCurr(t) = p.Results.sigma;
    end
    pChoices(t) = 1/(1+exp(-(c(t,2)-c(t,1))/sigmaCurr(t)));
    choices(t) = binornd(1,pChoices(t));
    if choices(t) > 0
        rewards(t) = binornd(1,rewardsP(t,2));
        if rewards(t) > 0
            c(t+1,2) = p.Results.aF*c(t,2) + p.Results.aP*(1 - c(t,2));
        else
            c(t+1,2) = p.Results.aF*c(t,2) + p.Results.aN*(0-c(t,2));
        end
        c(t+1,1) = p.Results.aF*c(t,1);
        choiceConf(t) = pChoices(t);
    else
        rewards(t) = binornd(1,rewardsP(t,1));
        if rewards(t) > 0
            c(t+1,1) = p.Results.aF*c(t,1) + p.Results.aP*(1 - c(t,1));
        else
            c(t+1,1) = p.Results.aF*c(t,1) + p.Results.aN*(0-c(t,1));
        end
        c(t+1,2) = p.Results.aF*c(t,2);
        choiceConf(t) = 1-pChoices(t);
    end
    % update for meta-learning
    confBar(t+1) = confBar(t) + p.Results.aConf*(choiceConf(t)-confBar(t));
    if t > 1
        if rewards(t) == 0 && choiceConf(t) - (1-confBar(t)) > 0
            confSup(t) = p.Results.aFConf*confSup(t-1) + p.Results.aSup*(choiceConf(t));
        else
            confSup(t) = p.Results.aFConf*confSup(t-1);
        end
    end
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
behSim.pChoices = pChoices;
behSim.c = c;
behSim.rewardsP = rewardsP;
behSim.confSup = confSup;
behSim.confBar = confBar;
behSim.sigmaCurr = sigmaCurr;

if plotFlag
    figure; 
    subplot(2,1,1); hold on;
    plot(c(:,1),'color','r','linewidth',1.5);
    plot(c(:,2),'color','g','linewidth',1.5);
    plot(pChoices, 'color','k','linewidth',1.5)
    legend({'red', 'green'})
    title('c')
    subplot(2,1,2); hold on;
    tickLen = 1;
    for t = 1:length(choices)
        line([t t], [0 sign(choices(t)-0.5)*(0.5*tickLen + 0.5*rewards(t)*tickLen)]);
    end
    plot(confSup, 'color', [1 0.3 0.3], 'linewidth', 1.5);
    plot(rewardsP(:,2), 'color', [0.5 0.5 0.5], 'linewidth', 1.5)
    plot(sigmaCurr, 'color', [0.3 0.3 1], 'linewidth', 1.5);
    plot(1-confBar, 'color', [0.3 1 0.3], 'linewidth', 1.5);
    %plot(pePe, 'color', [0.3 1 0.3], 'linewidth', 1.5);
    title([num2str(extension)]);
    
end
end