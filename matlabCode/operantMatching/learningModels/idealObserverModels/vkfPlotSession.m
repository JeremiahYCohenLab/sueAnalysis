function vkfPlotSession(sessionName)

[behSessionData, unCorrectedBlockSwitch, out] = loadBehavioralData([sessionName '.asc'], 0);
s = parseBehavioralData(behSessionData, unCorrectedBlockSwitch);
choices = s.allChoices;
outcomes = abs(s.allRewards);

lambda = 0.549245864306393;  % volatility update rate
v0 = 0.422311804752216;      % initial volatility
omega = 0.975697922267773;   % observation noise
beta = 1.341230225800368;    % inverse-temp parameter for softmax decision function

numT = length(choices);
w = zeros(numT+1, 2);
v = [v0 v0; zeros(numT, 2)];
m = zeros(numT+1, 2);

for t = 1:numT
    % Select action
    mDiff = m(t, 2) -  m(t, 1);
    pRight = logistic(beta*mDiff);
    choiceProb(t,:) = [1-pRight pRight];

    if choices(t) == -1 % left choice selected probabilistically
        o = outcomes(t);
        
        k = (w(t, 1) + v(t, 1)) / (w(t, 1) + v(t, 1) + omega);
        alpha = sqrt(w(t, 1) + v(t, 1));
        m(t+1, 1) = m(t, 1) + alpha * (o - logistic(m(t, 1)));
        w(t+1, 1) = (1 - k) * (w(t, 1) + v(t, 1));
        wcov = (1 - k) * w(t, 1);
        v(t+1, 1) = v(t, 1) + lambda * ( (m(t+1, 1) - m(t, 1))^2 + w(t, 1) + w(t+1, 1) - 2 * wcov - v(t, 1) );

        m(t+1, 2) = m(t, 2); 
        w(t+1, 2) = w(t, 2); 
        v(t+1, 2) = v(t, 2); 

    else
        o = outcomes(t);

        k = (w(t, 2) + v(t, 2)) / (w(t, 2) + v(t, 2) + omega);
        alpha = sqrt(w(t, 2) + v(t, 2));
        m(t+1, 2) = m(t, 2) + alpha * (o - logistic(m(t, 2)));
        w(t+1, 2) = (1 - k) * (w(t, 2) + v(t, 2));
        wcov = (1 - k)* w(t, 2);
        v(t+1, 2) = v(t, 2) + lambda * ( (m(t+1, 2) - m(t, 2))^2 + w(t, 2) + w(t+1, 2) - 2 * wcov - v(t, 2) );

        m(t+1, 1) = m(t, 1); 
        w(t+1, 1) = w(t, 1); 
        v(t+1, 1) = v(t, 1);
    end

end

figure; hold on

plot(v(:,1), '-c', 'linewidth', 1.5)
plot(v(:,2), '-m', 'linewidth', 1.5)
legend('volatility L', 'volatility R')
yU = max(max(v));
yL = min(min(v));

blockProbs = [behSessionData(unCorrectedBlockSwitch).rewardProbL; behSessionData(unCorrectedBlockSwitch).rewardProbR]';

for i = 1:length(s.blockSwitch)
    bs_loc = s.blockSwitch(i);
    plot([bs_loc bs_loc],[yL yU+0.5],'--k','linewidth',1)
    if rem(i,2) == 0
        labelOffset = yU + 0.42;
    else
        labelOffset = yU + 0.34;
    end
    b = num2str(blockProbs(i,1));
    c = '/';
    d = num2str(blockProbs(i,2));
    label = strcat(b,c,d);
    text(bs_loc,labelOffset,label);
    set(text,'FontSize',3);
end

text(0,yU+0.12,'L/R');

rMag = 0.3;
nrMag = rMag/2;

% trial plot
j = 1;
for i = 1:length(choices)
    if choices(i) == 1
        if outcomes(i) == 1 % R side rewarded
            plot([i i],[yU yU+rMag],'k')
        else
            plot([i i],[yU yU+nrMag],'color', [0.4 0.4 0.4]) % R side not rewarded
        end
    else
        if outcomes(i) == 1 % L side rewarded
            plot([i i],[-1*rMag + yL yL],'k')
        else
            plot([i i],[-1*nrMag + yL yL], 'color', [0.4 0.4 0.4])
        end
    end
end

xlim([0 numT])
ylim([-1*rMag + yL yU+0.5])
set(gca,'tickdir', 'out')
set(gcf, 'renderer', 'painters', 'position', [-1697 322 1461 531])
suptitle(sessionName)