function [succMin, succ] = vkf(startValues)

lambda = startValues(1);  % volatility update rate
v0 = startValues(2);      % initial volatility
omega = startValues(3);   % observation noise
beta = startValues(4);    % inverse-temp parameter for softmax decision function

numT = 1000;
numCorr = 0;
numRwds = 0;
totalT = numT * 10;
rSeed = rng;
rSeed = rSeed.Seed;

for r = 1:10
    rSeed = rSeed + 1;
    p = RestlessBanditDecoupled('RandomSeed',rSeed,'BlockLength', [20 35],'MaxTrials', numT,...
        'RewardProbabilities', [90 50 10]);
    
    w = zeros(numT+1, 2);
    v = [v0 v0; zeros(numT, 2)];
    m = zeros(numT+1, 2);

    for t = 1:numT
        % Select action
        mDiff = m(t, 2) -  m(t, 1);
        pRight = logistic(beta*mDiff);
        choiceProb(t,:) = [1-pRight pRight];

        if binornd(1, pRight) == 0 % left choice selected probabilistically
            p.inputChoice([1 0]);
            o = p.AllRewards(t, 1);

            k = (w(t, 1) + v(t, 1)) / (w(t, 1) + v(t, 1) + omega);
            alpha = sqrt(w(t, 1) + v(t, 1));
            m(t+1, 1) = m(t, 1) + alpha * (o - logistic(m(t, 1)));
            w(t+1, 1) = (1 - k) * (w(t, 1) + v(t, 1));
            wcov = (1 - k) * w(t, 1);
            v(t+1, 1) = v(t, 1) + lambda * ( (m(t+1, 1) - m(t, 1))^2 + w(t, 1) + w(t+1, 1) - 2 * wcov - v(t, 1) );

            m(t+1, 2) = m(t, 2); 
            w(t+1, 2) = w(t, 2); 
            v(t+1, 2) = v(t, 2); 

            if p.RewardProbabilities(1) >= p.RewardProbabilities(2)
                numCorr = numCorr + 1;
            end

        else
            p.inputChoice([0 1]);
            o = p.AllRewards(t, 2);

            k = (w(t, 2) + v(t, 2)) / (w(t, 2) + v(t, 2) + omega);
            alpha = sqrt(w(t, 2) + v(t, 2));
            m(t+1, 2) = m(t, 2) + alpha * (o - logistic(m(t, 2)));
            w(t+1, 2) = (1 - k) * (w(t, 2) + v(t, 2));
            wcov = (1 - k)* w(t, 2);
            v(t+1, 2) = v(t, 2) + lambda * ( (m(t+1, 2) - m(t, 2))^2 + w(t, 2) + w(t+1, 2) - 2 * wcov - v(t, 2) );

            m(t+1, 1) = m(t, 1); 
            w(t+1, 1) = w(t, 1); 
            v(t+1, 1) = v(t, 1);

            if p.RewardProbabilities(2) >= p.RewardProbabilities(1)
                numCorr = numCorr + 1;
            end
        end

    end
    numRwds = numRwds + sum(sum(p.AllRewards));
end


succ = numCorr/totalT + numRwds/totalT;
succMin = 2 - succ;


% figure; hold on
% smoothKern = ones(1,5);
% smoothKernSize = length(smoothKern);
% smoothKern = smoothKern/length(smoothKern);
% 
% allChoices = ones(1, numT);
% allChoices(p.AllChoices(:,1) == 1) = -1;
% allRewards = zeros(1, numT);
% allRewards(p.AllRewards(:,1) == 1) = -1;
% allRewards(p.AllRewards(:,2) == 1) = 1;
% 
% %smooth choices and outcomes
% smoothChoice = conv(allChoices, smoothKern);
% smoothRewards = conv(allRewards, smoothKern);
% 
% %plot smoothed curves
% plot(smoothChoice,'k','linewidth',2);
% plot(smoothRewards,'-','Color',[30 144 255]./255,'linewidth',2)
% xlabel('Trials')
% ylabel('<-- Left       Right -->')
% xlim([1 numT])
% 
% blockSwitch = unique(sort([p.BlockSwitchL p.BlockSwitchR]));
% blockSwitch = blockSwitch(blockSwitch < numT);
% 
% for i = 1:length(blockSwitch)
%     bs_loc = blockSwitch(i);
%     plot([bs_loc bs_loc],[-1 1.5],'--k','linewidth',1)
%     if rem(i,2) == 0
%         labelOffset = 1.42;
%     else
%         labelOffset = 1.34;
%     end
%     a = num2str(p.BlockProbs(i,1));
%     b = '/';
%     c = num2str(p.BlockProbs(i,2));
%     label = strcat(a,b,c);
%     text(bs_loc,labelOffset,label);
%     set(text,'FontSize',3);
% end
% 
% text(0,1.12,'L/R');
% 
% rMag = 0.3;
% nrMag = rMag/2;
% 
% % trial plot
% for i = 1:numT
%     if allChoices(i) == 1
%         if allRewards(i) == 1 % R side rewarded
%             plot([i i],[1 1+rMag],'k')
%         else
%             plot([i i],[1 1+nrMag],'color', [0.4 0.4 0.4]) % R side not rewarded
%         end
%     else
%         if allRewards(i) == -1 % L side rewarded
%             plot([i i],[-1*rMag - 1 -1],'k')
%         else
%             plot([i i],[-1*nrMag - 1 -1], 'color', [0.4 0.4 0.4])
%         end
%     end
% end
% 
% set(gca,'tickdir', 'out')
% set(gcf, 'renderer', 'painters', 'position', [-1919          41        1920         963])
% 
% 
