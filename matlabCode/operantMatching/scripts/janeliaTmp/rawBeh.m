% trial plot
figure; 
subplot(2,1,2); hold on

j = 1;
for i = 1:length(allChoices)
    if allChoices(i) == 1
        if allRewards(i) == 1 % R side rewarded
            plot([i i],[1.05 1.35],'b', 'linewidth', 1.4)
        else
            plot([i i],[1.05 1.2],'r', 'linewidth', 1.4) % R side not rewarded
        end
    else 
        if allRewards(i) == -1 % L side rewarded
            plot([i i],[-1.35 -1.05],'b', 'linewidth', 1.4)
        else
            plot([i i],[-1.2 -1.05],'r', 'linewidth', 1.4)
        end
    end
end
for i = 1:length(blockSwitch)
    bs_loc = blockSwitch(i) - 0.5;
    plot([bs_loc bs_loc],[-1 1],'--','linewidth',1.1, 'Color', [0.3 0.3 0.3])
    if rem(i,2) == 0
        labelOffset = 1.45;
    else
        labelOffset = 1.55;
    end
    a = num2str(allProbsL(blockSwitch(i)+1));
    b = '/';
    c = num2str(allProbsR(blockSwitch(i)+1));
    label = strcat(a,b,c);
    text(bs_loc,labelOffset,label);
    set(text,'FontSize',3);
end
plot(conv(allChoices,normKern)/max(conv(allChoices,normKern)),'k','linewidth',2);
plot(conv(allRewards,normKern)/max(conv(allRewards,normKern)),'--k','linewidth',2)
ylabel('Probability ofchoosing right spout')
xlim([0 462])
ylim([-1.35 1.35])
yticks([-1 0 1]); yticklabels([0 0.5 1]);
ax = gca;
ax.TickDir = 'out';

subplot(2,1,1); hold on; 
plot(allProbsL/100,'--','linewidth',1.4, 'color', [0.7 0 1])
plot(allProbsR/100,'--','linewidth',1.4, 'color', 'm')
xlim([0 462])
ylim([0 0.8])
ylabel('Probability of reward delivery')
yticks([0 0.1 0.4 0.7])
ax = gca;
ax.TickDir = 'out';
set(gcf, 'Position', [-1910 159 1905 750]);