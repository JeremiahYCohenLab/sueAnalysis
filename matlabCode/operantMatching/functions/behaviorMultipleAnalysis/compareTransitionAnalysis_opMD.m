function compareTransitionAnalysis_opMD(xlFile, sheet, pre, post, varargin)


p = inputParser;
% default parameters if none given
p.addParameter('rwdProbs', [90 50 10]);
p.addParameter('tranWin', 5);
p.parse(varargin{:});

range = 20;


[wslsPre, transLowPre, transHighPre] = transitionAnalysis_opMD(xlFile, sheet, pre, p.Results.rwdProbs, p.Results.tranWin);
[wslsPost, transLowPost, transHighPost] = transitionAnalysis_opMD(xlFile, sheet, post, p.Results.rwdProbs, p.Results.tranWin);


lowAvgPre = mean(transLowPre,1);
highAvgPre = mean(transHighPre,1);
lowAvgPost = mean(transLowPost,1);
highAvgPost = mean(transHighPost,1);

colors = cool(4);

figure; 
x = [-range+1:range];
subplot(4,2,1); hold on
plot(x, lowAvgPre, '-', 'Color', colors(1,:), 'linewidth', 1.3)
plot(x, highAvgPre, '-', 'Color', colors(3,:), 'linewidth', 1.3)
plot([-range range], [0 0], ':k');
ylim([-1 1])
linetype = 'k';
vline(0, linetype)
set(gca, 'tickdir', 'out')
legend('medium -> low', 'high -> low')
title('pre')

subplot(4,2,2); hold on
plot(x, lowAvgPost, '-', 'Color', colors(2,:), 'linewidth', 1.3)
plot(x, highAvgPost, '-', 'Color', colors(4,:), 'linewidth', 1.3)
plot([-range range], [0 0], ':k');
ylim([-1 1])
vline(0, linetype)
set(gca, 'tickdir', 'out')
legend('medium -> low', 'high -> low')
title('post')

subplot(4,2,3); hold on
plot(x,[highAvgPre - lowAvgPre], '-k', 'linewidth', 2)
plot([-range range], [0 0], ':k');
ylabel('Choice average difference')
legend('high - medium')
set(gca, 'tickdir', 'out')
y = get(gca, 'ylim');
plot([0 0], y, ':k')

subplot(4,2,4); hold on
plot(x,[highAvgPost - lowAvgPost], '-k', 'linewidth', 2)
plot([-range range], [0 0], ':k');
ylabel('Choice average difference')
set(gca, 'tickdir', 'out')
y = get(gca, 'ylim');
plot([0 0], y, ':k')


yws = [0.7 1]; yls = [0 0.5];

subplot(4,2,5); hold on
errorbar(x, wslsPre.wS_low, wslsPre.sem_wS_low, 'Color', colors(1,:), 'linewidth', 1.3)
errorbar(x, wslsPre.wS_high, wslsPre.sem_wS_high, 'Color', colors(3,:), 'linewidth', 1.3)
ylabel('Win-stay')
ylim([yws])
plot([0 0], yws, '--k')
set(gca, 'tickdir', 'out')

subplot(4,2,6); hold on
errorbar(x, wslsPost.wS_low, wslsPost.sem_wS_low, 'Color', colors(1,:), 'linewidth', 1.3)
errorbar(x, wslsPost.wS_high, wslsPost.sem_wS_high, 'Color', colors(3,:), 'linewidth', 1.3)
ylabel('Win-stay')
ylim([yws])
plot([0 0], yws, '--k')
set(gca, 'tickdir', 'out')


subplot(4,2,7); hold on
errorbar(x, wslsPre.lS_low, wslsPre.sem_lS_low, 'Color', colors(1,:), 'linewidth', 1.3)
errorbar(x, wslsPre.lS_high, wslsPre.sem_lS_low, 'Color', colors(3,:), 'linewidth', 1.3)
ylabel('Lose-shift')
xlabel('Trials from switch')
ylim([yls])
plot([0 0], yls, '--k')
set(gca, 'tickdir', 'out')


subplot(4,2,8); hold on
errorbar(x, wslsPost.lS_low, wslsPost.sem_lS_low, 'Color', colors(1,:), 'linewidth', 1.3)
errorbar(x, wslsPost.lS_high, wslsPost.sem_lS_low, 'Color', colors(3,:), 'linewidth', 1.3)
ylabel('Lose-shift')
xlabel('Trials from switch')
ylim([yls])
plot([0 0], yls, '--k')
set(gca, 'tickdir', 'out')

suptitle(['Choice at block transitions: ' sheet]); 
set(gcf, 'renderer', 'painters', 'position', [-1365 18 733 978])