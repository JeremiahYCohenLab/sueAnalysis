clear all;
load('F:\tmpData\catWithOutcomePe.mat');
catsPe = zeros(size(catsOutcome))';
catsPe(allSigsPe(:,1)<0.05&allTstatsPe(:,1)<0, :) = 1;
catsPe(allSigsPe(:,1)<0.05&allTstatsPe(:,1)>0, :) = 2;
%%
figure2; hold on
scatter(allTstats(:, 1), allTstats(:,3), 40, [0.7 0.7 0.7], 'filled');
scatter(allTstats(catsOutcome==1, 1), allTstats(catsOutcome==1,3), 40, [0 0.8 0.8], 'filled');
scatter(allTstats(catsOutcome==2, 1), allTstats(catsOutcome==2,3), 40, 'm', 'filled');
plot([0 0], [min(allTstats(:,3)) max(allTstats(:,3))], 'Color', [0.4 0.4 0.4], 'LineStyle', ':');
plot([min(allTstats(:,1)) max(allTstats(:,1))], [0 0], 'Color', [0.4 0.4 0.4], 'LineStyle', ':');
scatter(allTstats(catsOutcome==1, 1), allTstats(catsOutcome==1,3), 30, [0 0.8 0.8], 'filled');
scatter(allTstats(catsQ==0, 1), allTstats(catsQ==0,3), 12, [1 1 1], 'filled');
scatter(allTstats(catsPe==1,1), allTstats(catsPe==1,3), 50, [0 0.8 0.8],'d', 'LineWidth',2)
scatter(allTstats(catsPe==2,1), allTstats(catsPe==2,3), 50, 'm', 'd', 'LineWidth',2)
%%
set(gca, 'TickDir', 'Out')
set(gca, 'XTick', [-10:10:20], 'FontSize', 14)
set(gca, 'YTick', [-10:5:5], 'FontSize', 14)
set(gca, 'Box', 'off')
ylim([-10 5]); xlim([-15 20])
xlabel('Outcome t-stats', 'FontSize', 18)
ylabel('Expected value t-stats', 'FontSize', 18)
%%