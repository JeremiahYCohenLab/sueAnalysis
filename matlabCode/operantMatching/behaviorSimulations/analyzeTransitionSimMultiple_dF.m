function analyzeTransitionSimMultiple_dF

for i = 1:8
[zeroCross(i), zeroCross_actual(i), threshCrossLow(i), threshCrossHigh(i), threshCrossLow_actual(i), threshCrossHigh_actual(i)] = ...
    analyzeTransitionSim_opMD('goodBehDays.xlsx', ['DR0' num2str(i)], 'clean',  'clean', 'tranWin', 5, 'modelName', 'sevenParam_absPePeAN_scale_int_bias', 'runs', 500);
end

colors = cool(8);
figure; 
subplot(1,3,1); hold on;
scatter(zeroCross, zeroCross_actual, ones(1,8)*150, colors, 'filled')
xl = xlim; yl = ylim;
plot([0 max([xl yl])], [0 max([xl yl])], '--k') 
ylim([0 max([xl yl])]);
xlim([0 max([xl yl])]);
xlabel('zero cross sim')
ylabel('zero cross actual')
set(gca, 'tickdir', 'out')

subplot(1,3,2); hold on;
scatter(threshCrossLow, threshCrossLow_actual, ones(1,8)*150, colors, 'filled')
xl = xlim; yl = ylim;
plot([0 max([xl yl])], [0 max([xl yl])], '--k') 
ylim([0 max([xl yl])]);
xlim([0 max([xl yl])]);
xlabel('medium thresh cross sim')
ylabel('medium thresh cross actual')
set(gca, 'tickdir', 'out')

subplot(1,3,3); hold on;
scatter(threshCrossHigh, threshCrossHigh_actual, ones(1,8)*150, colors, 'filled')
xl = xlim; yl = ylim;
plot([0 max([xl yl])], [0 max([xl yl])], '--k') 
ylim([0 max([xl yl])]);
xlim([0 max([xl yl])]);
xlabel('high thresh cross sim')
ylabel('high thresh cross actual')
set(gca, 'tickdir', 'out')

