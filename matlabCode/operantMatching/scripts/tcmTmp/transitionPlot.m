figure; 
subplot(2,2,1); hold on
title('Acutal behavior')
x = [-range+1:range];
plot(x,lowAvg, '-r', 'linewidth', 1.5)
plot(x,highAvg, '-b', 'linewidth', 1.5)

xx = [1:range+1];
mdlFitLow = singleExpFitInt(xx,lowAvg(range:end));
mdlFitHigh = singleExpFitInt(xx,highAvg(range:end));
expConvLow = mdlFitLow.a*exp(-(mdlFitLow.b)*(1:range+1)) + mdlFitLow.c;
expConvHigh = mdlFitHigh.a*exp(-(mdlFitHigh.b)*(1:range+1)) + mdlFitHigh.c;
plot([0:range], expConvLow, '--r', 'linewidth', 1.5)
plot([0:range], expConvHigh,'--b', 'linewidth', 1.5)

linetype = 'k';
vline(0, linetype)
ylim([-1 1])
set(gca, 'tickdir', 'out')
ylabel('Choice average')
legend('Medium -> Low', 'High -> Low')


subplot(2,2,3); hold on
title('Simulated behavior')
plot(x,lowAvgSim, '-r', 'linewidth', 1.5)
plot(x,highAvgSim, '-b', 'linewidth', 1.5)

mdlFitLowSim = singleExpFitInt(xx,lowAvgSim(range:end));
mdlFitHighSim = singleExpFitInt(xx,highAvgSim(range:end));
expConvLowSim = mdlFitLowSim.a*exp(-(mdlFitLowSim.b)*(1:range+1)) + mdlFitLowSim.c;
expConvHighSim = mdlFitHighSim.a*exp(-(mdlFitHighSim.b)*(1:range+1)) + mdlFitHighSim.c;
plot([0:range], expConvLowSim, '--r', 'linewidth', 1.5)
plot([0:range], expConvHighSim,'--b', 'linewidth', 1.5)

linetype = 'k';
vline(0, linetype)
ylim([-1 1])
set(gca, 'tickdir', 'out')
ylabel('Choice average')
xlabel('Trials from switch')

subplot(2,2,2); hold on
taus = [1/mdlFitLow.b 1/mdlFitHigh.b 1/mdlFitLowSim.b 1/mdlFitHighSim.b];
scatter([1:4], taus, 'k', 'filled')
xticks([1:4])
xticklabels({'Med->Low', 'High->Low', 'Sim: Med->Low', 'Sim: High->Low'})
ylabel('\tau')
xlim([0 5])
set(gca, 'tickdir', 'out')

subplot(2,2,4); hold on
plot(x,[highAvg - lowAvg], '-k', 'linewidth', 1.5)
plot(x,[highAvgSim - lowAvgSim], '-', 'Color', [0.7 0.7 0.7], 'linewidth', 1.5)
plot(x, zeros(1,length(x)), ':k')
vline(0, linetype)
ylim([-.2 .4])

ylabel('Choice average difference')
xlabel('Trials from switch')
legend('Actual Behavior', 'Simulated')
set(gca, 'tickdir', 'out')
set(gcf, 'renderer', 'painters')


%% plot two models

figure; 
subplot(2,2,3); hold on
x = [-range+1:range];
xx = [1:range+1];
title('Volatility --> NPE learning rate')
plot(x,lowAvgSim, '-r', 'linewidth', 1.5)
plot(x,highAvgSim, '-b', 'linewidth', 1.5)

mdlFitLowSim = singleExpFitInt(xx,lowAvgSim(range:end));
mdlFitHighSim = singleExpFitInt(xx,highAvgSim(range:end));
expConvLowSim = mdlFitLowSim.a*exp(-(mdlFitLowSim.b)*(1:range+1)) + mdlFitLowSim.c;
expConvHighSim = mdlFitHighSim.a*exp(-(mdlFitHighSim.b)*(1:range+1)) + mdlFitHighSim.c;
plot([0:range], expConvLowSim, '--r', 'linewidth', 1.5)
plot([0:range], expConvHighSim,'--b', 'linewidth', 1.5)

linetype = 'k';
vline(0, linetype)
ylim([-1 1])
set(gca, 'tickdir', 'out')
ylabel('Choice average')
xlabel('Trials from switch')

subplot(2,2,4); hold on
plot(x,[highAvg - lowAvg], '-k', 'linewidth', 1.5)
plot(x,[highAvgSim - lowAvgSim], '-', 'Color', [0.7 0.7 0.7], 'linewidth', 1.5)
plot(x, zeros(1,length(x)), ':k')
vline(0, linetype)
ylim([-0.5 1])

ylabel('Choice average difference')
xlabel('Trials from switch')
legend('Actual Behavior', 'Simulated')
set(gca, 'tickdir', 'out')
set(gcf, 'renderer', 'painters')