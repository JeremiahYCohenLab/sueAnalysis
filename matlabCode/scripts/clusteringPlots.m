
color1 = [0, 0.8, 0.8];
color2 = [1, 0.3, 1];
colors = [color1; color2];
%%
load F:\tmpData\catWithOutcome.mat
ind = cats;
figure2; hold on;
scatter(score(:,1), score(:,2), 25, [0.7 0.7 0.7], 'filled');
for i = 1:2
    scatter(score(ind==i,1), score(ind==i,2), 25, colors(i,:), 'filled');
end
legend({'NA', 'neg', 'pos'}, 'FontSize', 12)
xlabel('PC1');
ylabel('PC2');
xlim([-17 20])
ylim([-20 12])
set(gca, 'YTick', [-10, 0, 10]);
set(gca, 'XTick', [-10, 0, 10, 20]);
set(gca, 'TickDir', 'out');
title('outcome', 'FontSize', 15)
%%
load F:\tmpData\catWithOutcome.mat
ind = cats;
figure2; hold on;
scatter(score(:,1), score(:,2), 25, [0.7 0.7 0.7], 'filled');
for i = 1:2
    scatter(score(ind==i,1), score(ind==i,2), 25, colors(i,:), 'filled');
end
legend({'NA', 'neg', 'pos'}, 'FontSize', 12)
xlabel('PC1');
ylabel('PC2');
xlim([-17 20])
ylim([-20 12])
set(gca, 'YTick', [-10, 0, 10]);
set(gca, 'XTick', [-10, 0, 10, 20]);
set(gca, 'TickDir', 'out');
title('outcome', 'FontSize', 15)
%% 
load F:\tmpData\catWithQchosen.mat

figure2; hold on;
scatter(score(:,1), score(:,2), 25, [0.7 0.7 0.7], 'filled');
for i = 1:2
    scatter(score(ind==i,1), score(ind==i,2), 25, colors(i,:), 'filled');
end
legend({'pos', 'NA', 'neg'}, 'FontSize', 12)
xlabel('PC1');
ylabel('PC2');
xlim([-17 20])
ylim([-20 12])
set(gca, 'YTick', [-10, 0, 10]);
set(gca, 'XTick', [-10, 0, 10, 20]);
set(gca, 'TickDir', 'out');
title('expected value', 'FontSize', 15)
%%
load F:\tmpData\catWithWFFeatures.mat

figure2; hold on;
for i = 1:2
    scatter(score(ind==i,1), score(ind==i,2), 25, colors(i,:), 'filled');
end
legend({'TypeI', 'TypeII'}, 'FontSize', 12)
xlabel('PC1');
ylabel('PC2');
xlim([-17 20])
ylim([-20 12])
set(gca, 'YTick', [-10, 0, 10]);
set(gca, 'XTick', [-10, 0, 10, 20]);
set(gca, 'TickDir', 'out');
title('waveform clustering', 'FontSize', 15)
%% waveform
figure2; hold on;
for i = 1:2
    mat = waveformsSession(ind==i,:);
    plotFilled(1:size(mat,2), mat, colors(i,:));
end
plot([10 26], [-0.4 -0.4], 'k', 'LineWidth', 2)
text(12, -0.5, '0.5 ms', 'FontSize', 15)
set(gca, 'XColor', 'none');
set(gca, 'YColor', 'none');
%% evoke
load F:\tmpData\excit.mat
figure2; hold on;
scatter(excit(:,1), excit(:,2), 18, [0.7 0.7 0.7], 'filled');
for i = 1:2
    scatter(excit(ind==i,1), excit(ind==i,2), 25, colors(i,:), 'filled');
end
plot([0 10], [0 10], 'LineStyle', '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 2); 
xlabel('baseline')
ylabel('evoke')
set(gca, 'YTick', [ 0, 10, 20]);
set(gca, 'XTick', [0, 5, 10]);
set(gca, 'TickDir', 'out');

figure2; hold on;
edges = linspace(min((excit(:,2)-excit(:,1))./excit(:,1))-0.00001, max((excit(:,2)-excit(:,1))./excit(:,1))+0.00001, 15);
for i = 1:2
    histogram((excit(ind==i,2)-excit(ind==i,1))./excit(ind==i,1), edges, 'FaceColor', colors(i,:), 'FaceAlpha', 0.4, 'Normalization', 'probability');
end
set(gca, 'XTick', [ 0, 2, 4]);
set(gca, 'YTick', [0, 0.2]);
set(gca, 'TickDir', 'out');
plot([0 0], [0 0.3], 'LineStyle', '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 2); 

xlabel('evoke ratio')
%%
