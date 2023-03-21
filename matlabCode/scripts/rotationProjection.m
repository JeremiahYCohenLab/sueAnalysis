% Tstats projection
clear all;
load('F:\tmpData\catWithOutcomePe.mat');
outcomeInd = cell2mat(cellfun(@(x) strcmp(x, 'outcome'), regressors, 'UniformOutput', false));
qInd = cell2mat(cellfun(@(x) strcmp(x, 'Qchosen'), regressors, 'UniformOutput', false));
%% clustering based on waveform
color1 = [0 0.8 0.8];
color2 = [1 0.2 1];
[coeff,scores,latent, ~, explained, mu] = pca(zscore(waveformsSession, [], 1));
indAll = {};
dis = {};
for a = 1:10
    [indAll{a}, ~, dis{a}] = kmeans([scores(:, 1:5)], 2);
end
[~,optiInds] = min(cellfun(@mean, dis));
cats = indAll{optiInds};
if sum(cats==1)<sum(cats==2)
    cats = 3-cats;
end

figure2;
subplot(1,2,1)
scatter(scores(cats==1, 1), scores(cats==1, 2), 12, color1, 'filled')
hold on;
scatter(scores(cats==2, 1), scores(cats==2, 2), 12, color2, 'filled')

subplot(1,2,2); hold on
plotFilled(1:size(waveformsSession, 2), waveformsSession(cats==1,:), color1);
plotFilled(1:size(waveformsSession, 2), waveformsSession(cats==2,:), color2);
%% scatter in feature space
ind = allSigs(:,outcomeInd)<=0.05;
ind = ind.*sign(allTstats(:,outcomeInd));
indQ = allSigs(:,qInd)<=0.05;
indQ = indQ.*sign(allTstats(:,qInd));
figure2; hold on;
scatter(scores(:,1), excitRatio, 50, [0.6 0.6 0.6], 'filled')
scatter(scores(ind == -1,1), excitRatio(ind == -1), 50, color1, 'filled')
scatter(scores(ind == 1,1), excitRatio(ind == 1), 50, color2, 'filled')

%% scatter in tstats space
figure2; hold on;
% scatter(allTstats(:, outcomeInd), allTstats(:,qInd), 30, [0.7 0.7 0.7], 'filled');
scatter(allTstats(cats==1, outcomeInd), allTstats(cats==1,qInd), 50, [0 0 0], 'o', 'LineWidth', 2.5);
scatter(allTstats(cats==2, outcomeInd), allTstats(cats==2,qInd), 50, [0.6 0.6 0.6],  'o', 'LineWidth', 2.5);
scatter(allTstats(ind==-1, outcomeInd), allTstats(ind==-1,qInd), 5, color1, 'o', 'LineWidth', 2, 'MarkerEdgeAlpha', 1);
scatter(allTstats(ind==1, outcomeInd), allTstats(ind==1,qInd), 5, color2,  'o', 'LineWidth', 2, 'MarkerEdgeAlpha', 1);
% scatter(allTstats(allSigs(:,qInd)>0.05, outcomeInd), allTstats(allSigs(:,qInd)>0.05,qInd), 60, 'r',  'o', 'LineWidth', 1, 'MarkerEdgeAlpha', 0.7);
plot([-15 15], [0 0], 'LineStyle', '--', 'LineWidth', 2, 'Color', [0.6 0.6 0.6]);
plot([0 0], [-10 5], 'LineStyle', '--', 'LineWidth', 2, 'Color', [0.6 0.6 0.6]);
%% scatter in coeff space
figure2; hold on;
% scatter(allTstats(:, outcomeInd), allTstats(:,qInd), 30, [0.7 0.7 0.7], 'filled');
scatter(allCoeffs(cats==1, outcomeInd), allCoeffs(cats==1,qInd), 50, [0 0 0], 'x', 'LineWidth', 2);
scatter(allCoeffs(cats==2, outcomeInd), allCoeffs(cats==2,qInd), 50, [0.6 0.6 0.6],  'x', 'LineWidth', 2);
scatter(allCoeffs(ind==-1, outcomeInd), allCoeffs(ind==-1,qInd), 50, color1, '.', 'LineWidth', 2);
scatter(allCoeffs(ind==1, outcomeInd), allCoeffs(ind==1,qInd), 50, color2,  '.', 'LineWidth', 2);
xlim([-1.5 1.5])
ylim([-1.5 1.5])
%% polar histogram
edges = linspace(-pi, pi, 30);
allVec = [allTstats(:, outcomeInd), allTstats(:,qInd)];
[theta1, rho] = cart2pol(allVec((ind~=0|indQ~=0)&cats==1,1), allVec((ind~=0|indQ~=0)&cats==1,2));
[theta2, rho] = cart2pol(allVec((ind~=0|indQ~=0)&cats==2,1), allVec((ind~=0|indQ~=0)&cats==2,2));
figure2; 
polarhistogram(theta1,edges, 'FaceColor', [0.1 0.1 0.1], 'FaceAlpha',.7, 'Normalization', 'Probability', 'EdgeColor', 'none');
hold on; 
polarhistogram(theta2,edges, 'FaceColor', [0.7 0.7 0.7], 'FaceAlpha',.6, 'Normalization', 'Probability', 'EdgeColor', 'none');
%% color gradient on outcome
% calculate boundary tstats
colors = zeros(length(allSessions),3);
colorsp = zeros(length(allSessions),3);
for i = 1:length(allSessions)
    s = behAnalysisNoPlot_opMD(allSessions{i}, 'simpleFlag', 1);
    tBound = zeros(1,2);
    tBound(1) = tinv(0.025, length(s.allChoices)-length(regressors));
    tBound(2) = tinv(0.975, length(s.allChoices)-length(regressors));
    p = tcdf(allTstats(i,qInd), length(s.allChoices)-length(regressors));
    
    if ind(i) == 1
        colors(i,:) = color2;
    else
        if ind(i) == -1
            colors(i,:) = color1;
        else
            ratioTmp = (abs(allTstats(i,outcomeInd)/tBound(1)))^0.5;
            if allTstats(i,outcomeInd)> 0 
                colors(i,:) = ratioTmp*color2 + (1-ratioTmp)*[1 1 1];       
            else
                colors(i,:) = ratioTmp*color1 + (1-ratioTmp)*[1 1 1];   
            end
        end
    end
    
    if p>0.5
        colorsp(i,:) = ((p-0.5)/0.5)*color2 + ((1-p)/0.5)*[1 1 1];
    else
        colorsp(i,:) = ((0.5-p)/0.5)*color1 + (p/0.5)*[1 1 1];
    end
end
%%
figure; hold on;
scatter(allTstats(:, outcomeInd), allTstats(:,qInd), 50, colorsp, 'o', 'filled');
scatter(allTstats(:, outcomeInd), allTstats(:,qInd), 60, [0.5 0.5 0.5], 'o');
xlim([-15 15])
ylim([-10 5])
%%
figure2; hold on;
scatter(scores(:,1), excitRatio, 50, colorsp, 'o', 'filled');
scatter(scores(:,1), excitRatio, 60, [0.5 0.5 0.5], 'o');
% xlim([-15 15])
% ylim([-10 5])
bins = 100;
myMap = [linspace(0, 1, bins)' linspace(0.8, 1, bins)' linspace(0.8, 1, bins)'; ...
         linspace(1, 1, bins)' linspace(1, 0.2, bins)' linspace(1, 1, bins)'];
colorbar;
colormap(myMap);
a = colorbar;
ylabel(a, 'alpha', 'FontSize', 18)
colorbar('Ticks', 0:0.5:1, 'FontSize', 12);
%% color gradient on outcome and Qchosen
% calculate boundary tstats
numBins = 1000;
colorsBi = zeros(length(allSessions),3);
colorRainbow = hsv(2*numBins - 1);
angleMap = [linspace(0, pi, numBins) linspace(-pi, 0, numBins)];
angleMap = [angleMap(1:numBins), angleMap(numBins+2:end)];
theta = zeros(size(allSessions));
for i = 1:length(allSessions)
    s = behAnalysisNoPlot_opMD(allSessions{i}, 'simpleFlag', 1);
    pO = tcdf(allTstats(i, outcomeInd), length(s.allChoices)-length(regressors));
    pQ = tcdf(allTstats(i, qInd), length(s.allChoices)-length(regressors));
    
    pO = 2*(pO-0.5);
    pQ = 2*(pQ-0.5);
    [theta(i), rho] = cart2pol(pO, pQ);

    rhoNew = max(abs([pO pQ]));
    
    [~,colorLoc] = min(abs(angleMap-theta(i)));
    colorFull = colorRainbow(colorLoc,:);
    colorsBi(i,:) = rhoNew*colorFull + (1-rhoNew)*[1 1 1];
end
%%
figure; hold on;
scatter(allTstats(:, outcomeInd), allTstats(:,qInd), 50, colorsBi, 'o', 'filled');
scatter(allTstats(:, outcomeInd), allTstats(:,qInd), 60, [0.5 0.5 0.5], 'o');
xlim([-15 15])  
ylim([-10 5])
%%
figure2; hold on;
scatter(scores(:,1), excitRatio, 50, colorBi, 'o', 'filled');
scatter(scores(:,1), excitRatio, 60, [0.5 0.5 0.5], 'o');
% xlim([-15 15])
% ylim([-10 5])
bins = 100;
myMap = [linspace(0, 1, bins)' linspace(0.8, 1, bins)' linspace(0.8, 1, bins)'; ...
         linspace(1, 1, bins)' linspace(1, 0.2, bins)' linspace(1, 1, bins)'];
colorbar;
colormap(myMap);
a = colorbar;
ylabel(a, 'alpha', 'FontSize', 18)
colorbar('Ticks', 0:0.5:1, 'FontSize', 12);
%% histogram
edges = linspace(min(allTstats(:,outcomeInd))-0.01, max(allTstats(:,outcomeInd))+0.01, 20);
figure2;
hold on;
histogram(allTstats(cats==1, outcomeInd), edges, 'FaceColor', [0 0 0], 'EdgeColor', 'none', 'EdgeAlpha', 0.8, 'Normalization', 'probability');
histogram(allTstats(cats==2, outcomeInd), edges, 'FaceColor', [0.7 0.7 0.7], 'EdgeColor', 'none', 'EdgeAlpha', 0.7, 'Normalization', 'probability');
title('p<0.001', 'FontSize', 14)
%%
edges = linspace(min(allTstats(:,qInd))-0.01, max(allTstats(:,qInd))+0.01, 20);
figure2;
hold on;
histogram(allTstats(cats==1, qInd), edges, 'FaceColor', [0 0 0], 'EdgeColor', 'none', 'EdgeAlpha', 0.8, 'Normalization', 'probability');
histogram(allTstats(cats==2, qInd), edges, 'FaceColor', [0.7 0.7 0.7], 'EdgeColor', 'none', 'EdgeAlpha', 0.7, 'Normalization', 'probability');
title('p<0.001', 'FontSize', 14)
%%

