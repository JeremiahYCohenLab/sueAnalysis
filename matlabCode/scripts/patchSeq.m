
%% preparation
% filtering
load('C:\Users\zhixi\Downloads\meanAPs.mat', 'mat')
meanSpikes = mat';
%%
% fs = 50000;
% fc2 = [300]; %cutoff frequency in Hz (anything higher than this f passes)
% [b2, a2] = butter(2,fc2/(fs/2),'low'); 
% for i = 1:size(meanSpikes, 2)
%     if sum(isnan(meanSpikes(:,i))) == 0
%         meanSpikes(:,i) = filtfilt(b2,a2,meanSpikes(:,i));
%     end
% end
% % meanSpikes = meanSpikes';
%%
baseline = mean(meanSpikes(:, 1:400), 2);
meanSpike_bl = meanSpikes - baseline;
peak = meanSpike_bl(:, 501);
meanSpike_nl = meanSpike_bl./peak;
%% first derivative
meanSpike_1D = diff(meanSpike_nl, 1, 2);
% filtering
for i = 1:size(meanSpike_1D, 2)
    if sum(isnan(meanSpike_1D(:,i))) == 0
        meanSpike_1D(:,i) = filtfilt(b2,a2,meanSpike_1D(:,i));
    end
end
%%
baseline_1D = mean(meanSpike_1D(:, 1:200), 2);
peak1D = zeros(size(meanSpike_1D,1), 1);
peakInd = zeros(size(meanSpike_1D,1), 1);
meanSpike_1D_realign = zeros(size(meanSpike_1D,1), 401);
meanSpike_1D_bl = meanSpike_1D - baseline_1D;
for i = 1:size(meanSpike_1D_bl,1)
    [peak1D(i), peakInd(i)] = max(meanSpike_1D_bl(i,:));
    meanSpike_1D_realign(i,:) = meanSpike_1D_bl(i, (peakInd(i)-120):(peakInd(i)+280));
end
meanSpike_1D_nl = meanSpike_1D_realign./peak1D;
%%
% zscore
meanSpike_1D_nlZS = NaN(size(meanSpike_1D_nl));
for i = 1:size(meanSpike_1D_nlZS, 2)
    tmp = meanSpike_1D_nl(:,i);
    tmp(~isnan(tmp)) = zscore(tmp(~isnan(tmp)));
    meanSpike_1D_nlZS(:,i) = tmp;
end
%%
meanSpike_nlZS = NaN(size(meanSpike_nl));
for i = 1:size(meanSpike_nlZS, 2)
    tmp = meanSpike_nl(:,i);
    tmp(~isnan(tmp)) = zscore(tmp(~isnan(tmp)));
    meanSpike_nlZS(:,i) = tmp;
end
%%
[coeff,score,latent,tsquared,explained,mu] = pca(meanSpike_nl(:, 500-120:500+280));
%%
allWFs = meanSpike_nl(:, 500-120:500+280);
Mdl = rica(allWFs, 5);
score = transform(Mdl, allWFs);
%%
figure;
scatter(score(:,1), score(:,2));
%% GMM
% initialization
scoreTmp = score(sum(isnan(score),2)==0,:);
[expectations, theta] = fitGaussMixture(scoreTmp(:, 1:5),2,'kmeans');
%% Kmeans
numCat = 2;
indAll = {};
dis = {};
for a = 1:10
    [indAll{a}, ~, dis{a}] = kmeans([score(:, 1:3)], numCat);
end
[~,optiInds] = min(cellfun(@mean, dis));
indRaw = indAll{optiInds};
%%
confThresh = 0.8;
[confG, indG] = max(expectations.z, [], 2);
indG(confG<confThresh) = NaN;
indRaw = NaN(1,size(score, 1));
indRaw(~isnan(sum(score, 2))) = indG;
%%
y = tsne(score(:,1:5));
colors = [0 0 0;...
          0.7 0.7 0.7];
figure;
gscatter(y(~isnan(indG),1), y(~isnan(indG),2), indRaw(~isnan(indG)));
figure;
gscatter(score(~isnan(indRaw),1), score(~isnan(indRaw),2), indRaw(~isnan(indRaw)), colors);
%%
figure;
time = [-120:280]/50000 * 1000;
plotFilled(time, meanSpike_nl(indRaw==1,500-120:500+280), 'k');
plotFilled(time, meanSpike_nl(indRaw==2,500-120:500+280), [0.7 0.7 0.7]);
%%
figure;
typeIInd = find(indRaw == 2);
for i = 1:length(typeIInd)
    subplot(5,8,i);
    plot(time, meanSpike_nl(typeIInd(i), 500-120:500+280),'k', 'LineWidth',2);
end
sgtitle('TypeII')
%%
figure;
pdInd = 1;
edges = linspace(min(score(:,pdInd)-0.001), max(score(:,pdInd)+0.001), 30);
hold on;
histogram(score(indRaw==1,pdInd),edges, 'FaceColor','k','Normalization','probability');
histogram(score(indRaw==2,pdInd),edges, 'FaceColor',[0.4, 0.4 0.4],'Normalization','probability');
%%
figure; 
subplot(1,2,1);
plot(time, meanSpike_nl(indRaw == 1, 500-120:500+280)', 'k')
title('Type I')
subplot(1,2,2);
plot(time, meanSpike_nl(indRaw == 2, 500-120:500+280)', 'Color', [0.7 0.7 0.7])
title('Type II')
%%