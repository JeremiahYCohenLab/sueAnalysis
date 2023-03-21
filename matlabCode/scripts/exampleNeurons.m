ind1 = populationSig(5,1,:)<0;
ind2 = populationSig(5,1,:)>0;
ind3 = populationSigSW(5,2,:)~=0;
figure2; hold on;
color1 = [0, 0.8, 0.8];
color2 = [1, 0.3, 1];
scatter3(populationTStats(5, 1, :), populationTStats(5,3,:), populationTStatsSW(5,2,:), 12, [0.7, 0.7, 0.7], 'filled')
% scatter3(populationTStats(5, 1, ind3), populationTStats(5,3,ind3), populationTStatsSW(5,2,ind3), 40, [1, 1, 0], 'filled');
scatter3(populationTStats(5, 1, ind1), populationTStats(5,3,ind1), populationTStatsSW(5,2,ind1), 12, [0, 0.8, 0.8], 'filled');
scatter3(populationTStats(5, 1, ind2), populationTStats(5,3,ind2), populationTStatsSW(5,2,ind2), 12, [1, 0.3, 1], 'filled');
% scatter3(populationTStats(5, 1, ind1&ind2), populationTStats(5,3,ind1&ind2), populationTStatsSW(5,2,ind1&ind2), 12, [0.5, 0.5, 0.9], 'filled');

% xlabel('outcome'); ylabel('Qchosen'); zlabel('switch')
xlim([-15 20]);
set(gca, 'XTick', [-20:10:20])
set(gca, 'YTick', [-10:5:5])
set(gca,'tickdir', 'out')
ylim([-10 5]);
%%
figure2;
hold on;
histogram(populationTStats(5,1,ind1), -15:1:20, 'FaceColor', [0 0.8 0.8]);
histogram(populationTStats(5,1,ind2), -15:1:20, 'FaceColor', [1, 0.3, 1]);
set(gca, 'YColor', 'none')
set(gca, 'XColor', 'none')
%%
figure2;
hold on;
histogram(populationTStats(5,3,ind1), -10:0.5:5, 'FaceColor', [0 0.8 0.8]);
histogram(populationTStats(5,3,ind2), -10:0.5:5, 'FaceColor', [1, 0.3, 1]);
set(gca, 'YColor', 'none')
set(gca, 'XColor', 'none')
%% rate - pe
load F:\tmpData\gridPE.mat
figure2;
currSpikes = [focusSpikes{1,1}; focusSpikes{2,1}; focusSpikes{3,1}];
currTarget = [allPe{1,1}; allPe{2,1}; allPe{3,1}];
spikeMeans = mean(currSpikes, 'omitnan');
spikeSems = sem(currSpikes);
targetMeans = mean(currTarget, 'omitnan');

subplot(2,2,1); hold on;
plot(targetMeans, spikeMeans, 'color', color1, 'linewidth', 2, 'Marker', 'none', 'LineStyle', '-');
fill([targetMeans fliplr(targetMeans)], [spikeMeans+spikeSems fliplr(spikeMeans-spikeSems)], color1, 'facealpha', 0.25, 'edgecolor', 'none')
set(gca,'tickdir', 'out')
set(gca, 'XTick', [])
set(gca, 'YTick', [-0.4:0.4:0.4])
set(gca, 'XColor', 'none')
ylim([-0.4 0.5])
xlim([-2.2 2.2])


currSpikes = [focusSpikes{1,3}; focusSpikes{2,3}; focusSpikes{3,3}];
currTarget = [allPe{1,3}; allPe{2,3}; allPe{3,3}];
spikeMeans = mean(currSpikes, 'omitnan');
spikeSems = sem(currSpikes);
targetMeans = mean(currTarget, 'omitnan');

subplot(2,2,3); hold on;
plot(targetMeans, spikeMeans, 'color', color2, 'linewidth', 2, 'Marker', 'none', 'LineStyle', '-');
fill([targetMeans fliplr(targetMeans)], [spikeMeans+spikeSems fliplr(spikeMeans-spikeSems)], color2, 'facealpha', 0.25, 'edgecolor', 'none')
set(gca,'tickdir', 'out')
set(gca, 'XTick', [-2:1:2])
set(gca, 'YTick', [-0.4:0.4:0.4])
ylim([-0.6 0.5])
xlim([-2.2 2.2])
%%
load F:\tmpData\gridPE.mat
load F:\tmpData\slideTimeForGrd.mat
colors = zeros(6,3);
for c = 1:3
    colors(c,:) = [1 1.2-c*2/5 1.2-c*2/5];
    colors(c+3,:) = [(c-1)*2/5 (c-1)*2/5 1];
end

figure2;
subplot(2,1,1); hold on;
currAllSpikes = cat(3,allSpikes{1,1}, allSpikes{2,1}, allSpikes{3,1});
for k = 1:size(currAllSpikes,1)
    spikesTemp = squeeze(currAllSpikes(k,:,:));% take all spikes from all neurons from first percental
    plotFilled(slideTime, spikesTemp', colors(k,:));
end
set(gca,'tickdir', 'out')
set(gca, 'XTick', [])
set(gca, 'YTick', [-0.5:0.5:1])
set(gca, 'XColor', 'none')
ylim([-0.5 1.0])
xlim([-1000 2000])
title('Group 1','FontSize', 12)
subplot(2,1,2); hold on;
currAllSpikes = cat(3,allSpikes{1,3}, allSpikes{2,3}, allSpikes{3,3});
for k = 1:size(currAllSpikes,1)
    spikesTemp = squeeze(currAllSpikes(k,:,:));% take all spikes from all neurons from first percental
    plotFilled(slideTime, spikesTemp', colors(k,:));
end
set(gca,'tickdir', 'out')
set(gca, 'XTick', [])
set(gca, 'YTick', [-0.5:0.5:1])
set(gca, 'XTick', [-1000:1000:2000])
ylim([-0.5 1.1])
xlim([-1000 2000])
title('Group 2','FontSize', 12)
xlabel('Time-choice(ms)', 'FontSize', 12)
%% pSwitch rate
load F:\tmpData\gridSW.mat

currSpikes = [focusSpikes{1,1}; focusSpikes{2,1}; focusSpikes{3,1}];
currTarget = [allPe{1,1}; allPe{2,1}; allPe{3,1}];
spikeMeans = mean(currSpikes, 'omitnan');
spikeSems = sem(currSpikes);
targetMeans = mean(currTarget, 'omitnan');

subplot(2,2,2); hold on;
plot(targetMeans, spikeMeans, 'color', color1, 'linewidth', 2, 'Marker', 'none', 'LineStyle', '-');
fill([targetMeans fliplr(targetMeans)], [spikeMeans+spikeSems fliplr(spikeMeans-spikeSems)], color1, 'facealpha', 0.25, 'edgecolor', 'none')
set(gca,'tickdir', 'out')
set(gca, 'XTick', [])
set(gca, 'YTick', [0:0.2:0.4])
set(gca, 'XColor', 'none')
ylim([0.1 0.4])
xlim([-2 4])


currSpikes = [focusSpikes{1,3}; focusSpikes{2,3}; focusSpikes{3,3}];
currTarget = [allPe{1,3}; allPe{2,3}; allPe{3,3}];
spikeMeans = mean(currSpikes, 'omitnan');
spikeSems = sem(currSpikes);
targetMeans = mean(currTarget, 'omitnan');

subplot(2,2,4); hold on;
plot(targetMeans, spikeMeans, 'color', color2, 'linewidth', 2, 'Marker', 'none', 'LineStyle', '-');
fill([targetMeans fliplr(targetMeans)], [spikeMeans+spikeSems fliplr(spikeMeans-spikeSems)], color2, 'facealpha', 0.25, 'edgecolor', 'none')
set(gca,'tickdir', 'out')
set(gca, 'XTick', [-2:2:4])
set(gca, 'YTick', [0:0.2:0.4])
ylim([0.1 0.4])
xlim([-2 4])
%% behavior
%% Plot smoothed rwd&choices
figure2;
subplot(1,2,1); hold on
normKern = normpdf(-15:15,0,4);
normKern = normKern / sum(normKern);
xVals = (1:(length(normKern) + length(allChoices) - 1));
plot(conv(allChoices,normKern)/max(abs(conv(allChoices,normKern))),xVals,'k','linewidth',2);
plot(conv(allRewards,normKern)/max(abs(conv(allRewards,normKern))),xVals,'--','Color',[0.2 0.6 1],'linewidth',2)
ylabel('Trials', 'FontSize', 12)
xlabel('<-- Left       Right -->', 'FontSize', 12)
ylim([1 length(allChoice_R)])
xlim([-1.2 1.2])
ylim([-3 length(behSessionData)])
set(gca, 'YTick', [0:100:100+length(behSessionData)])
set(gca, 'XTick', [])
set(gca,'tickdir', 'out')
legend({'Choices','Rewards'}, 'FontSize', 12)
title('Animal behavior', 'FontSize', 12)
%% raster
spikeRasters_dF_choiceAlign('mZS061d20210324','cellName','TT3_SS_01', 'saveFigFlag', 0);
figure2;
subplot(1,2,2);
plotSpikeRaster(allTrial_spike_choice(1,:),'PlotType','vertline'); hold on;
set(gca, 'YColor', 'none')
ylim([-3 450])
xlim([-1000 2000])
set(gca, 'YDir','reverse')
xlabel('Time-choice(ms)', 'FontSize', 12)
set(gca,'YtickLabel',flip(0:100:100+length(allTrial_spike_choice)))
set(gca, 'Box', 'off')
title('Example neuron', 'FontSize', 12)
%% rwd/norwd
figure2; 
%example 1
spikeRasters_dF_choiceAlign('mZS061d20210324','cellName','TT3_SS_01', 'saveFigFlag', 0);
subplot(3,1,1); hold on;
mySDF_rwd = allTrial_spikeMatx_slide(s.rwd_Inds,:);
mySDF_noRwd = allTrial_spikeMatx_slide(s.nrwd_Inds,:);
plotFilled(slideTime, mySDF_rwd,[0 0 1])
plotFilled(slideTime, mySDF_noRwd,[0.7 0 1])
line([s.rwdDelay s.rwdDelay], [0 maxFreq], 'color', [1 0.3 0.3], 'LineStyle', '--')
ylim([0 maxFreq]);
xlim([-1000 2000]);
set(gca, 'YTick', [0:5:10])
set(gca, 'XColor', 'none')
set(gca,'tickdir', 'out')
title('Example Unit 1', 'FontSize', 12)

%%
%example 2
subplot(3,1,2); hold on;
spikeRasters_dF_choiceAlign('mZS061d20210405','cellName','TT5_SS_01', 'saveFigFlag', 0);
mySDF_rwd = allTrial_spikeMatx_slide(s.rwd_Inds,:);
mySDF_noRwd = allTrial_spikeMatx_slide(s.nrwd_Inds,:);
plotFilled(slideTime, mySDF_rwd,[0 0 1])
plotFilled(slideTime, mySDF_noRwd,[0.7 0 1])
line([s.rwdDelay s.rwdDelay], [0 maxFreq], 'color', [1 0.3 0.3], 'LineStyle', '--')
ylim([0 maxFreq]);
xlim([-1000 2000]);
set(gca, 'YTick', [0:5:10])
set(gca, 'XColor', 'none')
set(gca,'tickdir', 'out')
title('Unit 2', 'FontSize', 12)
%%
% %example 3
% subplot(3,1,3); hold on;
% spikeRasters_dF_choiceAlign('mZS061d20210326','cellName','TT1_SS_02', 'saveFigFlag', 0);
% mySDF_rwd = allTrial_spikeMatx_slide(s.rwd_Inds,:);
% mySDF_noRwd = allTrial_spikeMatx_slide(s.nrwd_Inds,:);
% plotFilled(slideTime, mySDF_rwd,[0 0 1])
% plotFilled(slideTime, mySDF_noRwd,[0.7 0 1])
% line([s.rwdDelay s.rwdDelay], [0 5], 'color', [1 0.3 0.3], 'LineStyle', '--')
% ylim([0 5]);
% xlim([-1000 2000]);
% set(gca, 'YTick', [0:5:10])
% set(gca, 'XTick', [-1000:1000:2000])
% set(gca,'tickdir', 'out')
% title('Unit 3', 'FontSize', 12)
% xlabel('Time-choice(ms)','FontSize', 12)
%% rwd/norwd, PE PSTH& raster
figure2; 
%example 1
spikeRasters_dF_choiceAlign('mZS061d20210324','cellName','TT3_SS_01', 'saveFigFlag', 0);
subplot(3,1,1); hold on;
mySDF_rwd = allTrial_spikeMatx_slide(s.rwd_Inds,:);
mySDF_noRwd = allTrial_spikeMatx_slide(s.nrwd_Inds,:);
plotFilled(slideTime, mySDF_rwd,[0 0 1])
plotFilled(slideTime, mySDF_noRwd,[0.7 0 1])
line([s.rwdDelay s.rwdDelay], [0 maxFreq], 'color', [1 0.3 0.3], 'LineStyle', '--')
ylim([0 maxFreq]);
xlim([-1000 2000]);
set(gca, 'YTick', [0:5:10])
set(gca, 'XColor', 'none')
set(gca,'tickdir', 'out')
title('Example Unit 1', 'FontSize', 12)

%example 2
figure2; hold on;
spikeGLM_dF('mZS061d20210405', 'good', 'saveFigFlag', 0,'cellName', 'TT5_SS_01');
[~, ind] = sort(pe);
sortedSpikes = allTrial_spike_choice(ind);
plotSpikeRaster(sortedSpikes,'PlotType','vertline');
line([s.rwdDelay s.rwdDelay], [0 maxFreq], 'color', [1 0.3 0.3], 'LineStyle', '--')
xlim([-1000 2000]);
set(gca, 'YTick', [])
ylabel('sorted by pe 1 <----> -1')
title('Unit 2', 'FontSize', 12)
% psth sorted by PE
[~, ind] = sort(pe);
indPos = ismember(ind, find(pe>0));
indPos = ind(indPos);
indNeg = ismember(ind, find(pe<=0));
indNeg = ind(indNeg);
indI = indNeg(1:floor(1/2*length(indNeg)));
indII = indNeg(floor(1/2*length(indNeg))+1:end);
indIII = indPos(1:floor(1/2*length(indPos)));
indIV = indPos(floor(1/2*length(indPos))+1:end);

mySDF_I = allTrial_spikeMatx_slide(indI,:);
mySDF_II = allTrial_spikeMatx_slide(indII,:);
mySDF_III = allTrial_spikeMatx_slide(indIII,:);
mySDF_IV = allTrial_spikeMatx_slide(indIV,:);

figure2; hold on
plotFilled(slideTime, mySDF_I,[0 0 1])
plotFilled(slideTime, mySDF_II,[0.5 0.5 1])
plotFilled(slideTime, mySDF_III,[1 0.5 0.5])
plotFilled(slideTime, mySDF_IV,[1 0 0])

%%
%example 3
subplot(3,1,3); hold on;
spikeRasters_dF_choiceAlign('mZS061d20210326','cellName','TT1_SS_02', 'saveFigFlag', 0);
mySDF_rwd = allTrial_spikeMatx_slide(s.rwd_Inds,:);
mySDF_noRwd = allTrial_spikeMatx_slide(s.nrwd_Inds,:);
plotFilled(slideTime, mySDF_rwd,[0 0 1])
plotFilled(slideTime, mySDF_noRwd,[0.7 0 1])
line([s.rwdDelay s.rwdDelay], [0 5], 'color', [1 0.3 0.3], 'LineStyle', '--')
ylim([0 5]);
xlim([-1000 2000]);
set(gca, 'YTick', [0:5:10])
set(gca, 'XTick', [-1000:1000:2000])
set(gca,'tickdir', 'out')
title('Unit 3', 'FontSize', 12)
xlabel('Time-choice(ms)','FontSize', 12)
%% rwd norwd raster
spikeRasters_dF_choiceAlign('mZS061d20210324','cellName','TT3_SS_01', 'saveFigFlag', 0);

i = 1;
figure2;
plotSpikeRaster(allTrial_spike_choice(i,[s.rwd_Inds s.nrwd_Inds]),'PlotType','vertline'); hold on
plot([-5000 10000],[length(s.rwd_Inds) length(s.rwd_Inds)],'r')
xlim([-1000 2000])
line([s.rwdDelay s.rwdDelay], [0 length(s.responseRateInds)], 'color', 'r')
set(gca,'XTick', [-1000:1000:2000]);
xlabel('Time-choice(ms)')
ylabel('trials')
set(gca, 'YTick', []);
%% spikes-pe separated by side
%example 1
spikeGLM_dF('mZS061d20210324', 'good', 'saveFigFlag', 0, 'cellName','TT3_SS_01', 'numBins', 6);
figure2; 
subplot(1,2,1); hold on;
errorbar(targetMeans, spikeMeans, spikeSems,'linewidth', 2)
title('allTrials')
ylim([2 8])
xlim([-1 1])
set(gca, 'tickdir', 'out')
set(gca, 'YTick', 2:2:8)
set(gca, 'XTick', [-1 0 1])
ylabel('spikes/s')
xlabel('rpe')
subplot(1,2,2); hold on;
errorbar(targetMeansL, spikeMeansL, spikeSemsL,'linewidth', 2, 'color', colorL)
errorbar(targetMeansR, spikeMeansR, spikeSemsR,'linewidth', 2, 'color', colorR)
ylim([2 8])
xlim([-1 1])
set(gca, 'tickdir', 'out')
set(gca, 'YTick', 2:2:8)
set(gca, 'XTick', [-1 0 1])
legend({'L', 'R'})
title('split by side')
xlabel('rpe')

% example 2
spikeGLM_dF('mZS061d20210329', 'good', 'saveFigFlag', 0,'cellName', 'TT1_SS_03', 'numBins', 6);
figure2; 
subplot(1,2,1); hold on;
errorbar(targetMeans, spikeMeans, spikeSems,'linewidth', 2)
title('allTrials')
ylim([0.5 4])
xlim([-1 1])
set(gca, 'tickdir', 'out')
set(gca, 'YTick', 0:2:4)
set(gca, 'XTick', [-1 0 1])
ylabel('spikes/s')
xlabel('rpe')
subplot(1,2,2); hold on;
errorbar(targetMeansL, spikeMeansL, spikeSemsL,'linewidth', 2, 'color', colorL)
errorbar(targetMeansR, spikeMeansR, spikeSemsR,'linewidth', 2, 'color', colorR)
ylim([0.5 4])
xlim([-1 1])
set(gca, 'tickdir', 'out')
set(gca, 'YTick', 0:2:4)
set(gca, 'XTick', [-1 0 1])
legend({'L', 'R'})
title('split by side')
xlabel('rpe')
%% spikes/s rpe scatter
spikeGLM_dF('mZS061d20210324', 'good', 'cellName', 'TT3_SS_01', 'saveFigFlag', 0, 'numBins', 6);
focusWin = [300 1000];
focusSpikes = sum(allTrial_spikeMatx_choice(:, focusWin(1)+p.Results.tb*1000: focusWin(2)+p.Results.tb*1000), 2)/((focusWin(2) - focusWin(1))/1000);
%% rwd/norwd, PE PSTH& raster split by side 
figure2;  
%example 1
spikeRasters_dF_choiceAlign('mZS061d20210324','cellName','TT3_SS_01', 'saveFigFlag', 0, 'numBins', 6);
spikeGLM_dF('mZS061d20210324', 'good', 'cellName', 'TT3_SS_01', 'saveFigFlag', 0, 'numBins', 6);
figure2;
subplot(1,2,1); hold on;
title('Left')
mySDF_rwd = allTrial_spikeMatx_slide(intersect(os.rwd_Inds, os.lickL_Inds),:);
mySDF_noRwd = allTrial_spikeMatx_slide(intersect(os.nrwd_Inds, os.lickL_Inds),:);
plotFilled(slideTime, mySDF_rwd,[0 0 1])
plotFilled(slideTime, mySDF_noRwd,[0.7 0 1])
line([os.rwdDelay os.rwdDelay], [0 15], 'color', [1 0.3 0.3], 'LineStyle', '--')
ylim([0 15]);
xlim([-1000 2000]);
set(gca, 'YTick', [0:5:10])
set(gca, 'XTick', [-1000:1000:2000])
set(gca,'tickdir', 'out')
legend({'rwd', '', 'nrwd', ''})
subplot(1,2,2); hold on;
xlabel('time from respond /ms')
ylabel('spikes/s')
title('Right')
mySDF_rwd = allTrial_spikeMatx_slide(intersect(os.rwd_Inds, os.lickR_Inds),:);
mySDF_noRwd = allTrial_spikeMatx_slide(intersect(os.nrwd_Inds, os.lickR_Inds),:);
plotFilled(slideTime, mySDF_rwd,[0 0 1])
plotFilled(slideTime, mySDF_noRwd,[0.7 0 1])
line([os.rwdDelay os.rwdDelay], [0 15], 'color', [1 0.3 0.3], 'LineStyle', '--')
ylim([0 15]);
xlim([-1000 2000]);
set(gca, 'YTick', [0:5:10])
set(gca,'tickdir', 'out')
legend({'rwd', '', 'nrwd', ''})
xlabel('time from respond /ms')

%example 2
figure2; hold on;
spikeGLM_dF('mZS061d20210329', 'good', 'saveFigFlag', 0,'cellName', 'TT1_SS_03');
% both sides
[~, ind] = sort(pe);
subplot(1,2,1);
sortedSpikes = allTrial_spike_choice(ind);
plotSpikeRaster(sortedSpikes,'PlotType','vertline');
line([os.rwdDelay os.rwdDelay], [0 length(sortedSpikes)], 'color', [1 0.3 0.3], 'LineStyle', '--')
line([-1000 2000], [length(os.nrwd_Inds), length(os.nrwd_Inds)], 'LineStyle', '--', 'color', [1 0.3 0.3]);

xlim([-1000 2000]);
ylim([0 length(ind)])
set(gca, 'YTick', [])
xlabel('time from choice (ms)', 'FontSize', 18);
ylabel('sorted by pe 1 <----> -1', 'FontSize', 18);
text(2100, length(os.nrwd_Inds), 'rpe = 0')

% both sides psth
[~, ind] = sort(pe);
indPos = ismember(ind, find(pe>0));
indPos = ind(indPos);
indNeg = ismember(ind, find(pe<=0));
indNeg = ind(indNeg);
indI = indNeg(1:floor(1/2*length(indNeg)));
indII = indNeg(floor(1/2*length(indNeg))+1:end);
indIII = indPos(1:floor(1/2*length(indPos)));
indIV = indPos(floor(1/2*length(indPos))+1:end);

mySDF_I = allTrial_spikeMatx_slide(indI,:);
mySDF_II = allTrial_spikeMatx_slide(indII,:);
mySDF_III = allTrial_spikeMatx_slide(indIII,:);
mySDF_IV = allTrial_spikeMatx_slide(indIV,:);

figure;
hold on
title('all trials')
plotFilled(slideTime, mySDF_I,[0 0 1])
plotFilled(slideTime, mySDF_II,[0.5 0.5 1])
plotFilled(slideTime, mySDF_III,[1 0.5 0.5])
plotFilled(slideTime, mySDF_IV,[1 0 0])
line([os.rwdDelay os.rwdDelay], [0 7], 'color', [1 0.3 0.3], 'LineStyle', '--')
ylim([0 12]);
xlim([-1000 2000]);
set(gca, 'XTick', [-1000:1000:2000])
set(gca, 'YTick', [0:5:10])
set(gca,'tickdir', 'out')
xlabel('time from choice')
ylabel('spikes/s')



% two sides psth

[~, ind] = sort(pe);
indL = ismember(ind, os.lickL_Inds);
indL = ind(indL);
indR = ismember(ind, os.lickR_Inds);
indR = ind(indR);
figure2;
subplot(1,2,1);
title('L');
sortedSpikes = allTrial_spike_choice(indL);
plotSpikeRaster(sortedSpikes,'PlotType','vertline');
line([os.rwdDelay os.rwdDelay], [0 length(sortedSpikes)], 'color', [1 0.3 0.3], 'LineStyle', '--')
xlim([-1000 2000]);
ylim([0 max([length(indL) length(indR)])])
set(gca, 'YTick', [])
ylabel('sorted by pe 1 <----> -1')
subplot(1,2,2);
title('R');
sortedSpikes = allTrial_spike_choice(indR);
plotSpikeRaster(sortedSpikes,'PlotType','vertline');
line([os.rwdDelay os.rwdDelay], [0 length(sortedSpikes)], 'color', [1 0.3 0.3], 'LineStyle', '--')
xlim([-1000 2000]);
ylim([0 max([length(indL) length(indR)])])
set(gca, 'YTick', [])

% psth sorted by PE split by side
[~, ind] = sort(pe);
indL = ismember(ind, os.lickL_Inds);
indL = ind(indL);
indR = ismember(ind, os.lickR_Inds);
indR = ind(indR);
% L
indPos = ismember(indL, find(pe>0));
indPos = indL(indPos);
indNeg = ismember(indL, find(pe<=0));
indNeg = indL(indNeg);
indI = indNeg(1:floor(1/2*length(indNeg)));
indII = indNeg(floor(1/2*length(indNeg))+1:end);
indIII = indPos(1:floor(1/2*length(indPos)));
indIV = indPos(floor(1/2*length(indPos))+1:end);

mySDF_I = allTrial_spikeMatx_slide(indI,:);
mySDF_II = allTrial_spikeMatx_slide(indII,:);
mySDF_III = allTrial_spikeMatx_slide(indIII,:);
mySDF_IV = allTrial_spikeMatx_slide(indIV,:);

figure
subplot(1,2,1); hold on
title('L')
plotFilled(slideTime, mySDF_I,[0 0 1])
plotFilled(slideTime, mySDF_II,[0.5 0.5 1])
plotFilled(slideTime, mySDF_III,[1 0.5 0.5])
plotFilled(slideTime, mySDF_IV,[1 0 0])
line([os.rwdDelay os.rwdDelay], [0 7], 'color', [1 0.3 0.3], 'LineStyle', '--')
ylim([0 7]);
xlim([-1000 2000]);
set(gca, 'XTick', [-1000:1000:2000])
set(gca, 'YTick', [0:5:10])
set(gca,'tickdir', 'out')
xlabel('time from choice')
ylabel('spikes/s')


% R
indPos = ismember(indR, find(pe>0));
indPos = indR(indPos);
indNeg = ismember(indR, find(pe<=0));
indNeg = indR(indNeg);
indI = indNeg(1:floor(1/2*length(indNeg)));
indII = indNeg(floor(1/2*length(indNeg))+1:end);
indIII = indPos(1:floor(1/2*length(indPos)));
indIV = indPos(floor(1/2*length(indPos))+1:end);

mySDF_I = allTrial_spikeMatx_slide(indI,:);
mySDF_II = allTrial_spikeMatx_slide(indII,:);
mySDF_III = allTrial_spikeMatx_slide(indIII,:);
mySDF_IV = allTrial_spikeMatx_slide(indIV,:);
subplot(1,2,2); hold on
title('R')
plotFilled(slideTime, mySDF_I,[0 0 1])
plotFilled(slideTime, mySDF_II,[0.5 0.5 1])
plotFilled(slideTime, mySDF_III,[1 0.5 0.5])
plotFilled(slideTime, mySDF_IV,[1 0 0])
line([os.rwdDelay os.rwdDelay], [0 7], 'color', [1 0.3 0.3], 'LineStyle', '--')
ylim([0 7]);
xlim([-1000 2000]);
set(gca, 'YTick', [0:5:10])
set(gca,'tickdir', 'out')
set(gca, 'XTick', [-1000:1000:2000])
xlabel('time from choice')
%% pupil 
pupilSessionAnalysis('mZS070d20211021');

figure;hold on;
baseline = mean(sessionPupilChoice(:,1:round(2*FR)),2,'omitnan');
sessionPupilChange = sessionPupilChoice - baseline;
my_SDF_1 = sessionPupilChange;
plotFilled(time, my_SDF_1,[0.3 0.3 0.3]);

pupilSessionAnalysis('mZS070d20211022');
baseline = mean(sessionPupilChoice(:,1:round(2*FR)),2,'omitnan');
sessionPupilChange = sessionPupilChoice - baseline;
my_SDF_1 = sessionPupilChange;
plotFilled(time, my_SDF_1,[0 1 0.3]);

xlim([-1.5 6])
set(gca, 'XTick', [-1 0 2 4 6])
ylabel('pupil dilation/px',"FontSize",12)
set(gca, 'YTick', [-3:3:6]);
ylim([-3 7])
xlabel('Time-choice(s)',"FontSize",12)
%% simulation
tonic = 2;
excMax = 1;
inhiMax = -1;
bins = 11;
gain = 4;
gainNE = 6;
inhi = 1;
inhiNE = 1.1;
b = 0.8;

connectionR = [linspace(excMax,inhiMax,bins);linspace(inhiMax,excMax,bins)];
connectionL = [linspace(excMax,inhiMax,bins);linspace(inhiMax,excMax,bins)];
%%
resolutionBiasGain = zeros(bins, bins);
RR = zeros(bins, bins);
RN = zeros(bins, bins);
LR = zeros(bins, bins);
LN = zeros(bins, bins);
baseline = FIcurve(tonic,'a', gain, 'b', b, 'd', 8);
for i = 1:bins
    for j = 1:bins
        RR(i,j) = FIcurve(connectionR(1,i)+tonic,'a', gain, 'b', b, 'd', 8); 

        RN(i,j) = FIcurve(connectionR(2,i)+tonic,'a', gainNE, 'b', b, 'd', 8);

        LR(i,j) = FIcurve(connectionL(1,j)+tonic,'a', gain, 'b', b, 'd', 8); 

        LN(i,j) = FIcurve(connectionL(2,j)+tonic,'a', gainNE, 'b', b, 'd', 8);
        
%         resolutionBiasGain(i,j) = RN(i,j) - RR(i,j) + LN(i,j) - LR(i,j);
        resolutionBiasGain(i,j) = RN(i,j) + RR(i,j) - LN(i,j) - LR(i,j);
    end
end
figure2;
subplot(1,2,1);
stem([RR(1,bins)-baseline, RN(1,bins)-baseline, LR(1,bins)-baseline, LN(1,bins)-baseline]);
title('gain')
xlim([0 5]) 
subplot(1,2,2);
stem([RR(bins,1)-baseline, RN(bins,1)-baseline, LR(bins,1)-baseline, LN(bins,1)-baseline]);
xlim([0 5])
%%
resolutionBiasInhi = zeros(bins, bins);
RR = zeros(bins, bins);
RN = zeros(bins, bins);
LR = zeros(bins, bins);
LN = zeros(bins, bins);
baseline = FIcurve(tonic,'a', gain, 'b', b, 'd', 8);
for i = 1:bins
    for j = 1:bins
        RR(i,j) = FIcurve(connectionR(1,i)+tonic,'a', gain, 'b', b, 'd', 8); 
        if connectionR(2,i) < 0
            RN(i,j) = FIcurve(inhiNE * connectionR(2,i)+tonic,'a', gain, 'b', b, 'd', 8); 
        else
            RN(i,j) = FIcurve(inhi * connectionR(2,i)+tonic,'a', gain, 'b', b, 'd', 8);
        end
        LR(i,j) = FIcurve(connectionL(1,j)+tonic,'a', gain, 'b', b, 'd', 8); 
        if connectionL(2,j) < 0
            LN(i,j) = FIcurve(inhiNE * connectionL(2,j)+tonic,'a', gain, 'b', b, 'd', 8);
        else
            LN(i,j) = FIcurve(inhi * connectionL(2,j)+tonic,'a', gain, 'b', b, 'd', 8);
        end

%         resolutionBiasInhi(i,j) = RN(i,j) - RR(i,j) + LN(i,j) - LR(i,j);
        resolutionBiasInhi(i,j) = RN(i,j) + RR(i,j) - LN(i,j) - LR(i,j);
    end
end
figure2;
subplot(1,2,1);
stem([RR(1,bins)-baseline, RN(1,bins)-baseline, LR(1,bins)-baseline, LN(1,bins)-baseline]);
title('inhi')
xlim([0 5]) 
subplot(1,2,2);
stem([RR(bins,1)-baseline, RN(bins,1)-baseline, LR(bins,1)-baseline, LN(bins,1)-baseline]);
xlim([0 5])

%%
mycolormap = customcolormap([0 0.5 1], [0 0 1; 1 1 1; 1 0 0]);
% figure2;
% subplot(2,2,1)
% imagesc(RR);
% set(gca, 'YDir','normal');
% colorbar;
% colormap(mycolormap);
% title('RR')
% subplot(2,2,2)
% imagesc(RN);
% set(gca, 'YDir','normal');
% colorbar;
% colormap(mycolormap);
% title('RN')
% subplot(2,2,3)
% imagesc(LR);
% set(gca, 'YDir','normal');
% colorbar;
% colormap(mycolormap);
% title('LR')
% subplot(2,2,4)
% imagesc(LN);
% set(gca, 'YDir','normal');
% colorbar;
% colormap(mycolormap);
% title('RN')

% figure2;
% subplot(2,1,1)
% imagesc(connectionL);
% set(gca, 'YDir','normal');
% colorbar;
% colormap(mycolormap);
% title('connectionL')
% subplot(2,1,2)
% imagesc(connectionR);
% set(gca, 'YDir','normal');
% colorbar;
% colormap(mycolormap);
% title('connectionR')

figure;
subplot(2,1,1)
imagesc(resolutionBiasGain);
set(gca, 'YDir','normal');
colorbar;
colormap(mycolormap);
title('Gain')
subplot(2,1,2)
imagesc(resolutionBiasInhi);
set(gca, 'YDir','normal');
colorbar;
colormap(mycolormap);
title('Inhi')

% axis off;
%%
gainDiag = diag(flip(resolutionBiasGain));
%%
inhiDiag = diag(flip(resolutionBiasInhi));
%%
figure; hold on;
plot(gainDiag);
plot(inhiDiag);
%%
sizeX = 20;
allColors = 255 * ones(sizeX,sizeX,3);
midX = 1;
midY = 1;
colorA = [255, 0, 255];
colorB = [0, 0, 255];
for i = 1:sizeX
    for j = 1:sizeX
        if i^2+j^2 > 0 
            co = i/sqrt(i^2+j^2);
            d = acosd(co);
            if d <= 45.01
                oriColor = d/45 * colorA + (1-d/45) * colorB;
                dis = sqrt(i^2 + j^2);
                if dis > sizeX
                   allColors(i,j,:) = oriColor;
                else
                    allColors(i,j,:) = dis/sizeX * oriColor + (1-dis/sizeX) * [255 255 255];
                end 
                allColors(i,j,:) = abs(allColors(i,j,:));
                allColors(i,j,allColors(i,j,:)>255) = 255;
            end
        end 
    end
end
%% patch
figure2; hold on;
binSize = 10;
for i = 1:sizeX
    for j = 1:sizeX
        if sqrt(i^2 + j^2)<=20
            patch([binSize*(i-1), binSize*(i), binSize*(i), binSize*(i-1)]', [binSize*(j-1), binSize*(j-1), binSize*(j), binSize*(j)]', squeeze(allColors(i,j,:))'/255, 'EdgeColor','none')
        end
    end
end
%%
optoIDTT('mZS062d20210504', 'subfolder','10ms_1','PulseWidth', 10, 'Trains',10);
figure2; hold on;
plot(time{j}, rawTraces{j}, 'color', 'k');
for k = 1:p.Results.Pulses
    patch([x(k) xx(k) xx(k) x(k)], [max(rawTraces{j}) max(rawTraces{j}) max(rawTraces{j})+10 ma x(rawTraces{j})+10], 'b', 'FaceAlpha', 1, 'EdgeColor','none');
end
plot([-0.25 -0.25], [-150 -100],'color', 'k', 'lineWidth',2);
plot([-0.25 0.25], [-150 -150],'color', 'k', 'lineWidth',2);
text(-0.5,-125, '50 uV', 'HorizontalAlignment','left','FontSize',12)
text(-0.04,-160, '0.5 s', 'HorizontalAlignment','center', 'FontSize',12)
set(gca, 'XColor', 'none');
set(gca, 'YColor', 'none');
%% simutaneuos recordings
load('F:\tmpData\allUnitAUC.mat');
%% example 1
session = 'mZS061d20210328';
units = {'TT1_SS_03';'TT3_SS_01';'TT4_SS_01'}';
%% 2
session = 'mZS061d20210327';
units = {'TT1_SS_01';'TT3_SS_01'}';
%% 3
session = 'mZS061d20210323';
units = {'TT1_SS_02';'TT3_SS_01';'TT4_SS_01'}';
%% plot example
colors = [0.4 0.4 0.4;
          0.7 0.7 0.7];
modelName = '5params';
col = 'good';
numSamps = 2000;

numBins = 5;
[root, sep] = currComputer();
pd = parseSessionString_df(session, root, sep);
params = getStanModelParams_sampsOnly(pd.animalName, col, modelName, numSamps, 'sessionName', session, 'sessionParamsFlag', 1);
t = inferModelVar(session, params, modelName, 'perturb', 0);
target = t.pe;
figure2Wide;
for i = 1:length(units)
    id = strcmp(allSessions, session) & strcmp(allUnits, units{i});
    currSpikes = spikeFocus{id};
    edges = linspace(min(target)-0.0001, max(target)+0.0001, numBins+1);
%     edges = quantile(target, linspace(0, 1, numBins+1));
%     edges(1) = edges(1) - 0.00000001;
%     edges(end) = edges(end) + 0.00000001;
    meanPe = NaN(1, numBins);
    meanSpikes = NaN(1, numBins);
    semSpikes = NaN(1, numBins);
    for j = 1:numBins
        if ~isempty(find(target>edges(j) & target<=edges(j+1), 1))
            tmpInd = target>edges(j) & target<=edges(j+1); 
            meanPe(j) = mean(target(tmpInd));
            meanSpikes(j) = mean(currSpikes(tmpInd));
            semSpikes(j) = sem(currSpikes(tmpInd));
        end
    end
    meanSpikes = meanSpikes(~isnan(meanPe));
    semSpikes = semSpikes(~isnan(meanPe));
    meanPe = meanPe(~isnan(meanPe));
    subplot(1,length(units),i); hold on;
    errorbar(meanPe, meanSpikes, semSpikes, 'LineWidth', 2, 'Color', colors(ind(id),:));
   if i == 1
       ylabel('Spikes/s')
   end
end

sgtitle(session)
%%


