    
tMax = 10;
rwdMatx = [];
for i = 1:tMax
    rwdMatx(i,:) = [nan(1,i) allRewardsBin(1:end-i)];
end
rwdHx = nansum(rwdMatx, 1);
rwdHx(1) = 0;

meow = postCSspikeCount(s.responseInds);
for i = 1:11
    tmp = meow(rwdHx == (i-1));
    meanFR(i) = mean(tmp);
    semFR(i) = std(tmp)/sqrt(length(tmp));
end

figure; errorbar([0:9], meanFR(1:end-1), semFR(1:end-1), 'color', [0 0 0], 'linewidth', 2);
set(gca, 'tickdir', 'out', 'box', 'off')
xlim([-0.5 9.5])
ylabel('Post CS spike count')
xlabel('Number of reward in previous 10 trials')


for i = 1:11
    tmpInds = [tmpInds find(rwdHx == (i-1))];
end
figure;
plotSpikeRaster(allTrial_spike(1,tmpInds),'PlotType','vertline');
set(gca,'Xticklabel',[]);

figure;
LineFormat.LineWidth = 1;
plotSpikeRaster(allTrial_spike(i,rwdHx_Inds),'PlotType','vertline', 'LineFormat', LineFormat);



binInds = discretize(diffQ, 10);
woof = preCSspikeCount(s.responseInds);
for meow = 1:10
    tmp = woof(binInds == meow);
    meanFR(meow) = mean(tmp);
    semFR(meow) = std(tmp)/sqrt(length(tmp));
end

figure; errorbar(meanFR, semFR);
    


numBins = 8;
woof = ceil(length(s.responseInds)/numBins) * numBins + 1;
binEdges = [1:numBins:woof];
binnedRwds = histcounts([find(allRewardsBin == 1)], binEdges);
trialInds = discretize([1:length(s.responseInds)], binEdges);
trialInds(isnan(trialInds)) = max(trialInds) + 1;
meow = postCSspikeCount(s.responseInds);
binnedSpikes = [];
for tmp = 1:max(trialInds)
    binnedSpikes(tmp) = sum(meow(trialInds == tmp));
end

x = binEdges(1:end-1) + diff(binEdges)/2;
pp = spline(x,binnedRwds/8);
ppp = spline(x,binnedSpikes/numBins);


figure; hold on;
yyaxis left; plot(meow, '-', 'linewidth', 2, 'Color', [0.8 0.8 0.8]); fnplt(ppp, 'k');
ylabel('Post-CS spike count')
yyaxis right; fnplt(pp, 'c')
ylabel('Rewards per trial')
set(gca, 'tickdir', 'out')





numBins = 5;
woof = ceil(length(s.responseInds)/numBins) * numBins + 1;
binEdges = [1:numBins:woof];
trialInds = discretize([1:length(s.responseInds)], binEdges);
trialInds(isnan(trialInds)) = max(trialInds) + 1;
meow = zscore(maxFRcs(s.responseInds));
binnedSpikes = [];
binnedR = [];
for tmp = 1:max(trialInds)
    binnedSpikes(tmp) = sum(meow(trialInds == tmp));
    binnedR(tmp) = sum(R(trialInds == tmp));
end

x = binEdges(1:end-1) + diff(binEdges)/2;
pp = spline(x,binnedR/3);
ppp = spline(x,binnedSpikes/numBins);


figure; hold on;
yyaxis left; plot(meow, '-', 'linewidth', 2, 'Color', [0.8 0.8 0.8]); fnplt(ppp, 'k');
ylabel('Max firing rate during CS')
ylim([-5 4.2])
yyaxis right; plot(R, '-', 'Color', [0.7 0 1], 'linewidth', 2)
ylim([-0.2 1])
%fnplt(pp, 'c')
ylabel('rBar')
xlabel('Trials')
set(gca, 'tickdir', 'out')


