function photometryAnalysisSession(session)
% load data
[root, sep] = currComputer();
pd = parseSessionString_df(session, root, sep);
dataPath = [pd.sortedFolder session '_photometry.mat'];
category = 'good';
modelName = '5params';
numSamps = 2000;
psthBinNum = 4;
if ~exist(dataPath, 'file')
    fprintf('no photometry data')
end
load(dataPath);
myKernel = ones(2,1);
myKernel = myKernel/sum(myKernel);
LCmatChoiceFiltered = [];
mPFCmatChoiceFiltered = [];
for j = 1:size(LCmatChoice,1)
    temp = conv(LCmatChoice(j,:), myKernel);
    temp = temp((floor(0.5*length(myKernel))+1):(end-floor(0.5*length(myKernel))));
    LCmatChoiceFiltered(j,:) = temp;
    
    temp = conv(mPFCmatChoice(j,:), myKernel);
    temp = temp((floor(0.5*length(myKernel))+1):(end-floor(0.5*length(myKernel))));
    mPFCmatChoiceFiltered(j,:) = temp;        
end
midPointsFilter = midPoints + 0.5*length(myKernel)*stepSize;
midPointsFilter = midPointsFilter(1:length(temp));

os = behAnalysisNoPlot_opMD(session, 'simpleFlag', 1);
params = getStanModelParams_sampsOnly(pd.animalName, category, modelName, numSamps, 'sessionName', session);
t = inferModelVar(session, params, modelName);
pe = t.pe;
target = pe;
edges = [linspace(min(target)-0.01, 0, 0.5*psthBinNum+1) linspace(0, max(target)+0.01, 0.5*psthBinNum+1)];
edges = edges([1:0.5*psthBinNum 0.5*psthBinNum+2:end]);

focusWins = {[-2 -1] [0.3 1.8]};


lfig = figure;
sgtitle([session ' LC' ])
mfig = figure;
sgtitle([session ' mPFC' ])
figure(lfig);
subplot(4,6,1);
m1 = LCmat(os.responseInds(os.allRewards~=0),:);
m2 = LCmat(os.responseInds(os.allRewards==0),:);
plotFilled(midPoints, m1, 'b');
plotFilled(midPoints, m2, 'r');
legend({'rwd', '', 'norwd', ''})
figure(mfig);
subplot(4,6,1);
m1 = mPFCmat(os.responseInds(os.allRewards~=0),:);
m2 = mPFCmat(os.responseInds(os.allRewards==0),:);
plotFilled(midPoints, m1, 'b');
plotFilled(midPoints, m2, 'r');
legend({'rwd', '', 'norwd', ''})

figure(lfig);
subplot(4,6,2);
m1 = LCmat(os.responseInds(os.changeChoice_Inds),:);
m2 = LCmat(os.responseInds(os.stayChoice_Inds),:);
plotFilled(midPoints, m1, 'b');
plotFilled(midPoints, m2, 'r');
legend({'change', '', 'stay', ''})
figure(mfig);
subplot(4,6,2);
m1 = mPFCmat(os.responseInds(os.changeChoice_Inds),:);
m2 = mPFCmat(os.responseInds(os.stayChoice_Inds),:);
plotFilled(midPoints, m1, 'b');
plotFilled(midPoints, m2, 'r');
legend({'change', '', 'stay', ''})

figure(lfig);
subplot(4,6,3);
m1 = LCmat(os.CSplus,:);
m2 = LCmat(os.CSminus,:);
plotFilled(midPoints, m1, 'b');
plotFilled(midPoints, m2, 'r');
legend({'plus', '', 'minus', ''})
figure(mfig);
subplot(4,6,3);
m1 = mPFCmat(os.CSplus,:);
m2 = mPFCmat(os.CSminus,:);
plotFilled(midPoints, m1, 'b');
plotFilled(midPoints, m2, 'r');
legend({'plus', '', 'minus', ''})

figure(lfig);
subplot(4,6,4);
m1 = LCmat(os.responseInds(os.allChoices==1),:);
m2 = LCmat(os.responseInds(os.allChoices~=1),:);
plotFilled(midPoints, m1, 'b');
plotFilled(midPoints, m2, 'r');
legend({'right', '', 'left', ''})
figure(mfig);
subplot(4,6,4);
m1 = mPFCmat(os.responseInds(os.allChoices==1),:);
m2 = mPFCmat(os.responseInds(os.allChoices~=1),:);
plotFilled(midPoints, m1, 'b');
plotFilled(midPoints, m2, 'r');
legend({'right', '', 'left', ''})

figure(lfig);
subplot(4,6,5);
m1 = LCmat(os.responseInds(intersect(find(os.allChoices==1), os.rwd_Inds)),:);
m2 = LCmat(os.responseInds(intersect(find(os.allChoices==1), os.nrwd_Inds)),:);
plotFilled(midPoints, m1, 'b');
plotFilled(midPoints, m2, 'r');
legend({'Rrwd', '', 'Rnorwd', ''})
figure(mfig);
subplot(4,6,5);
m1 = mPFCmat(os.responseInds(intersect(find(os.allChoices==1), os.rwd_Inds)),:);
m2 = mPFCmat(os.responseInds(intersect(find(os.allChoices==1), os.nrwd_Inds)),:);
plotFilled(midPoints, m1, 'b');
plotFilled(midPoints, m2, 'r');
legend({'Rrwd', '', 'Rnorwd', ''})

figure(lfig);
subplot(4,6,6);
m1 = LCmat(os.responseInds(intersect(find(os.allChoices==-1), os.rwd_Inds)),:);
m2 = LCmat(os.responseInds(intersect(find(os.allChoices==-1), os.nrwd_Inds)),:);
plotFilled(midPoints, m1, 'b');
plotFilled(midPoints, m2, 'r');
legend({'Rrwd', '', 'Rnorwd', ''})
figure(mfig);
subplot(4,6,6);
m1 = mPFCmat(os.responseInds(intersect(find(os.allChoices==-1), os.rwd_Inds)),:);
m2 = mPFCmat(os.responseInds(intersect(find(os.allChoices==-1), os.nrwd_Inds)),:);
plotFilled(midPoints, m1, 'b');
plotFilled(midPoints, m2, 'r');
legend({'Lrwd', '', 'Lnorwd', ''})

figure(lfig);
subplot(4,3,[4,7,10]); hold on;
[~, sortedInd] = sort(pe);
imagesc(midPointsFilter, 1:size(LCmatChoiceFiltered,1), LCmatChoiceFiltered(sortedInd,:));
plot([0.001*os.rwdDelay 0.001*os.rwdDelay], [0 size(LCmatChoiceFiltered,1)],'Color', 'w', 'LineWidth', 2, 'LineStyle', '--');
xlim([-1 2.5])
myMap = [[linspace(0, 1, 50)', linspace(0, 1, 50)', linspace(1, 1, 50)'];
[linspace(1, 1, 100)', linspace(1, 0, 100)', linspace(1, 0, 100)']];
colormap(myMap)
colorbar
% maxScale = max(abs(waveformsSession), [], 'all');
caxis([-2 4])


currWin = focusWins{1};
bl = mean(LCmatChoice(:, midPoints>=currWin(1)&midPoints<currWin(2)), 2);
for i = 2:length(focusWins)
    currWin = focusWins{i};
    focusMean = mean(LCmatChoice(:, midPoints>=currWin(1)&midPoints<currWin(2)), 2) - bl;
    subplot(4,6,[11,17,23]); hold on;
    imagesc(currWin, [1 length(focusMean)], zscore(focusMean(sortedInd)));
    colormap(myMap);
    title('-bl')
    focusMean = mean(LCmatChoice(:, midPoints>=currWin(1)&midPoints<currWin(2)), 2);
    subplot(4,6,[12, 18, 24]); hold on;
    imagesc(currWin, [1 length(focusMean)], zscore(focusMean(sortedInd)));
    colormap(myMap);
end

figure(mfig);
subplot(4,3,[4,7,10]); hold on;
[~, sortedInd] = sort(pe);
imagesc(midPointsFilter, 1:size(mPFCmatChoiceFiltered,1), mPFCmatChoiceFiltered(sortedInd,:));
plot([0.001*os.rwdDelay 0.001*os.rwdDelay], [0 size(LCmatChoiceFiltered,1)],'Color', 'w', 'LineWidth', 2, 'LineStyle', '--');
xlim([-1 2.5])
myMap = [[linspace(0, 1, 50)', linspace(0, 1, 50)', linspace(1, 1, 50)'];
[linspace(1, 1, 100)', linspace(1, 0, 100)', linspace(1, 0, 100)']];
colormap(myMap)
colorbar
% maxScale = max(abs(waveformsSession), [], 'all');
caxis([-2 4])


currWin = focusWins{1};
bl = mean(mPFCmatChoice(:, midPoints>=currWin(1)&midPoints<currWin(2)), 2);
for i = 2:length(focusWins)
    currWin = focusWins{i};
    focusMean = mean(mPFCmatChoice(:, midPoints>=currWin(1)&midPoints<currWin(2)), 2) - bl;
    subplot(4,6,[11, 17, 23]); hold on;
    imagesc(currWin, [1 length(focusMean)], zscore(focusMean(sortedInd)));
    colormap(myMap);
    title('-bl')
    focusMean = mean(mPFCmatChoice(:, midPoints>=currWin(1)&midPoints<currWin(2)), 2);
    subplot(4,6,[12, 18, 24]); hold on;
    imagesc(currWin, [1 length(focusMean)], zscore(focusMean(sortedInd)));
    colormap(myMap);
end
colorPSTH = [1 0 0;
             1 0.4 0.4;
             0.4 0.4 1;
             0 0 1];
figure(lfig);
subplot(4,3, [5 8]); hold on;
for j = 1:psthBinNum
    plotFilled(midPointsFilter, LCmatChoiceFiltered(target>=edges(j)&target<edges(j+1),:), colorPSTH(j,:));
end

figure(mfig);
subplot(4,3, [5 8]); hold on;
for j = 1:psthBinNum
    plotFilled(midPointsFilter, mPFCmatChoiceFiltered(target>=edges(j)&target<edges(j+1),:), colorPSTH(j,:));
end

numBins = 4;
meanSignal = zeros(1, numBins);
meanPe = zeros(1, numBins);
semSignal = zeros(1, numBins);
target = t.pe;
edges = linspace(min(target)-0.01, max(target)+0.01, numBins+1);
rwdWin = focusWins{2};
signal = zscore(mean(LCmatChoice(:, midPoints>=rwdWin(1)&midPoints<rwdWin(2)), 2));
for i = 1:numBins
    meanSignal(i) = mean(signal(target>=edges(i) & target<edges(i+1)));
    meanPe(i) = mean(target(target>=edges(i) & target<edges(i+1)));
    semSignal(i) = sem(signal(target>=edges(i) & target<edges(i+1)));
end

figure(lfig)
subplot(4,3,11); hold on;
plot(meanPe, meanSignal, 'LineWidth', 2, 'Color', 'k');
patch([meanPe, flip(meanPe)], [meanSignal-semSignal flip(meanSignal+semSignal)], 'k', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
xlabel('RPE');
ylabel('dF/F zscored')
set(gca, 'TickDir', 'out');
set(gca, 'Box', 'off');
set(gca, 'XTick', -1:0.5:1);
set(gca, 'YTick', -0.5:0.5:1.0);
title('LC');

meanSignal = zeros(1, numBins);
semSignal = zeros(1, numBins);
target = t.pe;
edges = linspace(min(target)-0.01, max(target)+0.01, numBins+1);
signal = zscore(mean(mPFCmatChoice(:, midPoints>=rwdWin(1)&midPoints<rwdWin(2)), 2));
for i = 1:numBins
    meanSignal(i) = mean(signal(target>=edges(i) & target<edges(i+1)));
    semSignal(i) = sem(signal(target>=edges(i) & target<edges(i+1)));
end

figure(mfig)
subplot(4,3,11); hold on;
plot(meanPe, meanSignal, 'LineWidth', 2, 'Color', 'k');
patch([meanPe, flip(meanPe)], [meanSignal-semSignal flip(meanSignal+semSignal)], 'k', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
xlabel('RPE');
ylabel('dF/F zscored')
set(gca, 'TickDir', 'out');
set(gca, 'Box', 'off');
set(gca, 'XTick', -1:0.5:1);
set(gca, 'YTick', -0.5:0.5:1.0);
title('mPFC');

screen = get(0,'Screensize');
screen(4) = screen(4) - 100;
set(mfig, 'Position', screen)
set(lfig, 'Position', screen)
end