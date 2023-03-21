%% load data and parameters
load('F:\tmpData\allSessionUnit.mat');
load('F:\tmpData\waveformFiltered.mat');
tb = 1;
tf = 2;
stepSize = 200;
binSize = 1000;

transIdx = zeros(length(waveformsSession), 1);
for i = 1:length(transIdx)
    x = cellfun(@(x)strcmp(allSessionWF{i}, x), allSessions);
    y = cellfun(@(x)strcmp(allUnitWF{i}, x), allUnits);
    if ~isempty(find(x+y>=2,1))
        transIdx(i) = 1;
    end
end
waveformsSession = waveformsSession(logical(transIdx),:);
%% calculateAUC
allError = {};
for i = 1:length(allSessions)
    session = allSessions{i};
    unit = allUnits{i};
    [~, ~, matChoiceSlide, slideTime] = getUnitMatChoice(session, unit, tb, tf, stepSize, binSize);
    s = behAnalysisNoPlot_opMD(session, 'simpleFlag', 1);
    solenoidTime = 3100 - s.lickLat;
    soleFreeInd = repmat(slideTime+0.5*binSize, length(solenoidTime), 1) <= solenoidTime';
    matChoiceSlide(~soleFreeInd) = NaN; % for removing too late decision trials to clean up solenoid noise

    rwdAUC = zeros(1,size(matChoiceSlide,2));
    aucSig = zeros(1, size(matChoiceSlide,2));
    rwdTStats = zeros(1,size(matChoiceSlide,2)); 
    for j = 1:size(matChoiceSlide,2)
        spikeCounts = matChoiceSlide(:,j);
        outcome = abs(s.allRewards);
        outcome = outcome(~isnan(spikeCounts));
        rightSide = zeros(size(s.allChoices));
        rightSide(s.allChoices>0) = 1;
        rightSide = rightSide(~isnan(spikeCounts));
        spikeCounts = spikeCounts(~isnan(spikeCounts));
        glm = fitglm([spikeCounts], outcome, 'Distribution','binomial','Link','logit');

            
        score = glm.Fitted.Probability;
        if glm.Coefficients.tStat(2)>0
            [Xglm,Yglm,~,AUCglm] = perfcurve(outcome,score,1, 'NBoot',500, 'Alpha', 0.05);
            aucSig(j) = AUCglm(2)>0.5;
        else
            [Xglm,Yglm,~,AUCglm] = perfcurve(outcome,score,0, 'NBoot',500, 'Alpha', 0.05);
            aucSig(j) = AUCglm(3)<0.5;
        end
        
    %     if (AUCglm(2)-0.5)*(AUCglm(3)-0.5)>0
    %         AUCglmSig = 1;
    %     else
    %         AUCglmSig = 0;
    %     end
    %     if glm.Coefficients.Estimate(2) < 0
    %         AUCglm = 1-AUCglm;
    %     end
        rwdAUC(j) = AUCglm(1);
        rwdTStats(j) = glm.Coefficients.tStat(2); 
        allAUCZS(i,j) = auROCZS(spikeCounts(outcome==0), spikeCounts(outcome==1));
    end
    allAUC(i,:) = rwdAUC;
    allCoeff(i,:) = rwdTStats;
    allAUCsig(i,:) = aucSig;
end
%% clustering by waveform
color1 = [0.2 0.2 0.2];
color2 = [0.8 0.8 0.8];
% filtering
% fs = 32000;
% fc = [300 8000];
% [b,a] = butter(2,fc/(fs/2),'bandpass');
% 
% waveformsSession = filtfilt(b,a,waveformsSession')';
%% clustering
numCat = 2;
[coeff,score,latent, ~, explained, mu] = pca(zscore(waveformsSession, [], 1));
indAll = {};
dis = {};
for a = 1:10
    [indAll{a}, ~, dis{a}] = kmeans([score(:, 1:5)], numCat);
end
[~,optiInds] = min(cellfun(@mean, dis));
ind = indAll{optiInds};
figure2;
hold on;
plotFilled(1:size(waveformsSession,2), waveformsSession(ind==1,:), color1);
plotFilled(1:size(waveformsSession,2), waveformsSession(ind==2,:), color2);

plot([5 21], [-0.3 -0.3], 'LineWidth', 3,'Color', 'k');
text(10, -0.4, '0.5 ms', 'FontSize', 14)

set(gca, 'Box', 'off')
set(gca,'tickdir', 'out')
xlabel('expected value tstats', 'FontSize', 18)
set(gca, 'XTick', [-10:5:5], 'FontSize', 14)
set(gca, 'YTick', [0 0.1 0.2], 'FontSize', 14)
set(gca, 'XColor', 'none')
set(gca, 'YColor', 'none')
%% 
figure2;
maxScale = max(abs(waveformsSession), [], 'all');
myMap = [[linspace(0, 1, 100)', linspace(0, 1, 100)', linspace(1, 1, 100)'];
         [linspace(1, 1, 100)', linspace(1, 0, 100)', linspace(1, 0, 100)']];
     
subplot(1,2,1);
imagesc(1:size(waveformsSession,2), 1:sum(ind==1), waveformsSession(ind==1,:)); 
ylim([0 max([sum(ind==1) sum(ind==2)])])
colormap(myMap)
caxis([-maxScale maxScale])
title('Type II', 'FontSize', 18)
set(gca, 'Box', 'off')
set(gca,'tickdir', 'out')
xlabel('expected value tstats', 'FontSize', 18)
set(gca, 'XTick', [-10:5:5], 'FontSize', 14)
set(gca, 'YTick', [0 0.1 0.2], 'FontSize', 14)
set(gca, 'XColor', 'none')
set(gca, 'YColor', 'none')
colorbar; 
colorbar('Ticks', [-1.5 0 1.5])

subplot(1,2,2);
imagesc(1:size(waveformsSession,2), 1:sum(ind==2), waveformsSession(ind==2,:)); 
ylim([0 max([sum(ind==1) sum(ind==2)])])
colormap(myMap)
caxis([-maxScale maxScale])
title('Type I', 'FontSize', 18)
colorbar;
set(gca, 'Box', 'off')
set(gca,'tickdir', 'out')
xlabel('expected value tstats', 'FontSize', 18)
set(gca, 'XTick', [-10:5:5], 'FontSize', 14)
set(gca, 'YTick', [0 0.1 0.2], 'FontSize', 14)
set(gca, 'XColor', 'none')
set(gca, 'YColor', 'none')
colorbar('Ticks', [-1.5 0 1.5])
%% waveform 
figure2;
[~, sortInd] = sort(score(:, 1));
waveformsSessionSorted = waveformsSession(sortInd,:);
indSorted = ind(sortInd);
subplot(1,2,1); hold on;

plot(1:size(waveformsSession,2), waveformsSessionSorted(indSorted==1,:)+0.5*[1:sum(ind==1)]', 'k');
plot([10 26], [-5 -5], 'k', 'LineWidth', 2);
text(18, -7, '0.5 ms', 'HorizontalAlignment', 'center');
set(gca, 'XColor', 'none');
set(gca, 'YColor', 'none');
ylim([-10 0.5*sum(ind==1)+3]);
xlim([10 size(waveformsSession,2)-10])
title('Type II')

subplot(1,2,2);hold on;
plot(1:size(waveformsSession,2), waveformsSessionSorted(indSorted==2,:)+0.5*[1:sum(ind==2)]', 'k');
plot([10 26], [-5 -5], 'k', 'LineWidth', 2);
text(18, -7, '0.5 ms', 'HorizontalAlignment', 'center');
ylim([-10 0.5*sum(ind==1)+3]);
xlim([10 size(waveformsSession,2)-10])
title('Type I')
set(gca, 'XColor', 'none');
set(gca, 'YColor', 'none');
%%
% colorbar;
% caxis([-maxMap maxMap])
% set(gca, 'Box', 'off')
% set(gca, 'TickDir', 'Out')
% set(gca, 'XTick', [-500 0 500], 'FontSize', 12)
% set(gca, 'YColor', 'none')
% xlabel('Time from lick (ms)', 'FontSize', 15)
% set(gca, 'Box', 'off')
% ylim([-2 yH]); xlim([-500 1000])
% colorbar('Ticks', [-floor(maxMap) ceil(minMap) 0 floor(maxMap)], 'FontSize', 12)
%%
figure2;
myMap = [[linspace(0, 1, 100)', linspace(0, 1, 100)', linspace(1, 1, 100)'];
         [linspace(1, 1, 100)', linspace(1, 0, 100)', linspace(1, 0, 100)']];
subplot(1,4,1);
imagesc(slideTime, 1:size(allAUC, 1), allAUC)
title('glmauc')
colormap(myMap)
caxis([0 1])

subplot(1,4,2);
imagesc(slideTime, 1:size(allAUC, 1), allAUCZS)
title('auc')
colormap(myMap)
caxis([0 1])

subplot(1,4,3);
imagesc(slideTime, 1:size(allAUC, 1), allCoeff)
maxScale = max(abs(allCoeff), [], 'all');
title('glm')
colormap(myMap)
caxis([-maxScale maxScale])

subplot(1,4,4);
imagesc(slideTime, 1:size(allCoeff, 1), allCoeff(:,1))
title('lm')
maxScale = max(abs(allCoeff(:,1)), [], 'all');
colormap(myMap)
caxis([-maxScale maxScale])

%% sorting by maxWin
slideTimePostLick = slideTime(slideTime >= 0.5*binSize + s.rwdDelay & slideTime <= 2000 - 0.5*binSize);
allAUCPostLick = allAUC(:, slideTime >= 0.5*binSize + s.rwdDelay & slideTime <= 2000 - 0.5*binSize);
[~, maxBin] = max(abs(allAUCPostLick - 0.5),[],2);
loc = sub2ind(size(allAUCPostLick),1:length(allSessions),maxBin');
aucMax = allAUCPostLick(loc)'; 
aucMaxSig = allAUCsig(:,slideTime >= 0.5*binSize + s.rwdDelay & slideTime <= 2000 - 0.5*binSize);
aucMaxSig = aucMaxSig(loc)';
myMap = [[linspace(0, 1, 100)', linspace(0, 1, 100)', linspace(1, 1, 100)'];
         [linspace(1, 1, 100)', linspace(1, 0, 100)', linspace(1, 0, 100)']];
figure2;
subplot(1,2,1);
imagesc(aucMax(ind==1));
title('auc')
ylim([1 max([sum(ind==1) sum(ind==2)])])
colormap(myMap)
caxis([0 1])

subplot(1,2,2);
imagesc(aucMax(ind==2));
title('auc')
ylim([1 max([sum(ind==1) sum(ind==2)])])
colormap(myMap)
caxis([0 1])
%% compare coeff with maxAUC
figure2;

subplot(1,2,1);
imagesc(slideTime, 1:size(allAUC, 1), allCoeffs(:,1))
title('lm')
maxScale = max(abs(allCoeffs(:,1)), [], 'all');
colormap(myMap)
caxis([-maxScale maxScale])

subplot(1,2,2);
imagesc(slideTime, 1:size(allAUC, 1), aucMax)
title('maxAUC')
colormap(myMap)
caxis([0 1])
%% aucMax
figure2;
hold on;
bins = linspace(min(aucMax)-0.001, max(aucMax)+0.001, 30);
histogram(aucMax(ind==1), bins, 'FaceColor', color1, 'EdgeColor', 'none');
histogram(aucMax(ind==2), bins, 'FaceColor', color2, 'EdgeColor', 'none');
auc = auROCZS(aucMax(ind==2),aucMax(ind==1));
title('max auc', 'FontSize', 18)
set(gca,'tickdir', 'out')
set(gca, 'XTick', [0:0.25:1], 'FontSize', 14)
set(gca, 'YTick', [0 0.1 0.2], 'FontSize', 14)
set(gca,'tickdir', 'out')
set(gca, 'XTick', 0:0.5:1, 'FontSize', 14)
set(gca, 'YTick', 0:5:10, 'FontSize', 14)
xlim([0 1])
%% aucFirst
aucFirst = allAUC(:,slideTime>=300+0.5*binSize & slideTime<300+stepSize+0.5*binSize);

figure2;
hold on;
bins = linspace(min(aucMax)-0.001, max(aucMax)+0.001, 30);
histogram(aucFirst(ind==1), bins, 'FaceColor', color1,'EdgeColor', 'none');
histogram(aucFirst(ind==2), bins, 'FaceColor', color2, 'EdgeColor', 'none');
title('first auc', 'FontSize', 18)
set(gca,'tickdir', 'out')
set(gca, 'XTick', 0:0.5:1, 'FontSize', 14)
set(gca, 'YTick', 0:5:10, 'FontSize', 14)
xlim([0 1])
%% Tstats

figure2;
hold on;
bins = linspace(min(allTstats(:,1))-0.001, max(allTstats(:,1))+0.001, 30);
histogram(allTstats((ind==1),1), bins, 'FaceColor', color1, 'Normalization', 'probability', 'EdgeColor', 'none');
histogram(allTstats((ind==2),1), bins, 'FaceColor', color2, 'Normalization', 'probability', 'EdgeColor', 'none');
title('Tstats', 'FontSize', 18)
set(gca,'tickdir', 'out')
% set(gca, 'XTick', [0:0.25:1], 'FontSize', 14)
% set(gca, 'YTick', [0 0.1 0.2], 'FontSize', 14)
% xlim([0 1])


%% colormap
figure2;
myMap = [[linspace(0, 1, 100)', linspace(0, 1, 100)', linspace(1, 1, 100)'];
         [linspace(1, 1, 100)', linspace(1, 0, 100)', linspace(1, 0, 100)']];

subplot(1,2,1);
imagesc(slideTime, 1:sum(ind==1), allAUC(ind==1,:));
title('auc')
ylim([1 max([sum(ind==1) sum(ind==2)])])
colormap(myMap)
caxis([0 1])
title('Type II', 'FontSize', 18)
set(gca, 'Box', 'off')
set(gca,'tickdir', 'out')
xlabel('time from choice (ms)', 'FontSize', 18)
set(gca, 'XTick', [-500:500:1500])
set(gca, 'XColor', 'k')
set(gca, 'YColor', 'none')
colorbar; 
colorbar('Ticks', [0 0.5 1])

subplot(1,2,2);
imagesc(slideTime, 1:sum(ind==2), allAUC(ind==2,:));
title('auc')
ylim([1 max([sum(ind==1) sum(ind==2)])])
colormap(myMap)
caxis([0 1])
title('Type I', 'FontSize', 18)
set(gca, 'Box', 'off')
set(gca,'tickdir', 'out')
xlabel('time from choice (ms)', 'FontSize', 18)
set(gca, 'YColor', 'none')
set(gca, 'XTick', [-500:500:1500])
set(gca, 'XColor', 'k')
colorbar; 
colorbar('Ticks', [0 0.5 1])
%% compare glm and auc
figure2; hold on;
scatter(allTstats(:,1), aucFirst)
scatter(allTstats(allSigs(:,1)<0.05,1), aucFirst(allSigs(:,1)<0.05))
plot([0 0], [0 1],'Color','r', 'LineStyle', '--')
hold on; plot([0 0], [0 1],'Color','r', 'LineStyle', '--')
hold on; plot([-15 20], [0.5 0.5],'Color','r', 'LineStyle', '--')
title('First AUC vs First Tstats')
xlabel('Tstats'); ylabel('AUC')

figure2; hold on;
scatter(allTstats(:,1), aucMax)
scatter(allTstats(allSigs(:,1)<0.05,1), aucMax(allSigs(:,1)<0.05))
plot([0 0], [0 1],'Color','r', 'LineStyle', '--')
hold on; plot([0 0], [0 1],'Color','r', 'LineStyle', '--')
hold on; plot([-15 20], [0.5 0.5],'Color','r', 'LineStyle', '--')
title('Max AUC vs First Tstats')
xlabel('Tstats'); ylabel('AUC')
%%





%% use the new window to do new LM fitting
[root, sep] = currComputer();
col = 'good';
modelName = '5params';

numBins = 8;  
numBinsPSTH = 6;
numBinsSW = 4; 
tbPSTH = 1;
tfPSTH = 2.5;
binSizePSTH = 200;
stepSizePSTH = 100;

allPe = zeros(length(allUnits), numBins);
allSpikes = zeros(length(allUnits), numBins);
allPSTH = [];
allSW = zeros(length(allUnits), numBinsSW);
spikeFocus = cell(length(allSessions),1);
for i = 1:length(allSessions)
    session = allSessions{i};
    unit = allUnits{i};
    s = behAnalysisNoPlot_opMD(session, 'simpleFlag', 1);
    pd = parseSessionString_df(session, root, sep);
    sampFile = [pd.animalName col '_' modelName];
    path = [root pd.animalName sep pd.animalName 'sorted' sep 'stan' sep 'bernoulli' sep modelName sep col sep];
    [t,~,noSession] = getStanModelParams_samps(modelName, [path sampFile '.mat'], 4000, 'sessionName', session);
    
    Qchosen = NaN(size(s.allChoices));
    choice = s.allChoices;
    for j = 1:length(choice)
        if choice(j)>0
            Qchosen(j) = t.Q(j,2);
        else
            Qchosen(j) = t.Q(j,1);
        end
    end
    
    rightSide = 0.5*(s.allChoices + 3);
    
    [~, ~, matChoiceSlide, slideTime] = getUnitMatChoice(session, unit, tb, tf, stepSize, binSize);
    solenoidTime = 3100 - s.lickLat;
    soleFreeInd = repmat(slideTime+0.5*binSize, length(solenoidTime), 1) <= solenoidTime';
    matChoiceSlide(~soleFreeInd) = NaN; % for removing too late decision trials to clean up solenoid noise
    
    spikeCountsFocus = matChoiceSlide(:, sum(slideTime < 0.5*binSize + s.rwdDelay) + maxBin(i));
    spikeFocus{i} = spikeCountsFocus;
    outcome = abs(s.allRewards);
    outcome = outcome(~isnan(spikeCountsFocus));
    rightSide = rightSide(~isnan(spikeCountsFocus));
    Qchosen = Qchosen(~isnan(spikeCountsFocus));
    spikeCountsFocus = spikeCountsFocus(~isnan(spikeCountsFocus));
    
    auc(i) = auROCZS(spikeCountsFocus(outcome==0), spikeCountsFocus(outcome==1));
    lm = fitlm([outcome' Qchosen' rightSide'], spikeCountsFocus);
    
    lmMaxAll(i).coeffs = lm.Coefficients.Estimate(2:end);
    lmMaxAll(i).ps = lm.Coefficients.pValue(2:end);
    lmMaxAll(i).tStats = lm.Coefficients.tStat(2:end);
    % tuning curve and psth
    spikeCountsFocus = matChoiceSlide(:, sum(slideTime < 0.5*binSize + s.rwdDelay) + maxBin(i));
    spikeCountsFocus(~isnan(spikeCountsFocus)) = zscore(spikeCountsFocus(~isnan(spikeCountsFocus)));
    [~, ~, matChoiceSlide, slideTime] = getUnitMatChoice(session, unit, tbPSTH, tfPSTH, stepSizePSTH, binSizePSTH);
    matChoiceSlide = zscore(matChoiceSlide, 0, 'all');
    target = t.pe;
    tempPe = NaN(1,numBins);
    tempSpike = NaN(1,numBins); 
    tempPSTH = NaN(numBinsPSTH, length(slideTime));
    tempSW = NaN(1,numBinsSW);
    % tuning
    edges = [linspace(min(target)- 0.01, 0, 0.5*numBins+1) linspace(0, max(target)+ 0.01,0.5*numBins+1)];
    edges = [edges(1:0.5*numBins+1), edges(0.5*numBins+3:end)];
    for j = 1:numBins
        if ~isempty(find(target >= edges(j) & target < edges(j+1), 1))
            tempPe(j) = mean(target(target >= edges(j) & target < edges(j+1)));
            tempSpike(j) = mean(spikeCountsFocus(target >= edges(j) & target < edges(j+1)));
        end
    end
    allPe(i,:) = tempPe;
    allSpikes(i,:) = tempSpike;
   
    % PSTH
    edges = [linspace(min(target)- 0.01, 0, 0.5*numBinsPSTH+1) linspace(0, max(target)+ 0.01,0.5*numBinsPSTH+1)];
    edges = [edges(1:0.5*numBinsPSTH+1), edges(0.5*numBinsPSTH+3:end)];
    for j = 1:numBinsPSTH
        if ~isempty(find(target >= edges(j) & target < edges(j+1), 1))
            tempPSTH(j,:) = mean(matChoiceSlide(target >= edges(j) & target < edges(j+1),:));
        end
    end
    allPSTH = cat(3, allPSTH, tempPSTH);
    % pSwith
    target = spikeCountsFocus;
    edges = linspace(min(target)-0.01, max(target)+0.01, numBinsSW+1);
    svsNext = NaN(1,length(s.allChoices));
    svsNext(s.changeChoice_Inds-1) = 1;
    svsNext(s.stayChoice_Inds-1) = 0;
    for j = 1:numBinsSW
        tempSW(j) = mean(svsNext(target>=edges(j)&target<edges(j+1)), 'omitnan');
    end
    tempSW(~isnan(tempSW)) = zscore(tempSW(~isnan(tempSW)));
    allSW(i,:) = tempSW;
    
end
%%
coeffsMax = [lmMaxAll.coeffs]';
tStatsMax = [lmMaxAll.tStats]';
sigMax = [lmMaxAll.ps]'<0.05;
figure2;
hold on;
edges = linspace(min(tStatsMax(:,1))-0.01, max(tStatsMax(:,1))+0.01, 30);
histogram(tStatsMax((ind==1),1), edges, 'FaceColor', color1, 'EdgeColor', 'none');
histogram(tStatsMax((ind==2),1), edges, 'FaceColor', color2, 'EdgeColor', 'none');
[h, p] = ttest2(tStatsMax((ind==1),1), tStatsMax((ind==2),1));

title(['max Tstats Outcome' num2str(p)]);

set(gca, 'Box', 'off')
set(gca,'tickdir', 'out')
set(gca, 'XTick', [-15:5:15])
set(gca, 'YTick', [0:5:15])
set(gca, 'XColor', 'k')

figure2;
hold on;
edges = linspace(min(tStatsMax(:,2))-0.01, max(tStatsMax(:,2))+0.01, 30);
histogram(tStatsMax((ind==1),2), edges, 'FaceColor', color1, 'EdgeColor', 'none');
histogram(tStatsMax((ind==2),2), edges, 'FaceColor', color2, 'EdgeColor', 'none');
[h, p] = ttest2(tStatsMax((ind==1),2), tStatsMax((ind==2),2));
title(['max Tstats Qchosen' num2str(p)]);

set(gca, 'Box', 'off')
set(gca,'tickdir', 'out')
set(gca, 'XTick', [-15:5:15])
set(gca, 'YTick', [0:5:15])
set(gca, 'XColor', 'k')

%% see improvement on consistent and inconsistent results
% inconsistent ones
diffInd = sign(tStatsMax(:,1)) ~= sign(allTstats(:,1));
edges = linspace(min(allTstats(diffInd,1)), max(allTstats(diffInd,1)), 20);
figure2; subplot(2,1,1); hold on; 
histogram(allTstats(ind==1 & diffInd,1), edges, 'FaceColor', color1, 'Normalization', 'probability', 'EdgeColor', 'none')
histogram(allTstats(ind==2 & diffInd,1), edges, 'FaceColor', color2, 'Normalization', 'probability', 'EdgeColor', 'none')
auc = auROCZS(allTstats(ind==1 & diffInd,1),allTstats(ind==2 & diffInd,1));
title(['inconsistent on first window ' num2str(auc)])
subplot(2,1,2); hold on;
edges = linspace(min(tStatsMax(diffInd,1)), max(tStatsMax(diffInd,1)), 20);
histogram(tStatsMax(ind==1 & diffInd,1), edges, 'FaceColor', color1, 'Normalization', 'probability', 'EdgeColor', 'none')
histogram(tStatsMax(ind==2 & diffInd,1), edges, 'FaceColor', color2, 'Normalization', 'probability', 'EdgeColor', 'none')
auc = auROCZS(tStatsMax(ind==1 & diffInd,1),tStatsMax(ind==2 & diffInd,1));
title(['inconsistent on max window ' num2str(auc)])
%%
% consistent ones
diffInd = sign(tStatsMax(:,1)) == sign(allTstats(:,1));
edges = linspace(min(allTstats(diffInd,1)), max(allTstats(diffInd,1)), 20);
figure2; subplot(2,1,1); hold on; 
histogram(allTstats(ind==1 & diffInd,1), edges, 'FaceColor', color1, 'Normalization', 'probability', 'EdgeColor', 'none')
histogram(allTstats(ind==2 & diffInd,1), edges, 'FaceColor', color2, 'Normalization', 'probability', 'EdgeColor', 'none')
auc = auROCZS(allTstats(ind==1 & diffInd,1),allTstats(ind==2 & diffInd,1));
title(['consistent on first window ' num2str(auc)])
subplot(2,1,2); hold on;
edges = linspace(min(tStatsMax(diffInd,1)), max(tStatsMax(diffInd,1)), 20);
histogram(tStatsMax(ind==1 & diffInd,1), edges, 'FaceColor', color1, 'Normalization', 'probability', 'EdgeColor', 'none')
histogram(tStatsMax(ind==2 & diffInd,1), edges, 'FaceColor', color2, 'Normalization', 'probability', 'EdgeColor', 'none')
auc = auROCZS(tStatsMax(ind==1 & diffInd,1),tStatsMax(ind==2 & diffInd,1));
title(['consistent on max window ' num2str(auc)])
%%
%% scatter in tstats space
% maxWin
figure2; hold on;
outcomeInd = 1;
qInd = 2;
scatter(tStatsMax(ind==1, outcomeInd), tStatsMax(ind==1,qInd), 50, [0 0 0], 'o', 'LineWidth', 2.5);
scatter(tStatsMax(ind==2, outcomeInd), tStatsMax(ind==2,qInd), 50, [0.6 0.6 0.6],  'o', 'LineWidth', 2.5);
plot([-15 15], [0 0], 'LineStyle', '--', 'LineWidth', 2, 'Color', [0.6 0.6 0.6]);
plot([0 0], [-10 5], 'LineStyle', '--', 'LineWidth', 2, 'Color', [0.6 0.6 0.6]);
xlim([-18 20])
ylim([-8 6])
set(gca, 'Box', 'off')
set(gca,'tickdir', 'out')
set(gca, 'XTick', [-15:5:20])
set(gca, 'YTick', [-10:5:5])
set(gca, 'XColor', 'k')
xlabel('outcome', 'FontSize', 18)
ylabel('Qchosen', 'FontSize', 18)
% firstWin
figure2; hold on;
outcomeInd = 1;
qInd = 3;
scatter(allTstats(ind==1, outcomeInd), allTstats(ind==1,qInd), 50, [0 0 0], 'o', 'LineWidth', 2.5);
scatter(allTstats(ind==2, outcomeInd), allTstats(ind==2,qInd), 50, [0.6 0.6 0.6],  'o', 'LineWidth', 2.5);
plot([-15 15], [0 0], 'LineStyle', '--', 'LineWidth', 2, 'Color', [0.6 0.6 0.6]);
plot([0 0], [-10 5], 'LineStyle', '--', 'LineWidth', 2, 'Color', [0.6 0.6 0.6]);
xlim([-18 20])
ylim([-8 6])
set(gca, 'Box', 'off')
set(gca,'tickdir', 'out')
set(gca, 'XTick', [-15:5:20])
set(gca, 'YTick', [-10:5:5])
set(gca, 'XColor', 'k')
xlabel('outcome', 'FontSize', 18)
ylabel('Qchosen', 'FontSize', 18)
%% scatter in coeff space
figure2; hold on;
% scatter(allTstats(:, outcomeInd), allTstats(:,qInd), 30, [0.7 0.7 0.7], 'filled');
% scatter(allCoeffs(cats==1, outcomeInd), allCoeffs(cats==1,qInd), 50, [0 0 0], 'x', 'LineWidth', 2);
% scatter(allCoeffs(cats==2, outcomeInd), allCoeffs(cats==2,qInd), 50, [0.6 0.6 0.6],  'x', 'LineWidth', 2);
% scatter(allCoeffs(ind==-1, outcomeInd), allCoeffs(ind==-1,qInd), 50, color1, '.', 'LineWidth', 2);
% scatter(allCoeffs(ind==1, outcomeInd), allCoeffs(ind==1,qInd), 50, color2,  '.', 'LineWidth', 2);
% xlim([-1.5 1.5])
% ylim([-1.5 1.5])
%% polar histogram all
edges = linspace(-pi, pi, 20);
% max tStats
outcomeInd = 1;
qInd = 2;
sigO = sigMax(:, outcomeInd);
sigQ = sigMax(:, qInd);

allVec = [tStatsMax(:, outcomeInd), tStatsMax(:,qInd)];
allVec = [coeffsMax(:, outcomeInd), coeffsMax(:,qInd)];
[theta1, rho] = cart2pol(allVec(ind==1,1), allVec(ind==1,2));
[theta2, rho] = cart2pol(allVec(ind==2,1), allVec(ind==2,2));
figure2; 
polarhistogram(theta1,edges, 'FaceColor', [0.1 0.1 0.1], 'FaceAlpha',.7, 'EdgeColor', 'none', 'Normalization', 'Probability');
hold on; 
polarhistogram(theta2,edges, 'FaceColor', [0.7 0.7 0.7], 'FaceAlpha',.6, 'EdgeColor', 'none', 'Normalization', 'Probability');
title('maxWin, sig')


% frst tStats
outcomeInd = 1;
qInd = 3;
allVec = [allTstats(:, outcomeInd), allTstats(:,qInd)];
[theta1, rho] = cart2pol(allVec(ind==1,1), allVec(ind==1,2));
[theta2, rho] = cart2pol(allVec(ind==2,1), allVec(ind==2,2));
figure2; 
polarhistogram(theta1,edges, 'FaceColor', [0.1 0.1 0.1], 'FaceAlpha',.7, 'EdgeColor', 'none', 'Normalization', 'Probability');
hold on; 
polarhistogram(theta2,edges, 'FaceColor', [0.7 0.7 0.7], 'FaceAlpha',.6, 'EdgeColor', 'none', 'Normalization', 'Probability');
title('firstWin, all')
%% polar histogram sig ones
edges = linspace(-pi, pi, 20);
% max tStats
outcomeInd = 1;
qInd = 2;
sigO = sigMax(:, outcomeInd);
sigQ = sigMax(:, qInd);

allVec = [coeffsMax(:, outcomeInd), coeffsMax(:,qInd)];
[theta1, rho] = cart2pol(allVec((sigO)&ind==1,1), allVec((sigO)&ind==1,2));
[theta2, rho] = cart2pol(allVec((sigO)&ind==2,1), allVec((sigO)&ind==2,2));
figure2; 
polarhistogram(theta1,edges, 'FaceColor', [0.1 0.1 0.1], 'FaceAlpha',.7, 'EdgeColor', 'none', 'Normalization', 'Probability');
hold on; 
polarhistogram(theta2,edges, 'FaceColor', [0.7 0.7 0.7], 'FaceAlpha',.6, 'EdgeColor', 'none', 'Normalization', 'Probability');
title('maxWin, sig')
% frst tStats
outcomeInd = 1;
qInd = 3;
sigO = sigMax(:, outcomeInd);
sigQ = sigMax(:, qInd);

allVec = [allTstats(:, outcomeInd), allTstats(:,qInd)];
[theta1, rho] = cart2pol(allVec((sigO|sigQ)&ind==1,1), allVec((sigO|sigQ)&ind==1,2));
[theta2, rho] = cart2pol(allVec((sigO|sigQ)&ind==2,1), allVec((sigO|sigQ)&ind==2,2));
figure2; 
polarhistogram(theta1,edges, 'FaceColor', [0.1 0.1 0.1], 'FaceAlpha',.7, 'EdgeColor', 'none', 'Normalization', 'Probability');
hold on; 
polarhistogram(theta2,edges, 'FaceColor', [0.7 0.7 0.7], 'FaceAlpha',.6, 'EdgeColor', 'none', 'Normalization', 'Probability');
title('firstWin, sig')
%% unwrap
% maxWin
edges = linspace(0, 2*pi, 20);

outcomeInd = 1;
qInd = 2;
allVec = [tStatsMax(:, outcomeInd), tStatsMax(:,qInd)];
[theta1, rho] = cart2pol(allVec(ind==1,1), allVec(ind==1,2));
[theta2, rho] = cart2pol(allVec(ind==2,1), allVec(ind==2,2));
theta2Rotate = theta2;
theta2Rotate(theta2Rotate>1) = theta2Rotate(theta2Rotate>1) - 1.5;
theta2Rotate(theta2Rotate<1) = theta2Rotate(theta2Rotate<1) + 2*pi - 1.5;
theta1Rotate = theta1;
theta1Rotate(theta1Rotate>1) = theta1Rotate(theta1Rotate>1) - 1;
theta1Rotate(theta1Rotate<1) = theta1Rotate(theta1Rotate<1) + 2*pi - 1.5;

figure2; 
histogram(theta1Rotate,edges, 'FaceColor', [0.1 0.1 0.1], 'FaceAlpha',.7, 'EdgeColor', 'none', 'Normalization', 'Probability');
hold on; 
histogram(theta2Rotate,edges, 'FaceColor', [0.7 0.7 0.7], 'FaceAlpha',.6, 'EdgeColor', 'none', 'Normalization', 'Probability');
plot([pi-1 pi-1], [0 0.5], 'Color', 'k', 'LineStyle', '--')
plot([2*pi-1-pi/4 2*pi-1-pi/4], [0 0.5], 'Color', 'r', 'LineStyle', '--')

title('maxWin, all')
% first win
outcomeInd = 1;
qInd = 3;
allVec = [allTstats(:, outcomeInd), allTstats(:,qInd)];
[theta1, rho] = cart2pol(allVec(ind==1,1), allVec(ind==1,2));
[theta2, rho] = cart2pol(allVec(ind==2,1), allVec(ind==2,2));
theta2Rotate = theta2;
theta2Rotate(theta2Rotate>1) = theta2Rotate(theta2Rotate>1) - 1.5;
theta2Rotate(theta2Rotate<1) = theta2Rotate(theta2Rotate<1) + 2*pi - 1.5;
theta1Rotate = theta1;
theta1Rotate(theta1Rotate>1) = theta1Rotate(theta1Rotate>1) - 1;
theta1Rotate(theta1Rotate<1) = theta1Rotate(theta1Rotate<1) + 2*pi - 1.5;
figure2; 
histogram(theta1Rotate,edges, 'FaceColor', [0.1 0.1 0.1], 'FaceAlpha',.7, 'EdgeColor', 'none', 'Normalization', 'Probability');
hold on; 
histogram(theta2Rotate,edges, 'FaceColor', [0.7 0.7 0.7], 'FaceAlpha',.6, 'EdgeColor', 'none', 'Normalization', 'Probability');
plot([pi-1 pi-1], [0 0.5], 'Color', 'k', 'LineStyle', '--')
plot([2*pi-1-pi/4 2*pi-1-pi/4], [0 0.5], 'Color', 'r', 'LineStyle', '--')
title('firstWin, all')
%% populational tuning curve
figure2;
hold on;
plotFilled(mean(allPe(sigMax(:,1) & tStatsMax(:,1)>0 & ind==2,:), 'omitnan'), allSpikes(sigMax(:,1) & tStatsMax(:,1)>0 & ind==2,:), color2);
plotFilled(mean(allPe(sigMax(:,1) & tStatsMax(:,1)<0 & ind==1,:), 'omitnan'), allSpikes(sigMax(:,1) & tStatsMax(:,1)<0 & ind==1,:), color1);

set(gca, 'Box', 'off')
set(gca,'tickdir', 'out')
set(gca, 'XTick', [-1 0 1])
set(gca, 'YTick', [-0.5 0 0.5])
xlabel('rpe', 'FontSize', 18)
ylabel('spikes/s (zscored)', 'FontSize', 18)
%% pouplational P(SW|spikes)
figure2;
hold on;
plotFilled(1:numBinsSW, allSW(sigMax(:,1) & tStatsMax(:,1)>0 & ind==2,:), color2);
plotFilled(1:numBinsSW, allSW(sigMax(:,1) & tStatsMax(:,1)<0 & ind==1,:), color1);

set(gca, 'Box', 'off')
set(gca,'tickdir', 'out')
set(gca, 'XTick', [1:4])
set(gca, 'YTick', [-0.5 0 0.5])
xlabel('quantile', 'FontSize', 18)
ylabel('P(switch)', 'FontSize', 18)
%% PSTH

colors = [1 0 0;
          1 0.3 0.3;
          1 0.6 0.6;
%           0.6 0.6 1;
%           1 0.6 0.6;
          0.6 0.6 1;
          0.3 0.3 1;
          0 0 1];
      
figure2;
subplot(1,2,1);
for i = 1:numBinsPSTH
    plotFilled(slideTime, squeeze(allPSTH(i,:,sigMax(:,1) & tStatsMax(:,1)<0 & ind==1))', colors(i,:));
end
set(gca, 'Box', 'off')
set(gca,'tickdir', 'out')
set(gca, 'XTick', [-1000:1000:3000])
set(gca, 'YTick', [-0.4:0.4:1.2])
xlabel('time from choice (ms)', 'FontSize', 18)
ylabel('spikes/s (zscored)', 'FontSize', 18)

title('Type II', 'FontSize', 18)

subplot(1,2,2);
for i = 1:numBinsPSTH
    plotFilled(slideTime, squeeze(allPSTH(i,:,sigMax(:,1) & tStatsMax(:,1)>0 & ind==2))', colors(i,:));
end
set(gca, 'Box', 'off')
set(gca,'tickdir', 'out')
set(gca, 'XTick', [-1000:1000:3000])
set(gca, 'YTick', [-0.4:0.4:1.2])
xlabel('time from choice (ms)', 'FontSize', 18)
ylabel('spikes/s (zscored)', 'FontSize', 18)

title('Type I', 'FontSize', 18)
%%
colors = [1 0 0;
%           1 0.3 0.3;
          1 0.6 0.6;
%           0.6 0.6 1;
%           1 0.6 0.6;
          0.6 0.6 1;
%           0.3 0.3 1;
          0 0 1];
      
[root, sep] = currComputer();
col = 'good';
modelName = '5params';

numBinsPSTH = 4;
tbPSTH = 1;
tfPSTH = 2.5;
binSizePSTH = 200;
stepSizePSTH = 100;

allPe = zeros(length(allUnits), numBins);
allPSTH = [];
tempInds = find(sigMax(:,1)>0 & tStatsMax(:,1)>0 & ind==2);
for i = 1:length(tempInds)
    i = 25;
    session = allSessions{tempInds(i)};
    unit = allUnits{tempInds(i)};
    s = behAnalysisNoPlot_opMD(session, 'simpleFlag', 1);
    pd = parseSessionString_df(session, root, sep);
    sampFile = [pd.animalName col '_' modelName];
    path = [root pd.animalName sep pd.animalName 'sorted' sep 'stan' sep 'bernoulli' sep modelName sep col sep];
    [t,~,noSession] = getStanModelParams_samps(modelName, [path sampFile '.mat'], 4000, 'sessionName', session);
    
    % psth
    [cellChoice, ~, matChoiceSlide, slideTime] = getUnitMatChoice(session, unit, tbPSTH, tfPSTH, stepSizePSTH, binSizePSTH);
%     matChoiceSlide = zscore(matChoiceSlide, 0, 'all');
    target = t.pe;
    tempPSTH =[];
    edges = [linspace(min(target)- 0.01, 0, 0.5*numBinsPSTH+1) linspace(0, max(target)+ 0.01,0.5*numBinsPSTH+1)];
    edges = [edges(1:0.5*numBinsPSTH+1), edges(0.5*numBinsPSTH+3:end)];
    figure2; hold on;
    subplot(1,3,[1 2])
    for j = 1:numBinsPSTH
        if ~isempty(find(target >= edges(j) & target < edges(j+1), 1))
            tempPSTH = matChoiceSlide(target >= edges(j) & target < edges(j+1),:);
            plotFilled(slideTime, tempPSTH, colors(j,:));
        end
    end
    subplot(1,3,3);  hold on;
    [~, sortInd] = sort(t.pe);
    cellChoice(cellfun(@isempty,cellChoice)) = {zeros(1,0)}; 
    plotSpikeRaster(cellChoice(sortInd),'PlotType','vertline');
    hold on;
    plot([-1000 2000], [length(s.nrwd_Inds) length(s.nrwd_Inds)], 'LineStyle', '--', 'Color', [0.5 0.5 0.5])
    sgtitle([session unit], 'Interpreter', 'none');
end
%% focus on nrwd trials
numBinsSWN = 4;
allSWNrwd = NaN(length(allSessions), numBinsSWN);
for i = 1:length(allSessions)
    spikesCurr = spikeFocus{i};
    s = behAnalysisNoPlot_opMD(allSessions{i}, 'simpleFlag', 1);
    target = spikesCurr(s.rwd_Inds);
    edges = linspace(min(target)-0.1, max(target)+0.1, numBinsSWN+1);
    svsNext = NaN(size(s.allChoices));
    svsNext(s.changeChoice_Inds-1) = 1;
    svsNext(s.stayChoice_Inds-1) = 0;
    svsNextNrwd = svsNext(s.rwd_Inds);
    tempSW = NaN(1,numBinsSWN);
    for j = 1:numBinsSWN
        tempSW(j) = mean(svsNextNrwd(target>=edges(j)&target<edges(j+1)), 'omitnan');
    end
    tempSW(~isnan(tempSW)) = zscore(tempSW(~isnan(tempSW)));
    allSWNrwd(i,:) = zscore(tempSW);
end
%% check positive narrow spikes
coeffsMax = [lmMaxAll.coeffs]';
tStatsMax = [lmMaxAll.tStats]';
sigMax = [lmMaxAll.ps]'<0.05;
outcomeInd = 1;
qInd = 2;
indCheck = find(ind==1 & aucMax > 0.5 & sigMax(:,outcomeInd)==1);
waveformCheck = waveformsSession(indCheck,:);
scoreCheck = score(indCheck, 1);
[~, sortInd] = sort(scoreCheck);
figure2; hold on;
plot(1:size(waveformCheck,2), waveformCheck(sortInd,:)+0.5*[1:length(indCheck)]', 'k');
plot([10 26], [-5 -5], 'k', 'LineWidth', 2);
text(18, -7, '0.5 ms', 'HorizontalAlignment', 'center');
set(gca, 'XColor', 'none');
set(gca, 'YColor', 'none');
ylim([-10 0.5*length(indCheck)+3]);
xlim([10 size(waveformCheck,2)-10])
title('Type II, positive coeffs')
%%
