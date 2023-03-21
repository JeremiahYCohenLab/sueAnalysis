load('F:\tmpData\catWithOutcome.mat')
[uniqueSessions, ia, ic] = unique(allSessions');
[root, sep] = currComputer();
for i = 1:length(uniqueSessions)
    pd = parseSessionString_df(uniqueSessions{i}, root,sep);
    savePath = [root pd.animalName sep pd.animalName 'sorted' sep 'itiLick' sep];
    if ~exist(savePath, 'dir')
        mkdir(savePath);
    end
    lickExtract(uniqueSessions{i}, 0);
%     saveFigurePDF(gcf, [savePath uniqueSessions{i} 'itiLick.pdf']);
%     close(gcf)
end
%% waveform classification
[coeff,scores,latent, ~, explained, mu] = pca(zscore(waveformsSession, [], 1));
indAll = {};
dis = {};
for a = 1:10
    [indAll{a}, ~, dis{a}] = kmeans([scores(:, 1:5)], 2);
end
[~,optiInds] = min(cellfun(@mean, dis));
cats = indAll{optiInds};
figure2;
scatter(scores(cats==1, 1), scores(cats==1, 2), 12, 'm')
hold on;
scatter(scores(cats==2, 1), scores(cats==2, 2), 12, 'c')

%%
tMin = 500; % time after solenoid coming back
tb = 500; % time before lick
tf = 1000; % time after lick
binSize = 200;
stepSize = 50;
time = -tb:tf;
slideTime = [0.5*binSize:stepSize:(length(time)-0.5*binSize)] - tb;
pAll = cell(length(allSessions),1);
diffAll = zeros(length(allSessions),length(slideTime));
numAll = zeros(length(allSessions),1);
for i = 1:length(allSessions)
    session = allSessions{i};
    unit = allUnits{i};
    pd = parseSessionString_df(session, root, sep);
    load([pd.sortedFolder session '_sessionData_nL.mat']);
    savePath = [root pd.animalName sep pd.animalName 'sorted' sep 'itiLick' sep 'neurons' sep];
    if ~exist(savePath, 'dir')
        mkdir(savePath);
    end
    
    lickS = lickExtract(session, 0);
    sideITI = lickS.sideITI;
    sideITI(sideITI==3) = 0;
    numAll(i) = sum(sideITI~=0);
    spikeMat = zeros(sum(sideITI~=0), length(time));
    allSpikes = {};
    for j = 1:length(sideITI)
        if sideITI(j) ~= 0
            currSpikes = sessionData(j).(unit)- (sessionData(j).CSon + lickS.solenoidTime) - lickS.itiLatAll(j);
            currSpikes = currSpikes(currSpikes > -tb & currSpikes < tf);
            if isempty(currSpikes)
                currSpikes = zeros(1,0);
            end
            allSpikes{length(allSpikes)+1} = currSpikes;
            tempSpike = currSpikes + tb; 
            if any(tempSpike == 0)
                tempSpike(tempSpike == 0) = 1; % this avoids error in next line
            end
            spikeMat(length(allSpikes), tempSpike) = 1;
            spikeMat(length(allSpikes), tempSpike(currSpikes + lickS.itiLatAll(j) < tMin)) = NaN;
        end
    end
    allSpikes(cellfun(@isempty,allSpikes)) = {zeros(1,0)}; 
    
    spikeMatSlide = zeros(length(allSpikes), length(slideTime));
    for j = 1:length(slideTime)
        spikeMatSlide(:,j) = sum(spikeMat(:, time >= slideTime(j)-0.5*binSize & time <= slideTime(j)+0.5*binSize), 2)/(binSize/1000);
    end
    %% sham
    ind = randi(length(sessionData),2*sum(sideITI~=0),1);
    indSide = rand(length(ind),1);
    indL = ind(indSide <= sum(sideITI==-1)/sum(sideITI~=0));
    indR = ind(indSide > sum(sideITI==-1)/sum(sideITI~=0));
    sideITISham = zeros(1, length(sessionData));
    sideITISham(indL) = -1;
    sideITISham(indR) = 1;
    itiLatAll = lickS.itiLatAll;
    itiLatAll = itiLatAll(sideITI~=0);
    itiLatAllSham = itiLatAll(randi(length(itiLatAll), length(sideITISham),1));
    spikeMatSham = zeros(sum(sideITISham~=0), length(time));
    allSpikesSham = {};
    for j = 1:length(sideITISham)
        if sideITISham(j) ~= 0
            currSpikes = sessionData(j).(unit)- (sessionData(j).CSon + lickS.solenoidTime) - itiLatAllSham(j);
            currSpikes = currSpikes(currSpikes > -tb & currSpikes < tf);
            if isempty(currSpikes)
                currSpikes = zeros(1,0);
            end
            allSpikesSham{length(allSpikesSham)+1} = currSpikes;
            tempSpike = currSpikes + tb; 
            if any(tempSpike == 0)
                tempSpike(tempSpike == 0) = 1; % this avoids error in next line
            end
            spikeMatSham(length(allSpikesSham), tempSpike) = 1;
        end
    end
    allSpikesSham(cellfun(@isempty,allSpikesSham)) = {zeros(1,0)}; 
    
    spikeMatSlideSham = zeros(length(allSpikesSham), length(slideTime));
    for j = 1:length(slideTime)
        spikeMatSlideSham(:,j) = sum(spikeMatSham(:, time >= slideTime(j)-0.5*binSize & time <= slideTime(j)+0.5*binSize), 2)/(binSize/1000);
    end
    
    %% test
    pValues = zeros(1, length(slideTime));
    for j = 1:length(slideTime)
        [~, pValues(j)] = ttest2(spikeMatSlide(:,j), spikeMatSlideSham(:,j));
    end   
    pAll{i} = pValues;
    diffAll(i,:) = (mean(spikeMatSlide, 'omitnan') - mean(spikeMatSlideSham, 'omitnan'))./mean(spikeMatSlideSham, 'all');

        %% plot
%     figure;
%     subplot(1,3,1); hold on;
%     lineFormat.LineWidth = 2;
%     [~, ind] = sort(sideITI(sideITI~=0));
%     plotSpikeRaster(allSpikes(ind), 'PlotType','vertline', 'LineFormat', lineFormat);
%     hold on; plot([0 0], [0 length(allSpikes)], 'LineWidth', 1, 'Color', 'r');
%     xlabel('time from lick (ms)')
%     subplot(3,3,2); hold on;
%     histogram(lickS.itiLatAll(sideITI~=0), 0:500:max(lickS.itiLatAll(sideITI~=0))+1, 'Normalization', 'probability');
%     histogram(itiLatAllSham, 0:500:max(lickS.itiLatAll(sideITI~=0))+1, 'Normalization', 'probability');
%     plot([tMin tMin], [0 0.3], 'LineWidth', 1, 'Color', 'r');
%     legend({'real', 'sham'})
%     xlabel('time from solenoid back (ms)')
%     subplot(3,3,5); hold on;
%     plotFilled(slideTime, spikeMatSlide, 'b'); 
%     plotFilled(slideTime, spikeMatSlideSham, [0.6 0.6 0.6]); 
%     plot([0 0], [0 5], 'LineStyle', '--', 'LineWidth', 1, 'Color', 'r');
%     xlabel('time from lick (ms)')
%     scatter(slideTime(pValues<0.05), mean(spikeMatSlide(:,pValues<0.05), 'omitnan')+0.8*ones(1,sum(pValues<0.05)), 11, 'r', 'filled');
%     sgtitle([session unit], 'Interpreter', 'none');
% %     
% %     screen = get(0,'Screensize');
% %     screen(4) = screen(4) - 100;
%     set(gcf, 'Position', screen);
    
%     saveFigurePDF(gcf, [savePath session 'SpikesItiLick.pdf']);
%     close(gcf)
end
%% pca of ITI respond shape
diffAllClean = diffAll(numAll>=30,:);
catsClean = cats(numAll>=30);
[coeff, scoresITI,latent,tsquared,explained,mu] = pca(zscore(diffAllClean, 1));
[~,sorted] = sort(scoresITI(:,1));
figure2; hold on;
imagesc(slideTime, 1:size(diffAllClean,1), diffAllClean(sorted, :));
plot([0 0], [1 size(diffAllClean,1)], 'w');
title('PC sorted')

[~, maxLoc] = max(diffAllClean, [], 2);
[~, sorted] = sort(maxLoc);
figure2;hold on;
imagesc(slideTime, 1:size(diffAllClean,1), diffAllClean(sorted, :));
title('max time sorted')
plot([0 0], [1 size(diffAllClean,1)], 'w');

maxPeak = sum(diffAllClean(:,slideTime>=-300&slideTime<=0), 2);
[~, sorted] = sort(maxPeak);
figure2; hold on;
imagesc(slideTime, 1:size(diffAllClean,1), diffAllClean(sorted, :));
title('preLick sorted')
plot([0 0], [1 size(diffAllClean,1)], 'w');

diffPeak = mean(diffAllClean(:,slideTime>=-300&slideTime<=0), 2) - mean(diffAllClean(:,slideTime>0&slideTime<=500), 2);
[~, sorted] = sort(diffPeak);
figure2; hold on;
imagesc(slideTime, 1:size(diffAllClean,1), diffAllClean(sorted, :));
title('diffSpikes sorted')
plot([0 0], [1 size(diffAllClean,1)], 'w');
%% compare two clusters
diffAllClean = diffAll(numAll>=30,:);
catsClean = cats(numAll>=30);
[coeff, scoresITI,latent,tsquared,explained,mu] = pca(zscore(diffAllClean, 1));

maxMap = max(diffAllClean,[], 'all');
minMap = min(diffAllClean,[], 'all');
negNum = round(300*maxMap/(maxMap-minMap));
posNum = round(300*maxMap/(maxMap-minMap));
myMap = [[linspace(0, 1, negNum)' linspace(0,1,negNum)' linspace(1,1,negNum)']; ... 
    [linspace(1, 1, posNum)' linspace(1,0,posNum)' linspace(1,0,posNum)']];
yH = max([sum(catsClean==1), sum(catsClean==2)]);

figure2;
subplot(1,2,1);hold on;
[~, sorted] = sort(maxPeak(catsClean==1));
temp = diffAllClean(catsClean==1,:);
imagesc(slideTime, 1:size(temp,1), temp(sorted,:));
plot([0 0], [1 size(temp,1)], 'k');
colormap(myMap)
colorbar;
caxis([-maxMap maxMap])
set(gca, 'Box', 'off')
set(gca, 'TickDir', 'Out')
set(gca, 'XTick', [-500 0 500], 'FontSize', 12)
set(gca, 'YColor', 'none')
xlabel('Time from lick (ms)', 'FontSize', 15)
set(gca, 'Box', 'off')
ylim([-2 yH]); xlim([-500 1000])
colorbar('Ticks', [-floor(maxMap) ceil(minMap) 0 floor(maxMap)], 'FontSize', 12)

subplot(1,2,2);hold on;
[~, sorted] = sort(maxPeak(catsClean==2));
temp = diffAllClean(catsClean==2,:);
imagesc(slideTime, 1:size(temp,1), temp(sorted,:));
plot([0 0], [1 size(temp,1)], 'k');
colorbar;
caxis([-maxMap maxMap])
colormap(myMap)
set(gca, 'Box', 'off')
set(gca, 'TickDir', 'Out')
set(gca, 'XTick', [-500 0 500], 'FontSize', 12)
set(gca, 'YColor', 'none')
xlabel('Time from lick (ms)', 'FontSize', 15)
set(gca, 'Box', 'off')
ylim([-2 yH]); xlim([-500 1000])
colorbar('Ticks', [-floor(maxMap) ceil(minMap) 0 floor(maxMap)], 'FontSize', 12)
%%