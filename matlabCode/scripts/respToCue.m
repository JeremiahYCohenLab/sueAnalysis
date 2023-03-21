load('F:\tmpData\allUnitAUC.mat');
myMap = repmat(linspace(0.98, 0.5, 1000)', 1, 3);
numBins = 6;
color1 = [0.2 0.2 0.2];
color2 = [0.7 0.7 0.7];
%%
blBin = [1 0];
goBin = [0 0.3];
stepSize = 100;
binSize = 100;
allBaseline = cell(size(allSessions));
allGoRaw = cell(size(allSessions));
for i = 1:length(allSessions)
    [~,matCue] = getUnitMatCue(allSessions{i}, allUnits{i}, blBin(1), blBin(2), stepSize, binSize);
    allBaseline{i} = sum(matCue, 2)/abs(blBin(2)-blBin(1));
    [~,matCue] = getUnitMatCue(allSessions{i}, allUnits{i}, goBin(1), goBin(2), stepSize, binSize);
    allGoRaw{i} = sum(matCue, 2)/abs(goBin(2)-goBin(1));
end
%%
meanBl = cellfun(@mean, allBaseline);
meanGo = cellfun(@mean, allGoRaw);
meanEvoke = (meanGo-meanBl)./meanBl;
figure2; hold on;
scatter(meanBl(ind==1), meanGo(ind==1), 15, color1);
scatter(meanBl(ind==2), meanGo(ind==2), 15, color2);
legend({'typeII', 'typeI'})
figure2; hold on;
edges = linspace(min(meanBl)-0.01, max(meanBl)+0.01, 15);
histogram(meanBl(ind==1), edges, 'FaceColor', color1, 'EdgeColor', 'none');
histogram(meanBl(ind==2), edges, 'FaceColor', color2, 'EdgeColor', 'none');
legend({'typeII', 'typeI'})
title('baseline')

figure2; hold on;
edges = linspace(min(meanGo)-0.01, max(meanGo)+0.01, 15);
histogram(meanGo(ind==1), edges, 'FaceColor', color1, 'EdgeColor', 'none');
histogram(meanGo(ind==2), edges, 'FaceColor', color2, 'EdgeColor', 'none');
legend({'typeII', 'typeI'})
title('goCue')

figure2; hold on;
edges = linspace(min(meanEvoke)-0.01, max(meanEvoke)+0.01, 15);
histogram(meanEvoke(ind==1), edges, 'FaceColor', color1, 'EdgeColor','none');
histogram(meanEvoke(ind==2), edges, 'FaceColor', color2, 'EdgeColor','none');
legend({'typeII', 'typeI'})
title('evoke')
%%
for i = 1:length(allSessions)

    currBl = allBaseline{i};
    currGo = allGoRaw{i};
%         edges = linspace(min(currBl)-0.1, max(currBl)+ 0.1, numBins+1);
    edges = binEqualSize(currBl, numBins+1);
    meanGo = zeros(1, length(edges)-1);
    meanBl = zeros(1, length(edges)-1);
    semGo = zeros(1, length(edges)-1);

    for j = 1:length(edges)-1
        if j == length(edges)-1
            tmpInd = currBl >= edges(j) & currBl <= edges(j+1);
        else
            tmpInd = currBl >= edges(j) & currBl < edges(j+1);
        end
        meanBl(j) = mean(currBl(tmpInd), 'omitnan');
        meanGo(j) = mean(currGo(tmpInd), 'omitnan');
        semGo(j) = sem(currGo(tmpInd));
    end
%     figure; hold on;
%     hist3([currBl, currGo], [8, length(unique(currGo))],  'CdataMode', 'auto', 'edgeColor', 'none');
%     view(2)
%     colormap(myMap)
%     plot3(meanBl, meanGo, 300*ones(1, length(edges)-1), 'Color', [0 0.8 0.8], 'LineWidth', 2);
%     fill3([meanBl flip(meanBl)], [meanGo-semGo flip(meanGo+semGo)], 300*ones(1, 2*(length(edges)-1)), [0 0.8 0.8], 'EdgeColor', 'none', 'FaceAlpha', 0.2);
%     title([allSessions{i}, ' ' allUnits{i}], 'Interpreter', 'none')
end
%%
numBins = 5;
popMeans = zeros(length(allSessions), numBins);
popSems = zeros(length(allSessions), numBins);
popBls = zeros(length(allSessions), numBins);

for i = 1:length(allSessions)
    currBl = allBaseline{i};
    currGo = allGoRaw{i};
    edges = linspace(min(currBl)-0.1, max(currBl)+ 0.1, numBins+1);
    meanGo = zeros(1, length(edges)-1);
    meanBl = zeros(1, length(edges)-1);
    semGo = zeros(1, length(edges)-1);

    for j = 1:length(edges)-1
        if j == length(edges)-1
            tmpInd = currBl >= edges(j) & currBl <= edges(j+1);
        else
            tmpInd = currBl >= edges(j) & currBl < edges(j+1);
        end
        meanBl(j) = mean(currBl(tmpInd), 'omitnan');
        meanGo(j) = mean(currGo(tmpInd), 'omitnan');
        semGo(j) = sem(currGo(tmpInd));
    end    
    
    popMeans(i,:) = meanGo;
    popSems(i,:) = semGo;
    popBls(i,:) = meanBl;
end
%% plot filled line 
figure2; hold on;
semAll = zeros(1, numBins);
for j = 1:numBins
    semAll(j) = sem(popMeans(ind==1,j));
end
plot(mean(popBls(ind==1,:), 'omitnan'), mean(popMeans(ind==1,:), 'omitnan'), 'Color', [0 0.8 0.8], 'LineWidth', 2);
fill([mean(popBls(ind==1,:), 'omitnan') flip(mean(popBls(ind==1,:), 'omitnan'))], [mean(popMeans(ind==1,:), 'omitnan')-semAll flip(mean(popMeans(ind==1,:), 'omitnan')+semAll)], [0 0.8 0.8], 'EdgeColor', 'none', 'FaceAlpha', 0.25)

semAll = zeros(1, numBins);
for j = 1:numBins
    semAll(j) = sem(popMeans(ind==2,j));
end
plot(mean(popBls(ind==2,:), 'omitnan'), mean(popMeans(ind==2,:), 'omitnan'), 'Color', [1 0 1], 'LineWidth', 2);
fill([mean(popBls(ind==2,:), 'omitnan') flip(mean(popBls(ind==2,:), 'omitnan'))], [mean(popMeans(ind==2,:), 'omitnan')-semAll flip(mean(popMeans(ind==2,:), 'omitnan')+semAll)], [1 0 1], 'EdgeColor', 'none', 'FaceAlpha', 0.25)
%% calculate regression
tStatsMax = zeros(length(allSessions),1);
sigMax = zeros(length(allSessions),1);
coeffsMax = zeros(length(allSessions),1);
for i = 1:length(allSessions)
    currBl = allBaseline{i};
    currGo = allGoRaw{i};
    lm = fitlm(currBl, currGo);
    tStatsMax(i) = lm.Coefficients.tStat(2);
    sigMax(i) = lm.Coefficients.pValue(2);
    coeffsMax(i) = lm.Coefficients.Estimate(2);
end
%% coeffs
figure; hold on;
edges = linspace(min(coeffsMax(:,1))-0.01, max(coeffsMax(:,1))+0.01, 31);
barMeams = 0.5*(edges(1:end-1) + edges(2:end));
count1 = histcounts(coeffsMax(ind==1&sigMax(:,1)>=0.05,1), edges);
count2 = histcounts(coeffsMax(ind==1&sigMax(:,1)<0.05,1), edges);
count3 = count1 + count2;
bar(barMeams, count3, 1, 'FaceColor', [0.7 0.7 0.7], 'FaceAlpha', 1, 'EdgeColor', 'none');
bar(barMeams, count2, 1, 'FaceColor', [0 0.8 0.8], 'FaceAlpha', 1, 'EdgeColor', 'none');
title('typeI')

figure; hold on;
edges = linspace(min(coeffsMax(:,1))-0.01, max(coeffsMax(:,1))+0.01, 31);
barMeams = 0.5*(edges(1:end-1) + edges(2:end));
count1 = histcounts(coeffsMax(ind==2&sigMax(:,1)>=0.05,1), edges);
count2 = histcounts(coeffsMax(ind==2&sigMax(:,1)<0.05,1), edges);
count3 = count1 + count2;
bar(barMeams, count3, 1, 'FaceColor', [0.7 0.7 0.7], 'FaceAlpha', 1, 'EdgeColor', 'none');
bar(barMeams, count2, 1, 'FaceColor', [1 0 1], 'FaceAlpha', 1, 'EdgeColor', 'none');
title('typeII')
%% tStats
figure; hold on;
edges = linspace(min(tStatsMax(:,1))-0.01, max(tStatsMax(:,1))+0.01, 20);
barMeams = 0.5*(edges(1:end-1) + edges(2:end));
count1 = histcounts(tStatsMax(ind==1&sigMax(:,1)>=0.05,1), edges);
count2 = histcounts(tStatsMax(ind==1&sigMax(:,1)<0.05,1), edges);
count3 = count1 + count2;
bar(barMeams, count3, 1, 'FaceColor', [0.7 0.7 0.7], 'FaceAlpha', 1, 'EdgeColor', 'none');
bar(barMeams, count2, 1, 'FaceColor', [0 0.8 0.8], 'FaceAlpha', 1, 'EdgeColor', 'none');
title('typeI')

figure; hold on;
edges = linspace(min(tStatsMax(:,1))-0.01, max(tStatsMax(:,1))+0.01, 20);
barMeams = 0.5*(edges(1:end-1) + edges(2:end));
count1 = histcounts(tStatsMax(ind==2&sigMax(:,1)>=0.05,1), edges);
count2 = histcounts(tStatsMax(ind==2&sigMax(:,1)<0.05,1), edges);
count3 = count1 + count2;
bar(barMeams, count3, 1, 'FaceColor', [0.7 0.7 0.7], 'FaceAlpha', 1, 'EdgeColor', 'none');
bar(barMeams, count2, 1, 'FaceColor', [1 0 1], 'FaceAlpha', 1, 'EdgeColor', 'none');
title('typeII')
%% map Q coeffcient ind to units
qInd = NaN(length(allSessions),1);
for i = 1:length(allSessions)
    sesInd = cellfun(@(x) strcmp(x, allSessions{i}), allSessionsQ);
    unitInd = cellfun(@(x) strcmp(x, allUnits{i}), allUnitsQ);
    if sum(sesInd&unitInd)>0
        qInd(i) = cats(sesInd&unitInd);
    end
end
%% tStats
figure; hold on;
edges = linspace(min(tStatsMax(:,1))-0.01, max(tStatsMax(:,1))+0.01, 31);
barMeams = 0.5*(edges(1:end-1) + edges(2:end));
count1 = histcounts(tStatsMax(qInd==0&sigMax(:,1)>=0.05,1), edges);
count2 = histcounts(tStatsMax(qInd==0&sigMax(:,1)<0.05,1), edges);
count3 = count1 + count2;
bar(barMeams, count3, 1, 'FaceColor', [0.7 0.7 0.7], 'FaceAlpha', 1, 'EdgeColor', 'none');
bar(barMeams, count2, 1, 'FaceColor', [0 0.8 0.8], 'FaceAlpha', 1, 'EdgeColor', 'none');
title('typeI')

figure; hold on;
edges = linspace(min(tStatsMax(:,1))-0.01, max(tStatsMax(:,1))+0.01, 31);
barMeams = 0.5*(edges(1:end-1) + edges(2:end));
count1 = histcounts(tStatsMax(qInd==1&sigMax(:,1)>=0.05,1), edges);
count2 = histcounts(tStatsMax(qInd==1&sigMax(:,1)<0.05,1), edges);
count3 = count1 + count2;
bar(barMeams, count3, 1, 'FaceColor', [0.7 0.7 0.7], 'FaceAlpha', 1, 'EdgeColor', 'none');
bar(barMeams, count2, 1, 'FaceColor', [1 0 1], 'FaceAlpha', 1, 'EdgeColor', 'none');
title('typeII')
%%  map outcome coeffcient ind to units
oInd = NaN(length(allSessions),1);
for i = 1:length(allSessions)
    sesInd = cellfun(@(x) strcmp(x, allSessions{i}), allSessionsO);
    unitInd = cellfun(@(x) strcmp(x, allUnits{i}), allUnitsO);
    if sum(sesInd&unitInd)>0
        oInd(i) = ind(sesInd&unitInd);
    end
end
%% coeffs
figure; hold on;
edges = linspace(min(coeffsMax(:,1))-0.01, max(coeffsMax(:,1))+0.01, 31);
barMeams = 0.5*(edges(1:end-1) + edges(2:end));
count1 = histcounts(coeffsMax(oInd==1&sigMax(:,1)>=0.05,1), edges);
count2 = histcounts(coeffsMax(oInd==1&sigMax(:,1)<0.05,1), edges);
count3 = count1 + count2;
bar(barMeams, count3, 1, 'FaceColor', [0.7 0.7 0.7], 'FaceAlpha', 1, 'EdgeColor', 'none');
bar(barMeams, count2, 1, 'FaceColor', [0 0.8 0.8], 'FaceAlpha', 1, 'EdgeColor', 'none');
title('typeI')

figure; hold on;
edges = linspace(min(coeffsMax(:,1))-0.01, max(coeffsMax(:,1))+0.01, 31);
barMeams = 0.5*(edges(1:end-1) + edges(2:end));
count1 = histcounts(coeffsMax(oInd==2&sigMax(:,1)>=0.05,1), edges);
count2 = histcounts(coeffsMax(oInd==2&sigMax(:,1)<0.05,1), edges);
count3 = count1 + count2;
bar(barMeams, count3, 1, 'FaceColor', [0.7 0.7 0.7], 'FaceAlpha', 1, 'EdgeColor', 'none');
bar(barMeams, count2, 1, 'FaceColor', [1 0 1], 'FaceAlpha', 1, 'EdgeColor', 'none');
title('typeII')
%%

