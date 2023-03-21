animalNames = {'ZS059', 'ZS060', 'ZS061', 'ZS062'};
[root, sep] = currComputer();
col = 'good';
model = '5params';
%% get corresponding unit files
sessionList = {};
unitList = {};
optoUnitList = {};
subFolders = {};
maxLratio = 0.05;
depths = [];
depthsRank = [];
for a = 1:length(animalNames)
    animalName = animalNames{a};
    xlFile = [animalName '.xlsx'];
    savePath = [root animalName sep animalName 'sorted' sep 'optoFiles' sep];
    if ~exist(savePath, 'dir')
        mkdir(savePath)
    end

    [nums, unitsInfo,~] = xlsread([root xlFile], 'neurons');
    sessionListTmp = unitsInfo(2:end,1); 
    unitListTmp = unitsInfo(2:end,2);
    optoUnitListTmp = unitsInfo(2:end,8);
    subFoldersTmp = unitsInfo(2:end,9);
    depthTmp = nums(:, 9);

    quality = nums(:,1);
    sessionListTmp = sessionListTmp(quality<=maxLratio);
    unitListTmp = unitListTmp(quality<=maxLratio);
    optoUnitListTmp = optoUnitListTmp(quality<=maxLratio);
    subFoldersTmp = subFoldersTmp(quality<=maxLratio);
    depthTmp = depthTmp(quality<=maxLratio);
    
    sessionList = [sessionList; sessionListTmp];
    unitList = [unitList; unitListTmp];
    optoUnitList = [optoUnitList; optoUnitListTmp];
    subFolders = [subFolders; subFoldersTmp];
    depths = [depths; depthTmp];
    [~, ~, ranked] = unique(depthTmp);
    depthsRank = [depthsRank; ranked];
end
%%
[root, sep] = currComputer();
qualInd = zeros(size(sessionList));
waveformsSession = [];
for i = 1:length(sessionList)
    session = sessionList{i};
    unit = unitList{i};
    pd = parseSessionString_df(session, root, sep);
    neuralynxDataPath = [pd.sortedFolder session '_sessionData_nL.mat'];
    unitMetDir = [pd.nLynxFolderSession session '_' unit '_met.mat'];
    sortedFolderLocation = [pd.sortedFolder];
    load(unitMetDir);

    respInds = find(met.spikeProp>=0.8);
    if isempty(respInds)
        
        continue
    else 
        respLat = nanmin(met.spikeLat(respInds));
        if respLat > 15000 || met.isiV > 0.001 || met.distance>0.3 || isnan(met.distance)
            continue
        end
    end
    qualInd(i) = 1;
    [~,peakCh] = max(max(metSess.waveform));
    [~,peakSamp] = max(metSess.waveform(:,peakCh));
    tempWF = metSess.waveform(:,peakCh);
    % filter
    fc = 6000;
    fs = 32000;
    [b, a] = butter(2,fc/(fs/2),'low');
    tempWF = filtfilt(b, a, tempWF);
    tempWF = tempWF(peakSamp-30:peakSamp+40);
    peak = max(tempWF);
    waveformsSession = [waveformsSession; tempWF'/peak];  
    
end
%%
%% clustering
color1 = [0.2 0.2 0.2];
color2 = [0.8 0.8 0.8];
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
%% depth
depthsRank = depthsRank(logical(qualInd));
depths = depths(logical(qualInd));
%%
numBins = 10;
figure2;hold on;
edges = linspace(min(depthsRank)-0.1, max(depthsRank)+0.1, numBins);
histogram(depthsRank(ind==1), edges, 'FaceColor', color1, 'Normalization', 'probability', 'EdgeColor', 'none');
histogram(depthsRank(ind==2), edges, 'FaceColor', color2, 'Normalization', 'probability', 'EdgeColor', 'none');
ylabel('probability')
set(gca, 'YTick', [0:0.25:0.75]);
title('rank')
xlabel('depth rank')
set(gca, 'TickDir', 'out');
set(gca, 'XTick', [0:10:30]);
legend({'Type II', 'TypeI'})

figure2; hold on;
count1 = histcounts(depthsRank(ind==1), edges);
count2 = histcounts(depthsRank(ind==2), edges);
plot(0.5*(edges(1:end-1) + edges(2:end)), count1./(count1+count2), 'Color', color1, 'LineWidth',2, 'LineStyle', '-');
plot(0.5*(edges(1:end-1) + edges(2:end)), count2./(count1+count2), 'Color', color2, 'LineWidth',2, 'LineStyle', '-');
ylabel('proportion of two types')
set(gca, 'YTick', [0:0.25:1]);
ylim([-0.1 1.1])
title('rank')
xlabel('depth rank')
set(gca, 'TickDir', 'out');
set(gca, 'XTick', [0:10:30]);
legend({'Type II', 'TypeI'})

%%
numBins = 10;
figure2;hold on;
edges = linspace(min(depths)-0.1, max(depths)+0.1, numBins);
histogram(depths(ind==1), edges, 'FaceColor', color1, 'Normalization', 'probability', 'EdgeColor', 'none');
histogram(depths(ind==2), edges, 'FaceColor', color2, 'Normalization', 'probability', 'EdgeColor', 'none');
title('raw')
xlabel('depth')
ylabel('probability')
set(gca, 'YTick', [0:0.25:0.75]);
title('rank')
set(gca, 'TickDir', 'out');
legend({'Type II', 'TypeI'})

figure2; hold on;
count1 = histcounts(depths(ind==1), edges);
count2 = histcounts(depths(ind==2), edges);
plot(0.5*(edges(1:end-1) + edges(2:end)), count1./(count1+count2), 'Color', color1, 'LineWidth',2, 'LineStyle', '-');
plot(0.5*(edges(1:end-1) + edges(2:end)), count2./(count1+count2), 'Color', color2, 'LineWidth',2, 'LineStyle', '-');
ylabel('proportion of two types')
set(gca, 'YTick', [0:0.25:1]);
ylim([-0.1 1.1])
title('rank')
xlabel('depth (mm)')
set(gca, 'TickDir', 'out');
% set(gca, 'XTick', [0:10:30]);
legend({'Type II', 'TypeI'})
%%