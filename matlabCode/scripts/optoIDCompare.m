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

    quality = nums(:,1);
    sessionListTmp = sessionListTmp(quality<=maxLratio);
    unitListTmp = unitListTmp(quality<=maxLratio);
    optoUnitListTmp = optoUnitListTmp(quality<=maxLratio);
    subFoldersTmp = subFoldersTmp(quality<=maxLratio);
    
    sessionList = [sessionList; sessionListTmp];
    unitList = [unitList; unitListTmp];
    optoUnitList = [optoUnitList; optoUnitListTmp];
    subFolders = [subFolders; subFoldersTmp];

end
%%
[root, sep] = currComputer();
qualInd = zeros(size(sessionList));
waveformsSession = [];
respPs = [];
respLats = [];
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
    
    % allPs
    respPs = [respPs; met.spikeProp];
    % allLats
    respLats = [respLats; met.spikeLat];
    
end

respLats = respLats/1000;
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
set(gca, 'XColor', 'none')
set(gca, 'YColor', 'none')
%% calculate peak respP
bestP = max(respPs, [], 2);
edge = linspace(min(bestP)-0.001, max(bestP)+0.001, 4);
figure2; hold on;
histogram(bestP(ind==1), edge, 'FaceColor', color1, 'EdgeColor', 'none', 'Normalization', 'probability');
histogram(bestP(ind==2), edge, 'FaceColor', color2, 'EdgeColor', 'none', 'Normalization', 'probability');
legend({'Type II', 'Type I'})
xlabel('max(P(spike))')
title('max(P(spike))')
%% calculate mean respP
meanP = mean(respPs, 2);
edge = linspace(min(meanP)-0.001, max(meanP)+0.001, 8);
figure2; hold on;
histogram(meanP(ind==1), edge, 'FaceColor', color1, 'EdgeColor', 'none', 'Normalization', 'probability');
histogram(meanP(ind==2), edge, 'FaceColor', color2, 'EdgeColor', 'none', 'Normalization', 'probability');
legend({'Type II', 'Type I'})
xlabel('mean(P(spike))')
title('mean(P(spike))')
%% calculate min respLat
minLat = min(respLats/1000, [], 2);
edge = linspace(min(minLat)-0.001, max(minLat)+0.001, 8);
figure2; hold on;
histogram(minLat(ind==1), edge, 'FaceColor', color1, 'EdgeColor', 'none', 'Normalization', 'probability');
histogram(minLat(ind==2), edge, 'FaceColor', color2, 'EdgeColor', 'none', 'Normalization', 'probability');
legend({'Type II', 'Type I'})
xlabel('min(spikeLat)')
title('min(spikeLat)')
%% calculate mean respLat
meanLat = mean(respLats, 2, 'omitnan');
edge = linspace(min(meanLat)-0.001, max(meanLat)+0.001, 8);
figure2; hold on;
histogram(meanLat(ind==1), edge, 'FaceColor', color1, 'EdgeColor', 'none', 'Normalization', 'probability');
histogram(meanLat(ind==2), edge, 'FaceColor', color2, 'EdgeColor', 'none', 'Normalization', 'probability');
legend({'Type II', 'Type I'})
xlabel('mean(spikeLat)')
title('mean(spikeLat)')
%% scatters
figure2Wide; 
subplot(1,2,1); hold on;
scatter(bestP(ind==1), minLat(ind==1), 15, color1, 'LineWidth', 1.5)
scatter(bestP(ind==2), minLat(ind==2), 15, color2, 'LineWidth', 1.5)
legend({'Type II', 'Type I'})
xlabel('maxP')
ylabel('minLat')

subplot(1,2,2); hold on;
scatter(meanP(ind==1), meanLat(ind==1), 15, color1, 'LineWidth', 1.5)
scatter(meanP(ind==2), meanLat(ind==2), 15, color2, 'LineWidth', 1.5)
legend({'Type II', 'Type I'})
xlabel('meanP')
ylabel('meanLat')
%%



%% 'adaptation'
% Presponse
pChange = mean(respPs(:,6:end), 2, 'omitnan') - mean(respPs(:,1:5), 2, 'omitnan');
edge = linspace(min(pChange)-0.001, max(pChange)+0.001, 15);
figure2; hold on;
histogram(pChange(ind==1), edge, 'FaceColor', color1, 'EdgeColor', 'none', 'Normalization', 'probability');
histogram(pChange(ind==2), edge, 'FaceColor', color2, 'EdgeColor', 'none', 'Normalization', 'probability');
p = ranksum(pChange(ind==1), pChange(ind==2));
legend({'Type II', 'Type I'})
xlabel('lateP-earlyP')
title(['lateP-earlyP' ' ' num2str(p)])
%% respLat
latChange = respLats(:,end) - respLats(:,1);
edge = linspace(min(latChange)-0.001, max(latChange)+0.001, 15);
figure2; hold on;
histogram(latChange(ind==1), edge, 'FaceColor', color1, 'EdgeColor', 'none', 'Normalization', 'probability');
histogram(latChange(ind==2), edge, 'FaceColor', color2, 'EdgeColor', 'none', 'Normalization', 'probability');
legend({'Type II', 'Type I'})
xlabel('lastLat-firstLat')
title('lastLat-firstLat')

%% first resp
figure2; hold on;
edge = linspace(min(respLats(:,1)), max(respLats(:,1)), 8);
histogram(respLats(ind==1, 1), edge, 'FaceColor', color1, 'EdgeColor', 'none', 'Normalization', 'probability');
histogram(respLats(ind==2, 1), edge, 'FaceColor', color2, 'EdgeColor', 'none', 'Normalization', 'probability');
legend({'Type II', 'Type I'})
xlabel('firstLat')
title('firstLat')
%% first respP
figure2; hold on;
edge = linspace(min(respPs(:,1)), max(respPs(:,1)), 11);
histogram(respPs(ind==1, 1), edge, 'FaceColor', color1, 'EdgeColor', 'none', 'Normalization', 'probability');
histogram(respPs(ind==2, 1), edge, 'FaceColor', color2, 'EdgeColor', 'none', 'Normalization', 'probability');
legend({'Type II', 'Type I'})
xlabel('firstP')
title('firstP')
%% first respP
figure2; hold on;
edge = linspace(min(respPs(:,end)), max(respPs(:,end)), 11);
histogram(respPs(ind==1, end), edge, 'FaceColor', color1, 'EdgeColor', 'none', 'Normalization', 'probability');
histogram(respPs(ind==2, end), edge, 'FaceColor', color2, 'EdgeColor', 'none', 'Normalization', 'probability');
legend({'Type II', 'Type I'})
xlabel('lastP')
title('lastP')
%% scatter
figure2Wide; 
subplot(1,2,1); hold on;
scatter(mean(respPs(ind==1, 1:2), 2), respPs(ind==1, end), 15, 'm', 'LineWidth', 1.5, 'MarkerEdgeAlpha', 1)
scatter(mean(respPs(ind==2, 1:2), 2), respPs(ind==2, end), 15, 'c', 'LineWidth', 1.5, 'MarkerEdgeAlpha', 1)
legend({'Type II', 'Type I'})
xlabel('firsP')
ylabel('firstLat')
%%
subplot(1,2,2); hold on;
scatter(meanP(ind==1), meanLat(ind==1), 15, color1, 'LineWidth', 1.5)
scatter(meanP(ind==2), meanLat(ind==2), 15, color2, 'LineWidth', 1.5)
legend({'Type II', 'Type I'})
xlabel('meanP')
ylabel('meanLat')

%%
empInd = sum(isnan(respLats), 2);
respPsZS = respPs(empInd==0,:);
respLatsZS = respLats(empInd==0,:);
indZS = ind(empInd==0);
%%
[coeffs, scores] = pca(zscore(respPsZS));
[coeffsLat, scoresLat] = pca(zscore(respLatsZS));
%%
i = 1;
figure2;
edges = linspace(min(scores(:,i))-0.01, max(scores(:,i))+0.01, 15);
hold on; histogram(scores(indZS==1, i), edges, 'FaceColor', color1, 'Normalization', 'probability')
hold on; histogram(scores(indZS==2, i), edges, 'FaceColor', color2, 'Normalization', 'probability')
[h, p] = ttest2(scores(indZS==1, i), scores(indZS==2, i));
title([num2str(i) ' ' num2str(p)])
%%
%%
figure2;
scatter(pChange(ind==2), meanLat(ind==2), 15, color2, 'filled', 'MarkerEdgeColor', 'none')
hold on;
scatter(pChange(ind==1), meanLat(ind==1), 15, color1, 'MarkerEdgeAlpha', 1, 'LineWidth', 1.5)
xlabel('pChange(lateHalf-firstHalf)')
ylabel('meanSpikeLat (ms)')
legend({'Type I', 'Type II'})
%%
figure2;
scatter(pChange(ind==1), meanLat(ind==1), 15, color1, 'filled')
hold on;
scatter(pChange(ind==2), meanLat(ind==2), 15, color2, 'filled')

xlabel('pChange(lateHalf-firstHalf)')
ylabel('meanSpikeLat (ms)')
legend({'Type I', 'Type II'})
%%
figure2;
hold on;
histogram2(pChange(ind==1), meanLat(ind==1), 'FaceColor', color1, 'Normalization', 'probability', 'EdgeColor', 'none', 'FaceAlpha', 0.7)
histogram2(pChange(ind==2), meanLat(ind==2), 'FaceColor', color2, 'Normalization', 'probability', 'EdgeColor', 'none', 'FaceAlpha', 0.7)
legend({'Type II', 'Type I'})
%%