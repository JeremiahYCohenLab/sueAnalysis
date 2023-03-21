function psthPattern(animalNames, category, varargin)
%task and model parameters
p = inputParser;
% default parameters if none given
p.addParameter('numCat', 2);
p.addParameter('maxTrial', 1000); 
p.addParameter('binSize', 150)% in ms
p.addParameter('stepSize', 100)
p.addParameter('saveFigFlag', 1);
p.addParameter('tb', 1)% in s
p.addParameter('tf', 3)% in s
p.addParameter('tfCue',0.75)% in s, time after cue aligned time
p.addParameter('tfRwd', 1.5)% in s, time after rwd delivery time
p.addParameter('focusWin', [300 1000])
p.parse(varargin{:});
maxTrial = p.Results.maxTrial;
% basic info
[root, sep] = currComputer();
% time window
time = -1000*p.Results.tb:1000*p.Results.tf;
midPoints = (0.5*p.Results.binSize + 1):p.Results.stepSize:(length(time)-0.5*p.Results.binSize);
slideTime = midPoints - p.Results.tb*1000;
%% combined information
baselineFreq = [];
firstSpikeFreq = [];
firstSpikeLat = [];
width = [];
spikeLatChange = [];
spikeFreqChange = [];
waveforms = [];
waveformsSession = [];
allSessions = {};
allUnits = {};
psthComb = [];

populationDepths = []; % all estimated positions
rwdvsNrwd = [];
excit = [];
preSessBaseline = [];
myKernel = [1, 2, 4, 2, 1];
myKernel = myKernel/sum(myKernel);
cats = []; % 0 for nothing, 1 for negative, 2 for positive
%% animal loop
for ani = 1:length(animalNames)
    % load model fitting results
    animalName = animalNames{ani};
    % load sessionList and unitList
    xlFile = [animalName '.xlsx'];
    [nums, unitsInfo,~] = xlsread([root xlFile], 'neurons');
    sessionList = unitsInfo(2:end,1); 
    unitList = unitsInfo(2:end,2);
    optoUnitList = unitsInfo(2:end,8);
    subFolders = unitsInfo(2:end,9);
    depths = nums(:, 9);

    drift = unitsInfo(2:end,10);
    drift = contains(drift, 'drift')|contains(drift, 'Drift')|contains(drift, 'duplicate')|contains(drift, 'Duplicate');
    quality = nums(:,1);
    sessionList = sessionList(quality<=0.05 & ~drift);
    unitList = unitList(quality<=0.05 & ~drift);
    depths = depths(quality<=0.05 & ~drift);
    optoUnitList = optoUnitList(quality<=0.05 & ~drift);
    subFolders = subFolders(quality<=0.05 & ~drift);
    
    % load good days
    [~ , goodDayList, ~] = xlsread([root xlFile], animalName);
    [~,col] = find(~cellfun(@isempty,strfind(goodDayList, category)) == 1);
    goodDayList = goodDayList(2:end,col);
    loadedSession = [];
    
%% session and unit loop
    for ses = 1:length(sessionList)
        session = sessionList{ses};
        unit = unitList{ses};
        fprintf([session unit '\n']);
    % paths
        pd = parseSessionString_df(session, root, sep);
        neuralynxDataPath = [pd.sortedFolder session '_sessionData_nL.mat'];
        unitMetDir = [pd.nLynxFolderSession session '_' unit '_met.mat'];
        sortedFolderLocation = [pd.sortedFolder 'session' sep];
    % decide is good behavior
    if ~ismember(session, goodDayList)
        continue
    end
    % decide if a good unit
    if exist(unitMetDir,'file')
        load(unitMetDir);
    else
        met = getClusterMetric(session, unit, 0, 1);
    end
    respInds = find(met.spikeProp>=0.8);
    if isempty(respInds)
        continue
    else 
        respLat = nanmin(met.spikeLat(respInds));
        if respLat > 15000 || met.isiV > 0.001 || met.distance>0.3 || isnan(met.distance)
            continue
        end
    end
    % load behavior and neurons
    if exist(neuralynxDataPath,'file')
        load(neuralynxDataPath)
    else
        sessionData = generateSessionData_nL_operantMatching(session);
    end
    % find position
    populationDepths = [populationDepths depths(ses)];
    if ~strcmp(session,loadedSession) % avoiding recomputing model variables for different unit from same session
        %% behavior preparation 
        % parse behavior
        os = behAnalysisNoPlot_opMD(session);
        if length(os.CSplus) ~= length(sessionData)
            fprintf([session ' error\n']);
            return
        end
        choice = os.allChoices';
        choice(choice<0) = 0;
        outcome = abs(os.allRewards)';
        choice = choice(1:min(length(choice), maxTrial));
        outcome = outcome(1:length(choice));
        responseInds = os.responseInds(1:min(length(choice), maxTrial)); 
        preRwd = [NaN abs(os.allRewards(1:end-1))]';
        %% behavior
        % switch
        svs = zeros(length(os.responseInds),1);
        svs(os.changeChoice_Inds) = 1;
        svsNext = [svs(2:end); NaN];
        
        % choice kernel
        choiceKernel = conv(os.allChoices, myKernel);
        choiceKernel = abs(choiceKernel(0.5*length(myKernel)-0.5+1:end-(0.5*length(myKernel)-0.5)))'/sum(myKernel);
        

        hmm = double(os.hmmStates==1)';
        lickLat = os.lickLatLogZ';
        rightSide = zeros(size(outcome));
        rightSide(os.allChoices>0)=1;
        preITI = os.timeBtwn';
  
        tbl = table(outcome, preRwd, rightSide, lickLat, hmm, preITI, svsNext, choiceKernel);
        names = tbl.Properties.VariableNames;
        % zscore all regressors
        for cols = 1:length(names)
            tmp = tbl.(names{cols});
            if ~isempty(setdiff(tmp(~isnan(tmp)), [0 1 -1 NaN]))
                tmp(~isnan(tmp)) = zscore(tmp(~isnan(tmp)));
                tbl.(names{cols}) = tmp;
            end
        end
        %% table for gng
        % correct (hit CR) and wrong (miss, FA)
        correct = zeros(length(os.behSessionData),1);
        correct(setxor(find(os.CSplus>0), find(~isnan(os.lickSide)>0))) = 1;
        
        loadedSession = session;
    end
    %% 
    % load metric
    baselineFreq = [baselineFreq met.baseline];
    firstSpikeFreq = [firstSpikeFreq 1000*mean(met.spikeNum(:,1))/20];
    firstSpikeLat = [firstSpikeLat met.spikeLat(1)/1000];
    width = [width metSess.width];
    spikeLatChange = [spikeLatChange (met.spikeLat(10)-met.spikeLat(1))/1000];
    spikeFreqChange = [spikeFreqChange 1000*(mean(met.spikeNum(:,10))-mean(met.spikeNum(:,1)))/met.pulseWidth];
%             spikeFreqChange = [spikeFreqChange mean(met.spikeNum(:,10))-mean(met.spikeNum(:,1))];
    [peak,peakCh] = max(max(met.optoWaveform));
    [~,peakSamp] = max(met.optoWaveform(:,peakCh));
     tempWaveform = [zeros(1, 16-peakSamp), met.optoWaveform(max(1,peakSamp-15):min(peakSamp+25,size(met.optoWaveform,1)),peakCh)', zeros(1, peakSamp+25 - size(met.optoWaveform,1))];    
    waveforms = [waveforms; tempWaveform/peak];
    [peak,~] = max(max(metSess.waveform));
    [~,peakSamp] = max(metSess.waveform(:,peakCh));
    tempWF = metSess.waveform(:,peakCh);
    tempWF = tempWF(peakSamp-30:peakSamp+40);
    waveformsSession = [waveformsSession; tempWF'/peak];
    %% count sessions and units
    allSessions = [allSessions, session];
    allUnits = [allUnits, unit];
    
    %% create spike and lick cell aligned to first lick
    
    spikeFields = fields(sessionData);
    clust = find(contains(spikeFields,unit));
    allTrial_spike_choice = {};

    for k = 1:length(os.responseInds)
            if os.responseInds(k) == 1
                prevTrial_spike = [];
            else
                prevTrial_spikeInd = [sessionData(os.responseInds(k)-1).(spikeFields{clust})] > (sessionData(os.responseInds(k)).respondTime-p.Results.tb*1000);
                prevTrial_spike = sessionData(os.responseInds(k)-1).(spikeFields{clust})(prevTrial_spikeInd) - sessionData(os.responseInds(k)).respondTime;
            end

            currTrial_spikeInd = sessionData(os.responseInds(k)).(spikeFields{clust}) < sessionData(os.responseInds(k)).respondTime+p.Results.tf*1000 ... 
                & sessionData(os.responseInds(k)).(spikeFields{clust}) > sessionData(os.responseInds(k)).respondTime-p.Results.tb*1000;
            currTrial_spike = sessionData(os.responseInds(k)).(spikeFields{clust})(currTrial_spikeInd) - sessionData(os.responseInds(k)).respondTime;

            allTrial_spike_choice{k} = [prevTrial_spike currTrial_spike];
    end

    % sometimes no licks/spikes are considered 1x0 and sometimes they are []
    % plotSpikeRaster does not place nicely with [] so this converts all empty indices to 1x0
    allTrial_spike_choice(cellfun(@isempty,allTrial_spike_choice)) = {zeros(1,0)}; 
    % spike matrix choice
    allTrial_spikeMatx_choice = zeros(length(os.responseInds),length(time));
    trialDurDiff = zeros(1, length(os.responseInds));
    for j = 1:length(os.responseInds)
        trialDurDiff(j) = (sessionData(os.responseInds(j)).trialEnd - sessionData(os.responseInds(j)).CSon)- p.Results.tf*1000;
    end
   for j = 1:length(allTrial_spike_choice)
        tempSpike = allTrial_spike_choice{j};
        tempSpike = tempSpike + p.Results.tb*1000; % add this to pad time for SDF
        if any(tempSpike == 0)
            tempSpike(tempSpike == 0) = 1;
        end
        allTrial_spikeMatx_choice(j,tempSpike) = 1;
        if trialDurDiff(j) < 0
            allTrial_spikeMatx_choice(j, isnan(allTrial_spikeMatx_choice(j, 1:end+trialDurDiff(j)))) = 0;  %converts within trial duration NaNs to 0's
        else
            allTrial_spikeMatx_choice(j, isnan(allTrial_spikeMatx_choice(j,:))) = 0;
        end
   end
  % slide window
    allTrial_spikeMatx_slide = zeros(length(os.responseInds), length(midPoints));
    for w = 1:length(midPoints)
        allTrial_spikeMatx_slide(:,w) = ...
            nansum(allTrial_spikeMatx_choice(:,midPoints(w)-0.5*p.Results.binSize:midPoints(w)+0.5*p.Results.binSize-1),2)*1000/p.Results.binSize;
    end
    
%     allTrial_spikeMatx_slide = zscore(allTrial_spikeMatx_slide, 0, 'all');
    %%  create spike cell aligned to cue
    allTrial_spike_cue = {};

    for k = 1:length(os.behSessionData)
            if k == 1
                prevTrial_spike = [];
            else
                prevTrial_spikeInd = [sessionData(k-1).(spikeFields{clust})] > (sessionData(k).CSon-p.Results.tb*1000);
                prevTrial_spike = sessionData(k-1).(spikeFields{clust})(prevTrial_spikeInd) - sessionData(k).CSon;
            end

            currTrial_spikeInd = sessionData(k).(spikeFields{clust}) < sessionData(k).CSon+p.Results.tf*1000 ... 
                & sessionData(k).(spikeFields{clust}) > sessionData(k).CSon-p.Results.tb*1000;
            currTrial_spike = sessionData(k).(spikeFields{clust})(currTrial_spikeInd) - sessionData(k).CSon;

            allTrial_spike_cue{k} = [prevTrial_spike currTrial_spike];
    end

    % sometimes no licks/spikes are considered 1x0 and sometimes they are []
    % plotSpikeRaster does not place nicely with [] so this converts all empty indices to 1x0
    allTrial_spike_cue(cellfun(@isempty,allTrial_spike_cue)) = {zeros(1,0)}; 
    %  spike matrix cue 
    allTrial_spikeMatx_cue = zeros(length(os.behSessionData),length(time));
    trialDurDiff = zeros(1, length(os.behSessionData));
    for j = 1:length(os.behSessionData)
        trialDurDiff(j) = (sessionData(j).trialEnd - sessionData(j).CSon)- p.Results.tf*1000;
    end
   for j = 1:length(allTrial_spike_cue)
        tempSpike = allTrial_spike_cue{j};
        tempSpike = tempSpike + p.Results.tb*1000; % add this to pad time for SDF
        if any(tempSpike == 0)
            tempSpike(tempSpike == 0) = 1;
        end
        allTrial_spikeMatx_cue(j,tempSpike) = 1;
        if trialDurDiff(j) < 0
            allTrial_spikeMatx_cue(j, isnan(allTrial_spikeMatx_cue(j, 1:end+trialDurDiff(j)))) = 0;  %converts within trial duration NaNs to 0's
        else
            allTrial_spikeMatx_cue(j, isnan(allTrial_spikeMatx_cue(j,:))) = 0;
        end
   end
  % slide window
    allTrial_spikeMatx_cue = allTrial_spikeMatx_cue(os.responseInds,:);
    allTrial_spikeMatx_slide_cue = zeros(length(os.responseInds), length(midPoints));
    for w = 1:length(midPoints)
        allTrial_spikeMatx_slide_cue(:,w) = ...
            nansum(allTrial_spikeMatx_cue(:,midPoints(w)-0.5*p.Results.binSize:midPoints(w)+0.5*p.Results.binSize-1),2)*1000/p.Results.binSize;
    end
    
%     allTrial_spikeMatx_slide_cue = zscore(allTrial_spikeMatx_slide_cue, 0, 'all');
    %%
    % focus window
    spikeCounts = nansum(allTrial_spikeMatx_choice(:,p.Results.tb*1000+p.Results.focusWin(1):p.Results.tb*1000+p.Results.focusWin(2)),2);
    spikeCountsCue = nansum(allTrial_spikeMatx_cue(:,p.Results.tb*1000+p.Results.focusWin(1):p.Results.tb*1000+p.Results.focusWin(2)),2);

    % calculate and combine rwdvsNrwd
    rwdvsNrwd = [rwdvsNrwd; [mean(spikeCounts(os.rwd_Inds)), mean(spikeCounts(os.nrwd_Inds))]];
    % baseline 
    baseline = nansum(allTrial_spikeMatx_cue(:,1:p.Results.tb*1000),2)/p.Results.tb;
    % go-cue
    goSpikeCounts = nansum(allTrial_spikeMatx_cue(:,p.Results.tb*1000:p.Results.tb*1000+300),2)/0.3;
    % excitability
    excit = [excit; [mean(baseline), mean(goSpikeCounts)]];  
    % pre session baseline
    % find recording start time: 
    [timestamps, eventID, TTL, Evstring] = Nlx2MatEV([pd.nLynxFolderSession 'Events.nev'],[1 1 1 0 1], 0, 1);
    timestamps = round(timestamps/1000);
    clustAll = clust - find(contains(spikeFields,'allSpikes'));
    startTime = min(timestamps);
    currpreSbl = sessionData(clustAll).allSpikes;
    currpreSbl = currpreSbl(currpreSbl<sessionData(1).CSon & currpreSbl>(sessionData(1).CSon - 20000)); % within 20s before session started
    currpreSbl = length(currpreSbl)/min(20000, sessionData(1).CSon - startTime) * 1000;
    preSessBaseline = [preSessBaseline currpreSbl];
    %% put psths together
    % psth of all trials aligned to go cue and then rwd and nrwd after rwd
    meanCue = mean(allTrial_spikeMatx_slide_cue(:, slideTime >= (-1)*p.Results.tb*1000 & slideTime <= p.Results.tfCue*1000));
    meanRwd = mean(allTrial_spikeMatx_slide(os.rwd_Inds, slideTime>= os.rwdDelay & slideTime<=os.rwdDelay + 1000*p.Results.tfRwd));
    meanNrwd = mean(allTrial_spikeMatx_slide(os.nrwd_Inds, slideTime>= os.rwdDelay & slideTime<=os.rwdDelay + 1000*p.Results.tfRwd));
    
    meanComb = zeros(1, length(meanCue) + length(meanRwd) + length(meanNrwd));
    meanComb = [meanCue, meanRwd, meanNrwd];
    meanComb = meanComb/mean(meanCue(1:length(find((slideTime+0.5*p.Results.binSize)<=0))));
    psthComb = [psthComb; meanComb];
%%
    % zscore 
    spikeCounts = zscore(nansum(allTrial_spikeMatx_choice(:,p.Results.tb*1000+p.Results.focusWin(1):p.Results.tb*1000+p.Results.focusWin(2)),2));
    %%
    % calculate regression results.
    currTbl = addvars(tbl, spikeCounts);
    lm = fitlm(currTbl, ['spikeCounts~1+outcome+rightSide']);
    catRegressorInd = find(strcmp(lm.CoefficientNames, 'outcome'));    
    
    if lm.Coefficients.pValue(catRegressorInd)>0.05
       cats = [cats 0]; 
    else
        if lm.Coefficients.Estimate(catRegressorInd)>0
            cats = [cats 2];
%             spikeGLM_dF(session, 'good','cellName', unit,'regressors','1  + Qchosen + outcome + rightSide');
%             spikeGLM_dF(session, 'good','cellName', unit,'regressors','1  + Qchosen + outcome*rightSide');
        else
            cats = [cats 1];
        end
    end
     
    end 
end
% PCs and features
% m =  mean(waveforms,1);
% waveforms = waveforms - m;
numCat = p.Results.numCat;
[coeff,score,latent, ~, explained, mu] = pca(zscore(waveformsSession, [], 1));
excitRatio = (excit(:,2)-excit(:,1))./excit(:,1);
matrix = [zscore(excit(:,1)) zscore(excit(:,2)) width' excitRatio score(:,1) score(:,2)];
names = {'baselineFreq', 'respFreq', 'width', 'excitRat  io', 'PC1', 'PC2'};
matrix = zscore(matrix,0,1);
[coeffAll,scoreAll,latentAll, ~, explainedAll, muAll] = pca(matrix);
%% plot everything
figure2;
numCat = 2;
ind = cats;
colors = cool(numCat);
subplot(5, 4, [5 9]); hold on;
scatter3(scoreAll(:,1), scoreAll(:,2), scoreAll(:,3), 15, [0.7 0.7 0.7], 'filled');
for i = 1:numCat
    scatter3(scoreAll(ind==i,1), scoreAll(ind==i,2), scoreAll(ind==i,3), 12, colors(i,:), 'filled');
end
xlabel('PC1')
ylabel('PC2')
zlabel('PC3')
subplot(5,size(matrix,2)+2, 1:2); hold on;
for i = 1:numCat
    mat = waveformsSession(ind==i,:);
    plotFilled(1:size(mat,2), mat, colors(i,:));
end
title([num2str(sum(ind==1)) '/' num2str(sum(ind==2))]);

for i = 1:size(matrix,2)
    subplot(5, size(matrix,2)+2, i+2);
    hold on
    edges = linspace(min(matrix(:,i))-0.001, max(matrix(:,i)+0.001),10);
    for j = 1:numCat
        histogram(matrix(ind==j,i), edges, 'FaceColor', colors(j,:), 'Normalization', 'probability');
    end
    title(names{i})
end
% sp-pe
subplot(5, 4, [15, 19]); hold on;
for i = 1:numCat
    currSpikes = focusSpikes(ind==i,:);
    currTarget = allPe(ind==i, :);
    spikeMeans = mean(currSpikes, 'omitnan');
    spikeSems = sem(currSpikes);
    targetMeans = mean(currTarget, 'omitnan');

    plot(targetMeans, spikeMeans, 'color', colors(i,:), 'linewidth', 2, 'Marker', 'none', 'LineStyle', '-');
    fill([targetMeans fliplr(targetMeans)], [spikeMeans+spikeSems fliplr(spikeMeans-spikeSems)], colors(i,:), 'facealpha', 0.25, 'edgecolor', 'none')
end
xlabel('pe');
ylabel('zscored(spikes/s)')


colorsG = zeros(6,3);
for c = 1:3
    colorsG(c,:) = [1 1.2-c*2/5 1.2-c*2/5];
    colorsG(c+3,:) = [(c-1)*2/5 (c-1)*2/5 1];
end

subplot(5, 4, [16]); hold on;
currAllSpikes = allSpikes(:,:,ind==1);
for k = 1:size(currAllSpikes,1)
    spikesTemp = squeeze(currAllSpikes(k,:,:));% take all spikes from all neurons from first percental
    plotFilled(slideTime, spikesTemp', colorsG(k,:));
end
title('neg')
subplot(5, 4, [20]); hold on;
currAllSpikes = allSpikes(:,:,ind==2);
for k = 1:size(currAllSpikes,1)
    spikesTemp = squeeze(currAllSpikes(k,:,:));% take all spikes from all neurons from first percental
    plotFilled(slideTime, spikesTemp', colorsG(k,:));
end
xlabel('time from lick');
title('pos')
%% SW
subplot(5, 4, [7, 11]); hold on;
for i = 1:numCat
    currSpikes = focusSpikesBins(ind==i,:);
    currTarget = allSW(ind==i, :);
    spikeMeans = mean(currSpikes, 'omitnan');
    targetSems = sem(currSpikes);
    targetMeans = mean(currTarget, 'omitnan');

    plot(spikeMeans, targetMeans, 'color', colors(i,:), 'linewidth', 2, 'Marker', 'none', 'LineStyle', '-');
    fill([spikeMeans fliplr(spikeMeans)], [targetMeans+targetSems fliplr(targetMeans-targetSems)], colors(i,:), 'facealpha', 0.25, 'edgecolor', 'none')
end
ylabel('p(sw)');
xlabel('zscored(spikes/s)')
% switchy
subplot(5, 4, [8, 12]); hold on;
for i = 1:numCat
    currSpikes = focusSpikesBins(ind==i,:);
    currTarget = allKernel(ind==i, :);
    spikeMeans = mean(currSpikes, 'omitnan');
    targetSems = sem(currSpikes);
    targetMeans = mean(currTarget, 'omitnan');

    plot(spikeMeans, targetMeans, 'color', colors(i,:), 'linewidth', 2, 'Marker', 'none', 'LineStyle', '-');
    fill([spikeMeans fliplr(spikeMeans)], [targetMeans+targetSems fliplr(targetMeans-targetSems)], colors(i,:), 'facealpha', 0.25, 'edgecolor', 'none')
end
ylabel('choiceEnt');
xlabel('Spikes s^{-1} (zscored)')
%% postition
subplot(5, 4, [13, 17]); hold on;
edges = linspace(min(populationDepths)-0.000001, max(populationDepths)+0.000001, 15);
for i = 1:numCat
    histogram(populationDepths(ind==i), edges, 'FaceColor', colors(i,:), 'Normalization', 'probability');
end
xlabel('depth')

subplot(5, 4, [14, 18]); hold on;
edges = linspace(min(populationDepths)-0.000001, max(populationDepths)+0.000001, 15);
ratio = NaN(numCat, length(edges)-1);
meanDepth = NaN(numCat, length(edges)-1);
for i = 1:numCat
    tempDepths = populationDepths(ind==i);
    for j = 1:length(edges)-1
        meanDepth(i,j) = mean(edges(j:j+1));
        ratio(i, j) = sum(tempDepths>=edges(j)&tempDepths<edges(j+1))/sum(populationDepths>=edges(j)&populationDepths<edges(j+1));
    end
    plot(meanDepth(i,:), ratio(i, :), 'color', colors(i,:), 'linewidth', 2, 'Marker', 'none', 'LineStyle', '-');
end
xlabel('depth (mm)')
ylabel('ratio')
 
%% clusters
subplot(5, 4, [6, 10]); hold on;
scatter(score(:,1), score(:,2), 20, [0.7 0.7 0.7], 'filled');
for i = 1:numCat
    scatter(score(ind==i,1), score(ind==i,2), 15, colors(i,:), 'filled');
end

xlabel('PC1')
ylabel('width')

suptitle([p.Results.regressors '-' p.Results.catRegressor]);
screen = get(0,'Screensize');
screen(4) = screen(4) - 100;
set(gcf, 'Position', screen)

figure;
subplot(3, 4, 1); hold on;
scatter(excit(:,1), excit(:,2), 18, [0.7 0.7 0.7], 'filled');
for i = 1:numCat
    scatter(excit(ind==i,1), excit(ind==i,2), 15, colors(i,:), 'filled');
end
plot([0 0], [0 0.3], 'LineStyle', '--', 'Color', [0.2 0.2 0.2], 'LineWidth', 2);
xlabel('baseline')
ylabel('evoke')

subplot(3, 4, 2); hold on;
edges = linspace(min((excit(:,2)-excit(:,1))./excit(:,1))-0.00001, max((excit(:,2)-excit(:,1))./excit(:,1))+0.00001, 15);
for i = 1:numCat
    histogram((excit(ind==i,2)-excit(ind==i,1))./excit(ind==i,1), edges, 'FaceColor', colors(i,:), 'FaceAlpha', 0.4, 'Normalization', 'probability');
end
plot([0 0], [0 0.3], 'LineStyle', '--', 'Color', [0.2 0.2 0.2], 'LineWidth', 2);

xlabel('(evoke-baseline)/baseline')

subplot(3, 4, 3); hold on;
edges = linspace(min(excit(:,1))-0.00001, max(excit(:,1))+0.00001, 15);
for i = 1:numCat
    histogram(excit(ind==i,1), edges, 'FaceColor', colors(i,:), 'FaceAlpha', 0.4, 'Normalization', 'probability');
end
plot([0 0], [0 0.3], 'LineStyle', '--', 'Color', [0.2 0.2 0.2], 'LineWidth', 2);

xlabel('bl')

subplot(3, 4, 4); hold on;
edges = linspace(min(excit(:,2))-0.00001, max(excit(:,2))+0.00001, 15);
for i = 1:numCat
    histogram(excit(ind==i,2), edges, 'FaceColor', colors(i,:), 'FaceAlpha', 0.4, 'Normalization', 'probability');
end
plot([0 0], [0 0.3], 'LineStyle', '--', 'Color', [0.2 0.2 0.2], 'LineWidth', 2);

xlabel('evoked')

subplot(3,4,5); hold on;
for i = 1:numCat
    scatter(excit(ind==i, 1), preSessBaseline(ind==i), 12, colors(i,:), 'filled');
end
plot([0 8], [0 8], 'Color', [0.7 0.7 0.7], 'LineStyle', '--');
xlabel(['preCueBaseline'])
ylabel(['preSessionBaseline'])

subplot(3,4,6); hold on;
edges = linspace(min(preSessBaseline)-0.00001, max(preSessBaseline)+0.00001, 15);
for i = 1:numCat
    histogram(preSessBaseline(ind==i), edges, 'FaceColor', colors(i,:), 'FaceAlpha', 0.4, 'Normalization', 'probability');
end
xlabel('preSessBaseline')

screen = get(0,'Screensize');
screen(4) = screen(4) - 100;
set(gcf, 'Position', screen)

figure;
subplot(2,2,1);
hold on;
for i = 1:numCat
    currSpikes = focusSpikesRwdBins(ind==i,:);
    currTarget = allSWRwd(ind==i, :);
    spikeMeans = mean(currSpikes, 'omitnan');
    targetSems = sem(currSpikes);
    targetMeans = mean(currTarget, 'omitnan');

    plot(spikeMeans, targetMeans, 'color', colors(i,:), 'linewidth', 2, 'Marker', 'none', 'LineStyle', '-');
    fill([spikeMeans fliplr(spikeMeans)], [targetMeans+targetSems fliplr(targetMeans-targetSems)], colors(i,:), 'facealpha', 0.25, 'edgecolor', 'none')
end
ylabel('p(SW|rwd)');
xlabel('Spikes s^{-1} (zscored)')

title('P(SW|rwd)')

subplot(2,2,2);
hold on;
for i = 1:numCat
    currSpikes = focusSpikesRwdBins(ind==i,:);
    currTarget = allKernelRwd(ind==i, :);
    spikeMeans = mean(currSpikes, 'omitnan');
    targetSems = sem(currSpikes);
    targetMeans = mean(currTarget, 'omitnan');

    plot(spikeMeans, targetMeans, 'color', colors(i,:), 'linewidth', 2, 'Marker', 'none', 'LineStyle', '-');
    fill([spikeMeans fliplr(spikeMeans)], [targetMeans+targetSems fliplr(targetMeans-targetSems)], colors(i,:), 'facealpha', 0.25, 'edgecolor', 'none')
end
ylabel('choice entropy rwd');
xlabel('Spikes s^{-1} (zscored)')
title('choiceEnt|rwd')

subplot(2,2,3);
hold on;
for i = 1:numCat
    currSpikes = focusSpikesNrwdBins(ind==i,:);
    currTarget = allSWNrwd(ind==i, :);
    spikeMeans = mean(currSpikes, 'omitnan');
    targetSems = sem(currSpikes);
    targetMeans = mean(currTarget, 'omitnan');

    plot(spikeMeans, targetMeans, 'color', colors(i,:), 'linewidth', 2, 'Marker', 'none', 'LineStyle', '-');
    fill([spikeMeans fliplr(spikeMeans)], [targetMeans+targetSems fliplr(targetMeans-targetSems)], colors(i,:), 'facealpha', 0.25, 'edgecolor', 'none')
end
ylabel('p(SW|no rwd)');
xlabel('Spikes s^{-1} (zscored)')
title('P(SW)|no rwd')
subplot(2,2,4);
hold on;
for i = 1:numCat
    currSpikes = focusSpikesNrwdBins(ind==i,:);
    currTarget = allKernelNrwd(ind==i, :);
    spikeMeans = mean(currSpikes, 'omitnan');
    targetSems = sem(currSpikes);
    targetMeans = mean(currTarget, 'omitnan');

    plot(spikeMeans, targetMeans, 'color', colors(i,:), 'linewidth', 2, 'Marker', 'none', 'LineStyle', '-');
    fill([spikeMeans fliplr(spikeMeans)], [targetMeans+targetSems fliplr(targetMeans-targetSems)], colors(i,:), 'facealpha', 0.25, 'edgecolor', 'none')
end
ylabel('choice entropy nrwd');
xlabel('Spikes s^{-1} (zscored)')
title('choiceEnt no rwd')
screen = get(0,'Screensize');
screen(4) = screen(4) - 100;
set(gcf, 'Position', screen)


















