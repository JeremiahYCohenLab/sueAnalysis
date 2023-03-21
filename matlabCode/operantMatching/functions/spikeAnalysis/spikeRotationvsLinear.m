function spikeRotationvsLinear(animalNames, category, varargin)
%task and model parameters
p = inputParser;
% default parameters if none given
p.addParameter('numCat', 2);
p.addParameter('focusWin', [301 1300]); % in ms from respond/cue
p.addParameter('binSize', 100)% in ms
p.addParameter('tb', 2)% in s
p.addParameter('tf', 3)% in s
p.addParameter('stepSize', 100)
p.addParameter('sepOutcome', 1);
p.parse(varargin{:});
numBins = 4;
% basic info
[root, sep] = currComputer();
% time window
time = -1000*p.Results.tb:1000*p.Results.tf;
midPoints = (0.5*p.Results.binSize + 1):p.Results.stepSize:(length(time)-0.5*p.Results.binSize);
slideTime = midPoints - p.Results.tb*1000;
%% classification
width = [];
waveforms = [];
waveformsSession = [];
allSessions = {};
allUnits = {};
allGoRaw = {};
allBaseline = [];
%% combined information
populationDepths = []; % all estimated positions
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
        sortedFolderLocation = [pd.sortedFolder];
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
    %% 
    % load metric
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
    
    %% create spike and lick cell
    spikeFields = fields(sessionData);
    clust = find(contains(spikeFields,unit));
    allTrial_spike_choice = {};
    os = behAnalysisNoPlot_opMD(session, 'simpleFlag', 1);
    if length(os.CSplus)~=length(sessionData)
        return
    end
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
    
    allTrial_spikeMatx_slide = zscore(allTrial_spikeMatx_slide, 0, 'all');
    %%  create spike cell aligned to cue
    allTrial_spike_cue = {};

    for k = 1:length(sessionData)
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
    allTrial_spikeMatx_cue = zeros(length(sessionData),length(time));
    trialDurDiff = zeros(1, length(sessionData));
    for j = 1:length(sessionData)
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
    allTrial_spikeMatx_slide_cue = zeros(length(sessionData), length(midPoints));
    for w = 1:length(midPoints)
        allTrial_spikeMatx_slide_cue(:,w) = ...
            nansum(allTrial_spikeMatx_cue(:,midPoints(w)-0.5*p.Results.binSize:midPoints(w)+0.5*p.Results.binSize-1),2)*1000/p.Results.binSize;
    end
    
    allTrial_spikeMatx_slide_cue = zscore(allTrial_spikeMatx_slide_cue, 0, 'all');
    %%
    % focus window
    spikeCounts = nansum(allTrial_spikeMatx_choice(:,p.Results.tb*1000+p.Results.focusWin(1):p.Results.tb*1000+p.Results.focusWin(2)),2);
    spikeCountsCue = nansum(allTrial_spikeMatx_cue(:,p.Results.tb*1000+p.Results.focusWin(1):p.Results.tb*1000+p.Results.focusWin(2)),2);
    % baseline 
    baseline = nansum(allTrial_spikeMatx_cue(:,1:p.Results.tb*1000),2)/p.Results.tb;
    allBaseline{length(allSessions)} = baseline;
    % go-cue
    goSpikeCounts = nansum(allTrial_spikeMatx_cue(:,p.Results.tb*1000:p.Results.tb*1000+300),2)/0.3;
    allGoRaw{length(allSessions)} = goSpikeCounts;
    % excitability
    
    % fit two tables linear or rotation
    increase = goSpikeCounts - baseline;
    intercept = ones(size(increase));
    tbl = table(goSpikeCounts, increase, baseline, intercept);
    
    lmLinear = fitlm(tbl, 'increase ~ intercept','Intercept',false);
    lmRotation = fitlm(tbl, 'increase ~ baseline','Intercept',false);
    lmAll = fitlm(tbl,  ' goSpikeCounts ~ baseline + intercept','Intercept',false);
    
    aicL(length(allSessions),:) = [lmLinear.ModelCriterion.AIC];
    aicR(length(allSessions),:) = [lmRotation.ModelCriterion.AIC];
    aicAll(length(allSessions),:) = [lmAll.ModelCriterion.AIC];
    RsqL(length(allSessions),:) = [lmLinear.Rsquared];
    RsqR(length(allSessions),:) = [lmRotation.Rsquared];
    RsqAll(length(allSessions),:) = [lmAll.Rsquared];
    bicL(length(allSessions),:) = [lmLinear.ModelCriterion.BIC];
    bicR(length(allSessions),:) = [lmRotation.ModelCriterion.BIC];
    bicAll(length(allSessions),:) = [lmAll.ModelCriterion.BIC];
    coeffL(length(allSessions),:) = [lmLinear.Coefficients.Estimate];
    coeffR(length(allSessions),:) = [lmRotation.Coefficients.Estimate];
    coeffAll(length(allSessions),:) = [lmAll.Coefficients.Estimate];
    statsL(length(allSessions),:) = [lmLinear.Coefficients.tStat];
    statsR(length(allSessions),:) = [lmRotation.Coefficients.tStat];
    statsAll(length(allSessions),:) = [lmAll.Coefficients.tStat];
    rmseL(length(allSessions),:) = [lmLinear.RMSE];
    rmseR(length(allSessions),:) = [lmRotation.RMSE];
    rmseAll(length(allSessions),:) = [lmAll.RMSE];
    llL(length(allSessions),:) = [lmLinear.LogLikelihood]/length(sessionData);
    llR(length(allSessions),:) = [lmRotation.LogLikelihood]/length(sessionData);
    llAll(length(allSessions),:) = [lmAll.LogLikelihood]/length(sessionData);
    sigL(length(allSessions),:) = [lmLinear.Coefficients.pValue];
    sigR(length(allSessions),:) = [lmRotation.Coefficients.pValue];
    sigAll(length(allSessions),:) = [lmAll.Coefficients.pValue];
     
    % bin spikeNum by baseline

    edges = linspace(min(baseline)-0.001, max(baseline)+0.001, numBins+1);
    tempbl = zeros(1,numBins);
    tempGo = zeros(1,numBins);
    tempGoSem = zeros(1,numBins);
    for b = 1:numBins
        tempbl(b) = mean(baseline(baseline>=edges(b) & baseline<edges(b+1)));
        tempGo(b) = mean(goSpikeCounts(baseline>=edges(b) & baseline<edges(b+1)));
        tempGoSem(b) = sem(goSpikeCounts(baseline>=edges(b) & baseline<edges(b+1)));
    end
    allbl(length(allSessions),:) = tempbl;
    allGo(length(allSessions),:) = tempGo;
    allGoSem(length(allSessions),:) = tempGoSem;
    end 
end
% PCs and features
% m =  mean(waveforms,1);
% waveforms = waveforms - m;
numCat = p.Results.numCat;
[coeff,scores,latent, ~, explained, mu] = pca(zscore(waveformsSession, [], 1));
indAll = {};
dis = {};
for a = 1:10
    [indAll{a}, ~, dis{a}] = kmeans([scores(:, 1:5), width'], numCat);
end
[~,optiInds] = min(cellfun(@mean, dis));
ind = indAll{optiInds};
%% plot everything
% compare models
figure2;
subplot(1,2,1); hold on;
scatter(bicL, bicR);
plot([min(bicL) max(bicL)], [min(bicL) max(bicL)], 'LineStyle', '--', 'LineWidth',2,'color', 'r');
xlabel('Linear');
ylabel('Rotation');
subplot(1,2,2); hold on;
scatter(bicL, bicAll);
plot([min(bicL) max(bicL)], [min(bicL) max(bicL)], 'LineStyle', '--', 'LineWidth',2,'color', 'r');
xlabel('Linear');
ylabel('All');
sgtitle('BIC')

figure2;
subplot(1,2,1); hold on;
scatter(llL, llR);
plot([min(llL) max(llL)], [min(llL) max(llL)], 'LineStyle', '--', 'LineWidth',2,'color', 'r');
xlabel('Linear');
ylabel('Rotation');
subplot(1,2,2); hold on;
scatter(llL, llAll);
plot([min(llL) max(llL)], [min(llL) max(llL)], 'LineStyle', '--', 'LineWidth',2,'color', 'r');
xlabel('Linear');
ylabel('All');
sgtitle('log likelihood')
%% plot coeffs
figure2;
subplot(1,2,1); hold on;
histogram(statsAll(:,1));
title('baseline');
subplot(1,2,2); hold on;
histogram(statsAll(:,2));
title('intercept');
sgtitle('tstats')

color1 = [0 0.8 0.8];
color2 = [1 0 1];
figure2;
subplot(1,2,1); hold on;
histogram(coeffAll(ind == 1,1), linspace(min(coeffAll(:,1))-0.01, max(coeffAll(:,1))+0.01, 9), 'FaceColor', color1, 'Normalization', 'probability'); 
histogram(coeffAll(ind == 2,1), linspace(min(coeffAll(:,1))-0.01, max(coeffAll(:,1))+0.01, 9),'FaceColor', color2, 'Normalization', 'probability'); 
title('baseline');
subplot(1,2,2); hold on;
histogram(coeffAll(ind == 1,2), linspace(min(coeffAll(:,2))-0.01, max(coeffAll(:,2))+0.01, 9),'FaceColor', color1, 'Normalization', 'probability'); 
histogram(coeffAll(ind == 2,2), linspace(min(coeffAll(:,2))-0.01, max(coeffAll(:,2))+0.01, 9),'FaceColor', color2, 'Normalization', 'probability');  
title('intercept');
sgtitle('coeffs')

figure2;hold on;
plot(allbl', allGo', 'Color', [0.7 0.7 0.7]);
plotFilled(nanmean(allbl(ind==1,:)), allGo(ind==1,:), color1)
plotFilled(nanmean(allbl(ind==2,:)), allGo(ind==2,:), color2)
ylabel('resp to go cue zscored'); xlabel('baseline zscored')



















