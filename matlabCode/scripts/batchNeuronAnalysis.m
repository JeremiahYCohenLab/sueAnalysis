animalNames = {'ZS059','ZS060','ZS061','ZS062'};
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

    drift = unitsInfo(2:end,10);
    drift = contains(drift, 'drift')|contains(drift, 'Drift') | contains(drift, 'duplicate')|contains(drift, 'Duplicate');
    quality = nums(:,1);
    sessionListTmp = sessionListTmp(quality<=maxLratio & ~drift);
    unitListTmp = unitListTmp(quality<=maxLratio & ~drift);
    optoUnitListTmp = optoUnitListTmp(quality<=maxLratio & ~drift);
    subFoldersTmp = subFoldersTmp(quality<=maxLratio & ~drift);
    
    sessionList = [sessionList; sessionListTmp];
    unitList = [unitList; unitListTmp];
    optoUnitList = [optoUnitList; optoUnitListTmp];
    subFolders = [subFolders; subFoldersTmp];
end
%% get behavior & pupi list
animalNames = {'ZS059', 'ZS060', 'ZS061', 'ZS062'};
[root, sep] = currComputer();
category = 'good';
sessionList = {};
for a = 1:length(animalNames)
    animalName = animalNames{a};
    sheet = animalName;
    xlFile = [animalName '.xlsx'];
    [~, sessionInfo] = xlsread([root xlFile], sheet);
    [~,col] = find(contains(sessionInfo, category));
    sessionInfo = sessionInfo(2:end,col);
    endInd = find(cellfun(@isempty,sessionInfo),1);
    if ~isempty(endInd)
        sessionInfo = sessionInfo(1:endInd-1,:);
    end
    sessionList = [sessionList; sessionInfo];
end

%% quality check
for i = 1:length(sessionList)
    rasters = qualityCheck(sessionList{i},1,'unit',unitList{i});
    saveFigurePDF(rasters, [savePath sep sessionList{i} '_' unitList{i} '.pdf']);
end  
%% opto metrices
savePath = [root animalName sep animalName 'sorted' sep 'optoMetFiles' sep];
for i = 1:length(sessionList)
    getClusterMetric(sessionList{i}, unitList{i}, 0, 1);
%     saveFigurePDF(gcf, [savePath sep sessionList{i} '_' unitList{i} '.pdf']);
end 
%% find longer waveform in CSC file from session recordings
SamplingFreq = 32000;
for i = 132:length(sessionList)
    pd = parseSessionString_df(sessionList{i}, root, sep);
    sortedPath = [pd.nLynxFolderSession];
    load([pd.nLynxFolderSession sessionList{i} '_' unitList{i} '_met.mat']);
    [TTname, optoUnitNum] = strtok(unitList{i}, 'SS');
    TTname = TTname(1:end-1);
    ttNum = str2double(strtok(TTname,'TT'));
    chan = (ttNum*4 -3):ttNum*4;
    % get trace
    [ts, samp0] = Nlx2MatCSC([sortedPath 'CSC' num2str(chan(1)) '.ncs'], [1 0 0 0 1], 0, 1, []);
    [samp1] = Nlx2MatCSC([sortedPath 'CSC' num2str(chan(2)) '.ncs'], [0 0 0 0 1], 0, 1, []);
    [samp2] = Nlx2MatCSC([sortedPath 'CSC' num2str(chan(3)) '.ncs'], [0 0 0 0 1], 0, 1, []);
    [samp3] = Nlx2MatCSC([sortedPath 'CSC' num2str(chan(4)) '.ncs'], [0 0 0 0 1], 0, 1, []);
    samp = cat(3, samp0, samp1, samp2, samp3);
    samp = reshape(samp, [] ,4);
    % get spike times
    sortedFiles = dir(pd.nLynxFolderSession);
    ind = find(contains({sortedFiles.name}, [unitList{i} '.txt']));
    spikeTimes = load([pd.nLynxFolderSession sortedFiles(ind).name]); % in us
    % realign
    TTdir = fullfile(sortedPath,[TTname '.ntt']);
    [tt_ts, tt_sig] = Nlx2MatSpike(TTdir, [1 0 0 0 1], 0, 1, 1);
    allWaveForm = squeeze(tt_sig(:, :, ismember(tt_ts, spikeTimes)));
    meanWaveForm = mean(allWaveForm, 3);
    [~, peakChannel] = max(max(meanWaveForm,[], 1));
    [~, peakTimes] = max(allWaveForm(:,peakChannel,:), [], 1);
    tSamp = 1/SamplingFreq * 1e6; % time per sample in microseconds
    if length(spikeTimes) == length(squeeze(peakTimes))
     spikeTimes = spikeTimes + tSamp*(squeeze(peakTimes) - 11); % change depend on how different were the peak from 11th sample
    end
    
    spikeTimes = round(spikeTimes);
    % rescale time
    
    if length(unique(diff(ts))) == 1 % no pausing
        ts_interp = ts(1):tSamp:ts(1) + tSamp*(size(samp,1) - 1);
    else % if pausing/skip due to data loss, use the proper for loop
        ts_interp = NaN(1, size(samp,1)); 
        for j = 1:length(ts)
            ts_interp(512*(j - 1) + 1:512*(j)) = ts(j):tSamp:ts(j) + tSamp*511;
        end
    end
    
    % find signature channel
    [~, peakChannel] = max(max(met.waveform,[], 1));
    % 
    leg = -100:100; % 101 samples    
    peakInds = find(ismember(floor(ts_interp), spikeTimes)|ismember(ceil(ts_interp), spikeTimes));
    intactInds = (peakInds+leg(1))>0 & (peakInds+leg(end))<size(samp,1);
    peakInds = peakInds(intactInds);
    spikeTimes = spikeTimes(intactInds);
    waveforms = zeros(length(leg),4,length(peakInds));
    % waveform from CSC
    for w = 1:length(leg)
        waveforms(w,:,:) = samp(peakInds + leg(w), :)';
    end
    
    
    
    % find samples not affected by spikes as baseline trace
    
    nospikeSamp = [];
    nospikeSampLate = [];
    for w = 1:20
        tmp = samp(peakInds(2:end)- 80 - w, :)';
        nospikeSamp = [nospikeSamp, tmp];
        tmp = samp(peakInds(1:end-1)+ 80 + w, :)';
        nospikeSampLate = [nospikeSampLate, tmp];
    end
    peakLagAll = zeros(4,1);
    peakEndAll = zeros(4,1);
    tmpH = zeros(length(leg),4);
    tmpHLate = zeros(length(leg),4);

    baseline = mean(nospikeSamp,2);
    baselineLate = mean(nospikeSampLate,2);
    for w = 1:4
        for k = 1:length(leg)
            h = ttest2(squeeze(waveforms(k,w,:)), squeeze(nospikeSamp(w,:)),'Alpha',0.001);
            tmpH(k,w) = h;
            h = ttest2(squeeze(waveforms(k,w,:)), squeeze(nospikeSampLate(w,:)),'Alpha',0.001);
            tmpHLate(k,w) = h;
        end
        tmpHCov = conv(tmpH(:,w), ones(1,8));
        tmpHCov = tmpHCov(8:end);
        tmpHCovLate = conv(tmpHLate(:,w), ones(1,8));
        tmpHCovLate = tmpHCovLate(8:end);
        if ~isempty(find(tmpHCov(1:end-2) == 7 & tmpHCov(2:end-1)>5 & tmpHCov(3:end)>5))
          peakLagAll(w) = min(find(tmpHCov(1:end-2) == 7 & tmpHCov(2:end-1)>5 & tmpHCov(3:end)>5)); % find first continued 5 sig points
          peakEndAll(w) = max(find(tmpHCovLate == 5))+4; % find the end of the spike
        else
            peakLagAll(w) = 51;
            peakEndAll(w) = 51;
        end
        
    end
    
    width = peakEndAll(peakChannel) - peakLagAll(peakChannel);
    metSess = struct;
    metSess.waveform = mean(waveforms,3);
    metSess.width = width;
    save([pd.nLynxFolderSession sessionList{i} '_' unitList{i} '_met.mat'],'metSess','-append');
end
%% recalcualte spike width
SamplingFreq = 32000;
waveform = [];
maxInds = [];
minInds = [];
widthNew = [];
widthOld = [];

for i = 1:length(sessionList)
    pd = parseSessionString_df(sessionList{i}, root, sep);
    sortedPath = [pd.nLynxFolderSession];
    load([pd.nLynxFolderSession sessionList{i} '_' unitList{i} '_met.mat']);
    [~,peakCh] = max(max(metSess.waveform));
    [~,peakSamp] = max(metSess.waveform(:,peakCh));
    tempWF = metSess.waveform(:,peakCh);
    tempWF = tempWF(peakSamp-50:peakSamp+40);
    [~, maxInd] = max(tempWF);
    tempWFAP = tempWF(51:end);   
    [~, minInd] = min(tempWFAP);
    minInd = minInd + 50;
    if minInd < 50+30
        metSess.width = (minInd - maxInd)+1;
    else
        metSess.width = 0;
    end
    waveform = [waveform; tempWF'];
    maxInds = [maxInds, maxInd];
    minInds = [minInds, minInd];
    if isempty(met.width)
        fprintf([sessionList{i}, unitList{i}])
    end
    widthNew = [widthNew, metSess.width];
    widthOld = [widthOld, met.width];
  figure; hold on; plot(tempWF);scatter([maxInd], [tempWF(maxInd)]); scatter([minInd], [tempWF(minInd)]);
  title([sessionList{i} unitList{i} '  ' num2str(metSess.width)], 'Interpreter', 'none');
%      save([pd.nLynxFolderSession sessionList{i} '_' unitList{i} '_met.mat'],'metSess', 'met');
end
%% plot metrices
baselineFreq = [];
firstSpikeFreq = [];
firstSpikeLat = [];
width = [];
spikeLatChange = [];
spikeFreqChange = [];
waveforms = [];
waveformsSession = [];
maxInd = [];
widthSess = [];

for i = 1:length(sessionList)
    pd = parseSessionString_df(sessionList{i}, root, sep);
    if exist([pd.nLynxFolderSession sessionList{i} '_' unitList{i} '_met.mat'],'file')
        load([pd.nLynxFolderSession sessionList{i} '_' unitList{i} '_met.mat']);
    else
        met = getClusterMetric(sessionList{i}, unitList{i}, 0, 1);
    end
     respInds = find(met.spikeProp>=0.8);
    if length(respInds)>=1
        respLat = nanmin(met.spikeLat(respInds));
        if respLat <= 15000 && met.isiV <= 0.001 && met.distance<0.3
            baselineFreq = [baselineFreq met.baseline];
%             firstSpikeFreq = [firstSpikeFreq mean(met.spikeNum(:,1))];
            firstSpikeFreq = [firstSpikeFreq 1000*mean(met.spikeNum(:,1))/20];
            firstSpikeLat = [firstSpikeLat met.spikeLat(1)/1000];
            width = [width met.width];
            widthSess = [widthSess metSess.width];
            spikeLatChange = [spikeLatChange (met.spikeLat(10)-met.spikeLat(1))/1000];
            spikeFreqChange = [spikeFreqChange 1000*(mean(met.spikeNum(:,10))-mean(met.spikeNum(:,1)))/met.pulseWidth];
%             spikeFreqChange = [spikeFreqChange mean(met.spikeNum(:,10))-mean(met.spikeNum(:,1))];
            [peak,peakCh] = max(max(met.optoWaveform));
            [~,peakSamp] = max(met.optoWaveform(:,peakCh));
            tempWaveform = [zeros(1, 16-peakSamp), met.optoWaveform(max(1,peakSamp-15):min(peakSamp+25,size(met.optoWaveform,1)),peakCh)', zeros(1, peakSamp+25 - size(met.optoWaveform,1))];
            waveforms = [waveforms; tempWaveform/peak];
            [peak,~] = max(max(metSess.waveform));
            [~, minCh] = min(max(metSess.waveform));
            [~,peakSamp] = max(metSess.waveform(:,peakCh));
            tempWF = metSess.waveform(:,peakCh);
            tempWF = tempWF(peakSamp-50:peakSamp+50);
            maxInd = [maxInd peakSamp];
            waveformsSession = [waveformsSession; tempWF'/peak];
            if length(maxInd)==58
                fprintf([num2str(i), '\n'])
            end
        end
    end 
    
end
%%
% pca of spikeWaveform
[coeff,score,latent, ~, explained, mu] = pca(zscore(waveforms, 0,2));

[coeffSess,scoreSess,latentSess, ~, explainedSess, muSess] = pca(zscore(waveformsSession, 0,2));
%% plot all waveforms
[~, ind] = sort(scoreSess(:,1));
figure;
for i = 1:length(ind)
subplot(12, 13, i); hold on;
plot(waveformsSession(ind(i),5:100), 'LineWidth', 1.5);
title(num2str(score(ind(i),1))); set(gca, 'XColor', 'none', 'YColor', 'none')
end
%%
figure;
matrix = [baselineFreq' firstSpikeFreq' firstSpikeLat' spikeLatChange' spikeFreqChange' width' widthSess', score(:,1) score(:,2), scoreSess(:,1), scoreSess(:,2)];
names = {'baselineFreq', 'firstSpikeFreq', 'firstSpikeLat', 'spikeLatChange', 'spikeFreqChange', 'width', 'widthSess', 'PC1', 'PC2', 'PC1Sess', 'PC2Sess'};
scatterAll(matrix, names, 5,'c');
%% compare two waveforms
% pca of spikeWaveform

%%  lick distribution
animalName = 'ZS062';
col = 'good';
sampFile = [animalName col '_', model];
path = [root animalName sep animalName 'sorted' sep 'stan' sep 'bernoulli' sep model sep col sep];
load([path sampFile '.mat'], 'dayList');
for i = 1:length(dayList)
    lickCatAnalysis(dayList{i}, model, col);
end
%% spike raster & GLM

for i = 1:length(sessionList)
    pd = parseSessionString_df(sessionList{i}, root, sep);
    if exist([pd.nLynxFolderSession sessionList{i} '_' unitList{i} '_met.mat'],'file')
        load([pd.nLynxFolderSession sessionList{i} '_' unitList{i} '_met.mat']);
    else
        met = getClusterMetric(sessionList{i}, unitList{i}, 0, 1);
    end
    respInds = find(met.spikeProp>=0.8);
    if ~isempty(respInds)
        respLat = nanmin(met.spikeLat(respInds));
        if respLat <= 15000 && met.isiV <= 0.001
            [animalName, date] = strtok(sessionList{i}, 'd'); 
            animalName = animalName(2:end);
            savePath = [root animalName sep animalName 'sorted' sep 'rasters' sep 'FA' sep];
            savePathGLM = [root animalName sep animalName 'sorted' sep 'glm' sep];
            if ~exist(savePath, 'dir')
                mkdir(savePath)
            end
            if ~exist(savePathGLM, 'dir')
                mkdir(savePathGLM)
            end         
            os = behAnalysisNoPlot_opMD(sessionList{i});
            sessionData = os.behSessionData;
            lickInds  = [];
            for k = 1:length(os.behSessionData)
                if ~isnan(sessionData(k).rewardL)
                    lickInds = [lickInds k];
                elseif ~isnan(sessionData(k).rewardR)
                    lickInds = [lickInds k];
                elseif ismember(k, find(os.CSminus>0))
                    templick = [[sessionData(k).licksL], [sessionData(k).licksR]];
                    if min(templick)<sessionData(k).CSon + 1800 %1800 is resp window
                        lickInds = [lickInds k];
                    end
                end
            end
            if length(setdiff(lickInds,find(os.CSplus>0)))>=5 && length(setdiff(find(os.CSplus>0), lickInds))>=5           
                spikeRasters_dF_cueAlign(sessionList{i},'cellName',unitList{i});
                saveFigurePDF(gcf, [savePath sep sessionList{i} '_' unitList{i} 'cueAlign' '.pdf']);
            end
%             spikeGLM_dF_allRegressors(sessionList{i}, 'good','cellName', unitList{i},'regressors','1+pe+ outcome + rightSide*outcome');
%             saveFigurePDF(gcf, [savePathGLM sep sessionList{i} '_' unitList{i} 'GLMChoiceAlign_pePe' '.pdf']);
        end
    end  
end 
%% append pdfs
savePath = ['F:\allUnits\pupilCorr\'];
allFiles = dir(savePath);
allFiles = {allFiles([allFiles.bytes]>0).name}';
allFiles = strcat(savePath, allFiles);
append_pdfs([savePath 'combinePupilCorr.pdf'],allFiles{:});
%% pupil align
errors = ones(size(sessionList));
for i = 1:32
    error = timeAlign(sessionList{i},1);
    if ~isnan(error)
        errors(i) = error;
    end       
end
%% simultaneuos units
sessionListNoRepeat = unique(sessionList);
for i = 1:length(sessionListNoRepeat)
    sessionUnitsPE(sessionListNoRepeat{i}, 'good');
end
%% get clean list
sessionListGood = {};
unitListGood = {};

for i = 1:length(sessionList)
    pd = parseSessionString_df(sessionList{i}, root, sep);
    if exist([pd.nLynxFolderSession sessionList{i} '_' unitList{i} '_met.mat'],'file')
        load([pd.nLynxFolderSession sessionList{i} '_' unitList{i} '_met.mat']);
    else
        met = getClusterMetric(sessionList{i}, unitList{i}, 0, 1);
    end
    respInds = find(met.spikeProp>=0.8);
    if ~isempty(respInds)
        respLat = nanmin(met.spikeLat(respInds));
        if respLat <= 15000 && met.isiV <= 0.001 && met.distance <= 0.3 
            sessionListGood = [sessionListGood, sessionList{i}];
            unitListGood = [unitListGood, unitList{i}];
        end
    end  
end 
%%
savePath = 'F:\allUnits\pupilCorr\';
for i = 101:length(sessionListGood)
    pupillCorr(i) = unitCorrPupil(sessionListGood{i}, unitListGood{i}, 'lag', 3000, 'binSize', 200, 'binSizePost', 200, 'tf', 0.2, 'binSizePre', 200);
    saveFigurePDF(gcf, [savePath sessionListGood{i} unitListGood{i} '_pupilCorr.pdf']);
end
%%

