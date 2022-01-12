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
    drift = contains(drift, 'drift')|contains(drift, 'Drift');
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
savePath = [root 'allUnits\optoID\'];
for i = 1:length(sessionList)
    rasters = qualityCheck(sessionList{i},0,'unit',unitList{i});
    saveFigurePDF(rasters, [savePath sep sessionList{i} '_' unitList{i} '.pdf']);
end  
%% opto metrices
savePath = [root animalName sep animalName 'sorted' sep 'optoMetFiles' sep];
for i = 1:length(sessionList)
    getClusterMetric(sessionList{i}, unitList{i}, 0, 1);
%     saveFigurePDF(gcf, [savePath sep sessionList{i} '_' unitList{i} '.pdf']);
end 
%% plot metrices
baselineFreq = [];
firstSpikeFreq = [];
firstSpikeLat = [];
width = [];
spikeLatChange = [];
spikeFreqChange = [];
waveforms = [];

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
        if respLat <= 20000 && met.isiV <= 0.001 && met.distance<0.3
            baselineFreq = [baselineFreq met.baseline];
%             firstSpikeFreq = [firstSpikeFreq mean(met.spikeNum(:,1))];
            firstSpikeFreq = [firstSpikeFreq 1000*mean(met.spikeNum(:,1))/20];
            firstSpikeLat = [firstSpikeLat met.spikeLat(1)/1000];
            width = [width met.width];
            spikeLatChange = [spikeLatChange (met.spikeLat(10)-met.spikeLat(1))/1000];
            spikeFreqChange = [spikeFreqChange 1000*(mean(met.spikeNum(:,10))-mean(met.spikeNum(:,1)))/met.pulseWidth];
%             spikeFreqChange = [spikeFreqChange mean(met.spikeNum(:,10))-mean(met.spikeNum(:,1))];
            [peak,peakCh] = max(max(met.optoWaveform));
            [~,peakSamp] = max(met.optoWaveform(:,peakCh));
            tempWaveform = [zeros(1, 16-peakSamp), met.optoWaveform(max(1,peakSamp-15):min(peakSamp+25,size(met.optoWaveform,1)),peakCh)', zeros(1, peakSamp+25 - size(met.optoWaveform,1))];
            waveforms = [waveforms; tempWaveform/peak];
        end
    end 
    
end
%%
% pca of spikeWaveform
m = mean(waveforms,1);
waveforms = waveforms - mean(waveforms,1);
%%
[coeff,score,latent, ~, explained, mu] = pca(waveforms);
figure;
matrix = [baselineFreq' firstSpikeFreq' firstSpikeLat' spikeLatChange' spikeFreqChange' width' score(:,1) score(:,2)];
names = {'baselineFreq', 'firstSpikeFreq', 'firstSpikeLat', 'spikeLatChange', 'spikeFreqChange', 'width','PC1', 'PC2'};
scatterAll(matrix, names, 7,'c');

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
 savePath = ['F:\allUnits\spikeGLM\'];
for i = 1:length(sessionList)
    pd = parseSessionString_df(sessionList{i}, root, sep);
    if exist([pd.nLynxFolderSession sessionList{i} '_' unitList{i} '_met.mat'],'file')
        load([pd.nLynxFolderSession sessionList{i} '_' unitList{i} '_met.mat']);
    else
        met = getClusterMetric(sessionList{i}, unitList{i}, 0, 1);
    end
    respInds = find(met.spikeProp>=0.8);
    if length(respInds)>=2
        respLat = nanmin(met.spikeLat(respInds));
        if respLat <= 20000 && met.isiV <= 0.001 && met.distance < 0.3
            [animalName, date] = strtok(sessionList{i}, 'd'); 
            animalName = animalName(2:end);
%             savePath = [root animalName sep animalName 'sorted' sep 'rasters' sep];
%             savePathGLM = [root animalName sep animalName 'sorted' sep 'glm' sep];
            if ~exist(savePath, 'dir')
                mkdir(savePath)
            end
%             if ~exist(savePathGLM, 'dir')
%                 mkdir(savePathGLM)
%             end         
%             
%              spikeRasters_dF_cueAlign(sessionList{i},'cellName',unitList{i}, 'saveFigFlag', 0);
%              saveFigurePDF(gcf, [savePath sep sessionList{i} '_' unitList{i} 'choiceAlign' '.pdf']);
            spikeGLM_dF(sessionList{i}, 'good','cellName', unitList{i},'regressors','1+pe+ outcome + rightSide*outcome');
            saveFigurePDF(gcf, [savePath sep sessionList{i} '_' unitList{i} 'GLMChoiceAlign_pePe' '.pdf']);
        end
    end  
end 
%% append pdfs
savePath = ['F:\allUnits\spikeSimutaneous\'];
allFiles = dir(savePath);
allFiles = {allFiles([allFiles.bytes]>0).name}';
allFiles = strcat(savePath, allFiles);
append_pdfs([savePath 'combine.pdf'],allFiles{:});
%% pupil align
errors = ones(size(sessionList));
for i = 1:32
    error = timeAlign(sessionList{i},1);
    if ~isnan(error)
        errors(i) = error;
    end       
end
%% simultaneous units
savePath = ['F:\allUnits\spikeSimutaneous\'];
if ~exist(savePath, 'dir')
    mkdir(savePath)
end
sessionListNoRepeat = unique(sessionList);
for i = 1:length(sessionListNoRepeat)
    sessionUnitsPE(sessionListNoRepeat{i}, 'good');
    saveFigurePDF(gcf, [savePath sep sessionListNoRepeat{i} 'SimUnits' '.pdf']);
    close all;
end
%%

