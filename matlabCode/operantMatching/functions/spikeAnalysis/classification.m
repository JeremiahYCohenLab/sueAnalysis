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
        if respLat <= 15000 && met.isiV <= 0.001 && met.distance<0.3
            %% opto ID file
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
            %% session file
            baselineFreq = [baselineFreq ];
            
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
savePath = ['F:\allUnits\spikeRasters\missVSHit\'];
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
%% simultaneuos units
sessionListNoRepeat = unique(sessionList);
for i = 1:length(sessionListNoRepeat)
    sessionUnitsPE(sessionListNoRepeat{i}, 'good');
end
%%

