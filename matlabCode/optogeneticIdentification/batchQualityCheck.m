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
%% plot metrices
baselineFreq = [];
firstSpikeFreq = [];
firstSpikeLat = [];
width = [];
spikeLatChange = [];
spikeFreqChange = [];

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
        if respLat <= 15000 && met.isiV <= 0.001 && met.distance<0.3
            baselineFreq = [baselineFreq met.baseline];
            firstSpikeFreq = [firstSpikeFreq 1000*mean(met.spikeNum(:,1))/met.pulseWidth];
            firstSpikeLat = [firstSpikeLat met.spikeLat(1)/1000];
            width = [width met.width];
            spikeLatChange = [spikeLatChange (met.spikeLat(10)-met.spikeLat(1))/1000];
            spikeFreqChange = [spikeFreqChange 1000*(mean(met.spikeNum(:,10))-mean(met.spikeNum(:,1)))/met.pulseWidth];
        end
    end 
    
end

figure;
matrix = [baselineFreq' firstSpikeFreq' firstSpikeLat' spikeLatChange' spikeFreqChange' width'];
names = {'baselineFreq', 'firstSpikeFreq', 'firstSpikeLat', 'spikeLatChange', 'spikeFreqChange', 'width'};
scatterAll(matrix, names, 7);

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
            savePath = [root animalName sep animalName 'sorted' sep 'rasters' sep];
            savePathGLM = [root animalName sep animalName 'sorted' sep 'glm' sep];
            if ~exist(savePath, 'dir')
                mkdir(savePath)
            end
            if ~exist(savePathGLM, 'dir')
                mkdir(savePathGLM)
            end         
            
%             spikeRasters_dF_choiceAlign(sessionList{i},'cellName',unitList{i});
%             saveFigurePDF(gcf, [savePath sep sessionList{i} '_' unitList{i} 'choiceAlign' '.pdf']);
            spikeGLM_dF_allRegressors(sessionList{i}, 'good','cellName', unitList{i},'regressors','1+pe+ outcome + rightSide*outcome');
%             saveFigurePDF(gcf, [savePathGLM sep sessionList{i} '_' unitList{i} 'GLMChoiceAlign_pePe' '.pdf']);
        end
    end  
end 
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

