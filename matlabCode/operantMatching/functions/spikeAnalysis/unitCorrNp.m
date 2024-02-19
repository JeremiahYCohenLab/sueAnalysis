function s = unitCorrNp(session, unit1, unit2, varargin)
%task and model parameters
p = inputParser;
% default parameters if none given
p.addParameter('plotFlag', 1);
p.addParameter('maxTrial', 1000);
p.addParameter('binSize', 200)% in ms  
p.addParameter('binSizePost', 50)% in ms
% p.addParameter('stepSize', 2)  
p.addParameter('tb', 4)% in s
p.addParameter('tf', 5)% in s
p.addParameter('saveFigFlag', 0);
p.parse(varargin{:});

binSize = p.Results.binSize;
stepSize = binSize;
binSizePost = p.Results.binSizePost;
stepSizePost = binSizePost;
% basic info
[root, sep] = currComputer();
[animalName, date] = strtok(session, 'd'); 
animalName = animalName(2:end); 
date = date(1:9);
sessionFolder = ['m' animalName date];  
% paths
if isstrprop(session(end), 'alpha')
    neuralynxDataPath = [root animalName sep sessionFolder sep 'sorted' sep 'session ' session(end) sep session '_sessionData_nL.mat'];
    savepath = [root animalName sep session(1:end-1) sep  'figures' sep 'session ' sessionName(end) sep];
    unitPath = [root animalName sep session(1:end-1) sep 'neuralynx' sep 'session' sep];
else
    neuralynxDataPath = [root animalName sep sessionFolder sep 'sorted' sep 'session' sep session '_sessionData_nL.mat'];
    savepath = [root animalName sep session sep  'figures' sep];   
    unitPath = [root animalName sep session sep 'neuralynx' sep 'session' sep];
end
% load behavior & neuron
if exist(neuralynxDataPath,'file')
    load(neuralynxDataPath)
else
    sessionData = generateSessionData_nL_operantMatching(session);
end
% time window
time = -1000*p.Results.tb:1000*p.Results.tf;
% save path
if isempty(dir(savepath))
    mkdir(savepath)
end
%% find identified units, create spike and lick cell

spikeFields = fields(sessionData);
xlFile = [animalName '.xlsx'];
[nums, unitsInfo,~] = xlsread([root xlFile], 'neurons');
sessionList = unitsInfo(2:end,1); 
unitList = unitsInfo(2:end,2);
optoUnitList = unitsInfo(2:end,8);
subFolders = unitsInfo(2:end,9);

drift = unitsInfo(2:end,10);
drift = contains(drift, 'drift')|contains(drift, 'Drift');
quality = nums(:,1);
maxLratio = 0.06;
sessionList = sessionList(quality<=maxLratio & ~drift);
unitList = unitList(quality<=maxLratio & ~drift);
optoUnitList = optoUnitList(quality<=maxLratio & ~drift);
subFolders = subFolders(quality<=maxLratio & ~drift);

sessionID = contains(sessionList,session);
sessionUnits = unitList(sessionID);
sessionSubFolders = subFolders(sessionID);
unitsQual = zeros(1,length(sessionUnits));
for a = 1:length(sessionUnits)
    if exist([unitPath session '_' sessionUnits{a} '_met.mat'],'file')
        load([unitPath session '_' sessionUnits{a} '_met.mat']);
    else
        met = getClusterMetric(session, sessionUnits{a}, 0, 1);
    end
     respInds = find(met.spikeProp>=0.8);
    if ~isempty(respInds)
        respLat = nanmin(met.spikeLat(respInds));
        if respLat <= 15000 && met.isiV <= 0.001 && met.distance<0.3
            unitsQual(a) = 1;
        end
    end 
end

sessionUnits = sessionUnits(unitsQual>0);
clust = zeros(size(sessionUnits));
for i = 1:length(clust)
    clust(i) = find(contains(spikeFields,sessionUnits{i}));
end
allClustsInds = find(contains(spikeFields, 'TT'));
allSpikesInds = find(contains(spikeFields, 'allSpikes'));

% focus unit1
focusInd1 = clust(contains(sessionUnits, unit1));
if isempty(focusInd1)
    fprintf([session unit1 ' not good' '\n'])
    return
end

% focus unit2
focusInd2 = clust(contains(sessionUnits, unit2));
if isempty(focusInd2)
    fprintf([session unit2 ' not good' '\n'])
    return
end
clust = [focusInd1; focusInd2];

%% put spikes by trials into cells, aligned by cue
allTrial_spike_cue = {};
for k = 1:length(sessionData)
    for i = 1:length(clust)
        if k == 1
            prevTrial_spike = [];
        else
            prevTrial_spikeInd = [sessionData(k-1).(spikeFields{clust(i)})] > (sessionData(k).CSon-p.Results.tb*1000);
            prevTrial_spike = sessionData(k-1).(spikeFields{clust(i)})(prevTrial_spikeInd) - sessionData(k).CSon;
        end
        
        currTrial_spikeInd = sessionData(k).(spikeFields{clust(i)}) < sessionData(k).CSon+p.Results.tf*1000 ... 
            & sessionData(k).(spikeFields{clust(i)}) > sessionData(k).CSon-p.Results.tb*1000;
        currTrial_spike = sessionData(k).(spikeFields{clust(i)})(currTrial_spikeInd) - sessionData(k).CSon;
        
        allTrial_spike_cue{i,k} = [prevTrial_spike currTrial_spike];

    end
end

% sometimes no licks/spikes are considered 1x0 and sometimes they are []
% plotSpikeRaster does not place nicely with [] so this converts all empty indices to 1x0
allTrial_spike_cue(cellfun(@isempty,allTrial_spike_cue)) = {zeros(1,0)}; 

%% loop through neurons
for j = 1:length(sessionData)
    trialDurDiff(j) = (sessionData(j).trialEnd - sessionData(j).CSon)- p.Results.tf*1000;
end
trialDurDiff(end) = 0; 

allSpikes = [];
allSpikesPre = [];
allSpikesAllTime = [];
allSpikesPreSession = [];
allSpikesPreCue = [];

for i = 1:length(clust)
    % Initialize matrices for SDF
    allTrial_spikeMatx_cue = zeros(length(sessionData),length(time));
   for j = 1:length(allTrial_spike_cue)
        tempSpike = allTrial_spike_cue{i,j};
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
    % slide window for pre
    midPoints = (0.5*binSize + 1):stepSize:(length(time(time<0))-0.5*binSize);
    allTrial_spikeMatx_slide = zeros(length(sessionData), length(midPoints));
    for w = 1:length(midPoints)
        allTrial_spikeMatx_slide(:,w) = ...
            nansum(allTrial_spikeMatx_cue(:,midPoints(w)-0.5*binSize:midPoints(w)+0.5*binSize-1),2)*1000/binSize;
    end
    % slide window for in trial
    midPointsPost = 1000*p.Results.tb+(0.5*binSizePost + 1):stepSizePost:(length(time)-0.5*binSizePost);
    allTrial_spikeMatx_slidePost = zeros(length(sessionData), length(midPointsPost));
    for w = 1:length(midPointsPost)
        allTrial_spikeMatx_slidePost(:,w) = ...
            nansum(allTrial_spikeMatx_cue(:,midPointsPost(w)-0.5*binSizePost:midPointsPost(w)+ 0.5*binSizePost-1),2)*1000/binSizePost;
    end    
    % currAllSpikes into bins
    currAllSpikes = sessionData(clust(i)- min(allClustsInds)+1).allSpikes;
    currAllSpikes = currAllSpikes - sessionData(1).CSon  + 1;
    currAllSpikes = currAllSpikes(currAllSpikes>0);
    if i==1
        binStarts = 1:stepSize:(max(currAllSpikes)+stepSize);
    end
%     for a = 1:length(binStarts)
%         currAllSpikesBin(a) = sum(currAllSpikes>=binStarts(a) & currAllSpikes<binStarts(a)+p.Results.binSize);
%     end
    currAllSpikesBin = histcounts(currAllSpikes,binStarts)*1000/binSize;
    % current pre session spikes
    % find recording start time: 
    pd = parseSessionString_df(session, root, sep);
    [timestamps, eventID, TTL, Evstring] = Nlx2MatEV([pd.nLynxFolderSession 'Events.nev'],[1 1 1 0 1], 0, 1);
    timestamps = round(timestamps/1000);
    startTime = min(timestamps);
    currAllSpikesPreSess = sessionData(clust(i)-min(allClustsInds)+1).allSpikes;
    currAllSpikesPreSess = currAllSpikesPreSess(currAllSpikesPreSess<sessionData(1).CSon & currAllSpikesPreSess>sessionData(1).CSon-20000);
    if i==1
        binStartsPre = max(sessionData(1).CSon-20000, startTime):stepSize:sessionData(1).CSon;
    end
    currAllSpikesPreSessBin = histcounts(currAllSpikesPreSess,binStartsPre)*1000/binSize;
    % zscore spike matrix
    allSpikes = [allSpikes, reshape((allTrial_spikeMatx_slidePost)', [], 1)];
    allSpikesPre = [allSpikesPre, reshape((allTrial_spikeMatx_slide)', [], 1)];
    allSpikesAllTime = [allSpikesAllTime; currAllSpikesBin];
    allSpikesPreSession = [allSpikesPreSession; currAllSpikesPreSessBin];
    allSpikesPreCue = [allSpikesPreCue; mean(allTrial_spikeMatx_slide, 2, 'omitnan')'];
end
%% calculate CorrCoeff
[hIn, pIn] = corrcoef(allSpikes);
[hPre, pPre] = corrcoef(allSpikesPre);
[hAll, pAll] = corrcoef(allSpikesAllTime');
[hPreSess, pPreSess] = corrcoef(allSpikesPreSession');
[hPreCue, pPreCue] = corrcoef(allSpikesPreCue');
s.hIn = hIn(1,2);
s.pIn = pIn(1,2);
s.hPre = hPre(1,2);
s.pPre = pPre(1,2);
s.hAll = hAll(1,2);
s.pAll = pAll(1,2);
s.hPreSess = hPreSess(1,2);
s.pPreSess = pPreSess(1,2);
s.pPreCue = pPreCue(1,2);
s.hPreCue = hPreCue(1,2);


   
