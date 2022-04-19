function s = unitCorrBaseline(session, unit1, unit2, varargin)
%task and model parameters
p = inputParser;
% default parameters if none given
p.addParameter('plotFlag', 1);
p.addParameter('maxTrial', 1000);
p.addParameter('binSize', 0.2)% in ms  
p.addParameter('lag', 5)% in ms
p.addParameter('saveFigFlag', 0);
p.addParameter('numShuffle', 1000);
p.addParameter('jitter', 5); % in ms
p.parse(varargin{:});

binSize = p.Results.binSize;
stepSize = binSize;
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
    sortedDir = [root animalName sep session(1:end-1) sep 'neuralynx' sep 'session ' sessionName(end) sep];
else
    neuralynxDataPath = [root animalName sep sessionFolder sep 'sorted' sep 'session' sep session '_sessionData_nL.mat'];
    savepath = [root animalName sep session sep  'figures' sep];   
    unitPath = [root animalName sep session sep 'neuralynx' sep 'session' sep];
    sortedDir = [root animalName sep session sep 'neuralynx' sep 'session' sep];
end
% load behavior & neuron
if exist(neuralynxDataPath,'file')
    load(neuralynxDataPath)
else
    sessionData = generateSessionData_nL_operantMatching(session);
end
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


%% loop through neurons

allSpikesPreSession = [];
% Clustered data
sortedFiles = dir(sortedDir);
allUnits = [unit1, unit2]; 
for i = 1:length(clust)
    % currAllSpikes into bins
    % find recording start time: 
    pd = parseSessionString_df(session, root, sep);
    [timestamps, eventID, TTL, Evstring] = Nlx2MatEV([pd.nLynxFolderSession 'Events.nev'],[1 1 1 0 1], 0, 1);
    startTime = min(timestamps);
    spikeInds = contains({sortedFiles.name},allUnits{i}) & contains({sortedFiles.name},'txt');
    currAllSpikesPreSess = load([sortedDir sortedFiles(spikeInds).name]);

    currAllSpikesPreSess = currAllSpikesPreSess(currAllSpikesPreSess<sessionData(1).CSon*1000 & currAllSpikesPreSess>sessionData(1).CSon*1000-20000000); % 20s before first trial
    if i==1
        binStartsPre = max(sessionData(1).CSon*1000-20000000, startTime):stepSize*1000:sessionData(1).CSon*1000; % 20s before first trial
    end
    currAllSpikesPreSessBin = histcounts(currAllSpikesPreSess,binStartsPre)*1000/(binSize*1000);
    % zscore spike matrix
    allSpikesPreSession = [allSpikesPreSession; currAllSpikesPreSessBin];
end
% calculate crossCorr
lag = round(p.Results.lag/p.Results.binSize);
xPreSess = allSpikesPreSession{1};
yPreSess = allSpikesPreSession{2};
corrCoPreSess = xcorr(xPreSess,yPreSess,lag,'normalized');
corrCoPreSess(0.5*(length(corrCoPreSess)+1)) = NaN;
s.corrCoPreSess = corrCoPreSess;
%% generating surrogate
% find recording start time: 
allPreSession = NaN(p.Results.numShuffle, length(s.corrCoPreSess));
names = {sortedFiles.name};
parfor a = 1:p.Results.numShuffle
    allSpikesPre = [];
    allSpikesPreSession = [];
    for i = 1:length(clust)
    spikeInds = contains(names,allUnits{i}) & contains(names,'txt');
    currAllSpikesPreSess = load([sortedDir names{spikeInds}]);
    currAllSpikesPreSess = currAllSpikesPreSess(currAllSpikesPreSess<sessionData(1).CSon*1000 & currAllSpikesPreSess>sessionData(1).CSon*1000-20000000); % 20s before first trial
    
        % surrogate
        currAllSpikesPreSess = currAllSpikesPreSess + p.Results.jitter*(2*rand(1, length(currAllSpikesPreSess)) - 1)*1000;
        currAllSpikesPreSessBin = histcounts(currAllSpikesPreSess,binStartsPre)*1000000/binSize; %% binSize, ms
        % zscore spike matrix
        allSpikesPreSession(i,:) = currAllSpikesPreSessBin;
    end

    xPreSess = allSpikesPreSession(1,:);
    yPreSess = allSpikesPreSession(2,:);
    

    corrCoPreSess = xcorr(xPreSess,yPreSess,lag,'normalized');
    corrCoPreSess(0.5*(length(corrCoPreSess)+1)) = NaN;

    allPreSession(a,:) = corrCoPreSess;
end
% find alpha boundary
  
CIPreSess = NaN(2,length(s.corrCo));
for i = 1:length(s.corrCo)
    if i ~= 0.5*(length(s.corrCo)+1)
        tmp = sort(allPreTrial(:,i));
        CIPreTrial(1,i) = tmp(floor(0.5*p.Results.alpha*p.Results.numShuffle));
        CIPreTrial(2,i) = tmp(ceil((1-0.5*p.Results.alpha)*p.Results.numShuffle));
        tmp = sort(allPreSession(:,i));
        CIPreSess(1,i) = tmp(floor(0.5*p.Results.alpha*p.Results.numShuffle));
        CIPreSess(2,i) = tmp(ceil((1-0.5*p.Results.alpha)*p.Results.numShuffle));
    end
end

s.CIPreTrial = CIPreTrial;
s.CIPreSess = CIPreSess;
s.meanPreTrial = mean(allPreTrial, 1, 'omitnan');
s.menPreSession = mean(allPreSession, 1, 'omitnan');


   
