function s = unitCorrPupil(session, unit, varargin)
%task and model parameters
p = inputParser;
% default parameters if none given
p.addParameter('plotFlag', 0);
p.addParameter('saveFlag', 0);
p.addParameter('maxTrial', 1000);
p.addParameter('binSize', 200)% in ms  
p.addParameter('binSizePost', 500)
p.addParameter('binSizePre', 300)% in ms  
p.addParameter('lag', 2000)% in ms
p.addParameter('lagPre', 2000)% in ms
p.addParameter('tb', 2); % in s, cannot be longer than 2s
p.addParameter('tf', 0.5);
p.parse(varargin{:});

binSize = p.Results.binSize;
stepSize = binSize;
binSizePost = p.Results.binSizePost;
stepSizePost = binSizePost;
binSizePre = p.Results.binSizePre;
stepSizePre = binSizePre;
% basic info
[root, sep] = currComputer();
[animalName, date] = strtok(session, 'd'); 
animalName = animalName(2:end); 
date = date(1:9);
sessionFolder = ['m' animalName date];  
% paths
pd = parseSessionString_df(session, root, sep);
if isstrprop(session(end), 'alpha')
    neuralynxDataPath = [root animalName sep sessionFolder sep 'sorted' sep 'session ' session(end) sep session '_sessionData_nL.mat'];
    savepath = [root animalName sep session(1:end-1) sep  'figures' sep 'session ' sessionName(end) sep];
    unitPath = [root animalName sep session(1:end-1) sep 'neuralynx' sep 'session' sep];
    pupilAlignPath = [pd.sortedFolder sep 'session' sep session '_pupil.mat'];
else
    neuralynxDataPath = [root animalName sep sessionFolder sep 'sorted' sep 'session' sep session '_sessionData_nL.mat'];
    savepath = [root animalName sep session sep  'figures' sep];   
    unitPath = [root animalName sep session sep 'neuralynx' sep 'session' sep];
    pupilAlignPath = [pd.sortedFolder sep 'session' sep session '_pupil.mat'];
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

% focus unit
focusInd = clust(contains(sessionUnits, unit));
if isempty(focusInd)
    fprintf([session unit ' not good' '\n'])
    return
end

%% put spikes by trials into cells, aligned by cue
allTrial_spike_cue = {};
for k = 1:length(sessionData)
    if k == 1
        prevTrial_spike = [];
    else
        prevTrial_spikeInd = [sessionData(k-1).(spikeFields{focusInd})] > (sessionData(k).CSon-p.Results.tb*1000);
        prevTrial_spike = sessionData(k-1).(spikeFields{focusInd})(prevTrial_spikeInd) - sessionData(k).CSon;
    end

    currTrial_spikeInd = sessionData(k).(spikeFields{focusInd}) < sessionData(k).CSon+p.Results.tf*1000 ... 
        & sessionData(k).(spikeFields{focusInd}) > sessionData(k).CSon-p.Results.tb*1000;
    currTrial_spike = sessionData(k).(spikeFields{focusInd})(currTrial_spikeInd) - sessionData(k).CSon;

    allTrial_spike_cue{1,k} = [prevTrial_spike currTrial_spike];
end

% sometimes no licks/spikes are considered 1x0 and sometimes they are []
% plotSpikeRaster does not place nicely with [] so this converts all empty indices to 1x0
allTrial_spike_cue(cellfun(@isempty,allTrial_spike_cue)) = {zeros(1,0)}; 

for j = 1:length(sessionData)
    trialDurDiff(j) = (sessionData(j).trialEnd - sessionData(j).CSon)- p.Results.tf*1000;
end
trialDurDiff(end) = 0; 

allSpikes = [];
allSpikesPre = [];
allSpikesAllTime = [];
allSpikesPreSession = [];


% Initialize matrices for SDF
allTrial_spikeMatx_cue = zeros(length(sessionData),length(time));
for j = 1:length(allTrial_spike_cue)
    tempSpike = allTrial_spike_cue{1,j};
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
midPoints = (0.5*binSizePre + 1):stepSizePre:(length(time(time<0))-0.5*binSizePre);
allTrial_spikeMatx_slide = zeros(length(sessionData), length(midPoints));
for w = 1:length(midPoints)
    allTrial_spikeMatx_slide(:,w) = ...
        nansum(allTrial_spikeMatx_cue(:,midPoints(w)-0.5*binSizePre:midPoints(w)+0.5*binSizePre-1),2)*1000/binSizePre;
end
% slide window for in trial
midPointsPost = 1000*p.Results.tb+(0.5*binSizePost + 1):stepSizePost:(length(time)-0.5*binSizePost);
allTrial_spikeMatx_slidePost = zeros(length(sessionData), length(midPointsPost));
for w = 1:length(midPointsPost)
    allTrial_spikeMatx_slidePost(:,w) = ...
        nansum(allTrial_spikeMatx_cue(:,midPointsPost(w)-0.5*binSizePost:midPointsPost(w)+ 0.5*binSizePost-1),2)*1000/binSizePost;
end    
% currAllSpikes into bins
currAllSpikes = sessionData(focusInd- min(allClustsInds)+1).allSpikes;
currAllSpikes = currAllSpikes - sessionData(1).CSon  + 1;
currAllSpikes = currAllSpikes(currAllSpikes>0);

binStarts = 1:stepSize:(max(currAllSpikes)+stepSize);

%     for a = 1:length(binStarts)
%         currAllSpikesBin(a) = sum(currAllSpikes>=binStarts(a) & currAllSpikes<binStarts(a)+p.Results.binSize);
%     end
currAllSpikesBin = histcounts(currAllSpikes,binStarts)*1000/binSize;
% current pre session spikes
% find recording start time: 
[timestamps, eventID, TTL, Evstring] = Nlx2MatEV([pd.nLynxFolderSession 'Events.nev'],[1 1 1 0 1], 0, 1);
timestamps = round(timestamps/1000);
startTime = min(timestamps);
currAllSpikesPreSess = sessionData(focusInd-min(allClustsInds)+1).allSpikes;
currAllSpikesPreSess = currAllSpikesPreSess(currAllSpikesPreSess<sessionData(1).CSon & currAllSpikesPreSess>sessionData(1).CSon-20000);

binStartsPre = max(sessionData(1).CSon-20000, startTime):stepSize:sessionData(1).CSon;

currAllSpikesPreSessBin = histcounts(currAllSpikesPreSess,binStartsPre)*1000/binSize;
% spikes time series
allSpikes = reshape((allTrial_spikeMatx_slidePost)', [], 1);
allSpikesPre = reshape((allTrial_spikeMatx_slide)', [], 1);
allSpikesAllTime = currAllSpikesBin;
allSpikesPreSession = currAllSpikesPreSessBin;

%% calculate all binned mean pupil diameter
% load pupil 
if exist(pupilAlignPath, 'file')
    load(pupilAlignPath);
    if errorProp >= 0.5
        fprintf([session ' badly aligned pupil \n'])
        s.allCoeffs = [];
        s.allPs = [];
        s.evokeCoeffs = [];
        s.evokePs = [];
        s.baselineCoeffs = [];
        s.baselinePs = [];        
    end
else
    fprintf([session ' no aligned pupil \n'])
    s.allCoeffs = [];
    s.allPs = [];
    s.evokeCoeffs = [];
    s.evokePs = [];
    s.baselineCoeffs = [];
    s.baselinePs = [];
    return
end
numLag = floor(p.Results.lag/p.Results.binSize); 
pupilLaggedAll = NaN(numLag, length(currAllSpikesBin));
cueFT(cueFT<=0) = 1;
for i = 1:numLag
    frameStartsTemp = cueFT(1)  + round(FR*p.Results.binSize/1000 * (i - 1))+ round(FR/1000*(binStarts-binStarts(1)));
    frameStartsTemp = frameStartsTemp(frameStartsTemp<=length(diaRealign));
    diaTemp = diaRealign(frameStartsTemp(1):frameStartsTemp(end));
    cellSizes = diff(frameStartsTemp);
    cellSizes(1) = cellSizes(1) + 1; % include the first bin into first package
    diaLag = mat2cell(diaTemp,1,cellSizes);
    diaLag = cellfun(@(x) mean(x, 'omitnan'), diaLag);
    pupilLaggedAll(i,1:length(diaLag)) = diaLag;
end
% calculate CorrCoeff
allCoeffs = NaN(1,numLag);
allPs = NaN(1,numLag);
for i = 1:numLag
    [hTmp, pTmp, RLTmp, RUTmp] = corrcoef(allSpikesAllTime, pupilLaggedAll(i,:), 'Rows', 'complete');
    allCoeffs(i) = hTmp(1,2);
    allPs(i) = pTmp(1,2);
end

s.allCoeffs = allCoeffs;
s.allPs = allPs;
%% calculate coeffs of evoked pupil
% find pupil after cue
tbPupil = round(FR*2)+1;
numLagPost = floor(p.Results.lag/binSizePost); 
pupilLaggedAfterCue = NaN(numLagPost, length(allSpikes));
baselineDia = mean(sessionPupilCue(:,1:(tbPupil-1)), 2, 'omitnan');
for i = 1:numLagPost
    pupilLaggedAfterCue(i,:) = mean(sessionPupilCue(:,(tbPupil:(tbPupil+round(FR*binSizePost/1000)))+(i-1)*round(binSizePost/1000*FR)), 2, 'omitnan')'./baselineDia';
%     pupilLaggedAfterCue(i,:) = mean(sessionPupilCue(:,(tbPupil:(tbPupil+round(FR*binSizePost/1000)))+(i-1)*round(binSizePost/1000*FR)), 2, 'omitnan')';

    pupilLaggedAfterCue(i,~qualInd) = NaN;
end

% calculate CorrCoeff
evokeCoeffs = NaN(1,numLagPost);
evokePs = NaN(1,numLagPost);
for i = 1:numLagPost
    [hTmp, pTmp] = corrcoef(allSpikes, pupilLaggedAfterCue(i,:)', 'Rows', 'complete');
    evokeCoeffs(i) = hTmp(1,2);
    evokePs(i) = pTmp(1,2);
end

s.evokeCoeffs = evokeCoeffs;
s.evokePs = evokePs;
%% calculate coeff of pre cue pupil
allSpikesPre = (allTrial_spikeMatx_slide(:,1)); % (only use the fist one)
numLagPre = floor(p.Results.lagPre/binSizePre); % counts from -2s
pupilLaggedPreCue = NaN(numLagPre, length(allSpikesPre));
for i = 1:numLagPre
    pupilLaggedPreCue(i,:) = mean(sessionPupilCue(:,(1:(1+round(FR*binSizePre/1000)))+(i-1)*round(binSizePre/1000*FR)), 2, 'omitnan')';
    pupilLaggedPreCue(i,~qualInd) = NaN;
end

% calculate CorrCoeff
baselineCoeffs = NaN(1,numLagPre);
baselinePs = NaN(1,numLagPre);
for i = 1:numLagPre
    [hTmp, pTmp] = corrcoef(allSpikesPre, pupilLaggedPreCue(i,:)', 'Rows', 'complete');
    baselineCoeffs(i) = hTmp(1,2);
    baselinePs(i) = pTmp(1,2);
end

s.baselineCoeffs = baselineCoeffs;
s.baselinePs = baselinePs;
currAllSpikes = sessionData(focusInd- min(allClustsInds)+1).allSpikes;
sessionTime = 1:(currAllSpikes(end)-currAllSpikes(1) + 1);
session_spikeArray = zeros(1,length(sessionTime));      
session_spikeArray(currAllSpikes - currAllSpikes(1) + 1) = 1000;

boxKern = ones(1,5000) / 5000;                                      smoothKern = ones(1,5000) / 5000;
sessionSpikeSDF = conv(session_spikeArray, boxKern);                sessionSpikeSDFsmooth = conv(session_spikeArray, smoothKern);
sessionSpikeSDF = sessionSpikeSDF(1:(end-length(boxKern)-1));       sessionSpikeSDFsmooth = sessionSpikeSDFsmooth(1:(end-length(smoothKern)+1));

realTime = 1/1000*(sessionTime - (sessionData(1).CSon - currAllSpikes(1)));

time = 1/FR*((1:length(diaRealign)) - cueFT(1));
smoothedDia = smoothdata(diaRealign, 10, 'omitnan');
newFrame = cueFT(1) + round(FR*([sessionData.CSon]-[sessionData(1).CSon])/1000);
%% 
if p.Results.plotFlag
    figure2;
    subplot(1,3,1); hold on
    timeLag = stepSize*(0:numLag-1);
    plot(timeLag, allCoeffs, 'LineWidth', 2);
    scatter(timeLag(allPs<0.05), allCoeffs(allPs<0.05), 12, 'r', 'filled');
    plot(timeLag, zeros(size(allCoeffs)), 'LineWidth', 2, 'LineStyle', '--', 'Color', [0.7 0.7 0.7]);
    title('across all time');
 
    subplot(2,3,2); hold on
    timeLag = stepSizePost*(0:numLagPost-1);
    plot(timeLag, evokeCoeffs, 'LineWidth', 2);
    scatter(timeLag(evokePs<0.05), evokeCoeffs(evokePs<0.05), 12, 'r', 'filled');
    plot(timeLag, zeros(size(evokeCoeffs)), 'LineWidth', 2, 'LineStyle', '--', 'Color', [0.7 0.7 0.7]);
    title('evoked');
    
    subplot(2,3,5); hold on
    timeP = 1000/FR * (1:size(sessionPupilCue,2)) - 2000;
    plotFilled(timeP,sessionPupilCue, 'm');
    title('pupil dynamics')
    
    subplot(1,3,3); hold on
    timeLag = stepSizePre*(0:numLagPre-1);
    plot(timeLag, baselineCoeffs, 'LineWidth', 2);
    scatter(timeLag(baselinePs<0.05), baselineCoeffs(baselinePs<0.05), 12, 'r', 'filled');
    plot(timeLag, zeros(size(baselineCoeffs)), 'LineWidth', 2, 'LineStyle', '--', 'Color', [0.7 0.7 0.7]);
    title('baseline');
    
    figure2;
    subplot(2,1,1); hold on;
    plot(realTime, sessionSpikeSDFsmooth, 'k', 'LineWidth',2);
    scatter([[sessionData.CSon]-sessionData(1).CSon]/1000, 3.2*ones(size([sessionData.CSon])), 18, 'r', 'filled');
    xlim([100 200])
    subplot(2,1,2); hold on;
    plot(time, smoothedDia, 'b', 'LineWidth',2);
    scatter([time(newFrame)], 120*ones(size(newFrame)), 18, 'r', 'Filled');
    xlim([100 200])
    
    sgtitle([session unit '_pupilCorrelation'],'Interpreter', 'none')
    screen = get(0,'Screensize');
    screen(4) = screen(4) - 100;
    set(gcf, 'Position', screen)
end



if p.Results.saveFlag
    if ~exist([savepath sep 'pupil' sep], 'dir')
        mkdir([savepath sep 'pupil' sep]);
    end
    saveFigurePDF(gcf, [savepath sep 'pupil' sep session unit '_pupilCorrelation.pdf']);

end

% edges = binEqualSize(allSpikes, 4);
% indSort1 = find(allSpikes==uniq(1));
% indSort2 = find(allSpikes==uniq(2));
% indSort3 = find(allSpikes==uniq(3));
% indSort4 = find(allSpikes==uniq(4)|allSpikes==uniq(5));
% 
% colors = zeros(4,3);
% for i = 1:4
%     colors(i,:) = [1, 1-(1/4)*i, 1-(1/4)*i];
% end
% figure2; hold on;
% time = linspace(-2,10,size(sessionPupilCue,2));
% currMat = sessionPupilCue(intersect(indSort1, find(qualInd>0)),:);
% plotFilled(time, currMat, colors(1,:));
% currMat = sessionPupilCue(intersect(indSort2, find(qualInd>0)),:);
% plotFilled(time, currMat, colors(2,:));
% currMat = sessionPupilCue(intersect(indSort3, find(qualInd>0)),:);
% plotFilled(time, currMat, colors(3,:));
% currMat = sessionPupilCue(intersect(indSort4, find(qualInd>0)),:);
% plotFilled(time, currMat, colors(4,:));
% 
% figure2; hold on;
% time = linspace(-2,10,size(sessionPupilCue,2));
% currMat = sessionPupilCue(intersect(indSort1, find(qualInd>0)),:);
% plot(time, mean(currMat, 1, 'omitnan'), 'Color', colors(1,:), 'LineWidth', 2);
% currMat = sessionPupilCue(intersect(indSort2, find(qualInd>0)),:);
% plot(time, mean(currMat, 1, 'omitnan'), 'Color', colors(2,:), 'LineWidth', 2);
% currMat = sessionPupilCue(intersect(indSort3, find(qualInd>0)),:);
% plot(time, mean(currMat, 1, 'omitnan'), 'Color', colors(3,:), 'LineWidth', 2);
% currMat = sessionPupilCue(intersect(indSort4, find(qualInd>0)),:);
% plot(time, mean(currMat, 1, 'omitnan'), 'Color', colors(4,:), 'LineWidth', 2);





   
