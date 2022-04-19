function s = unitAutoCorr(session, unit1, unit2, varargin)
%task and model parameters
p = inputParser;
% default parameters if none given
p.addParameter('plotFlag', 1);
p.addParameter('maxTrial', 1000);
p.addParameter('binSize', 2)% in ms  
% p.addParameter('stepSize', 2)  
p.addParameter('tb', 4)% in s
p.addParameter('tf', 5)% in s
p.addParameter('saveFigFlag', 0);
p.addParameter('lag', 1000);% in ms
p.addParameter('jitter', 5); % in ms
p.addParameter('numShuffle', 500); % pairs of surrogates to generate
p.addParameter('alpha', 0.01);
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
%% find identified units, create spike cell

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
midPoints = (0.5*binSize + 1):stepSize:(length(time)-0.5*binSize);
slideTime = midPoints - p.Results.tb*1000;
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
allSpikesAllTime = {};
allSpikesPreSession = {};

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
    % slide window
    allTrial_spikeMatx_slide = zeros(length(sessionData), length(midPoints));
    for w = 1:length(midPoints)
        allTrial_spikeMatx_slide(:,w) = ...
            nansum(allTrial_spikeMatx_cue(:,midPoints(w)-0.5*binSize:midPoints(w)+0.5*p.Results.binSize-1),2)*1000/p.Results.binSize;
    end
    % currAllSpikes into bins
    currAllSpikes = sessionData(clust(i)-min(allClustsInds)+1).allSpikes;
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
    allSpikes = cat(3, allSpikes, (allTrial_spikeMatx_slide(:,find(slideTime>0,1):end)));
    allSpikesPre = cat(3, allSpikesPre, (allTrial_spikeMatx_slide(:,1:find(slideTime<0,1,'last'))));
    allSpikesAllTime = [allSpikesAllTime; currAllSpikesBin];
    allSpikesPreSession = [allSpikesPreSession; currAllSpikesPreSessBin];
end
%% calculate autoCorr

    lag = round(p.Results.lag/p.Results.binSize);
    x = allSpikes(:,:,1);
    xPre = allSpikesPre(:,:,1);
    xAll = allSpikesAllTime{1};
    xPreSess = allSpikesPreSession{1};

    y = allSpikes(:,:,2);
    yPre = allSpikesPre(:,:,2);
    yAll = allSpikesAllTime{2};
    yPreSess = allSpikesPreSession{2};
      % sep in trial and inter trial
    corrCo = [];
    for t = 1:size(x,1)
      mTemp = xcorr(x(t,:),y(t,:),lag,'normalized');
      mTemp(0.5*(length(mTemp)+1)) = NaN;
      corrCo = [corrCo; mTemp];
    end     
    corrCoPre = [];
    for t = 1:size(xPre,1)
      mTemp = xcorr(xPre(t,:),yPre(t,:),lag,'normalized');
      mTemp(0.5*(length(mTemp)+1)) = NaN;
      corrCoPre = [corrCoPre; mTemp];
    end     
    corrCoAll = xcorr(xAll,yAll,lag,'normalized');
    corrCoAll(0.5*(length(corrCoAll)+1)) = NaN;

    % combine all spikes preSession
    corrCoPreSess = xcorr(xPreSess,yPreSess,lag,'normalized');
    corrCoPreSess(0.5*(length(corrCoPreSess)+1)) = NaN;
    
    % consolidate into struct
    s.corrCo = mean(corrCo, 'omitnan');
    s.corrCoPre = mean(corrCoPre, 'omitnan');
    s.corrCoAll = corrCoAll;
    s.corrCoPreSess = corrCoPreSess;
    s.lag = (-lag:lag)*stepSize;
    
    
    %% generating surrogate
    % find recording start time: 
    pd = parseSessionString_df(session, root, sep);
    [timestamps, eventID, TTL, Evstring] = Nlx2MatEV([pd.nLynxFolderSession 'Events.nev'],[1 1 1 0 1], 0, 1);
    timestamps = round(timestamps/1000);
    startTime = min(timestamps);
    
    allPreTrial = NaN(p.Results.numShuffle, length(s.corrCo));
    allPreSession = NaN(p.Results.numShuffle, length(s.corrCo));
    binStartsPre = max(sessionData(1).CSon-20000, startTime):stepSize:sessionData(1).CSon;
    
    parfor a = 1:p.Results.numShuffle
        allSpikesPre = [];
        allSpikesPreSession = [];
        for i = 1:length(clust)
            % Initialize matrices for SDF
            allTrial_spikeMatx_cue = zeros(length(sessionData),length(time));
           for j = 1:length(allTrial_spike_cue)
                tempSpike = allTrial_spike_cue{i,j};
                tempSpike = tempSpike + p.Results.tb*1000; % add this to pad time for SDF
                % surrogate
                tempSpike = tempSpike + randi(2*p.Results.jitter+1, size(tempSpike,1), size(tempSpike,2)) - p.Results.jitter -1;
                tempSpike = tempSpike(tempSpike > 0 & tempSpike < length(time));
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
            allTrial_spikeMatx_slide = zeros(length(sessionData), length(midPoints));
            for w = 1:length(midPoints)
                allTrial_spikeMatx_slide(:,w) = ...
                    nansum(allTrial_spikeMatx_cue(:,midPoints(w)-0.5*binSize:midPoints(w)+0.5*p.Results.binSize-1),2)*1000/p.Results.binSize;
            end
            currAllSpikesPreSess = sessionData(clust(i)-min(allClustsInds)+1).allSpikes;
            currAllSpikesPreSess = currAllSpikesPreSess(currAllSpikesPreSess<sessionData(1).CSon & currAllSpikesPreSess>sessionData(1).CSon-20000);
            % surrogate
            currAllSpikesPreSess = currAllSpikesPreSess + randi(2*p.Results.jitter+1, size(currAllSpikesPreSess,1), size(currAllSpikesPreSess,2)) - p.Results.jitter - 1;
            currAllSpikesPreSessBin = histcounts(currAllSpikesPreSess,binStartsPre)*1000/binSize;
            % zscore spike matrix
            allSpikesPre(:,:,i) =  allTrial_spikeMatx_slide(:,1:find(slideTime<0,1,'last'));
            allSpikesPreSession(i,:) = currAllSpikesPreSessBin;
        end
        
        xPre = allSpikesPre(:,:,1);
        xPreSess = allSpikesPreSession(1,:);

        yPre = allSpikesPre(:,:,2);
        yPreSess = allSpikesPreSession(2,:);
   
        corrCoPre = [];
        for t = 1:size(xPre,1)
          mTemp = xcorr(xPre(t,:),yPre(t,:),lag,'normalized');
          mTemp(0.5*(length(mTemp)+1)) = NaN;
          corrCoPre = [corrCoPre; mTemp];
        end     
        
        corrCoPreSess = xcorr(xPreSess,yPreSess,lag,'normalized');
        corrCoPreSess(0.5*(length(corrCoPreSess)+1)) = NaN;
        
        allPreTrial(a,:) = mean(corrCoPre, 'omitnan');
        allPreSession(a,:) = corrCoPreSess;
    end
    % find alpha boundary
    CIPreTrial = NaN(2,length(s.corrCo));
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
%% plotting
if p.Results.plotFlag
    figure2;
    subplot(1,3,1); hold on;
    bar((-lag:lag)*stepSize, corrCoAll,'k');
    legend('allSpikes');
    
    subplot(1,3,2); hold on;
    bar((-lag:lag)*stepSize,mean(corrCo, 'omitnan'),'m', 'FaceAlpha', 0.7);
    bar((-lag:lag)*stepSize,mean(corrCoPre, 'omitnan'),'c', 'FaceAlpha', 0.7);
    plot((-lag:lag)*stepSize, mean(allPreTrial, 1, 'omitnan'), 'Color', 'c', 'LineWidth', 2);
    plot(repmat((-lag:lag)'*stepSize, 1, 2), CIPreTrial', 'Color', 'c', 'LineWidth', 2, 'LineStyle', '--');
    legend({'in trial', 'pre trial'})
    title([unit1 unit2 ' crossCorr'], 'Interpreter','none');

    subplot(1,3,3); hold on;
    bar((-lag:lag)*stepSize,corrCoPreSess,'k');
    plot((-lag:lag)*stepSize, mean(allPreSession, 1, 'omitnan'), 'Color', [0.7 0.7 0.7], 'LineWidth', 2);
    plot(repmat((-lag:lag)'*stepSize', 1, 2), CIPreSess', 'Color', [0.7 0.7 0.7], 'LineWidth', 2, 'LineStyle', '--');
    legend('preSession');
    % subplot(2, length(clust)+1, 1); hold on;
    % spikeTimes = sessionData(clust(1)-min(allClustsInds)+1).allSpikes;
    % histogram(diff(spikeTimes), 'Normalization', 'probability');
    % title('ISI')
    % xlabel('ms')
    % set(gca, 'XScale', 'log')
    % plot([2 2], [0 0.2], 'LineStyle', '--', 'Color', 'r');
    % xlim([0 max(diff(sessionData(clust(1)-min(allClustsInds)+1).allSpikes))])

    sgtitle([session ' ' unit1 unit2], 'Interpreter', 'none');
    screen = get(0,'Screensize');
    screen(4) = screen(4) - 100;
    set(gcf, 'Position', screen)
end
   
