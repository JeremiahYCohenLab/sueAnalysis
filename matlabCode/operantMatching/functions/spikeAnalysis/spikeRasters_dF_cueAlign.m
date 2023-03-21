  function spikeRasters_dF_cueAlign(sessionName, varargin)

p = inputParser;
% default parameters if none given
p.addParameter('cellName', 'all');
p.addParameter('saveFigFlag', 1)
p.addParameter('intanFlag',0)
p.addParameter('revForFlag',0)
p.addParameter('modelsFlag',0)
p.addParameter('timeMax', 121000)
p.addParameter('timeBins', 12)
p.addParameter('tb', 1.5);
p.addParameter('tf', 5); % in s
p.addParameter('binSize', 150); %in ms
p.addParameter('stepSize', 100); % in ms
p.parse(varargin{:});

cellName = p.Results.cellName;

% Path
[root,sep] = currComputer();

[animalName] = strtok(sessionName, 'd');
animalName = animalName(2:end);

if isstrprop(sessionName(end), 'alpha')
    sortedFolderLocation = [root animalName sep sessionName(1:end-1) sep 'sorted' sep 'session ' sessionName(end) sep];
    savepath = [root animalName sep sessionName(1:end-1) sep  'figures' sep 'session ' sessionName(end) sep];
else
    sortedFolderLocation = [root animalName sep sessionName sep 'sorted' sep 'session' sep];
    savepath = [root animalName sep sessionName sep  'figures' sep];
end

    

if p.Results.intanFlag == 1
    if exist([sortedFolderLocation sessionName '_sessionData_intan.mat'],'file')
        load([sortedFolderLocation sessionName '_sessionData_intan.mat'])
    else
        [sessionData] = generateSessionData_nL_operantMatching(sessionName);
    end
else
    if exist([sortedFolderLocation sessionName '_sessionData_nL.mat'],'file')
        load([sortedFolderLocation sessionName '_sessionData_nL.mat'])
    else
        [sessionData] = generateSessionData_nL_operantMatching(sessionName);
    end
end

[s] = behAnalysisNoPlot_opMD(sessionName, 'revForFlag', p.Results.revForFlag);

if isempty(dir(savepath))
    mkdir(savepath)
end


time = -1000*p.Results.tb:1000*p.Results.tf;

omitInds = isnan([sessionData.rewardTime]);
tempBlockSwitch = s.blockSwitch;
blockSwitch_CSminCorrected = s.blockSwitch;
for i = 2:length(blockSwitch_CSminCorrected)
    subVal = sum(omitInds(tempBlockSwitch(i-1):tempBlockSwitch(i)));
    blockSwitch_CSminCorrected(i:end) = blockSwitch_CSminCorrected(i:end) - subVal;
end

CSminus_Inds = find(strcmp({sessionData.trialType},'CSminus')>0);
CSplus_Inds = find(strcmp({sessionData.trialType},'CSplus')>0);
%% Sort all spikes into a raster-able matrix
spikeFields = fields(sessionData);
if iscell(cellName)
    for i = 1:length(cellName)
        clust(i) = find(~cellfun(@isempty,strfind(spikeFields,cellName{i})));
    end
elseif regexp(cellName, 'all')
        clust = find(~cellfun(@isempty,strfind(spikeFields,'C_')) | ~cellfun(@isempty,strfind(spikeFields,'TT')));
else
    clust = find(~cellfun(@isempty,strfind(spikeFields,cellName)));
end
    
allTrial_spike_choice = {};
lickInds = [];
currTrial_lick = [];
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
        
        allTrial_spike_choice{i,k} = [prevTrial_spike currTrial_spike];

    end
    
    if ~isnan(sessionData(k).rewardL)
        currTrial_lickInd = [sessionData(k).licksL] < (sessionData(k).CSon + p.Results.tf*1000);
        currTrial_lick = sessionData(k).licksL(currTrial_lickInd) - sessionData(k).CSon;
        lickInds = [lickInds k];
    elseif ~isnan(sessionData(k).rewardR)
        currTrial_lickInd = [sessionData(k).licksR] < (sessionData(k).CSon + p.Results.tf*1000);
        currTrial_lick = sessionData(k).licksR(currTrial_lickInd) - sessionData(k).CSon;  
        lickInds = [lickInds k];
    elseif ismember(k, CSminus_Inds)
        templick = [[sessionData(k).licksL], [sessionData(k).licksR]];
        if min(templick)<sessionData(k).CSon + 1800 %1800 is resp window
            lickInds = [lickInds k];
            if min([sessionData(k).licksL]) < min([sessionData(k).licksR])
                templick = [sessionData(k).licksL];
            else
                templick = [sessionData(k).licksR];
            end
        currTrial_lickInd = templick < (sessionData(k).CSon + p.Results.tf*1000);
        currTrial_lick = templick(currTrial_lickInd) - sessionData(k).CSon;  
        end
    else
        currTrial_lick = 0;
    end
    allTrial_lick{k} = [currTrial_lick];
end

% sometimes no licks/spikes are considered 1x0 and sometimes they are []
% plotSpikeRaster does not place nicely with [] so this converts all empty indices to 1x0
allTrial_spike_choice(cellfun(@isempty,allTrial_spike_choice)) = {zeros(1,0)}; 



%% Generate indices for block switches
blockS = [];
blockL = [];
blockR = [];
s.blockSwitch = s.blockSwitch + 1;
s.blockSwitch(1) = 1; %% this is to avoid an error on the first iteration of the loop because it's trying to index at 0

if p.Results.revForFlag
    for i = 1:length(s.blockSwitch)
        [rewardProbL, rewardProbR] = strtok(s.blockProbs(i), '/');
        rewardProbL = str2double(rewardProbL); rewardProbR = str2double(rewardProbR{1}(2:end));
        if rewardProbL > rewardProbR
            if i ~= length(s.blockSwitch)
                blockL = [blockL s.blockSwitch(i):s.blockSwitch(i+1)-1];
            else
                blockL = [blockL s.blockSwitch(i):length(s.responseInds)];
            end
        else
            if i ~= length(s.blockSwitch)
                blockR = [blockR s.blockSwitch(i):s.blockSwitch(i+1)-1];
            else
                blockR = [blockR s.blockSwitch(i):length(s.responseInds)];
            end
        end
    end
else
    for i = 1:length(s.blockSwitch)
        if s.behSessionData(s.blockSwitch(i)).rewardProbL == s.behSessionData(s.blockSwitch(i)).rewardProbR   
            if i ~= length(s.blockSwitch)
                blockS = [blockS s.blockSwitch(i):s.blockSwitch(i+1)-1];
            else
                blockS = [blockS s.blockSwitch(i):length(s.responseInds)];
            end
        elseif s.behSessionData(s.blockSwitch(i)).rewardProbL > s.behSessionData(s.blockSwitch(i)).rewardProbR
            if i ~= length(s.blockSwitch)
                blockL = [blockL s.blockSwitch(i):s.blockSwitch(i+1)-1];
            else
                blockL = [blockL s.blockSwitch(i):length(s.responseInds)];
            end
        else
            if i ~= length(s.blockSwitch)
                blockR = [blockR s.blockSwitch(i):s.blockSwitch(i+1)-1];
            else
                blockR = [blockR s.blockSwitch(i):length(s.responseInds)];
            end
        end
    end
end

%set max trial number for plotting
yLimMax = max([length(blockR) length(blockL) length(s.lickR_Inds) length(s.lickL_Inds)]);

%% smooth rewards over time
sessionTime = [sessionData(1).CSon:sessionData(end).CSon + 3000] - sessionData(1).CSon;       %pad time for reward on last trial

sessionRwds = [sessionData.rewardTime] - sessionData(1).CSon;     %baseline to start time
session_rwdsArray = zeros(1,length(sessionTime));
sessionRwds = sessionRwds(s.responseInds(s.rwd_Inds));
session_rwdsArray(round(sessionRwds)) = 1;

boxKern = ones(1,60000);                                       %smooth rewards over time
sessionRwdsSmooth = conv(session_rwdsArray, boxKern);
sessionRwdsSmooth = sessionRwdsSmooth(1:(end-(length(boxKern)-1)));


%% Generate smoothed choice-history values       


%% Generate smoothed reward-history values over trials

[~,rwdHx_Inds] = sort(s.rwdHx);

%outcome indices for rwd hist
[rwdHxRwd_Inds] = ismember(rwdHx_Inds, s.rwd_Inds); 
[rwdHxNoRwd_Inds] = ismember(rwdHx_Inds, s.nrwd_Inds); 

%for tercile analysis
tercile = floor(length(rwdHx_Inds)/3);
rwdHxI_Inds = rwdHx_Inds(1:tercile);
rwdHxII_Inds = rwdHx_Inds(tercile+1:tercile*2);
rwdHxIII_Inds = rwdHx_Inds(tercile*2+1:end);


%outcome indices for rwd hist divisions
rwdHxIrwd_Inds = s.responseInds(intersect(rwdHxI_Inds, s.rwd_Inds));          rwdHxInoRwd_Inds = s.responseInds(intersect(rwdHxI_Inds, s.nrwd_Inds));
rwdHxIIrwd_Inds = s.responseInds(intersect(rwdHxII_Inds, s.rwd_Inds));          rwdHxIInoRwd_Inds = s.responseInds(intersect(rwdHxII_Inds, s.nrwd_Inds));
rwdHxIIIrwd_Inds = s.responseInds(intersect(rwdHxIII_Inds, s.rwd_Inds));          rwdHxIIInoRwd_Inds = s.responseInds(intersect(rwdHxIII_Inds, s.nrwd_Inds));

%choice indices for rwd hist divisions
% [~,rwdHxIR_Inds,~] = intersect(rwdHxI_Inds, lickR_Inds);          [~,rwdHxIL_Inds,~] = intersect(rwdHxI_Inds, lickL_Inds);
% [~,rwdHxIIR_Inds,~] = intersect(rwdHxII_Inds, lickR_Inds);          [~,rwdHxIIL_Inds,~] = intersect(rwdHxII_Inds, lickL_Inds);
% [~,rwdHxIIIR_Inds,~] = intersect(rwdHxIII_Inds, lickR_Inds);          [~,rwdHxIIIL_Inds,~] = intersect(rwdHxIII_Inds, lickL_Inds);

% divide lick latency into three percentile
[~,lickLat_Inds] = sort(s.lickLatZ);

latI_Inds = lickLat_Inds(1:tercile);
latII_Inds = lickLat_Inds(tercile+1:tercile*2);
latIII_Inds = lickLat_Inds(tercile*2+1:end);

latI_Inds = s.responseInds(latI_Inds); 
latII_Inds = s.responseInds(latII_Inds);
latIII_Inds = s.responseInds(latIII_Inds);
%% generate indeces for Q learning model values

if p.Results.modelsFlag
    Qmdls = qLearning_fit2LR([sessionName '.asc']);
    pe = [Qmdls.fourParams_twoLearnRates_tForget.pe];
    peSD = movstd(pe, 10);
    %rBar = [Qmdls.fiveParams_opponency.rBar];

    [pe_Sorted,pe_Inds] = sort(pe);
    [peSD_Sorted,peSD_Inds] = sort(peSD(1:end-1));

    peSD_Inds = peSD_Inds + 1;
    pe_Lim = find(pe_Sorted > 0, 1);
    
    
    [rBar_sorted rBar_Inds] = sort(rBar);
    tercile = floor(length(rBar_Inds)/3);
    rBarI_Inds = rBar_Inds(1:tercile);
    rBarII_Inds = rBar_Inds(tercile+1:tercile*2);
    rBarIII_Inds = rBar_Inds(tercile*2+1:end);
end





%% sort licking for sdfs

for j = 1:length(sessionData)
    trialDurDiff(j) = (sessionData(j).trialEnd - sessionData(j).CSon)- p.Results.tf*1000;
end
trialDurDiff(end) = 0;  %to account for no trialEnd timestamp on last trial

%initialize lick matrices
allTrial_lickMatx = NaN(length(sessionData),length(time)); 

for j = 1:length(allTrial_lick)
    tempLick = allTrial_lick{1,j};
    tempLick = tempLick + p.Results.tb*1000; % add this to pad time for SDF
    allTrial_lickMatx(j,tempLick) = 1;
    if trialDurDiff(j) < 0
        allTrial_lickMatx(j, isnan(allTrial_lickMatx(j, 1:end+trialDurDiff(j)))) = 0;  %converts within trial duration NaNs to 0's
    else
        allTrial_lickMatx(j, isnan(allTrial_lickMatx(j,:))) = 0;
    end
    if sum(allTrial_lickMatx(j,:)) == 0     %if there is no spike data for this trial, don't count it
        allTrial_lickMatx(j,:) = NaN;
    end
end

% calculate slide window lickRate
midPoints = (0.5*p.Results.binSize + 1):p.Results.stepSize:(length(time)-0.5*p.Results.binSize);
slideTime = midPoints - p.Results.tb*1000;
allTrial_lickMatx_slide = zeros(length(sessionData), length(midPoints));
for w = 1:length(midPoints)
    allTrial_lickMatx_slide(:,w) = ...
        sum(allTrial_lickMatx(:,midPoints(w)-0.5*p.Results.binSize:midPoints(w)+0.5*p.Results.binSize-1),2)*1000/p.Results.binSize;
    
end

%sort by receipt of reward
    rwd_lickMatx = allTrial_lickMatx_slide(s.responseInds(s.rwd_Inds),:);
    noRwd_lickMatx = allTrial_lickMatx_slide(s.responseInds(s.nrwd_Inds),:);


%sort by lick choice
    choiceR_lickMatx = allTrial_lickMatx_slide(s.responseInds(s.lickR_Inds),:);
    choiceL_lickMatx = allTrial_lickMatx_slide(s.responseInds(s.lickL_Inds),:);


%sort licks by reward history
    rwdHxI_lickMatx = allTrial_lickMatx_slide(s.responseInds(rwdHxI_Inds),:);
    rwdHxII_lickMatx = allTrial_lickMatx_slide(s.responseInds(rwdHxII_Inds),:);
    rwdHxIII_lickMatx = allTrial_lickMatx_slide(s.responseInds(rwdHxIII_Inds),:);



%% sort spikes for trial averages
for i = 1:length(clust)
    % Initialize matrices for SDF
    
    allTrial_spikeMatx_choice = zeros(length(sessionData),length(time));         
  
   for j = 1:length(allTrial_spike_choice)
        tempSpike = allTrial_spike_choice{i,j};
        tempSpike = tempSpike + p.Results.tb*1000; % add this to pad time for SDF
        if any(tempSpike == 0)
            tempSpike(tempSpike == 0) = 1; % this avoids error in next line
        end
        allTrial_spikeMatx_choice(j,tempSpike) = 1;
        if trialDurDiff(j) < 0
            allTrial_spikeMatx_choice(j, isnan(allTrial_spikeMatx_choice(j, 1:end+trialDurDiff(j)))) = 0;  %converts within trial duration NaNs to 0's
        else
            allTrial_spikeMatx_choice(j, isnan(allTrial_spikeMatx_choice(j,:))) = 0;
        end
    end
    % slide window
    allTrial_spikeMatx_slide = zeros(length(sessionData), length(midPoints));
    for w = 1:length(midPoints)
        allTrial_spikeMatx_slide(:,w) = ...
            nansum(allTrial_spikeMatx_choice(:,midPoints(w)-0.5*p.Results.binSize:midPoints(w)+0.5*p.Results.binSize-1),2)*1000/p.Results.binSize;
    end
    
    rwd_spikeMatx = allTrial_spikeMatx_slide(s.responseInds(s.rwd_Inds),:);
    noRwd_spikeMatx = allTrial_spikeMatx_slide(s.responseInds(s.nrwd_Inds),:);
    rwdPrev_spikeMatx = allTrial_spikeMatx_slide(s.responseInds(s.rwd_Inds(1:end-1)+1),:);
    noRwdPrev_spikeMatx = allTrial_spikeMatx_slide(s.responseInds(s.nrwd_Inds(1:end-1)+1),:);
    changeChoice_spikeMatx = allTrial_spikeMatx_slide(s.responseInds(s.changeChoice_Inds),:);
    stayChoice_spikeMatx = allTrial_spikeMatx_slide(s.responseInds(s.stayChoice_Inds),:);
  
   
    %sort spikes by reward history
    rwdHxI_spikeMatx = allTrial_spikeMatx_slide(s.responseInds(rwdHxI_Inds),:);
    rwdHxII_spikeMatx = allTrial_spikeMatx_slide(s.responseInds(rwdHxII_Inds),:);
    rwdHxIII_spikeMatx = allTrial_spikeMatx_slide(s.responseInds(rwdHxIII_Inds),:);
    % Downsample raster if necessary
    tempSpike = [sessionData.(spikeFields{clust(i)})];
    avgFiringRate = length(tempSpike)/((tempSpike(end)-tempSpike(1))/1000);

    downsampFlag = false;
    if avgFiringRate >= 12
        downsampFactor = floor(avgFiringRate/6);
        downsampFlag = true;
        for j = 1:size(allTrial_spike_choice,2)
            if ~isempty(allTrial_spike_choice{i,j})
              allTrial_spike_choice{i,j} = downsample(allTrial_spike_choice{i,j}, downsampFactor);
            end
        end   
    end
 
    
    %% lick latency and recent rwd hist comparison
    
    smoothWin = 250;
    trialBeg = p.Results.tb*1000;
    trialEnd = p.Results.tb*1000 + 1500;
    
    %% Plotting
    
    maxFreq = 1.2*max(nanmean(allTrial_spikeMatx_slide), [], 'all');
    maxLick = 1.2*max(nanmean(allTrial_lickMatx_slide), [], 'all');
    
    screenSize = get(0,'Screensize');
    screenSize(4) = screenSize(4) - 100;
    rasters = figure; hold on
    axisColor = [0 0 0];
    set(rasters,'defaultAxesColorOrder',[axisColor; axisColor]);
    set(rasters, 'Position', screenSize)

    % All trials
    r(1) = subplot(8,7,[1 8 15 22]); t(1) = title('All Trials');
    plotSpikeRaster(allTrial_spike_choice(i,:),'PlotType','vertline'); hold on
    plot(repmat([-5000 10000],length(s.blockSwitch),1)', [s.blockSwitch, s.blockSwitch]','r');
    line([0 0], [0 length(sessionData)], 'color', 'r')
    
    
    % Block_R trials
    r(4) = subplot(8,7,[2 9]); t(4) = title('higher prob on R');
    if ~all(cellfun(@isempty,allTrial_spike_choice(i,s.responseInds(blockR))))
        plotSpikeRaster(allTrial_spike_choice(i,s.responseInds(blockR)),'PlotType','vertline'); hold on
        plot([-5000 5000],[length(blockR) length(blockR)],'Color',[192 192 192]/255);
        bR_switch = find(diff(blockR) > 1);
        line([0 0], [0 length(blockR)], 'color', 'r');
        plot(repmat([-5000 10000],length(bR_switch),1)', [bR_switch; bR_switch],'r');
        set(gca,'Xticklabel',[]);
    end
    
    
    % Block_L trials
    r(5) = subplot(8,7,[16 23]); t(5) = title('higher prob on L');
    if ~all(cellfun(@isempty,allTrial_spike_choice(i,s.responseInds(blockL))))
        plotSpikeRaster(allTrial_spike_choice(i,s.responseInds(blockL)),'PlotType','vertline'); hold on
        plot([-5000 5000],[length(blockL) length(blockL)],'Color',[192 192 192]/255);
        bL_switch = find(diff(blockL) > 1);
        plot(repmat([-5000 10000],length(bL_switch),1)', [bL_switch; bL_switch],'r');
        line([0 0], [0 length(blockL)], 'color', 'r');
        set(gca,'Xticklabel',[]);
    end
    
    % R lick broken by rwd vs no rwd
    r(6) = subplot(8,7,[3 10]); t(6) = title('R lick; no rwd v. rwd');
    R_rwd = s.responseInds(intersect(s.lickR_Inds,s.rwd_Inds));
    R_norwd = s.responseInds(intersect(s.lickR_Inds,s.nrwd_Inds));
    if ~all(cellfun(@isempty,allTrial_spike_choice(i,[R_norwd R_rwd])))
        plotSpikeRaster(allTrial_spike_choice(i,[R_norwd R_rwd]),'PlotType','vertline'); hold on
        plot([-5000 10000],[length(R_norwd) length(R_norwd)],'r');
        plot([-5000 10000],[length([R_norwd R_rwd]) length([R_norwd R_rwd])],'Color',[192 192 192]/255);
        line([0 0], [0 length(s.lickR_Inds)], 'color', 'r');
        set(gca,'Xticklabel',[]);
    end
    
    % L lick broken by rwd vs no rwd
    r(7) = subplot(8,7,[17 24]); t(7) = title('L lick; no rwd v. rwd');
    L_rwd = s.responseInds(intersect(s.lickL_Inds,s.rwd_Inds));
    L_norwd = s.responseInds(intersect(s.lickL_Inds,s.nrwd_Inds));
    if ~all(cellfun(@isempty,allTrial_spike_choice(i,[L_norwd L_rwd])))
        plotSpikeRaster(allTrial_spike_choice(i,[L_norwd L_rwd]),'PlotType','vertline'); hold on
        plot([-5000 10000],[length(L_norwd) length(L_norwd)],'r');
        plot([-5000 10000],[length([L_norwd L_rwd]) length([L_norwd L_rwd])],'Color',[192 192 192]/255);
        line([0 0], [0 length(s.lickL_Inds)], 'color', 'r');
        set(gca,'Xticklabel',[]);
    end
    
    % R lick sorted by latency
    r(8) = subplot(8,7,[4 11]); t(8) = title('R licks');
    if ~all(cellfun(@isempty,allTrial_spike_choice(i,s.responseInds(s.lickR_Inds))))
        plotSpikeRaster(allTrial_spike_choice(i,s.responseInds(s.lickR_Inds)),'PlotType','vertline'); hold on
        plot([-5000 10000],[length(s.lickR_Inds) length(s.lickR_Inds)],'Color',[192 192 192]/255)
        line([0 0], [0 length(s.lickR_Inds)], 'color', 'r')
        set(gca,'Xticklabel',[]);
    end
        
    % L lick sorted by latency
    r(9) = subplot(8,7,[18 25]); t(9) = title('L licks');
    if ~all(cellfun(@isempty,allTrial_spike_choice(i,s.responseInds(s.lickL_Inds))))
        plotSpikeRaster(allTrial_spike_choice(i,s.responseInds(s.lickL_Inds)),'PlotType','vertline'); hold on  
        plot([-5000 10000],[length(s.lickL_Inds) length(s.lickL_Inds)],'Color',[192 192 192]/255)
        line([0 0], [0 length(s.lickL_Inds)], 'color', 'r')
        set(gca,'Xticklabel',[]);
    end
    
    % raster  switch vs stay
    r(10) = subplot(8,7,[5 12]); t(10) = title('switch vs stay');
    if ~all(cellfun(@isempty,allTrial_spike_choice(i,s.changeChoice_Inds)))
        plotSpikeRaster(allTrial_spike_choice(i,s.responseInds([s.changeChoice_Inds s.stayChoice_Inds])),'PlotType','vertline'); hold on
        plot([-5000 10000],[length(s.changeChoice_Inds) length(s.changeChoice_Inds)],'color', 'r')
        line([0 0], [0 length(s.responseInds)], 'color', 'r')
        title('switch vs stay');
        set(gca,'Xticklabel',[]);
    end
    
    % raster explore vs exploit
    r(11) = subplot(8,7,[19 26]); t(11) = title('explore vs exploit');
    if ~all(cellfun(@isempty,allTrial_spike_choice(i,s.hmmStates==1)))
        plotSpikeRaster(allTrial_spike_choice(i,s.responseInds([find(s.hmmStates==1) find(s.hmmStates~=1)])),'PlotType','vertline'); hold on
        plot([-5000 10000],[sum(s.hmmStates==1) sum(s.hmmStates==1)],'color', 'r')
        line([0 0], [0 length(s.responseInds)], 'color', 'r')
        set(gca,'Xticklabel',[]);
    end

    % rwd-history by outcome
    r(13) = subplot(8,7,[6 13]); t(13) = title('rwd hist L->H; no rwd v. rwd');
    rwdHxOutcome = [allTrial_spike_choice(i,s.responseInds(rwdHx_Inds(rwdHxNoRwd_Inds))) allTrial_spike_choice(i,s.responseInds(rwdHx_Inds(rwdHxRwd_Inds)))];
    plotSpikeRaster(rwdHxOutcome,'PlotType','vertline'); hold on
    plot([-5000 10000],[length(rwdHxNoRwd_Inds) length(rwdHxNoRwd_Inds)],'-r')
    line([0 0], [0 length(s.responseInds)], 'color', 'r')
    set(gca,'Xticklabel',[]);

    
%     % sorted by prediction error
%     r(15) = subplot(8,7,[48 55]); t(15) = title('pe L -> H');
%     plotSpikeRaster(allTrial_spike(i,pe_Inds),'PlotType','vertline'); hold on
%     plot(repmat([-5000 10000],length(pe_Lim),1)', [pe_Lim; pe_Lim],'r')
%     set(gca,'Xticklabel',[]);
%     
%     % sorted by prediction error variance
%     r(16) = subplot(8,7,[49 56]); t(16) = title('pe sd L -> H');
%     plotSpikeRaster(allTrial_spike(i,peSD_Inds),'PlotType','vertline'); hold on
%     set(gca,'Xticklabel',[]);
    
    % rwd/noRwd and lick rate SDFs
    z(2) = subplot(8,7,[29 36]); hold on
    mySDF_rwd = allTrial_spikeMatx_slide(s.responseInds(s.rwd_Inds),:);
    mySDF_noRwd = allTrial_spikeMatx_slide(s.responseInds(s.nrwd_Inds),:);
    plotFilled(slideTime, mySDF_rwd,[0 0 1])
    plotFilled(slideTime, mySDF_noRwd,[0.7 0 1])
    line([0 0], [0 maxFreq], 'color', 'r')
    ylim([0 maxFreq]);
    legend({'rwd','','no Rwd',''},'FontSize',6,'Location','northeast')
    
    % rdwHx terciles SDFs
    z(3) = subplot(8,7,[31 38]); t(5) = title('rwd hist'); hold off
    mySDF_rwdHxIlicks = rwdHxI_lickMatx; hold on;
    mySDF_rwdHxIIlicks = rwdHxII_lickMatx;
    mySDF_rwdHxIIIlicks = rwdHxIII_lickMatx;
    mySDF_rwdHxI = rwdHxI_spikeMatx;
    mySDF_rwdHxII = rwdHxII_spikeMatx;
    mySDF_rwdHxIII = rwdHxIII_spikeMatx;
    yyaxis left
    plotFilled(slideTime, mySDF_rwdHxIII, [0 0 1])
    plotFilled(slideTime, mySDF_rwdHxII,[0.4 0.4 1])
    plotFilled(slideTime, mySDF_rwdHxI, [0.8 0.8 1])
    line([s.rwdDelay s.rwdDelay], [0 maxFreq], 'color', 'r')
    ylim([0 maxFreq]);
    yyaxis right
    plotFilled(slideTime, mySDF_rwdHxIIIlicks,[0 0 0])
    plotFilled(slideTime, mySDF_rwdHxIIlicks, [0.4 0.4 0.4])
    plotFilled(slideTime, mySDF_rwdHxIlicks, [0.8 0.8 0.8])
    ylim([0 5*maxLick]);
    legend({'High','','Middle','','Low',''},'FontSize',5,'Location','northeast')
    set(gca,'Xticklabel',[]);
    
    
    
    
    %sdf by reward hist, rewarded trials
    z(4) = subplot(8,7,32); t(5) = title('rwd hist - rwd'); hold on
    mySDF_rwdHxIrwd = allTrial_spikeMatx_slide(rwdHxIrwd_Inds,:);
    mySDF_rwdHxIIrwd = allTrial_spikeMatx_slide(rwdHxIIrwd_Inds,:);
    mySDF_rwdHxIIIrwd = allTrial_spikeMatx_slide(rwdHxIIIrwd_Inds,:);
    plotFilled(slideTime, mySDF_rwdHxIIIrwd, [0 0 1])
    plotFilled(slideTime, mySDF_rwdHxIIrwd, [0.4 0.4 1])
    plotFilled(slideTime, mySDF_rwdHxIrwd, [0.8 0.8 1])
    line([0 0], [0 maxFreq], 'color', 'r')
    ylim([0 maxFreq]);
    
    %sdf by reward hist, NON-rewarded trials
    z(5) = subplot(8,7,39); t(6) = title('rwd hist - no rwd'); hold on
    mySDF_rwdHxInoRwd = allTrial_spikeMatx_slide(rwdHxInoRwd_Inds,:);
    mySDF_rwdHxIInoRwd = allTrial_spikeMatx_slide(rwdHxIInoRwd_Inds,:);
    mySDF_rwdHxIIInoRwd = allTrial_spikeMatx_slide(rwdHxIIInoRwd_Inds,:);
%    mySDF_rwdHxIVnoRwd = fastsmooth(nanmean(rwdHxIV_spikeMatx(rwdHxIVnoRwd_Inds,:), 1)*1000, 250);
%    plot(time, mySDF_rwdHxIVnoRwd(1:length(time)),'g','LineWidth',2)
    plotFilled(slideTime, mySDF_rwdHxIIInoRwd, [1 0 0])
    plotFilled(slideTime, mySDF_rwdHxIInoRwd, [1 0.4 0.4])
    plotFilled(slideTime, mySDF_rwdHxInoRwd, [1 0.8 0.8])
    line([0 0], [0 maxFreq], 'color', 'r')
    ylim([0 maxFreq]);
 
    % sdf sep by lick latency
    subplot(8,7,[33, 40]); t(6) = title('lick latency'); hold on
    mySDF_lickLatI = allTrial_spikeMatx_slide(latI_Inds,:);
    mySDF_lickLatII = allTrial_spikeMatx_slide(latII_Inds,:);
    mySDF_lickLatIII = allTrial_spikeMatx_slide(latIII_Inds,:);
%    mySDF_rwdHxIVnoRwd = fastsmooth(nanmean(rwdHxIV_spikeMatx(rwdHxIVnoRwd_Inds,:), 1)*1000, 250);
%    plot(time, mySDF_rwdHxIVnoRwd(1:length(time)),'g','LineWidth',2)
    plotFilled(slideTime, mySDF_lickLatI, [0 0 1])
    plotFilled(slideTime, mySDF_lickLatII, [0.4 0.4 1])
    plotFilled(slideTime, mySDF_lickLatIII, [0.8 0.8 1])
    line([0 0], [0 maxFreq], 'color', 'r')
    ylim([0 maxFreq]);
    legend('short','median','long');
    
    
    
    
    z(8) = subplot(8,7,[45 52]); hold on;
    mySDF_stay = stayChoice_spikeMatx;
    mySDF_change = changeChoice_spikeMatx;
    plotFilled(slideTime, mySDF_change,[0 0 1])
    plotFilled(slideTime, mySDF_stay,[0.7 0 1])
    line([0 0], [0 maxFreq], 'color', 'r')
    legend('switch', '', 'stay', '')
    title('stay vs switch')
    ylim([0 maxFreq]);


    z(9) = subplot(8,7,[46 53]); title('switch/stay rwd'); hold on;
    mySDF_stayRwd = allTrial_spikeMatx_slide(s.responseInds(intersect(s.rwd_Inds, s.stayChoice_Inds)),:);
    mySDF_switchRwd = allTrial_spikeMatx_slide(s.responseInds(intersect(s.rwd_Inds, s.changeChoice_Inds)),:);
    plotFilled(slideTime, mySDF_switchRwd,[0 0 1])
    plotFilled(slideTime, mySDF_stayRwd,[0.7 0 1])    
    line([0 0], [0 maxFreq], 'color', 'r')    
    legend({'switch-rwd','','stay-rwd',''},'FontSize',5)
    ylim([0 maxFreq]);
    
    z(10) = subplot(8,7,[47 54]); title('switch/stay noRwd'); hold on;
    mySDF_stayNoRwd = allTrial_spikeMatx_slide(s.responseInds(intersect(s.nrwd_Inds, s.stayChoice_Inds)),:);
    mySDF_switchNoRwd = allTrial_spikeMatx_slide(s.responseInds(intersect(s.nrwd_Inds, s.changeChoice_Inds)),:);
    plotFilled(slideTime, mySDF_switchNoRwd,[0 0 1]) 
    plotFilled(slideTime, mySDF_stayNoRwd,[0.7 0 1]) 
    line([0 0], [0 maxFreq], 'color', 'r')    
    legend({'switch-nrwd','','stay-nrwd',''},'FontSize',5)
    ylim([0 maxFreq]);
    
    z(11) = subplot(8,7,[48 55]); title('switch pre-rwd/noRwd'); hold on;
    mySDF_switchRwd = allTrial_spikeMatx_slide(s.responseInds(intersect(s.rwd_Inds + 1, s.changeChoice_Inds)),:);
    mySDF_switchNoRwd = allTrial_spikeMatx_slide(s.responseInds(intersect(s.nrwd_Inds + 1, s.changeChoice_Inds)),:);
    if size(mySDF_switchRwd,1) > 2
        plotFilled(slideTime, mySDF_switchRwd,[0 0 1])
    else
        if size(mySDF_switchRwd,1) == 1
            plot(slideTime, mySDF_switchRwd)
        end        
    end
    plotFilled(slideTime, mySDF_switchNoRwd,[0.8 0 1]) 
    line([0 0], [0 maxFreq], 'color', 'r')  
    if size(mySDF_switchRwd,1) > 0
        legend({'switch-prerwd','','switch-prenrwd',''},'FontSize',5)
    else
        legend({'switch-prenrwd',''},'FontSize',5)
    end
    ylim([0 maxFreq]);
    
    z(12) = subplot(8,7,[49 56]); title('stay pre-rwd/noRwd'); hold on;
    mySDF_stayRwd = allTrial_spikeMatx_slide(s.responseInds(intersect(s.rwd_Inds + 1, s.stayChoice_Inds)),:);
    mySDF_stayNoRwd = allTrial_spikeMatx_slide(s.responseInds(intersect(s.nrwd_Inds + 1, s.stayChoice_Inds)),:);
    plotFilled(slideTime, mySDF_stayRwd,[0 0 1])   
    plotFilled(slideTime, mySDF_stayNoRwd,[0.9 0 1]) 
    line([0 0], [0 maxFreq], 'color', 'r')    
    legend({'stay-prerwd','','stay-prenrwd',''},'FontSize',5)
    ylim([0 maxFreq]);
   
    
    z(13) = subplot(8,7,[34 41]); title('explore vs exploit'); hold on;
    mySDF_ore = allTrial_spikeMatx_slide(s.responseInds(s.hmmStates==1),:);
    mySDF_oit = allTrial_spikeMatx_slide(s.responseInds(s.hmmStates~=1),:);
    plotFilled(slideTime, mySDF_ore,[0 0 1])
    plotFilled(slideTime, mySDF_oit,[0.7 0 1])
    line([0 0], [0 maxFreq], 'color', 'r')
    legend('explore', '', 'exploit', '')
    ylim([0 maxFreq]);


    z(14) = subplot(8,7,[35 42]); title('explore rwd/nrwd'); hold on;
    mySDF_exploreRwd = allTrial_spikeMatx_slide(s.responseInds(intersect(s.rwd_Inds, find(s.hmmStates==1))),:);
    mySDF_exploreNoRwd = allTrial_spikeMatx_slide(s.responseInds(intersect(s.nrwd_Inds, find(s.hmmStates==1))),:);
    plotFilled(slideTime, mySDF_exploreRwd,[0.7 0 1])   
    plotFilled(slideTime, mySDF_exploreNoRwd,[0 0 1]) 
    line([0 0], [0 maxFreq], 'color', 'r')    
    legend({'explore-rwd','','explore-nrwd',''},'FontSize',5)
    ylim([0 maxFreq]);
    
    z(15) = subplot(8,7,[21 28]); title('exploit rwd/nrwd'); hold on;
    mySDF_exploitRwd = allTrial_spikeMatx_slide(s.responseInds(intersect(s.rwd_Inds, find(s.hmmStates~=1))),:);
    mySDF_exploitNoRwd = allTrial_spikeMatx_slide(s.responseInds(intersect(s.nrwd_Inds, find(s.hmmStates~=1))),:);
    plotFilled(slideTime, mySDF_exploitRwd,[0.7 0 1])
    plotFilled(slideTime, mySDF_exploitNoRwd,[0 0 1])
    line([0 0], [0 maxFreq], 'color', 'r')    
    legend({'exploit-rwd','','exploit-nrwd',''},'FontSize',5)
    ylim([0 maxFreq]);
    
    z(16) = subplot(8,7,[7 14]); title('exlpore/exploit rwd'); hold on;
    mySDF_exploreRwd = allTrial_spikeMatx_slide(s.responseInds(intersect(s.rwd_Inds, find(s.hmmStates==1))),:);
    mySDF_exploitRwd = allTrial_spikeMatx_slide(s.responseInds(intersect(s.rwd_Inds, find(s.hmmStates~=1))),:);
    plotFilled(slideTime, mySDF_exploreRwd,[0 0 1])
    plotFilled(slideTime, mySDF_exploitRwd,[0.7 0 1])    
    line([0 0], [0 maxFreq], 'color', 'r')    
    legend({'explore-rwd','','exploit-rwd',''},'FontSize',5)
    ylim([0 maxFreq]);
    
    z(17) = subplot(8,7,[20 27]); title('exlpore/exploit noRwd'); hold on;
    mySDF_exploreNoRwd = allTrial_spikeMatx_slide(s.responseInds(intersect(s.nrwd_Inds, find(s.hmmStates==1))),:);
    mySDF_exploitNoRwd = allTrial_spikeMatx_slide(s.responseInds(intersect(s.nrwd_Inds, find(s.hmmStates~=1))),:);
    plotFilled(slideTime, mySDF_exploreNoRwd,[0 0 1]) 
    plotFilled(slideTime, mySDF_exploitNoRwd,[0.7 0 1]) 
    line([0 0], [0 maxFreq], 'color', 'r')    
    legend({'explore-nrwd','','exploit-nrwd',''},'FontSize',5)
    ylim([0 maxFreq]);
    
    
    
    if ismember(animalName, {'ZS059','ZS060','ZS062'})
%         z(18) = subplot(8,7,44); hold on;

%         mySDF_shortPreRwd_spike = allTrial_spikeMatx_slide(s.responseInds(intersect(find(s.lickInds==0), s.rwd_Inds+1)),:);
%         mySDF_longPreRwd_spike = allTrial_spikeMatx_slide(s.responseInds(intersect(find(s.lickInds==1), s.rwd_Inds+1)),:);
%         plotFilled(slideTime, mySDF_shortPreRwd_spike,[0 0 1])
%         plotFilled(slideTime, mySDF_longPreRwd_spike,[0.7 0 1])
%         line([0 0], [0 maxLick], 'color', 'r')
%         title('shortAfterRwd vs longAfterRwd')
%         ylim([0 maxFreq]);
% 
%          subplot(8,7,51); hold on;
% 
%         mySDF_shortPreNoRwd_spike = allTrial_spikeMatx_slide(s.responseInds(intersect(find(s.lickInds==0), s.nrwd_Inds+1)),:);
%         mySDF_longPreNoRwd_spike = allTrial_spikeMatx_slide(s.responseInds(intersect(find(s.lickInds==1), s.nrwd_Inds+1)),:);
%         plotFilled(slideTime, mySDF_shortPreNoRwd_spike,[0 0 1])
%         plotFilled(slideTime, mySDF_longPreNoRwd_spike,[0.7 0 1])
%         line([0 0], [0 maxLick], 'color', 'r')
%         title('shortAfterNoRwd vs longAfterNoRwd')
%         ylim([0 maxFreq]);


        z(20) = subplot(8,7,[43 50]); hold on;
        mySDF_one_spike = allTrial_spikeMatx_slide(s.responseInds(s.lickInds==0),:);
        mySDF_two_spike = allTrial_spikeMatx_slide(s.responseInds(s.lickInds==1),:);
        plotFilled(slideTime, mySDF_one_spike,[0 0 1])
        plotFilled(slideTime, mySDF_two_spike,[0.7 0 1])
        line([0 0], [0 maxLick], 'color', 'r')
        legend('shortLicks', '', 'longLicks', '')
        title('short vs long')
        ylim([0 maxFreq]);
    end
    
    % CSplus vs CSminus
%     z(19) = subplot(8,7,[30 37]); hold on;
%     mySDF_CSplus_spike = allTrial_spikeMatx_slide(s.responseInds,:);
%     mySDF_CSminus_spike = allTrial_spikeMatx_slide(CSminus_Inds,:);
%     plotFilled(slideTime, mySDF_CSplus_spike,[0 0 1])
%     plotFilled(slideTime, mySDF_CSminus_spike,[0.7 0 1])
%     line([0 0], [0 maxLick], 'color', 'r')
%     legend('CSplus', '', 'CSminus', '')
%     title('CSplus vs CSminus')
%     ylim([0 maxFreq]);
    mySDF_CR_spike = [];
    mySDF_FA_spike = [];
    mySDF_miss_spike = [];
    mySDF_hit_spike = [];
    % FA vs CR
    subplot(8,7,[44]); hold on;
    if ~isempty(setdiff(CSplus_Inds,lickInds))
        mySDF_miss_spike = allTrial_spikeMatx_slide(setdiff(CSplus_Inds,lickInds),:);
        mySDF_CR_spike = allTrial_spikeMatx_slide(setdiff(CSminus_Inds,lickInds),:);
        if size(mySDF_CR_spike,1)>1 
            plotFilled(slideTime, mySDF_CR_spike,[0 0 1])
        else
            if size(mySDF_CR_spike,1)>0
                plot(slideTime, mySDF_CR_spike, 'color', [0 0 1])
            end
        end
            
        if size(setdiff(CSplus_Inds,lickInds),2)>1 
            plotFilled(slideTime, mySDF_miss_spike,[0.7 0 1])
            line([0 0], [0 maxLick], 'color', 'r')
            legend('CR', '', 'miss', '')
        else
            if size(setdiff(CSplus_Inds,lickInds))>0
               plot(slideTime, mySDF_miss_spike, 'color', [0.7 0 1])
            end
            line([0 0], [0 maxLick], 'color', 'r')
            legend('CR', '', 'miss')            
        end
        title('CR vs miss')
        ylim([0 maxFreq]); 
    end
    
    %% hit vs FA
    subplot(8,7,[51]); hold on;
    if ~isempty(setdiff(lickInds, CSplus_Inds))
        mySDF_FA_spike = allTrial_spikeMatx_slide(setdiff(lickInds, CSplus_Inds),:);
        mySDF_hit_spike = allTrial_spikeMatx_slide(s.responseInds,:);
        plotFilled(slideTime, mySDF_hit_spike,[0.7 0 1])
        if size(mySDF_FA_spike,1)>1 
            plotFilled(slideTime, mySDF_FA_spike,[0 0 1])
            
        else
            if size(mySDF_FA_spike,1)>0
                plot(slideTime, mySDF_FA_spike, 'color', [0 0 1])
            end
        end
        
        line([0 0], [0 maxLick], 'color', 'r')
        legend('hit', '', 'FA')
        
        title('FA vs hit')
        ylim([0 maxFreq]); 
    end
    %% hit vs miss
    subplot(8,7,[37]); hold on;
    if ~isempty(setdiff(CSplus_Inds, lickInds))
        plotFilled(slideTime, mySDF_hit_spike,[0.7 0 1])
        if size(mySDF_miss_spike,1)>1 
            plotFilled(slideTime, mySDF_miss_spike,[0 0 1])
            
        else
            if size(mySDF_FA_spike,1)>0
                plot(slideTime, mySDF_miss_spike, 'color', [0 0 1])
            end
        end
        
        line([0 0], [0 maxLick], 'color', 'r')
        legend('hit', '', 'miss')
        
        title('hit vs miss')
        ylim([0 maxFreq]); 
    end    
  
    %% CR vs FA
    subplot(8,7,30);
    if ~isempty(setdiff(lickInds, CSplus_Inds)) % if there's false alarm
        if size(mySDF_CR_spike,1)>1 
            plotFilled(slideTime, mySDF_CR_spike,[0 0 1])
        else
            if size(mySDF_CR_spike,1)>0
                plot(slideTime, mySDF_CR_spike, 'color', [0 0 1])
            end
        end

        
     
        if size(mySDF_FA_spike,1)>1 
            plotFilled(slideTime, mySDF_FA_spike,[0.7 0 1])
            line([0 0], [0 maxLick], 'color', 'r')
        else
            if size(mySDF_FA_spike,1)>0
               plot(slideTime, mySDF_FA_spike, 'color', [0.7 0 1])
            end
            line([0 0], [0 maxLick], 'color', 'r')          
        end

        
        if size(mySDF_CR_spike,1)>1
            legend('CR', '', 'FA')  
        else
            if size(mySDF_CR_spike,1)==1
                legend('CR', 'FA')
            else
                legend('FA')
            end
        end
        title('CR vs FA')
        ylim([0 maxFreq]); 
    end  
    
%%    
%     mySDF_rwd_lick = allTrial_lickMatx_slide(s.responseInds(s.rwd_Inds),:);
%     mySDF_noRwd_lick = allTrial_lickMatx_slide(s.responseInds(s.nrwd_Inds),:);
%     plotFilled(slideTime, mySDF_rwd_lick,[0 0 1])
%     plotFilled(slideTime, mySDF_noRwd_lick,[0.7 0 1])
%     line([0 0], [0 maxLick], 'color', 'r')
%     legend('rwd', '', 'noRwd', '')
%     title('lick: rwd vs noRwd')
%     ylim([0 maxLick]);
    
    
    %ISI histogram
    %reward and choice behavior
%       slopeSpace = subplot(8,7,[1 8 15]); t(16) = title('Behavior'); hol  d on
%     set(slopeSpace,'YDir','reverse')
%     plot(s.choiceSlope,1:length(s.choiceSlope),'b','linewidth',1.5)
%     plot(s.rwdSlope,1:length(s.rwdSlope),'k','linewidth',1.5)
%     xlim([0 90])
%     ylim([0 sum(s.CSplus_Inds)])
%     xlabel('<-- Left Choice (Slope) Right Choice -->')
%     ylabel('Trials')
%     plot(repmat([-5000 10000],length(blockSwitch_CSminCorrected),1)', [blockSwitch_CSminCorrected; blockSwitch_CSminCorrected],'r')
%     legend({'C','R'},'FontSize',5,'location','best')
%     
%     %rwd hist and firing rate across session
%     subplot(8,7,[43 44; 50 51])
%     yyaxis left
%     plot(sessionSpikeSDF,'Color', [0.5 0.5 0.5], 'LineWidth',2); hold on;
%     plot(sessionSpikeSDFsmooth,'-','Color', [0.5 0 0.8], 'LineWidth',2); hold on;
%     ylabel('Spikes/s')
%     xlim([60000 length(sessionTime)])
%     yyaxis right
%     plot(sessionRwdsSmooth, 'b', 'LineWidth', 2);
%     xlabel('Time')
%     set(gca,'Xtick',[]);
%     xlim([60000 length(sessionTime)])
%     rMin = ylim;
%     rMag = rMin(2)*0.9;
%     for j = 1:length(s.behSessionData)
%         if s.behSessionData(j).rewardL == 1 || s.behSessionData(j).rewardR == 1
%             xTemp = s.behSessionData(j).rewardTime - s.behSessionData(1).CSon;
%             plot([xTemp xTemp],[rMag rMin(2)],'b')
%         end
%     end
%     
%     %determine minimum and maximum y values for SDFs
%     minY = min([mySDF_rwdHxI(trialBeg:trialEnd) mySDF_rwdHxII(trialBeg:trialEnd) mySDF_rwdHxIII(trialBeg:trialEnd)]);
%     maxY = max([mySDF_rwdHxI(trialBeg:trialEnd) mySDF_rwdHxII(trialBeg:trialEnd) mySDF_rwdHxIII(trialBeg:trialEnd)]);
%     if minY - 0.07*minY > 0
%         SDFyLimMin = minY - 0.07*minY;
%     else
%         SDFyLimMin = 0;
%     end
%     SDFyLimMax = maxY + 0.07*maxY;
%     
%     linkaxes([r z],'x')
%     xlim(r(1),[min(time)+500 max(time)-500])
%     ylim(r(5),[0 yLimMax])
%     linkaxes(r(5:14),'y')
%     ylim(z(5),[SDFyLimMin SDFyLimMax])
%     linkaxes(z(1:9),'y')
%     legend('fr', 'fr smooth', 'rwd hist', 'location','best')
%     
%     
    
    axes( 'Position', [0, 0.95, 1, 0.05] ) ; % set axes for the 'text' call below
    set( gca, 'Color', 'None', 'XColor', 'None', 'YColor', 'None' ) ;
    if downsampFlag == true
        text( 0.5, 0, [sessionName ': ' spikeFields{clust(i)} ' Aligned to cue' '. Downsampled By ' num2str(downsampFactor)], 'FontSize', 14', 'FontWeight', 'Bold', ...
          'HorizontalAlignment', 'Center', 'VerticalAlignment', 'Bottom','Interpreter','none') ;
    else
        text( 0.5, 0, [sessionName ': ' spikeFields{clust(i)} ' Aligned to cue'], 'FontSize', 14', 'FontWeight', 'Bold', ...
          'HorizontalAlignment', 'Center', 'VerticalAlignment', 'Bottom','Interpreter','none') ;
    end
    
    
    if p.Results.saveFigFlag == 1 
        saveFigurePDF(rasters,[savepath sep sessionName '_' spikeFields{clust(i)} 'cueAlinged'])
    end
end