function sessionTrialTypePCsPE(session, col, varargin)
%task and model parameters
p = inputParser;
% default parameters if none given
p.addParameter('plotFlag', 1);
p.addParameter('maxTrial', 1000);
p.addParameter('modelNameOld', '5params');
p.addParameter('modelName', '5params');
p.addParameter('binSize', 200)% in ms  
p.addParameter('stepSize', 50)  
p.addParameter('tb', 1)% in s
p.addParameter('tf', 2)% in s
p.addParameter('saveFigFlag', 1);
p.parse(varargin{:});

paramNames = getParamNames_dF(p.Results.modelName,1);
maxTrial = p.Results.maxTrial;
% basic info
[root, sep] = currComputer();
[animalName, date] = strtok(session, 'd'); 
animalName = animalName(2:end); 
date = date(1:9);
sessionFolder = ['m' animalName date];  
% paths
if isstrprop(session(end), 'alpha')
    neuralynxDataPath = [root animalName sep sessionFolder sep 'sorted' sep 'session ' session(end) sep session '_sessionData_nL.mat'];
    sortedFolderLocation = [root animalName sep session(1:end-1) sep 'sorted' sep 'session ' sessionName(end) sep];
    savepath = [root animalName sep session(1:end-1) sep  'figures' sep 'session ' sessionName(end) sep];
    unitPath = [root animalName sep session(1:end-1) sep 'neuralynx' sep 'session' sep];
else
    neuralynxDataPath = [root animalName sep sessionFolder sep 'sorted' sep 'session' sep session '_sessionData_nL.mat'];
    sortedFolderLocation = [root animalName sep session sep 'sorted' sep 'session' sep];
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
if length(sessionUnits)<2
    fprintf(['no multiple units in ' session '\n'])
    return
end
sessionSubFolders = subFolders(sessionID);
unitsQual = zeros(1,length(sessionUnits));
for a = 1:length(sessionUnits)
    if exist([unitPath session '_' sessionUnits{a} '_met.mat'],'file')
        load([unitPath session '_' sessionUnits{a} '_met.mat']);
    else
        met = getClusterMetric(session, sessionUnits{a}, 0, 1);
    end
     respInds = find(met.spikeProp>=0.7);
    if ~isempty(respInds)
        respLat = nanmin(met.spikeLat(respInds));
        if respLat <= 15000 && met.isiV <= 0.001 && met.distance<0.3
            unitsQual(a) = 1;
        end
    end 
end

sessionUnits = sessionUnits(unitsQual>0);
if length(sessionUnits)<3
    fprintf(['no multiple units in ' session '\n'])
    return
end
clust = zeros(size(sessionUnits));
for i = 1:length(clust)
    clust(i) = find(contains(spikeFields,sessionUnits{i}));
end
%% load model fitting results and calculate DVs
sampFile = [animalName col '_', p.Results.modelNameOld];
path = [root animalName sep animalName 'sorted' sep 'stan' sep 'bernoulli' sep p.Results.modelNameOld sep col sep];
load([path sampFile '.mat'], 'dayList');
samples = load([path sampFile '.mat'], sampFile);
samples = samples.(sampFile);

%% behavior preparation 
% parse behavior
os = behAnalysisNoPlot_opMD(session);
choice = os.allChoices';
choice(choice<0) = 0;
outcome = abs(os.allRewards);
choice = choice(1:min(length(choice), maxTrial));
outcome = outcome(1:length(choice));
responseInds = os.responseInds(1:min(length(choice), maxTrial)); 
ITI = os.timeBtwn(1:length(choice)); 
preRwd = [NaN abs(os.allRewards(1:end-1))]';
%% behavior
% switch
svsTemp = find(choice(2:end) ~= choice(1:end-1)) + 1;
svs = zeros(1,length(responseInds));
svs(svsTemp) = 1;
% right choice
rightSide = zeros(size(svs));
rightSide(os.lickR_Inds) = 1;
% model 
id = find(strcmp(dayList,session),1);
if isempty(id)
    fprintf([session ' not ' col ' behavior \n']);
    return
end
%generate best estimates of parameters
paramEsts = [];
allSamples = [];
edges = cell(1,length(paramNames));
for a = 1:length(paramNames)
    tmp = samples.(paramNames{a})(:,id);
    allSamples = [allSamples tmp];
    edges{a} = linspace(min(tmp), max(tmp),50);
end
n = histcnd(allSamples,edges); %bin samples by multiple dimensions
[~, inds] = myMaxAll(n); %find the bin with max num in bin
for a = 1:length(paramNames) %use median in bin as best estimate
    tmp = allSamples(:,a);
    edgeTmp = edges{a};
    if inds(a) < 50
        paramEsts(a) = median(tmp(tmp >= edgeTmp(inds(a)) & tmp < edgeTmp(inds(a)+1)));
    else
        paramEsts(a) = edgeTmp(inds(a));
    end
end

% decide if input includes time forget
if contains(p.Results.modelName,'tF')
    input = 'choice, outcome, ITI)';
else
    input = 'choice, outcome)';
end
eval(['[LL,probC,Q,pe] = qLearningModel_' p.Results.modelName '(paramEsts,' input ';'])
% diff value
Qdiff = abs(Q(1:end-1,2)-Q(1:end-1,1));
% total value
Qsum = sum(Q(1:end-1,:),2);
% prepe
prePe = [NaN; pe(1:end-1)];

%time in session
timeTemp = [sessionData(responseInds).CSon] - sessionData(responseInds(1)).CSon;

%% 
biasSide = zeros(size(responseInds));
biasInd = ismember('bias', paramNames);
if paramEsts(biasInd)>0
    biasSide(os.lickR_Inds)=1;
else
    biasSide(os.lickL_Inds)=1;
end


%% put spikes by trials into cells, aligned by cue
allTrial_spike_choice = {};
for k = 1:length(os.responseInds)
    for i = 1:length(clust)
        if os.responseInds(k) == 1
            prevTrial_spike = [];
        else
            prevTrial_spikeInd = [sessionData(os.responseInds(k)-1).(spikeFields{clust(i)})] > (sessionData(os.responseInds(k)).CSon-p.Results.tb*1000);
            prevTrial_spike = sessionData(os.responseInds(k)-1).(spikeFields{clust(i)})(prevTrial_spikeInd) - sessionData(os.responseInds(k)).CSon;
        end
        
        currTrial_spikeInd = sessionData(os.responseInds(k)).(spikeFields{clust(i)}) < sessionData(os.responseInds(k)).CSon+p.Results.tf*1000 ... 
            & sessionData(os.responseInds(k)).(spikeFields{clust(i)}) > sessionData(os.responseInds(k)).CSon-p.Results.tb*1000;
        currTrial_spike = sessionData(os.responseInds(k)).(spikeFields{clust(i)})(currTrial_spikeInd) - sessionData(os.responseInds(k)).CSon;
        
        allTrial_spike_choice{i,k} = [prevTrial_spike currTrial_spike];

    end
    
    if ~isnan(sessionData(os.responseInds(k)).rewardL)
        currTrial_lickInd = [sessionData(os.responseInds(k)).licksL] < (sessionData(os.responseInds(k)).CSon + p.Results.tf*1000);
        currTrial_lick = sessionData(os.responseInds(k)).licksL(currTrial_lickInd) - sessionData(os.responseInds(k)).CSon;
    elseif ~isnan(sessionData(os.responseInds(k)).rewardR)
        currTrial_lickInd = [sessionData(os.responseInds(k)).licksR] < (sessionData(os.responseInds(k)).CSon + p.Results.tf*1000);
        currTrial_lick = sessionData(os.responseInds(k)).licksR(currTrial_lickInd) - sessionData(os.responseInds(k)).CSon;  
    else
        currTrial_lick = 0;
    end
    allTrial_lick{k} = [currTrial_lick];
end

% sometimes no licks/spikes are considered 1x0 and sometimes they are []
% plotSpikeRaster does not place nicely with [] so this converts all empty indices to 1x0
allTrial_spike_choice(cellfun(@isempty,allTrial_spike_choice)) = {zeros(1,0)}; 
%% initialize lick matrices
allTrial_lickMatx = NaN(length(os.responseInds),length(time)); 
for j = 1:length(os.responseInds)
    trialDurDiff(j) = (sessionData(os.responseInds(j)).trialEnd - sessionData(os.responseInds(j)).CSon)- p.Results.tf*1000;
end
trialDurDiff(end) = 0; 
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
allTrial_lickMatx_slide = zeros(length(os.responseInds), length(midPoints));
for w = 1:length(midPoints)
    allTrial_lickMatx_slide(:,w) = ...
        sum(allTrial_lickMatx(:,midPoints(w)-0.5*p.Results.binSize:midPoints(w)+0.5*p.Results.binSize-1),2)*1000/p.Results.binSize;   
end
%% zscore all varialbles
% pe = outcome - probC';
hmm = double(os.hmmStates==1)';
combineMat = [outcome', rightSide',os.lickInds];
regressors = {'outcome', 'rightSide','lick'};
for k = 1:size(combineMat,2)
    combineMat(~isnan(combineMat(:,k)),k) = zscore(combineMat(~isnan(combineMat(:,k)),k),0,1);
end

%% loop through neurons


GLM = figure; 
screen = get(0,'Screensize');
screen(4) = screen(4) - 100;
set(GLM, 'Position', screen)
suptitle([session 'allNEunits---pe'])
colors = cool(size(combineMat,2));
unitColors = cool(length(clust));

allSpikes = [];
allCOSpikes = [];
for i = 1:length(clust)
    % Initialize matrices for SDF
    allTrial_spikeMatx_choice = zeros(length(os.responseInds),length(time));         
   for j = 1:length(allTrial_spike_choice)
        tempSpike = allTrial_spike_choice{i,j};
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
    
    % zscore spike matrix
    allTrial_spikeMatx_slide = zscore(allTrial_spikeMatx_slide,0,'all');
    LR = mean(allTrial_spikeMatx_slide(intersect(os.rwd_Inds,os.lickL_Inds),:));
    LN = mean(allTrial_spikeMatx_slide(intersect(os.nrwd_Inds,os.lickL_Inds),:));
    RR = mean(allTrial_spikeMatx_slide(intersect(os.rwd_Inds,os.lickR_Inds),:));
    RN = mean(allTrial_spikeMatx_slide(intersect(os.nrwd_Inds,os.lickR_Inds),:));
    
    COspikesTemp = [LR, LN, RR, RN]';
    allCOSpikes = [allCOSpikes COspikesTemp];
    allSpikes = [allSpikes reshape(allTrial_spikeMatx_slide', [numel(allTrial_spikeMatx_slide),1])];
    %% zscore and fit model with licks
% 
%     allTrial_lickMatx_slideZS = zscore(allTrial_lickMatx_slide,0,1);
%     coeff = zeros(length(midPoints),size(combineMat,2)+2,3);
%     rsq = zeros(size(midPoints));
%     for k = 1:length(midPoints)
%         lm = fitglm([combineMat,allTrial_lickMatx_slideZS(:,k)],allTrial_spikeMatx_slide(:,k),'linear','Distribution','poisson');
%         for j = 1:size(combineMat,2)+2
%         coeff(k,j,1) = lm.Coefficients.Estimate(j);
%         ci = coefCI(lm);
%         coeff(k,j,2:3) = ci(j,:);
%         end
%         rsq(k)=lm.Rsquared.Adjusted;
%     end
    %% zscore and fit model without licks
end
    % r
    % allSpikes = allSpikes - mean(allSpikes,2);
    [PCs,scores,~,tsquared,explained, mu] = pca(allCOSpikes);
    allScores = (allSpikes-mu)*PCs;
for i = 1:length(clust) %all PCs
    coeff = zeros(length(midPoints),size(combineMat,2)+1,3);
    rsq = zeros(size(midPoints));
    sigs = zeros(length(midPoints),size(combineMat,2)+1);
    allTrial_spikeMatx_slide = reshape(allScores(:,i), [length(midPoints), length(os.responseInds)])';
    
    
    for k = 1:length(midPoints)
        lm = fitlm(combineMat,allTrial_spikeMatx_slide(:,k));
        for j = 1:size(combineMat,2)+1
        coeff(k,j,1) = lm.Coefficients.Estimate(j);
        ci = coefCI(lm);
        coeff(k,j,2:3) = ci(j,:);
        sigs(k,j) = double(lm.Coefficients.pValue(j)<0.05);
        end
        rsq(k)=lm.Rsquared.Adjusted;
    end
    %% bin no.spike by pe
    numBins = 10; % bins in pe
    width = 1000; %in ms
    startTime = 400; % in ms
    endTime = startTime + width; 
%     % decide start of the time windows
%     solenoidTimeEnd = max([sessionData(responseInds).CSon] + 3100 - [sessionData(responseInds).respondTime]); % 3600 ms is 500ms after sol coming back after 3100ms
%     solenoidTimeStart = min([sessionData(responseInds).CSon] + 3100 - [sessionData(responseInds).respondTime]);
%     
%     peInd = find(contains(regressors,'pe'));
%     if ~isempty(peInd) % if pe is regressor
%         sigIndEarly = find(sum(sigs(1:end-1,peInd+1),2)>0 & sum(sigs(2:end,peInd+1),2)>0 & slideTime(1:end-1)'<solenoidTimeStart & slideTime(1:end-1)'>os.rwdDelay); 
%         if ~isempty(sigIndEarly) % if any significance
%             peEffect = sum(squeeze(abs(coeff(sigIndEarly,peInd+1,1))),2);
%             [~,maxInd] = max(peEffect);
%             startTime = max(midPoints(sigIndEarly(maxInd)) - 0.5*width, 1000*p.Results.tb + os.rwdDelay);
%         else 
%             startTime = 1000*p.Results.tb + os.rwdDelay;
%         end
%         endTime  = min(startTime+width, solenoidTimeStart+1000*p.Results.tb);
% 
%         sigIndLate = find(sum(sigs(:,peInd+1),2)>0 & slideTime' > solenoidTimeEnd & slideTime'<1000*p.Results.tf);
%         if ~isempty(sigIndLate)
%             peLateEffect = sum(squeeze(abs(coeff(sigIndLate,peInd+1,1))),2);
%             [~,maxInd] = max(peLateEffect);
%             startTimeLate = max(midPoints(sigIndLate(maxInd)) - 0.5*width, solenoidTimeEnd + 1000*p.Results.tb);
%         else 
%             startTimeLate = solenoidTimeEnd + 1000*p.Results.tb;
%         end
%         endTimeLate  = min(startTimeLate+width, 1000*(p.Results.tb+p.Results.tf));
%     else
%         startTime = 1000*p.Results.tb + os.rwdDelay;
%         endTime  = min(startTime+width, solenoidTimeStart+1000*p.Results.tb);
%         startTimeLate = solenoidTimeEnd + 1000*p.Results.tb;
%         endTimeLate  = min(startTimeLate+width, 1000*(p.Results.tb+p.Results.tf));
%     end
    spikeCounts = zscore(mean(allTrial_spikeMatx_slide(:,find(slideTime > startTime,1):find(slideTime<endTime,1,'last')),2,'omitnan'));
    %% caculate on both sides
    edges = binEqualSize(pe, numBins);
    spikeMeans = zeros(numBins,1);
    spikeSems = zeros(numBins,1);
    peMeans = zeros(numBins,1);
%     spikeMeansLate = zeros(numBins,1);
%     spikeSemsLate = zeros(numBins,1);
    for k = 1:numBins
        if k < numBins
            spikeNumsTemp = spikeCounts(pe >= edges(k) & pe < edges(k+1));
%             spikeNumsTempLate = spikeCountsLate(pe >= edges(k) & pe < edges(k+1));
            peMeans(k) = mean(pe(pe >= edges(k) & pe < edges(k+1)));
        else
            spikeNumsTemp = spikeCounts(pe >= edges(k) & pe <= edges(k+1));
%             spikeNumsTempLate = spikeCountsLate(pe >= edges(k) & pe <= edges(k+1));
            peMeans(k) = mean(pe(pe >= edges(k) & pe <= edges(k+1)));
        end
        spikeMeans(k) = mean(spikeNumsTemp);
%         spikeMeansLate(k) = mean(spikeNumsTempLate);
        spikeSems(k) = sem(spikeNumsTemp);
%         spikeSemsLate(k) = sem(spikeNumsTempLate);
    end
    
    %% separate left
    peL = pe(os.lickL_Inds);
    spikeCountsL = spikeCounts(os.lickL_Inds);
%     spikeCountsLLate = spikeCountsLate(os.lickL_Inds);
    edgesL = binEqualSize(peL, numBins);
    spikeMeansL = zeros(numBins,1);
    spikeSemsL = zeros(numBins,1);
%     spikeMeansLLate = zeros(numBins,1);
%     spikeSemsLLate = zeros(numBins,1);
    peMeansL = zeros(numBins,1);
    for k = 1:numBins
        if k < numBins
            spikeNumsTemp = spikeCountsL(peL >= edgesL(k) & peL < edgesL(k+1));
%             spikeNumsTempLate = spikeCountsLLate(peL >= edgesL(k) & peL < edgesL(k+1));
            peMeansL(k) = mean(peL(peL >= edgesL(k) & peL < edgesL(k+1)));
        else
            spikeNumsTemp = spikeCountsL(peL >= edgesL(k) & peL <= edgesL(k+1));
%             spikeNumsTempLate = spikeCountsLLate(peL >= edgesL(k) & peL <= edgesL(k+1));
            peMeansL(k) = mean(peL(peL >= edgesL(k) & peL <= edgesL(k+1)));
        end
        spikeMeansL(k) = mean(spikeNumsTemp);
%         spikeMeansLLate(k) = mean(spikeNumsTempLate);
        spikeSemsL(k) = sem(spikeNumsTemp);
%         spikeSemsLLate(k) = sem(spikeNumsTempLate);
    end
    
    %% separate right
    peR = pe(os.lickR_Inds);
    spikeCountsR = spikeCounts(os.lickR_Inds);
%     spikeCountsRLate = spikeCountsLate(os.lickR_Inds);
    edgesR = binEqualSize(peR, numBins);
    spikeMeansR = zeros(numBins,1);
    spikeSemsR = zeros(numBins,1);
    peMeansR = zeros(numBins,1);
%     spikeMeansRLate = zeros(numBins,1);
%     spikeSemsRLate = zeros(numBins,1);
    for k = 1:numBins
        if k < numBins
            spikeNumsTemp = spikeCountsR(peR >= edgesR(k) & peR < edgesR(k+1));
%             spikeNumsTempLate = spikeCountsRLate(peR >= edgesR(k) & peR < edgesR(k+1));
            peMeansR(k) = mean(peR(peR >= edgesR(k) & peR < edgesR(k+1)));
        else
            spikeNumsTemp = spikeCountsR(peR >= edgesR(k) & peR <= edgesR(k+1));
%             spikeNumsTempLate = spikeCountsRLate(peR >= edgesR(k) & peR <= edgesR(k+1));
            peMeansR(k) = mean(peR(peR >= edgesR(k) & peR <= edgesR(k+1)));
        end
        spikeMeansR(k) = mean(spikeNumsTemp);
%         spikeMeansRLate(k) = mean(spikeNumsTempLate);
        spikeSemsR(k) = sem(spikeNumsTemp);
%         spikeSemsRLate(k) = sem(spikeNumsTempLate);        
    end
    
    %% plots 
    
    if ismember('bias', paramNames)
        biasInd = contains(paramNames, 'bias');
        if paramEsts(biasInd) > 0
            colorR = [1 0 1];
            colorL = [0 1 1];
        else
            colorR = [0 1 1];
            colorL = [1 0 1];    
        end
    else
        colorL = [1 0 1];
        colorR = [0 1 1];
    end
        
    subplot(3,length(clust)+1,i+1); hold on;
    x = slideTime;
    yyaxis right
    plot(x, coeff(:,1,1), 'Color', [0.5 0.5 0.5], 'LineStyle','-', 'Marker','none', 'linewidth', 2);
    ylim([1.2*min(coeff(:,:,2),[],'all') 1.2*max(coeff(:,:,3),[],'all')])
    ylabel('intercept')

    yyaxis left 
    for j = 1:size(combineMat,2)
        plot(x, coeff(:,j+1,1), 'Color', colors(j,:), 'LineStyle','-', 'Marker','none', 'linewidth', 2); 
    end

    for j = 1:size(combineMat,2)
       fill([x fliplr(x)], [coeff(:,j+1,2)' fliplr(coeff(:,j+1,3)')], colors(j,:), 'facealpha', 0.25, 'edgecolor', 'none')
    end

%     yyaxis right
%     fill([x fliplr(x)], [coeff(:,1,2)' fliplr(coeff(:,1,3)')], [0.5 0.5 0.5], 'facealpha', 0.25, 'edgecolor', 'none')
%     plt = gca;
%     plt.YAxis(2).Color = [0.2 0.2 0.2];

    yyaxis left
    line(minmax(x), [0 0],'Color', [0.2 0.2 0.2], 'LineStyle','--')
    line([0 0], 1.2*[min(coeff(:,:,2),[],'all'),max(coeff(:,:,3),[],'all')],'Color', [0.2 0.2 0.2], 'LineStyle','--')
    line([os.rwdDelay os.rwdDelay], 1.2*[min(coeff(:,:,2),[],'all'),max(coeff(:,:,3),[],'all')],'Color', [1 0 0], 'LineStyle','--')
    line([startTime endTime], [max(coeff(:,2:end,3),[],'all'),max(coeff(:,2:end,3),[],'all')],'Color', [1 0 0], 'LineWidth',4)
%     line([startTimeLate-1000*p.Results.tb endTimeLate-1000*p.Results.tb],[max(coeff(:,2:end,3),[],'all'),max(coeff(:,2:end,3),[],'all')],'Color', [0 0 1], 'LineWidth',4)
    ylim([1.2*min(coeff(:,2:end,2),[],'all') 1.2*max(coeff(:,2:end,3),[],'all')]);
    xlim(minmax(x))
    ylabel('\beta coefficients')
    xlabel('respond time')
    if i == 1
        legend(regressors, 'Location', 'best')
    end
    title(['PC' num2str(i)], 'interpreter', 'none');

    subplot(3,length(clust)+1,i+1+length(clust)+1); hold on;
    errorbar(peMeans, spikeMeans, spikeSems,'linewidth', 2)
    title('PC change with pe')
    
    subplot(3,length(clust)+1,i+1+2*(length(clust)+1)); hold on;
    errorbar(peMeansL, spikeMeansL, spikeSemsL,'linewidth', 2, 'color', 'c')
    errorbar(peMeansR, spikeMeansR, spikeSemsR,'linewidth', 2, 'color', 'm')
    legend({'L','R'})
    title('PC change with left/right pe')
    
    subplot(3,length(clust)+1,1); hold on;
    errorbar(peMeans, spikeMeans, spikeSems,'linewidth', 2, 'color', unitColors(i,:))   

    subplot(3,length(clust)+1,1+(length(clust)+1)); hold on;
    errorbar(peMeansL, spikeMeansL, spikeSemsL,'linewidth', 2, 'color', unitColors(i,:))
    
    subplot(3,length(clust)+1,1+2*(length(clust)+1)); hold on;
    errorbar(peMeansR, spikeMeansR, spikeSemsR,'linewidth', 2, 'color', unitColors(i,:))
    
end
   
    subplot(3,length(clust)+1,1); hold on;
    title('PC change with pe')  

    subplot(3,length(clust)+1,1+(length(clust)+1)); hold on;
    title(['PC change with left pe' num2str(length(os.lickL_Inds))])
    
    subplot(3,length(clust)+1,1+2*(length(clust)+1)); hold on;
    title(['PC change with right pe' num2str(length(os.lickR_Inds))])
    
    if p.Results.saveFigFlag == 1 
        saveFigurePDF(GLM,[savepath sep session '_' 'allPe'])
    end
    
figure;
subplot(1,3,1);
plot(explained);
title('varExplained');

subplot(1,3,2);
imagesc(PCs);
ylabel('Units');
xlabel('PCs');
title('coeffs')

subplot(1,3,3); hold on;
colors = cool(4);
for i = 1:4
    plot3(scores(1+(i-1)*length(midPoints):i*length(midPoints),1),scores(1+(i-1)*length(midPoints):i*length(midPoints),2),scores(1+(i-1)*length(midPoints):i*length(midPoints),3),'linewidth', 2, 'color',colors(i,:))
end
 rotate3d on
legend({'LR','LN','RR','RN'})