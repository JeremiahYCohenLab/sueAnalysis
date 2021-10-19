function sessionPspike(session, col, varargin)
%task and model parameters
p = inputParser;
% default parameters if none given
p.addParameter('plotFlag', 1);
p.addParameter('maxTrial', 1000);
p.addParameter('modelNameOld', '5params');
p.addParameter('modelName', '5params');
p.addParameter('binSize', 10)% in ms  
p.addParameter('stepSize', 10)  
p.addParameter('tb', 2)% in s
p.addParameter('tf', 4)% in s
p.addParameter('saveFigFlag', 0);
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
suptitle([session 'allNEunitsCorr'])
colors = cool(size(combineMat,2));
unitColors = cool(length(clust));

allSpikes = [];
allSpikesPre = [];

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
    allTrial_spikeMatx_slide = allTrial_spikeMatx_slide;
    allSpikes = cat(3, allSpikes, zscore(allTrial_spikeMatx_slide(:,find(slideTime>0,1):end),0,'all'));
    allSpikesPre = cat(3, allSpikesPre, zscore(allTrial_spikeMatx_slide(:,1:find(slideTime<0,1,'last')),0,'all'));
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
lag = round(1000/p.Results.binSize);
for i = 1:length(clust) %all PCs
  for j = 1:length(clust)
      if j>=i
          x = allSpikes(:,:,i);
          y = allSpikes(:,:,j);
          corrCo = [];
          for t = 1:size(x,1)
              mTemp = xcorr(x(t,:),y(t,:),lag,'normalized');
              mTemp(0.5*(length(mTemp)+1)) = NaN;
              corrCo = [corrCo; mTemp];
          end     
          x = allSpikesPre(:,:,i);
          y = allSpikesPre(:,:,j);
          corrCoPre = [];
          for t = 1:size(x,1)
              mTemp = xcorr(x(t,:),y(t,:),lag,'normalized');
              mTemp(0.5*(length(mTemp)+1)) = NaN;
              corrCoPre = [corrCoPre; mTemp];
          end               
          subplot(length(clust),length(clust),(i-1)*length(clust)+j); hold on;
          plot(-lag:lag,mean(corrCo),'linewidth',2,'color','m');
          plot(-lag:lag,mean(corrCoPre),'linewidth',2,'color','c');
          if i==1 && j == 1
              legend({'In Trial', 'Inter Trial'});
          end
          
          if i == j
              xlabel(sessionUnits{j},'Interpreter', 'none')
              ylabel(sessionUnits{i},'Interpreter', 'none')
          end
      end
  end
end
   
