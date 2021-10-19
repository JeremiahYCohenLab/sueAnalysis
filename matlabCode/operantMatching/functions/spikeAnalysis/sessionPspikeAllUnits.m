function sessionPspikeAllUnits(session, varargin)
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
pd = parseSessionString_df(session, root, sep);
sortedPath = [pd.nLynxFolder 'session' sep];
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
% all suspect units: 
suspectID = contains(sessionList,session);
suspectUnits = unitList(suspectID);
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
        if respLat <= 15000 && met.isiV <= 0.001 && met.distance<0.35
            unitsQual(a) = 1;
        end
    end 
end
% all units
allClust = find(contains(spikeFields,'TT'));
% suspected ones 
suspectClust = zeros(size(suspectUnits));
for i = 1:length(suspectUnits)
    suspectClust(i) = find(contains(spikeFields,suspectUnits{i}));
end
% potential other units
otherClust = setdiff(allClust, suspectClust);
% IDed units
sessionUnits = sessionUnits(unitsQual>0);
clust = zeros(size(sessionUnits));
for i = 1:length(clust)
    clust(i) = find(contains(spikeFields,sessionUnits{i}));
end
allClust = clust;
for i = 1:length(otherClust)
    spikeTimes = [load(strcat(sortedPath, spikeFields{otherClust(i)}, '.txt'))]';
    if ~isempty(spikeTimes)
        isi = 100*sum(diff(spikeTimes)<2000)/(length(spikeTimes)-1);
        meanRate = 1000000*length(spikeTimes)/(max(spikeTimes)-min(spikeTimes));
        if isi <0.1 && meanRate > 4
            allClust = [allClust; otherClust(i)];
        end
    end
end

unitNames = spikeFields(allClust);
%% put spikes by trials into cells, aligned by cue
allTrial_spike_choice = {};
os = behAnalysisNoPlot_opMD(session);
for j = 1:length(os.responseInds)
    trialDurDiff(j) = (sessionData(os.responseInds(j)).trialEnd - sessionData(os.responseInds(j)).CSon)- p.Results.tf*1000;
end
trialDurDiff(end) = 0; 
for k = 1:length(os.responseInds)
    for i = 1:length(allClust)
        if os.responseInds(k) == 1
            prevTrial_spike = [];
        else
            prevTrial_spikeInd = [sessionData(os.responseInds(k)-1).(spikeFields{allClust(i)})] > (sessionData(os.responseInds(k)).CSon-p.Results.tb*1000);
            prevTrial_spike = sessionData(os.responseInds(k)-1).(spikeFields{allClust(i)})(prevTrial_spikeInd) - sessionData(os.responseInds(k)).CSon;
        end
        
        currTrial_spikeInd = sessionData(os.responseInds(k)).(spikeFields{allClust(i)}) < sessionData(os.responseInds(k)).CSon+p.Results.tf*1000 ... 
            & sessionData(os.responseInds(k)).(spikeFields{allClust(i)}) > sessionData(os.responseInds(k)).CSon-p.Results.tb*1000;
        currTrial_spike = sessionData(os.responseInds(k)).(spikeFields{allClust(i)})(currTrial_spikeInd) - sessionData(os.responseInds(k)).CSon;
        
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
midPoints = (0.5*p.Results.binSize + 1):p.Results.stepSize:(length(time)-0.5*p.Results.binSize);
slideTime = midPoints - p.Results.tb*1000;
% sometimes no licks/spikes are considered 1x0 and sometimes they are []
% plotSpikeRaster does not place nicely with [] so this converts all empty indices to 1x0
allTrial_spike_choice(cellfun(@isempty,allTrial_spike_choice)) = {zeros(1,0)}; 

%% loop through neurons


Corr = figure; 
screen = get(0,'Screensize');
screen(4) = screen(4) - 100;
set(Corr, 'Position', screen)
suptitle([session 'allNEunitsCorr'])

allSpikes = [];
allSpikesPre = [];

for i = 1:length(allClust)
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
lag = round(500/p.Results.binSize);

for i = 1:length(allClust) %all PCs
  for j = 1:length(allClust)
      if j >= i
          x = allSpikes(:,:,i);
          y = allSpikes(:,:,j);
          varx = sum(x.*x,'all');
          vary = sum(y.*y,'all');
          corrCo = [];
          corrCoShuffle = [];
          corrCoPreShuffle = [];
          for t = 1:size(x,1)
              mTemp = xcorr(x(t,:),y(t,:),lag);
              mTemp(0.5*(length(mTemp)+1)) = NaN;
              corrCo = [corrCo; mTemp];
          end    
          corrCo = sum(corrCo)/sqrt(varx*vary);
          
          % shuffled data
          x = [x(round(0.5*length(os.responseInds)):end,:); x(1:round(0.5*length(os.responseInds))-1,:)]; %half shuffling
          for t = 1:size(x,1)
              mTemp = xcorr(x(t,:),y(t,:),lag);
              mTemp(0.5*(length(mTemp)+1)) = NaN;
              corrCoShuffle = [corrCoShuffle; mTemp];
          end    
          corrCoShuffle = sum(corrCoShuffle)/sqrt(varx*vary);
          % ITI
          x = allSpikesPre(:,:,i);
          y = allSpikesPre(:,:,j);
          varx = sum(x.*x,'all');
          vary = sum(y.*y,'all');
          corrCoPre = [];
          for t = 1:size(x,1)
              mTemp = xcorr(x(t,:),y(t,:),lag);
              mTemp(0.5*(length(mTemp)+1)) = NaN;
              corrCoPre = [corrCoPre; mTemp];
          end   
          corrCoPre = sum(corrCoPre)/sqrt(varx*vary);
          % shuffled data
          x = [x(round(0.5*length(os.responseInds)):end,:); x(1:round(0.5*length(os.responseInds))-1,:)]; %half shuffling
          for t = 1:size(x,1)
              mTemp = xcorr(x(t,:),y(t,:),lag);
              mTemp(0.5*(length(mTemp)+1)) = NaN;
              corrCoPreShuffle = [corrCoPreShuffle; mTemp];
          end    
          corrCoPreShuffle = sum(corrCoPreShuffle)/sqrt(varx*vary);          
          
          subplot(length(allClust),length(allClust),(i-1)*length(allClust)+j); hold on;
          plot((-lag:lag)*p.Results.binSize,corrCo,'linewidth',2,'color','m');
          plot((-lag:lag)*p.Results.binSize,corrCoShuffle,'linewidth',1,'color',[1 0.5 1]);
          
          plot((-lag:lag)*p.Results.binSize,corrCoPre,'linewidth',2,'color','c');
          plot((-lag:lag)*p.Results.binSize,corrCoPreShuffle,'linewidth',1,'color',[0.5 1 1]);
          if i==1 && j == 1
              legend({'In Trial','In Trial Shuffle', 'Inter Trial', 'Trial Shuffle'});
          end
          
          if j == i
              ylabel(unitNames{i},'Interpreter', 'none')
              xlabel(unitNames{j},'Interpreter', 'none')
          end
          
          
          if i <= length(clust) && j <= length(clust)
              set(gca,'Color',[1, 0.9 0.9]);
          else
              if i > length(clust) && j > length(clust)
                  set(gca,'Color',[0.9 1 0.9]);
              end
          end
      end
  end
end
   
