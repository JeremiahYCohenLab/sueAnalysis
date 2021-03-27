    function spikeGLM_dF(session, col, varargin)
%task and model parameters
p = inputParser;
% default parameters if none given
p.addParameter('cellName', ['all']);
p.addParameter('plotFlag', 1);
p.addParameter('maxTrial', 1000);
p.addParameter('modelNameOld', '4params');
p.addParameter('modelName', '4params');
p.addParameter('paramNames', {'aN', 'aP', 'aF', 'beta'});
p.addParameter('binSize', 300)% in ms
p.addParameter('stepSize', 100)
p.addParameter('tb', 1)% in s
p.addParameter('tf', 2)% in s
p.addParameter('saveFigFlag', 1);
p.parse(varargin{:});

paramNames = p.Results.paramNames;
maxTrial = p.Results.maxTrial;
cellName = p.Results.cellName;
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
else
    neuralynxDataPath = [root animalName sep sessionFolder sep 'sorted' sep 'session' sep session '_sessionData_nL.mat'];
    sortedFolderLocation = [root animalName sep session sep 'sorted' sep 'session' sep];
    savepath = [root animalName sep session sep  'figures' sep];   
end
% behavior & neuron
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
%% load model fitting results and calculate DVs
sampFile = [animalName col '_', p.Results.modelNameOld];
path = [root animalName sep animalName 'sorted' sep 'stan' sep 'bernoulli' sep p.Results.modelNameOld sep];
load([path sampFile '.mat'], 'dayList');
samples = load([path sampFile '.mat'], sampFile);
samples = samples.(sampFile);

%% behavior preparation 
% parse behavior
os = behAnalysisNoPlot_opMD(session);
choice = os.allChoices;
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
% model 
id = find(strcmp(dayList,session),1);
%generate best estimates of parameters
paramEsts = [];
allSamples = [];
edges = cell(1,length(paramNames));
for i = 1:length(paramNames)
    tmp = samples.(p.Results.paramNames{i})(:,id);
    allSamples = [allSamples tmp];
    edges{i} = linspace(min(tmp), max(tmp),50);
end
n = histcnd(allSamples,edges); %bin samples by multiple dimensions
[~, inds] = myMaxAll(n); %find the bin with max num in bin
for i = 1:length(paramNames) %use median in bin as best estimate
    tmp = allSamples(:,i);
    edgeTmp = edges{i};
    if inds(i) < 50
        paramEsts(i) = median(tmp(tmp >= edgeTmp(inds(i)) & tmp < edgeTmp(inds(i)+1)));
    else
        paramEsts(i) = edgeTmp(inds(i));
    end
end

% decide if input includes time forget
if contains(p.Results.modelName,'tF')
    input = 'choice, outcome, ITI)';
else
    input = 'choice, outcome)';
end
eval(['[~,~,Q,pe,cQ] = qLearningModel_' p.Results.modelName '(paramEsts,' input ';'])
% diff value
Qdiff = abs(Q(1:end-1,2)-Q(1:end-1,1));
% total value
Qsum = sum(Q(1:end-1,:),2);
% prepe
prePe = [NaN pe(1:end-1)'];

%time in session
timeTemp = [sessionData(responseInds).CSon] - sessionData(responseInds(1)).CSon;


%% create spike and lick cell

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
for k = 1:length(os.responseInds)
    for i = 1:length(clust)
        if os.responseInds(k) == 1
            prevTrial_spike = [];
        else
            prevTrial_spikeInd = [sessionData(os.responseInds(k)-1).(spikeFields{clust(i)})] > (sessionData(os.responseInds(k)).respondTime-p.Results.tb*1000);
            prevTrial_spike = sessionData(os.responseInds(k)-1).(spikeFields{clust(i)})(prevTrial_spikeInd) - sessionData(os.responseInds(k)).respondTime;
        end
        
        currTrial_spikeInd = sessionData(os.responseInds(k)).(spikeFields{clust(i)}) < sessionData(os.responseInds(k)).respondTime+p.Results.tf*1000 ... 
            & sessionData(os.responseInds(k)).(spikeFields{clust(i)}) > sessionData(os.responseInds(k)).respondTime-p.Results.tb*1000;
        currTrial_spike = sessionData(os.responseInds(k)).(spikeFields{clust(i)})(currTrial_spikeInd) - sessionData(os.responseInds(k)).respondTime;
        
        allTrial_spike_choice{i,k} = [prevTrial_spike currTrial_spike];

    end
    
    if ~isnan(sessionData(os.responseInds(k)).rewardL)
        currTrial_lickInd = [sessionData(os.responseInds(k)).licksL] < (sessionData(os.responseInds(k)).respondTime + p.Results.tf*1000);
        currTrial_lick = sessionData(os.responseInds(k)).licksL(currTrial_lickInd) - sessionData(os.responseInds(k)).respondTime;
    elseif ~isnan(sessionData(os.responseInds(k)).rewardR)
        currTrial_lickInd = [sessionData(os.responseInds(k)).licksR] < (sessionData(os.responseInds(k)).respondTime + p.Results.tf*1000);
        currTrial_lick = sessionData(os.responseInds(k)).licksR(currTrial_lickInd) - sessionData(os.responseInds(k)).respondTime;  
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
hmm = double(os.hmmStates==1)';
combineMat = [zscore(Qdiff).*os.timeBtwn', os.timeBtwn', Qdiff];
for k = 1:size(combineMat,2)
    combineMat(~isnan(combineMat(:,k)),k) = zscore(combineMat(~isnan(combineMat(:,k)),k),0,1);
end

%% loop through neurons

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

    coeff = zeros(length(midPoints),size(combineMat,2)+1,3);
    rsq = zeros(size(midPoints));
    for k = 1:length(midPoints)
        lm = fitglm(combineMat,allTrial_spikeMatx_slide(:,k),'linear','Distribution','poisson');
        for j = 1:size(combineMat,2)+1
        coeff(k,j,1) = lm.Coefficients.Estimate(j);
        ci = coefCI(lm);
        coeff(k,j,2:3) = ci(j,:);
        end
        rsq(k)=lm.Rsquared.Adjusted;
    end
    %% bin no.spike by pe
    startTrial = 5;
    numBins = 10;
    width = 500; %in ms
    startTime = os.rwdDelay + p.Results.tb*1000;
    endTime = startTime + width;
    peTemp = pe(startTrial:end);
    [~,edges] = histcounts(peTemp, numBins);
    spikeCounts = nansum(allTrial_spikeMatx_choice(startTrial:end,startTime:endTime),2)*1000/(width);
    spikeMeans = zeros(numBins,1);
    spikeSems = zeros(numBins,1);
    peMeans = zeros(numBins,1);
    for k = 1:numBins
        if k < numBins
            spikeNumsTemp = spikeCounts(peTemp >= edges(k) & peTemp < edges(k+1));
            peMeans(k) = mean(peTemp(peTemp >= edges(k) & peTemp < edges(k+1)));
        else
            spikeNumsTemp = spikeCounts(peTemp >= edges(k) & peTemp <= edges(k+1));
            peMeans(k) = mean(peTemp(peTemp >= edges(k) & peTemp <= edges(k+1)));
        end
        spikeMeans(k) = mean(spikeNumsTemp);
        spikeSems(k) = sem(spikeNumsTemp);
    end
    
    %% plots 
    GLM=figure2('position',[0 0 800 1500]);
    colors = jet(size(combineMat,2));
    subplot(4,1,[1 2]); hold on;
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

    yyaxis right
    fill([x fliplr(x)], [coeff(:,1,2)' fliplr(coeff(:,1,3)')], [0.5 0.5 0.5], 'facealpha', 0.25, 'edgecolor', 'none')
    plt = gca;
    plt.YAxis(2).Color = [0.2 0.2 0.2];

    yyaxis left
    line(minmax(x), [0 0],'Color', [0.2 0.2 0.2], 'LineStyle','--')
    line([0 0], 1.2*[min(coeff(:,:,2),[],'all'),max(coeff(:,:,3),[],'all')],'Color', [0.2 0.2 0.2], 'LineStyle','--')
    line([os.rwdDelay os.rwdDelay], 1.2*[min(coeff(:,:,2),[],'all'),max(coeff(:,:,3),[],'all')],'Color', [1 0 0], 'LineStyle','--')
    ylim([1.2*min(coeff(:,2:end,2),[],'all') 1.2*max(coeff(:,2:end,3),[],'all')]);
    xlim(minmax(x))
    ylabel('\beta coefficients')
    xlabel('respond time')
    legend('inter','preITI', 'Qdiff')
    titileText = [session ': ' spikeFields{clust(i)} ' Aligned to choice'];
    title(titileText,'interpreter','none');

    subplot(4,1,3); hold on;
    plot(x, rsq, 'linewidth', 2)
    ylim([0 1.2*max(rsq)])
    title('adjusted Rsquare')
    
    subplot(4,1,4); hold on;
    errorbar(peMeans, spikeMeans, spikeSems,'linewidth', 2)
    title('Spike rate change with pe')
    
    if p.Results.saveFigFlag == 1 
        saveFigurePDF(GLM,[savepath sep session '_' spikeFields{clust(i)} 'GLM'])
    end
    
end



%% analysis

% load behavior



