function spikeLvsR_dF(session, col, varargin)
%task and model parameters
p = inputParser;
% default parameters if none given
p.addParameter('cellName', ['all']);
p.addParameter('plotFlag', 1);
p.addParameter('maxTrial', 1000);
p.addParameter('modelName', '5params');
p.addParameter('regressors', '1+pe+biasSide+pe:biasSide')
p.addParameter('binSize', 300)% in ms
p.addParameter('stepSize', 100)
p.addParameter('tb', 1)% in s
p.addParameter('tf', 2)% in s
p.addParameter('saveFigFlag', 1);
p.addParameter('plotAgainst', 'pe')
p.addParameter('focusWindow', [300 1300]);
p.addParameter('numBins', 6);
p.addParameter('width', 1000);
p.parse(varargin{:});

numBins = p.Results.numBins; % bins in pe
width = p.Results.width; %in ms
paramNames = getParamNames_dF(p.Results.modelName, 1);
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
sampFile = [animalName col '_', p.Results.modelName];
path = [root animalName sep animalName 'sorted' sep 'stan' sep 'bernoulli' sep p.Results.modelName sep col sep];

%% behavior preparation 
% parse behavior
os = behAnalysisNoPlot_opMD(session);
choice = os.allChoices';
choice(choice<0) = 0;
outcome = abs(os.allRewards)';
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
% generate best estimates of parameters
[t,~,badSessionFlag] = getStanModelParams_samps(p.Results.modelName, [path sampFile '.mat'], 2000, 'sessionName', session);
if badSessionFlag
    fprintf(['not good behavior in ' session '\n'])
    return
end
% diff value
Qdiff = abs(t.Q(:,2)-t.Q(:,1));
% total value
Qsum = sum(t.Q,2);
% prepe
prePe = [NaN; t.pe(1:end-1)];
% prePe = [NaN; pe(1:end-1)'];
% chosen valie
Qchosen  = zeros(length(choice),1);
for j = 1:length(choice)
    if choice(j)>0
        Qchosen(j) = t.Q(j,2);
    else
        Qchosen(j) = t.Q(j,1);
    end
end
%time in session
timeInSession = [sessionData(responseInds).CSon]' - sessionData(responseInds(1)).CSon;

%% 
biasSide = zeros(size(responseInds))';
biasInd = contains( paramNames,'bias');
if mean(t.params(:,biasInd))>0
    biasSide(os.lickR_Inds)=1;
else
    biasSide(os.lickL_Inds)=1;
end
rightSide = zeros(length(choice),1);
rightSide(os.lickR_Inds) = 1;
choiceConf = 2*t.probChoice - 1;
pe = t.pe;
hmm = double(os.hmmStates==1)';
lickLat = os.lickLatLogZ';
if contains(p.Results.modelName, '7params_absPePeAN_scale_int_bias_ord')
    aN = t.aN;
    peBar = t.peBar;
    pePe = t.pePe;
    tbl = table(outcome, pe, prePe, preRwd, Qsum, Qdiff, choiceConf, biasSide, rightSide, timeInSession, lickLat, hmm, Qchosen, aN, peBar, pePe);
else
    tbl = table(outcome, pe, prePe, preRwd, Qsum, Qdiff, choiceConf, biasSide, rightSide, timeInSession, lickLat, hmm, Qchosen);
end
names = tbl.Properties.VariableNames;
% zscore all regressors
for cols = 1:length(names)
    tmp = tbl.(names{cols});
    tmp(~isnan(tmp)) = zscore(tmp(~isnan(tmp)));
    tbl.(names{cols}) = tmp;
end

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

    coeff = [];
    rsq = zeros(size(midPoints));
    sigs = [];
    for k = 1:length(midPoints)
        currSpikes = allTrial_spikeMatx_slide(:,k);
        currTbl = addvars(tbl, currSpikes);
        lm = fitlm(currTbl, ['currSpikes~' p.Results.regressors]);
        for j = 1:length(lm.CoefficientNames)
            coeff(k,j,1) = lm.Coefficients.Estimate(j);
            ci = coefCI(lm);
            coeff(k,j,2:3) = ci(j,:);
            sigs(k,j) = double(lm.Coefficients.pValue(j)<0.05);
        end
        rsq(k)=lm.Rsquared.Adjusted;
    end
    regressors = lm.CoefficientNames(2:end);
    %% bin no.spike by target regressor and plot

        eval(['target =' p.Results.plotAgainst ';']);

        % decide start of the time windows
        solenoidTimeEnd = max([sessionData(responseInds).CSon] + 3100 - [sessionData(responseInds).respondTime]); % 3600 ms is 500ms after sol coming back after 3100ms
        solenoidTimeStart = min([sessionData(responseInds).CSon] + 3100 - [sessionData(responseInds).respondTime]);

        peInd = find(contains(regressors,p.Results.plotAgainst));
        if ~isempty(peInd) % if pe is regressor
            sigIndEarly = find(sum(sigs(1:end-1,peInd+1),2)>0 & sum(sigs(2:end,peInd+1),2)>0 & slideTime(1:end-1)'<solenoidTimeStart & slideTime(1:end-1)'>os.rwdDelay); 
            if ~isempty(sigIndEarly) % if any significance
                peEffect = sum(squeeze(abs(coeff(sigIndEarly,peInd+1,1))),2);
                [~,maxInd] = max(peEffect);
%                 startTime = max(midPoints(sigIndEarly(maxInd)) - 0.5*width, 1000*p.Results.tb + os.rwdDelay);
                startTime = 1000*p.Results.tb + os.rwdDelay;
            else 
                startTime = 1000*p.Results.tb + os.rwdDelay;
            end
            endTime  = min(startTime+width, solenoidTimeStart+1000*p.Results.tb);

            sigIndLate = find(sum(sigs(:,peInd+1),2)>0 & slideTime' > solenoidTimeEnd & slideTime'<1000*p.Results.tf);
            if ~isempty(sigIndLate)
                peLateEffect = sum(squeeze(abs(coeff(sigIndLate,peInd+1,1))),2);
                [~,maxInd] = max(peLateEffect);
                startTimeLate = max(midPoints(sigIndLate(maxInd)) - 0.5*width, solenoidTimeEnd + 1000*p.Results.tb);
            else 
                startTimeLate = solenoidTimeEnd + 1000*p.Results.tb;
            end
            endTimeLate  = min(startTimeLate+width, 1000*(p.Results.tb+p.Results.tf));
        else
            startTime = 1000*p.Results.tb + os.rwdDelay;
            endTime  = min(startTime+width, solenoidTimeStart+1000*p.Results.tb);
            startTimeLate = solenoidTimeEnd + 1000*p.Results.tb;
            endTimeLate  = min(startTimeLate+width, 1000*(p.Results.tb+p.Results.tf));
        end
        startTime = 1000*p.Results.tb + os.rwdDelay;
        endTime  = min(startTime+width, solenoidTimeStart+1000*p.Results.tb);
        spikeCounts = nansum(allTrial_spikeMatx_choice(:,startTime:endTime),2)*1000/(width);
        spikeCountsLate = nansum(allTrial_spikeMatx_choice(:,startTimeLate:endTimeLate),2)*1000/(width);
        %% caculate on both sides
        edges = binEqualSize(target, numBins);
        edges = linspace(min(target)-0.001, max(target)+0.001, numBins+1);
        edges(0.5*numBins+1) = 0;
        spikeMeans = zeros(numBins,1);
        spikeSems = zeros(numBins,1);
        targetMeans = zeros(numBins,1);
        spikeMeansLate = zeros(numBins,1);
        spikeSemsLate = zeros(numBins,1);
        for k = 1:numBins
            if k < numBins
                spikeNumsTemp = spikeCounts(target >= edges(k) & target < edges(k+1));
                spikeNumsTempLate = spikeCountsLate(target >= edges(k) & target < edges(k+1));
                targetMeans(k) = mean(target(target >= edges(k) & target < edges(k+1)));
            else
                spikeNumsTemp = spikeCounts(target >= edges(k) & target <= edges(k+1));
                spikeNumsTempLate = spikeCountsLate(target >= edges(k) & target <= edges(k+1));
                targetMeans(k) = mean(target(target >= edges(k) & target <= edges(k+1)));
            end
            spikeMeans(k) = mean(spikeNumsTemp);
            spikeMeansLate(k) = mean(spikeNumsTempLate);
            spikeSems(k) = sem(spikeNumsTemp);
            spikeSemsLate(k) = sem(spikeNumsTempLate);
        end

        %% separate left
        targetL = target(os.lickL_Inds);
        spikeCountsL = spikeCounts(os.lickL_Inds);
        spikeCountsLLate = spikeCountsLate(os.lickL_Inds);
        edgesL = binEqualSize(targetL, numBins);
        edgesL = linspace(min(targetL)-0.001, max(targetL)+0.001, numBins+1);
        edgesL(0.5*numBins+1) = 0;
        spikeMeansL = zeros(numBins,1);
        spikeSemsL = zeros(numBins,1);
        spikeMeansLLate = zeros(numBins,1);
        spikeSemsLLate = zeros(numBins,1);
        targetMeansL = zeros(numBins,1);
        for k = 1:numBins
            if k < numBins
                spikeNumsTemp = spikeCountsL(targetL >= edgesL(k) & targetL < edgesL(k+1));
                spikeNumsTempLate = spikeCountsLLate(targetL >= edgesL(k) & targetL < edgesL(k+1));
                targetMeansL(k) = mean(targetL(targetL >= edgesL(k) & targetL < edgesL(k+1)));
            else
                spikeNumsTemp = spikeCountsL(targetL >= edgesL(k) & targetL <= edgesL(k+1));
                spikeNumsTempLate = spikeCountsLLate(targetL >= edgesL(k) & targetL <= edgesL(k+1));
                targetMeansL(k) = mean(targetL(targetL >= edgesL(k) & targetL <= edgesL(k+1)));
            end
            spikeMeansL(k) = mean(spikeNumsTemp);
            spikeMeansLLate(k) = mean(spikeNumsTempLate);
            spikeSemsL(k) = sem(spikeNumsTemp);
            spikeSemsLLate(k) = sem(spikeNumsTempLate);
        end

        %% separate right
        targetR = target(os.lickR_Inds);
        spikeCountsR = spikeCounts(os.lickR_Inds);
        spikeCountsRLate = spikeCountsLate(os.lickR_Inds);
        edgesR = binEqualSize(targetR, numBins);
        edgesR = linspace(min(targetR)-0.001, max(targetR)+0.001, numBins+1);
        edgesR(0.5*numBins+1) = 0;
        spikeMeansR = zeros(numBins,1);
        spikeSemsR = zeros(numBins,1);
        targetMeansR = zeros(numBins,1);
        spikeMeansRLate = zeros(numBins,1);
        spikeSemsRLate = zeros(numBins,1);
        for k = 1:numBins
            if k < numBins
                spikeNumsTemp = spikeCountsR(targetR >= edgesR(k) & targetR < edgesR(k+1));
                spikeNumsTempLate = spikeCountsRLate(targetR >= edgesR(k) & targetR < edgesR(k+1));
                targetMeansR(k) = mean(targetR(targetR >= edgesR(k) & targetR < edgesR(k+1)));
            else
                spikeNumsTemp = spikeCountsR(targetR >= edgesR(k) & targetR <= edgesR(k+1));
                spikeNumsTempLate = spikeCountsRLate(targetR >= edgesR(k) & targetR <= edgesR(k+1));
                targetMeansR(k) = mean(targetR(targetR >= edgesR(k) & targetR <= edgesR(k+1)));
            end
            spikeMeansR(k) = mean(spikeNumsTemp);
            spikeMeansRLate(k) = mean(spikeNumsTempLate);
            spikeSemsR(k) = sem(spikeNumsTemp);
            spikeSemsRLate(k) = sem(spikeNumsTempLate);        
        end

            
    end
    %% plots 
    GLM=figure2('position',[0 0 800 1600]);
    colors = cool(length(regressors));
    colorR = [1 0 1];
    colorL = [0 1 1];
    
    subplot(4,2,[1,2]); hold on;
    x = slideTime;
    yyaxis right
    plot(x, coeff(:,1,1), 'Color', [0.5 0.5 0.5], 'LineStyle','-', 'Marker','none', 'linewidth', 2);
    ylim([1.2*min(coeff(:,:,2),[],'all') 1.2*max(coeff(:,:,3),[],'all')])
    ylabel('intercept')

    yyaxis left 
    for j = 1:length(regressors)
        plot(x, coeff(:,j+1,1), 'Color', colors(j,:), 'LineStyle','-', 'Marker','none', 'linewidth', 2); 
    end

    for j = 1:length(regressors)
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
    line([startTime-1000*p.Results.tb endTime-1000*p.Results.tb], [max(coeff(:,2:end,3),[],'all'),max(coeff(:,2:end,3),[],'all')],'Color', [1 0 0], 'LineWidth',4)
    line([startTimeLate-1000*p.Results.tb endTimeLate-1000*p.Results.tb],[max(coeff(:,2:end,3),[],'all'),max(coeff(:,2:end,3),[],'all')],'Color', [0 0 1], 'LineWidth',4)
    ylim([1.2*min(coeff(:,2:end,2),[],'all') 1.2*max(coeff(:,2:end,3),[],'all')]);
    xlim(minmax(x))
    ylabel('\beta coefficients')
    xlabel('respond time')
    legend(regressors, 'Location', 'best')
    titileText = [session ': ' spikeFields{clust(i)} ' Aligned to choice'];
    title(titileText,'interpreter','none');


    subplot(4,2,3); hold on;
    errorbar(targetMeans, spikeMeans, spikeSems,'linewidth', 2)
    title('Spike rate change with pe')

    subplot(4,2,4); hold on;
    errorbar(targetMeans, spikeMeansLate, spikeSemsLate,'linewidth', 2)
    title('Spike rate after solenoid change with pe')

    subplot(4,2,5); hold on;
    errorbar(targetMeansL, spikeMeansL, spikeSemsL,'linewidth', 2, 'color', colorL)
    title(['Spike rate change with left pe' sprintf(('left = %d'), length(os.lickL_Inds))])

    subplot(4,2,6); hold on;
    errorbar(targetMeansL, spikeMeansLLate, spikeSemsLLate,'linewidth', 2, 'color', colorL);
    title('Spike rate after solenoid change with left pe')

    subplot(4,2,7); hold on;
    errorbar(targetMeansR, spikeMeansR, spikeSemsR,'linewidth', 2, 'color', colorR)
    title(['Spike rate change with right pe' sprintf(('right = %d'), length(os.lickR_Inds))])

    subplot(4,2,8); hold on;
    errorbar(targetMeansR, spikeMeansRLate, spikeSemsRLate,'linewidth', 2, 'color', colorR)
    title('Spike rate after solenoid change with right pe')

   
    if p.Results.saveFigFlag == 1 
        saveFigurePDF(GLM,[savepath sep session '_' spikeFields{clust(i)} 'GLM'])
    end
    close(GLM)
    split = figure;
    maxFreq = 1.5*max(mean(allTrial_spikeMatx_slide));
    % sorted by PE raster PSTH
    [~, ind] = sort(pe);
    indL = ismember(ind, os.lickL_Inds);
    indL = ind(indL);
    indR = ismember(ind, os.lickR_Inds);
    indR = ind(indR);
    subplot(3,5,[1 6]);
    title('L');
    sortedSpikes = allTrial_spike_choice(indL);
    LineFormat.LineWidth = 1.5;
    plotSpikeRaster(sortedSpikes,'PlotType','vertline', 'LineFormat', LineFormat);
    line([os.rwdDelay os.rwdDelay], [0 length(sortedSpikes)], 'color', [1 0.3 0.3], 'LineStyle', '--')
    xlim([-1000 2000]);
    ylim([0 max([length(indL) length(indR)])])
    set(gca, 'YTick', [])
    ylabel('sorted by pe 1 <----> -1')
    subplot(3,5,[2 7]);
    title('R');
    sortedSpikes = allTrial_spike_choice(indR); 
    plotSpikeRaster(sortedSpikes,'PlotType','vertline', 'LineFormat', LineFormat);
   
    line([os.rwdDelay os.rwdDelay], [0 length(sortedSpikes)], 'color', [1 0.3 0.3], 'LineStyle', '--')
    xlim([-1000 2000]);
    ylim([0 max([length(indL) length(indR)])])
    set(gca, 'YTick', [])
    % PSTH - rwd/nrwd[~, ind] = sort(pe);
    indL = ismember(ind, os.lickL_Inds);
    indL = ind(indL);
    peL = pe(indL);
    indR = ismember(ind, os.lickR_Inds);
    indR = ind(indR);
    peR = pe(indR);
    % L
    edges = [linspace(min(peL)-0.001, 0, 3) linspace(0, max(peL)+0.001, 3)];
    edges = unique(edges);
    indI = indL((peL>=edges(1) & peL<=edges(2)));
    indII = indL((peL>edges(2) & peL<=edges(3)));
    indIII = indL((peL>edges(3) & peL<=edges(4)));
    indIV = indL((peL>edges(4) & peL<=edges(5)));

    mySDF_I = allTrial_spikeMatx_slide(indI,:);
    mySDF_II = allTrial_spikeMatx_slide(indII,:);
    mySDF_III = allTrial_spikeMatx_slide(indIII,:);
    mySDF_IV = allTrial_spikeMatx_slide(indIV,:);

    subplot(3,5,13); hold on
    title('L')
    plotFilled(slideTime, mySDF_I,[1 0 0])
    plotFilled(slideTime, mySDF_II,[1 0.5 0.5])
    plotFilled(slideTime, mySDF_III,[0.5 0.5 1])
    plotFilled(slideTime, mySDF_IV,[0 0 1])
    line([os.rwdDelay os.rwdDelay], [0 7], 'color', [1 0.3 0.3], 'LineStyle', '--')
    ylim([0 maxFreq]);
    xlim([-1000 2000]);
    set(gca, 'XTick', [-1000:1000:2000])
    set(gca, 'YTick', [0:5:10])
    set(gca,'tickdir', 'out')
    xlabel('time from choice')
    ylabel('spikes/s')


    % R
    edges = [linspace(min(peR)-0.001, 0, 3) linspace(0, max(peR)+0.001, 3)];
    edges = unique(edges);
    indI = indR((peR>=edges(1) & peR<=edges(2)));
    indII = indR((peR>edges(2) & peR<=edges(3)));
    indIII = indR((peR>edges(3) & peR<=edges(4)));
    indIV = indR((peR>edges(4) & peR<=edges(5)));

    mySDF_I = allTrial_spikeMatx_slide(indI,:);
    mySDF_II = allTrial_spikeMatx_slide(indII,:);
    mySDF_III = allTrial_spikeMatx_slide(indIII,:);
    mySDF_IV = allTrial_spikeMatx_slide(indIV,:);
    
    subplot(3,5,14); hold on
    title('R')
    plotFilled(slideTime, mySDF_I,[1 0 0])
    plotFilled(slideTime, mySDF_II,[1 0.5 0.5])
    plotFilled(slideTime, mySDF_III,[0.5 0.5 1])
    plotFilled(slideTime, mySDF_IV,[0 0 1])
    line([os.rwdDelay os.rwdDelay], [0 7], 'color', [1 0.3 0.3], 'LineStyle', '--')
    ylim([0 maxFreq]);
    xlim([-1000 2000]);
    set(gca, 'YTick', [0:5:10])
    set(gca,'tickdir', 'out')
    set(gca, 'XTick', [-1000:1000:2000])
    xlabel('time from choice')
    
    % spikes/s - PE
    subplot(3,5,3); hold on;
    errorbar(targetMeans, spikeMeans, spikeSems,'linewidth', 2)
    title('allTrials')
    xlim([-1 1])
    set(gca, 'tickdir', 'out')
    set(gca, 'XTick', [-1 0 1])
    ylabel('spikes/s')
    xlabel('rpe')
    subplot(3,5,4); hold on;
    errorbar(targetMeansL, spikeMeansL, spikeSemsL,'linewidth', 2, 'color', colorL)
    errorbar(targetMeansR, spikeMeansR, spikeSemsR,'linewidth', 2, 'color', colorR)
    xlim([-1 1])
    set(gca, 'tickdir', 'out')
    set(gca, 'XTick', [-1 0 1])
    legend({'L', 'R'})
    title('split by side')
    xlabel('rpe') 
    
    % separate by rwd nrwd
    subplot(3,5,11); hold on;
    title('Left')
    mySDF_rwd = allTrial_spikeMatx_slide(intersect(os.rwd_Inds, os.lickL_Inds),:);
    mySDF_noRwd = allTrial_spikeMatx_slide(intersect(os.nrwd_Inds, os.lickL_Inds),:);
    plotFilled(slideTime, mySDF_rwd,[0 0 1])
    plotFilled(slideTime, mySDF_noRwd,[1 0 0])
    line([os.rwdDelay os.rwdDelay], [0 15], 'color', [1 0.3 0.3], 'LineStyle', '--')
    ylim([0 maxFreq]);
    xlim([-1000 2000]);
    set(gca, 'YTick', [0:5:10])
    set(gca, 'XTick', [-1000:1000:2000])
    set(gca,'tickdir', 'out')
    legend({'rwd', '', 'nrwd', ''})
    subplot(3,5,12); hold on;
    xlabel('time from respond /ms')
    ylabel('spikes/s')
    title('Right')
    mySDF_rwd = allTrial_spikeMatx_slide(intersect(os.rwd_Inds, os.lickR_Inds),:);
    mySDF_noRwd = allTrial_spikeMatx_slide(intersect(os.nrwd_Inds, os.lickR_Inds),:);
    plotFilled(slideTime, mySDF_rwd,[0 0 1])
    plotFilled(slideTime, mySDF_noRwd,[1 0 0])
    line([os.rwdDelay os.rwdDelay], [0 15], 'color', [1 0.3 0.3], 'LineStyle', '--')
    ylim([0 maxFreq]);
    xlim([-1000 2000]);
    set(gca, 'YTick', [0:5:10])
    set(gca,'tickdir', 'out')
    legend({'rwd', '', 'nrwd', ''})
    xlabel('time from respond /ms')

    screen = get(0,'Screensize');
    screen(4) = screen(4) - 100;
    set(gcf, 'Position', screen)
    
    
end
