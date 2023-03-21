% load data and params
load('F:\tmpData\catWithOutcome.mat')
tb = 2;
tf = 3;
focusWin = [300 1000];
modelName = '5params';
maxTrial = 1000;
binSize = 200;
stepSize = 100;
numBins = 8;
color1 = [0 0.8 0.8];
color2 = [1 0 1];
%% preparation
[root, sep] = currComputer();
paramNames = getParamNames_dF(modelName, 1);
time = -1000*tb:1000*tf;
midPoints = (0.5*binSize + 1):stepSize:(length(time)-0.5*binSize);
slideTime = midPoints - tb*1000;
allSigs = zeros(size(allSessions));
allTstats = zeros(size(allSessions));
allCoeffs = zeros(size(allSessions));
allPe = zeros(length(allSessions), numBins);
focusSpikes = zeros(length(allSessions), numBins);
allPeL = zeros(length(allSessions), numBins);
focusSpikesL = zeros(length(allSessions), numBins);
allPeR = zeros(length(allSessions), numBins);
focusSpikesR = zeros(length(allSessions), numBins);
allFocusSpikes = cell(size(allSessions));
allChoices = cell(size(allSessions));
allPeRaw = cell(size(allSessions));
%% loop through with side interation
loadedSession = [];
category = 'good';
for ses = 1:length(allSessions)
    session = allSessions{ses};
    unit = allUnits{ses};
    fprintf([session unit '\n']);
    % paths
    pd = parseSessionString_df(session, root, sep);
    neuralynxDataPath = [pd.sortedFolder session '_sessionData_nL.mat'];
    unitMetDir = [pd.nLynxFolderSession session '_' unit '_met.mat'];
    sortedFolderLocation = [pd.sortedFolder 'session' sep];
    sampFile = [pd.aniName category '_', modelName];
    path = [root pd.aniName sep pd.aniName 'sorted' sep 'stan' sep 'bernoulli' sep modelName sep category sep];

    % decide is good behavior
    % load good days
    xlFile = [pd.animalName '.xlsx'];
    [~ , goodDayList, ~] = xlsread([root xlFile], pd.aniName);
    [~,col] = find(~cellfun(@isempty,strfind(goodDayList, category)) == 1);
    goodDayList = goodDayList(2:end,col);
    
    if ~ismember(session, goodDayList)
        continue
    end
    % decide if a good unit
    if exist(unitMetDir,'file')
        load(unitMetDir);
    else
        met = getClusterMetric(session, unit, 0, 1);
    end
    respInds = find(met.spikeProp>=0.8);
    if isempty(respInds)
        continue
    else 
        respLat = nanmin(met.spikeLat(respInds));
        if respLat > 15000 || met.isiV > 0.001 || met.distance>0.3 || isnan(met.distance)
            continue
        end
    end
    % load behavior and neurons
    if exist(neuralynxDataPath,'file')
        load(neuralynxDataPath)
    else
        sessionData = generateSessionData_nL_operantMatching(session);
    end

    if ~strcmp(session,loadedSession) % avoiding recomputing model variables for different unit from same session
        %% behavior preparation 
        % parse behavior
        os = behAnalysisNoPlot_opMD(session);
        if length(os.CSplus) ~= length(sessionData)
            fprintf([session ' error\n']);
            return
        end
        choice = os.allChoices';
        choice(choice<0) = 0;
        outcome = abs(os.allRewards)';
        choice = choice(1:min(length(choice), maxTrial));
        outcome = outcome(1:length(choice));
        responseInds = os.responseInds(1:min(length(choice), maxTrial)); 
        preRwd = [NaN abs(os.allRewards(1:end-1))]';
        %% behavior
        % switch
        svs = zeros(length(os.responseInds),1);
        svs(os.changeChoice_Inds) = 1;
        svsNext = [svs(2:end); NaN];
        [t,~,noSession] = getStanModelParams_samps(modelName, [path sampFile '.mat'], 2000, 'sessionName', session);
        if noSession
            fprintf(['no good behavior in ' session '\n']);
            populationDepths = populationDepths(1:end-1);
            continue
        end

        % diff value
        Qdiff = abs(t.Q(:,2)-t.Q(:,1));
        % total value
        Qsum = sum(t.Q,2);
        % prepe
        prePe = [NaN; t.pe(1:end-1)];
        % pe
        pe = t.pe;
        % dawExp
        dawExp = double(t.probChoice <= 0.5);
        % confidence
        choiceConf = 2.*t.probChoice - 1;
        %time in session
        timeInSession = [sessionData(responseInds).CSon]' - sessionData(responseInds(1)).CSon;
        % chosen valie
        Qchosen  = zeros(length(choice),1);
        Qunchosen  = zeros(length(choice),1);
        QchosenUpdate = NaN(length(choice),1);
        for j = 1:length(choice)
            if j < length(choice)
                if choice(j)>0
                    Qchosen(j) = t.Q(j,2);
                    Qunchosen(j) = t.Q(j,1);
                    QchosenUpdate(j) = t.Q(j+2);
                else
                    Qchosen(j) = t.Q(j,1);
                    Qunchosen(j) = t.Q(j,2);
                    QchosenUpdate(j) = t.Q(j+1);
                end
            else                
                if choice(j)>0
                    Qchosen(j) = t.Q(j,2);
                    Qunchosen(j) = t.Q(j,1);
                else
                    Qchosen(j) = t.Q(j,1);
                    Qunchosen(j) = t.Q(j,2);
                end
                
            end
        end
        % bias side
        biasSide = zeros(size(responseInds))';
        biasInd = contains(paramNames, 'bias');
        if mean(t.params(:,biasInd))>0
            biasSide(os.lickR_Inds)=1;
        else
            biasSide(os.lickL_Inds)=1;
        end
        hmm = double(os.hmmStates==1)';
        lickLat = os.lickLatLogZ';
        rightSide = zeros(size(pe));
        rightSide(os.allChoices>0)=1;
        rightSide(os.allChoices<0)=-1;
        preITI = os.timeBtwn';

        tbl = table(outcome, pe, prePe, preRwd, Qsum, Qdiff, choiceConf, biasSide, rightSide, timeInSession, lickLat, hmm, Qchosen, Qunchosen, QchosenUpdate, preITI, dawExp, svsNext);

        names = tbl.Properties.VariableNames;
        % zscore all regressors
        for cols = 1:length(names)
            tmp = tbl.(names{cols});
            if ~isempty(setdiff(tmp(~isnan(tmp)), [0 1 -1 NaN]))
                tmp(~isnan(tmp)) = zscore(tmp(~isnan(tmp)));
                tbl.(names{cols}) = tmp;
            end
        end
        
        loadedSession = session;
    end
  
    %% create spike
    spikeFields = fields(sessionData);
    clust = find(contains(spikeFields,unit));
    allTrial_spike_choice = {};

    for k = 1:length(os.responseInds)
            if os.responseInds(k) == 1
                prevTrial_spike = [];
            else
                prevTrial_spikeInd = [sessionData(os.responseInds(k)-1).(spikeFields{clust})] > (sessionData(os.responseInds(k)).respondTime-tb*1000);
                prevTrial_spike = sessionData(os.responseInds(k)-1).(spikeFields{clust})(prevTrial_spikeInd) - sessionData(os.responseInds(k)).respondTime;
            end

            currTrial_spikeInd = sessionData(os.responseInds(k)).(spikeFields{clust}) < sessionData(os.responseInds(k)).respondTime+tf*1000 ... 
                & sessionData(os.responseInds(k)).(spikeFields{clust}) > sessionData(os.responseInds(k)).respondTime-tb*1000;
            currTrial_spike = sessionData(os.responseInds(k)).(spikeFields{clust})(currTrial_spikeInd) - sessionData(os.responseInds(k)).respondTime;

            allTrial_spike_choice{k} = [prevTrial_spike currTrial_spike];
    end

    % sometimes no licks/spikes are considered 1x0 and sometimes they are []
    % plotSpikeRaster does not place nicely with [] so this converts all empty indices to 1x0
    allTrial_spike_choice(cellfun(@isempty,allTrial_spike_choice)) = {zeros(1,0)}; 
    % spike matrix choice
    allTrial_spikeMatx_choice = zeros(length(os.responseInds),length(time));
    trialDurDiff = zeros(1, length(os.responseInds));
    for j = 1:length(os.responseInds)
        trialDurDiff(j) = (sessionData(os.responseInds(j)).trialEnd - sessionData(os.responseInds(j)).CSon)- tf*1000;
    end
   for j = 1:length(allTrial_spike_choice)
        tempSpike = allTrial_spike_choice{j};
        tempSpike = tempSpike + tb*1000; % add this to pad time for SDF
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
            nansum(allTrial_spikeMatx_choice(:,midPoints(w)-0.5*binSize:midPoints(w)+0.5*binSize-1),2)*1000/binSize;
    end
    
    allTrial_spikeMatx_slide = zscore(allTrial_spikeMatx_slide, 0, 'all');
    %%  create spike cell aligned to cue
    allTrial_spike_cue = {};

    for k = 1:length(os.behSessionData)
            if k == 1
                prevTrial_spike = [];
            else
                prevTrial_spikeInd = [sessionData(k-1).(spikeFields{clust})] > (sessionData(k).CSon-tb*1000);
                prevTrial_spike = sessionData(k-1).(spikeFields{clust})(prevTrial_spikeInd) - sessionData(k).CSon;
            end

            currTrial_spikeInd = sessionData(k).(spikeFields{clust}) < sessionData(k).CSon+tf*1000 ... 
                & sessionData(k).(spikeFields{clust}) > sessionData(k).CSon-tb*1000;
            currTrial_spike = sessionData(k).(spikeFields{clust})(currTrial_spikeInd) - sessionData(k).CSon;

            allTrial_spike_cue{k} = [prevTrial_spike currTrial_spike];
    end

    % sometimes no licks/spikes are considered 1x0 and sometimes they are []
    % plotSpikeRaster does not place nicely with [] so this converts all empty indices to 1x0
    allTrial_spike_cue(cellfun(@isempty,allTrial_spike_cue)) = {zeros(1,0)}; 
    
    %  spike matrix cue 
    allTrial_spikeMatx_cue = zeros(length(os.behSessionData),length(time));
    trialDurDiff = zeros(1, length(os.behSessionData));
    for j = 1:length(os.behSessionData)
        trialDurDiff(j) = (sessionData(j).trialEnd - sessionData(j).CSon)- tf*1000;
    end
   for j = 1:length(allTrial_spike_cue)
        tempSpike = allTrial_spike_cue{j};
        tempSpike = tempSpike + tb*1000; % add this to pad time for SDF
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
    allTrial_spikeMatx_slide_cue = zeros(length(os.behSessionData), length(midPoints));
    for w = 1:length(midPoints)
        allTrial_spikeMatx_slide_cue(:,w) = ...
            nansum(allTrial_spikeMatx_cue(:,midPoints(w)-0.5*binSize:midPoints(w)+0.5*binSize-1),2)*1000/binSize;
    end
    
    allTrial_spikeMatx_slide_cue = zscore(allTrial_spikeMatx_slide_cue, 0, 'all');
    %%
    % focus window zscore 
    spikeCounts = zscore(nansum(allTrial_spikeMatx_choice(:,tb*1000+focusWin(1):tb*1000+focusWin(2)),2));
    

    %% calculate and combine pe
    % bin spikes into pe bins and calcualted bin 
    numBins = 8;
    currTarget = tbl.pe;
    spikeMeans = zeros(1,numBins);
    targetMeans = zeros(1,numBins);
    spikeMeansR = zeros(1,numBins);
    targetMeansR = zeros(1,numBins);
    spikeMeansL = zeros(1,numBins);
    targetMeansL = zeros(1,numBins);
    % all
    edges = linspace(min(currTarget)-0.001, max(currTarget)+0.001, numBins+1);
    for k = 1:numBins
        spikeNumsTemp = spikeCounts(currTarget >= edges(k) & currTarget < edges(k+1));
        targetMeans(k) = mean(currTarget(currTarget >= edges(k) & currTarget < edges(k+1)), 'omitnan');     
        spikeMeans(k) = mean(spikeNumsTemp, 'omitnan');
    end
    % R
    currTarget = tbl.pe;
    currTarget = currTarget(os.allChoices==1);
    spikeCountsR = spikeCounts(os.allChoices==1);
    edges = linspace(min(currTarget)-0.001, max(currTarget)+0.001, numBins+1);
    for k = 1:numBins
        spikeNumsTemp = spikeCountsR(currTarget >= edges(k) & currTarget < edges(k+1));
        targetMeansR(k) = mean(currTarget(currTarget >= edges(k) & currTarget < edges(k+1)), 'omitnan');     
        spikeMeansR(k) = mean(spikeNumsTemp, 'omitnan');
    end
    % L
    currTarget = tbl.pe;
    currTarget = currTarget(os.allChoices==-1);
    spikeCountsL = spikeCounts(os.allChoices==-1);
    edges = linspace(min(currTarget)-0.001, max(currTarget)+0.001, numBins+1);
    for k = 1:numBins
        spikeNumsTemp = spikeCountsL(currTarget >= edges(k) & currTarget < edges(k+1));
        targetMeansL(k) = mean(currTarget(currTarget >= edges(k) & currTarget < edges(k+1)), 'omitnan');     
        spikeMeansL(k) = mean(spikeNumsTemp, 'omitnan');
    end
    
    % collecting neuron data
    allPe(ses,:) = targetMeans;
    focusSpikes(ses,:) = spikeMeans;
    allPeR(ses,:) = targetMeansR;
    focusSpikesR(ses,:) = spikeMeansR;
    allPeL(ses,:) = targetMeansL;
    focusSpikesL(ses,:) = spikeMeansL;
    %% calculate regression results with interations
    currTbl = addvars(tbl, spikeCounts);
    lm = fitlm(currTbl, ['spikeCounts~ 1+ rightSide * outcome + Qchosen']);
    catRegressorInd = find(strcmp(lm.CoefficientNames, 'outcome:rightSide'));    
    allCoeffs(ses) = lm.Coefficients.Estimate(catRegressorInd);
    allTstats(ses) = lm.Coefficients.tStat(catRegressorInd);
    if lm.Coefficients.pValue(catRegressorInd)>0.05
       allSigs(ses) = 0;
    else
        if lm.Coefficients.Estimate(catRegressorInd)>0
            allSigs(ses) = 1;
        else
            allSigs(ses) = -1;
        end
    end
    allChoices{ses} = os.allChoices;
    allFocusSpikes{ses} = spikeCounts;
    allPeRaw{ses} = pe;
end 
%% waveform analysis
numCat = 2;
[coeff,score,latent, ~, explained, mu] = pca(zscore(waveformsSession, [], 1));
indAll = {};
dis = {};
for a = 1:10
    [indAll{a}, ~, dis{a}] = kmeans([score(:, 1:5)], numCat);
end
[~,optiInds] = min(cellfun(@mean, dis));
ind = indAll{optiInds};
%%
edges = linspace(min(score(:,1))-0.01, max(score(:,1))-0.01, 30);
figure2;hold on;
histogram(score(ind==2,1), edges, "FaceColor", [1 0.3 1], 'EdgeColor', 'none');
histogram(score(ind==1,1), edges, "FaceColor", [0 0.8 0.8], 'EdgeColor', 'none')
%% tstats distribution
figure2;hold on;
histogram(allTstats(cats == 1), linspace(min(allTstats)-0.001,max(allTstats)+0.001,11), 'FaceColor', [0.6 0.6 0.6]);
histogram(allTstats(cats == 1& allSigs~=0), linspace(min(allTstats)-0.001,max(allTstats)+0.001,11), 'FaceColor', color1);
title('type I')

figure2;hold on;
histogram(allTstats(cats == 2), linspace(min(allTstats)-0.001,max(allTstats)+0.001,11), 'FaceColor', [0.6 0.6 0.6]);
histogram(allTstats(cats == 2& allSigs~=0), linspace(min(allTstats)-0.001,max(allTstats)+0.001,11), 'FaceColor', color2);
title('type II')

figure2;hold on; 
histogram(allTstats(cats == 0), linspace(min(allTstats)-0.001,max(allTstats)+0.001,11), 'FaceColor', [0.6 0.6 0.6]);
histogram(allTstats(cats == 0& allSigs~=0), linspace(min(allTstats)-0.001,max(allTstats)+0.001,11), 'FaceColor', [1 0.4 0.4]);
title('no type')

figure2; hold on;
histogram(allTstats(allSigs==0), linspace(min(allTstats)-0.001,max(allTstats)+0.001,18), 'FaceColor', [0.6 0.6 0.6], 'EdgeColor', 'none', 'FaceAlpha', 0.4);
histogram(allTstats(allSigs>0), linspace(min(allTstats)-0.001,max(allTstats)+0.001,18), 'FaceColor', [1 0.5 0], 'EdgeColor', 'none', 'FaceAlpha', 0.4);
histogram(allTstats(allSigs<0), linspace(min(allTstats)-0.001,max(allTstats)+0.001,18), 'FaceColor', [0 0.5 1], 'EdgeColor', 'none', 'FaceAlpha', 0.4);
title('all')
%% spike-pe trend
figure2; hold on;
currSpikes = focusSpikesR(cats==1,:);
currTarget = allPeR(cats==1, :);
spikeMeans = mean(currSpikes, 'omitnan');
spikeSems = sem(currSpikes);
targetMeans = mean(currTarget, 'omitnan');

plot(targetMeans, spikeMeans, 'color', 'r', 'linewidth', 2, 'Marker', 'none', 'LineStyle', '-');
fill([targetMeans fliplr(targetMeans)], [spikeMeans+spikeSems fliplr(spikeMeans-spikeSems)], 'r', 'facealpha', 0.25, 'edgecolor', 'none')

currSpikes = focusSpikesL(cats==1,:);
currTarget = allPeL(cats==1, :);
spikeMeans = mean(currSpikes, 'omitnan');
spikeSems = sem(currSpikes);
targetMeans = mean(currTarget, 'omitnan');

plot(targetMeans, spikeMeans, 'color', 'g', 'linewidth', 2, 'Marker', 'none', 'LineStyle', '-');
fill([targetMeans fliplr(targetMeans)], [spikeMeans+spikeSems fliplr(spikeMeans-spikeSems)], 'g', 'facealpha', 0.25, 'edgecolor', 'none')
%% type I
figure2; hold on;
currSpikes = focusSpikesR(cats==1,:);
currTarget = allPeR(cats==1, :);
spikeMeans = mean(currSpikes, 'omitnan');
spikeSems = sem(currSpikes);
targetMeans = mean(currTarget, 'omitnan');

plot(targetMeans, spikeMeans, 'color', 'r', 'linewidth', 2, 'Marker', 'none', 'LineStyle', '-');
fill([targetMeans fliplr(targetMeans)], [spikeMeans+spikeSems fliplr(spikeMeans-spikeSems)], 'r', 'facealpha', 0.25, 'edgecolor', 'none')

currSpikes = focusSpikesL(cats==1,:);
currTarget = allPeL(cats==1, :);
spikeMeans = mean(currSpikes, 'omitnan');
spikeSems = sem(currSpikes);
targetMeans = mean(currTarget, 'omitnan');

plot(targetMeans, spikeMeans, 'color', 'g', 'linewidth', 2, 'Marker', 'none', 'LineStyle', '-');
fill([targetMeans fliplr(targetMeans)], [spikeMeans+spikeSems fliplr(spikeMeans-spikeSems)], 'g', 'facealpha', 0.25, 'edgecolor', 'none')
title('type I')
%%
% type II
figure2; hold on;
currSpikes = focusSpikesR(cats==2,:);
currTarget = allPeR(cats==2, :);
spikeMeans = mean(currSpikes, 'omitnan');
spikeSems = sem(currSpikes);
targetMeans = mean(currTarget, 'omitnan');

plot(targetMeans, spikeMeans, 'color', 'r', 'linewidth', 2, 'Marker', 'none', 'LineStyle', '-');
fill([targetMeans fliplr(targetMeans)], [spikeMeans+spikeSems fliplr(spikeMeans-spikeSems)], 'r', 'facealpha', 0.25, 'edgecolor', 'none')

currSpikes = focusSpikesL(cats==2,:);
currTarget = allPeL(cats==2, :);
spikeMeans = mean(currSpikes, 'omitnan');
spikeSems = sem(currSpikes);
targetMeans = mean(currTarget, 'omitnan');

plot(targetMeans, spikeMeans, 'color', 'g', 'linewidth', 2, 'Marker', 'none', 'LineStyle', '-');
fill([targetMeans fliplr(targetMeans)], [spikeMeans+spikeSems fliplr(spikeMeans-spikeSems)], 'g', 'facealpha', 0.25, 'edgecolor', 'none')
title('type II')
legend({'R', '', 'L', ''})

figure2;

%%
% type NA
figure2; hold on;
currSpikes = focusSpikesR(cats==0,:);
currTarget = allPeR(cats==0, :);
spikeMeans = mean(currSpikes, 'omitnan');
spikeSems = sem(currSpikes);
targetMeans = mean(currTarget, 'omitnan');

plot(targetMeans, spikeMeans, 'color', 'r', 'linewidth', 2, 'Marker', 'none', 'LineStyle', '-');
fill([targetMeans fliplr(targetMeans)], [spikeMeans+spikeSems fliplr(spikeMeans-spikeSems)], 'r', 'facealpha', 0.25, 'edgecolor', 'none')

currSpikes = focusSpikesL(cats==0,:);
currTarget = allPeL(cats==0, :);
spikeMeans = mean(currSpikes, 'omitnan');
spikeSems = sem(currSpikes);
targetMeans = mean(currTarget, 'omitnan');

plot(targetMeans, spikeMeans, 'color', 'g', 'linewidth', 2, 'Marker', 'none', 'LineStyle', '-');
fill([targetMeans fliplr(targetMeans)], [spikeMeans+spikeSems fliplr(spikeMeans-spikeSems)], 'g', 'facealpha', 0.25, 'edgecolor', 'none')
title('type NA')
%% coeff scatter 
%% stats scatter
figure2; hold on
scatter(populationCoeffsL(5,:,:), populationCoeffsR(5,:,:), 30, [0.7 0.7 0.7], 'filled');
scatter(populationCoeffsL(5, 1, cats==1), populationCoeffsR(5, 1, cats==1), 30, [0 0.8 0.8], 'filled');
scatter(populationCoeffsL(5, 1, cats==2), populationCoeffsR(5, 1, cats==2), 30, [1 0.2 1], 'filled');
scatter(populationCoeffsL(5, 1, allSigs>0), populationCoeffsR(5, 1, allSigs>0), 50, [1 0.5 0.5], 'd', 'LineWidth',1.5);
scatter(populationCoeffsL(5, 1, allSigs<0), populationCoeffsR(5, 1, allSigs<0), 50, [0.5 0.5 1], 'd', 'LineWidth',1.5);
plot([min(populationCoeffsL(5,:,:)) max(populationCoeffsL(5,:,:))], [0 0], 'Color', [0.5 0.5 0.5], 'LineStyle', ':');   
plot([0 0], [min(populationCoeffsR(5,:,:)) max(populationCoeffsR(5,:,:))],'Color', [0.5 0.5 0.5], 'LineStyle', ':'); 
plot([min(populationCoeffsR(5,:,:)) max(populationCoeffsL(5,:,:))], [min(populationCoeffsR(5,:,:)) max(populationCoeffsL(5,:,:))],'Color', [0.5 0.5 0.5], 'LineStyle', ':'); 
%% stats scatter
figure2; hold on
scatter(populationTStatsL(5,:,:), populationTStatsR(5,:,:), 30, [0.7 0.7 0.7], 'filled');
scatter(populationTStatsL(5, 1, cats==1), populationTStatsR(5, 1, cats==1), 30, [0 0.8 0.8], 'filled');
scatter(populationTStatsL(5, 1, cats==2), populationTStatsR(5, 1, cats==2), 30, [1 0.2 1], 'filled');
scatter(populationTStatsL(5, 1, allSigs>0), populationTStatsR(5, 1, allSigs>0), 50, [1 0.3 0.3], 'd', 'LineWidth',1.5);
scatter(populationTStatsL(5, 1, allSigs<0), populationTStatsR(5, 1, allSigs<0), 50, [0.3 0.3 1], 'd', 'LineWidth',1.5);
plot([min(populationTStatsL(5,:,:)) max(populationTStatsL(5,:,:))], [0 0], 'Color', [0.5 0.5 0.5], 'LineStyle', ':');   
plot([0 0], [min(populationTStatsR(5,:,:)) max(populationTStatsR(5,:,:))],'Color', [0.5 0.5 0.5], 'LineStyle', ':'); 
plot([min(populationTStatsR(5,:,:)) max(populationTStatsL(5,:,:))], [min(populationTStatsR(5,:,:)) max(populationTStatsL(5,:,:))],'Color', [0.5 0.5 0.5], 'LineStyle', ':'); 

%%
set(gca, 'TickDir', 'Out')
set(gca, 'XTick', [-10:10:10], 'FontSize', 14)
set(gca, 'YTick', [0 10 20], 'FontSize', 14)
set(gca, 'Box', 'off')
ylim([0 25]); xlim([-15 15])
xlabel('PC1', 'FontSize', 18)
ylabel('count', 'FontSize', 18)
%%