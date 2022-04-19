function spikeFeatures_population(animalNames, category, varargin)
%task and model parameters
p = inputParser;
% default parameters if none given
p.addParameter('numCat', 2);
p.addParameter('maxTrial', 1000);
% p.addParameter('modelName','7params_absPePeAN_scale_int_bias_ord')
p.addParameter('modelName','5params')
p.addParameter('regressors', '1+Qchosen+outcome*rightSide')
p.addParameter('plotRegressor', 'pe')
p.addParameter('focusWin', [300 1300]); % in ms from respond
p.addParameter('binSize', 100)% in ms
p.addParameter('saveFigFlag', 1);
p.addParameter('tb', 3)% in s
p.addParameter('tf', 3)% in s
p.addParameter('stepSize', 100)
p.addParameter('sepOutcome', 1);
p.parse(varargin{:});

paramNames = getParamNames_dF(p.Results.modelName, 1);
maxTrial = p.Results.maxTrial;
% basic info
[root, sep] = currComputer();
% time window
time = -1000*p.Results.tb:1000*p.Results.tf;
midPoints = (0.5*p.Results.binSize + 1):p.Results.stepSize:(length(time)-0.5*p.Results.binSize);
slideTime = midPoints - p.Results.tb*1000;
%% classification
baselineFreq = [];
firstSpikeFreq = [];
firstSpikeLat = [];
width = [];
spikeLatChange = [];
spikeFreqChange = [];
waveforms = [];
waveformsSession = [];
rwdAUC = [];
rwdAUCSig  = [];
allSessions = {};
allUnits = {};
%% combined information
focusSpikes = []; % zscored within Session; vector 1*allTrials; 
allPe = []; % zscored within session; vector 1*allTrials;
allSpikes = [];
populationDepths = []; % all estimated positions
rwdvsNrwd = [];
focusSpikesBins = [];
allSW = [];
allKernel = [];
excit = []; % baseline vs peak
myKernel = [1, 2, 4, 2, 1];
myKernel = myKernel/sum(myKernel);
%% animal loop
for ani = 1:length(animalNames)
    % load model fitting results
    animalName = animalNames{ani};
    sampFile = [animalName category '_', p.Results.modelName];
    path = [root animalName sep animalName 'sorted' sep 'stan' sep 'bernoulli' sep p.Results.modelName sep category sep];
    % load sessionList and unitList
    xlFile = [animalName '.xlsx'];
    [nums, unitsInfo,~] = xlsread([root xlFile], 'neurons');
    sessionList = unitsInfo(2:end,1); 
    unitList = unitsInfo(2:end,2);
    optoUnitList = unitsInfo(2:end,8);
    subFolders = unitsInfo(2:end,9);
    depths = nums(:, 9);

    drift = unitsInfo(2:end,10);
    drift = contains(drift, 'drift')|contains(drift, 'Drift')|contains(drift, 'duplicate')|contains(drift, 'Duplicate');
    quality = nums(:,1);
    sessionList = sessionList(quality<=0.05 & ~drift);
    unitList = unitList(quality<=0.05 & ~drift);
    depths = depths(quality<=0.05 & ~drift);
    optoUnitList = optoUnitList(quality<=0.05 & ~drift);
    subFolders = subFolders(quality<=0.05 & ~drift);
    
    % load good days
    [~ , goodDayList, ~] = xlsread([root xlFile], animalName);
    [~,col] = find(~cellfun(@isempty,strfind(goodDayList, category)) == 1);
    goodDayList = goodDayList(2:end,col);
    loadedSession = [];
    
%% session and unit loop
    for ses = 1:length(sessionList)
        session = sessionList{ses};
        unit = unitList{ses};
        fprintf([session unit '\n']);
    % paths
        pd = parseSessionString_df(session, root, sep);
        neuralynxDataPath = [pd.sortedFolder 'session' sep session '_sessionData_nL.mat'];
        unitMetDir = [pd.nLynxFolderSession session '_' unit '_met.mat'];
        sortedFolderLocation = [pd.sortedFolder 'session' sep];
    % decide is good behavior
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
    % find position
    populationDepths = [populationDepths depths(ses)];
    if ~strcmp(session,loadedSession) % avoiding recomputing model variables for different unit from same session
        %% behavior preparation 
        % parse behavior
        os = behAnalysisNoPlot_opMD(session);
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
        [t,~,noSession] = getStanModelParams_samps(p.Results.modelName, [path sampFile '.mat'], 2000, 'sessionName', session);
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
        % choice kernel
        choiceKernel = conv(os.allChoices, myKernel);
        choiceKernel = abs(choiceKernel(0.5*length(myKernel)-0.5+1:end-(0.5*length(myKernel)-0.5)))'/sum(myKernel);
        
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
        preITI = os.timeBtwn';
        if contains(p.Results.modelName, '7params_absPePeAN_scale_int_bias_ord')
            aN = t.aN;
            peBar = t.peBar;
            pePe = t.pePe;
            scPe = pe.*(1-peBar);
            tbl = table(outcome, pe, prePe, preRwd, Qsum, Qdiff, choiceConf, biasSide, rightSide, timeInSession, lickLat, hmm, Qchosen, Qunchosen, QchosenUpdate, preITI, dawExp, svsNext, choiceKernel, scPe, aN, peBar, pePe);
        else
            tbl = table(outcome, pe, prePe, preRwd, Qsum, Qdiff, choiceConf, biasSide, rightSide, timeInSession, lickLat, hmm, Qchosen, Qunchosen, QchosenUpdate, preITI, dawExp, svsNext, choiceKernel);
        end
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
    %% 
    % load metric
    baselineFreq = [baselineFreq met.baseline];
    firstSpikeFreq = [firstSpikeFreq 1000*mean(met.spikeNum(:,1))/20];
    firstSpikeLat = [firstSpikeLat met.spikeLat(1)/1000];
    width = [width metSess.width];
    spikeLatChange = [spikeLatChange (met.spikeLat(10)-met.spikeLat(1))/1000];
    spikeFreqChange = [spikeFreqChange 1000*(mean(met.spikeNum(:,10))-mean(met.spikeNum(:,1)))/met.pulseWidth];
%             spikeFreqChange = [spikeFreqChange mean(met.spikeNum(:,10))-mean(met.spikeNum(:,1))];
    [peak,peakCh] = max(max(met.optoWaveform));
    [~,peakSamp] = max(met.optoWaveform(:,peakCh));
    tempWaveform = [zeros(1, 16-peakSamp), met.optoWaveform(max(1,peakSamp-15):min(peakSamp+25,size(met.optoWaveform,1)),peakCh)', zeros(1, peakSamp+25 - size(met.optoWaveform,1))];
    waveforms = [waveforms; tempWaveform/peak]; 
    [peak,~] = max(max(metSess.waveform));
    [~,peakSamp] = max(metSess.waveform(:,peakCh));
    tempWF = metSess.waveform(:,peakCh);
    tempWF = tempWF(peakSamp-50:peakSamp+50);
    waveformsSession = [waveformsSession; tempWF'/peak];
    %% create spike and lick cell
    spikeFields = fields(sessionData);
    clust = find(contains(spikeFields,unit));
    allTrial_spike_choice = {};
    allTrial_lick = {};

    for k = 1:length(os.responseInds)
            if os.responseInds(k) == 1
                prevTrial_spike = [];
            else
                prevTrial_spikeInd = [sessionData(os.responseInds(k)-1).(spikeFields{clust})] > (sessionData(os.responseInds(k)).respondTime-p.Results.tb*1000);
                prevTrial_spike = sessionData(os.responseInds(k)-1).(spikeFields{clust})(prevTrial_spikeInd) - sessionData(os.responseInds(k)).respondTime;
            end

            currTrial_spikeInd = sessionData(os.responseInds(k)).(spikeFields{clust}) < sessionData(os.responseInds(k)).respondTime+p.Results.tf*1000 ... 
                & sessionData(os.responseInds(k)).(spikeFields{clust}) > sessionData(os.responseInds(k)).respondTime-p.Results.tb*1000;
            currTrial_spike = sessionData(os.responseInds(k)).(spikeFields{clust})(currTrial_spikeInd) - sessionData(os.responseInds(k)).respondTime;

            allTrial_spike_choice{k} = [prevTrial_spike currTrial_spike];

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

    % spike matrix
    allTrial_spikeMatx_choice = zeros(length(os.responseInds),length(time));
    trialDurDiff = zeros(1, length(os.responseInds));
    for j = 1:length(os.responseInds)
        trialDurDiff(j) = (sessionData(os.responseInds(j)).trialEnd - sessionData(os.responseInds(j)).CSon)- p.Results.tf*1000;
    end
   for j = 1:length(allTrial_spike_choice)
        tempSpike = allTrial_spike_choice{j};
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
    
    allTrial_spikeMatx_slide = zscore(allTrial_spikeMatx_slide, 0, 'all');
    % focus window
    spikeCounts = nansum(allTrial_spikeMatx_choice(:,p.Results.tb*1000+p.Results.focusWin(1):p.Results.tb*1000+p.Results.focusWin(2)),2);
    % baseline 
    baseline = nansum(allTrial_spikeMatx_choice(:,1:2000),2);
    % go-cue
    goSpikeCounts = nansum(allTrial_spikeMatx_choice(:,p.Results.tb*1000-500:p.Results.tb*1000),2);
    % excitability
    excit = [excit; [mean(baseline)/2, mean(goSpikeCounts)/0.5]];   
    % calculate and combine rwdvsNrwd
    rwdvsNrwd = [rwdvsNrwd; [mean(spikeCounts(os.rwd_Inds)), mean(spikeCounts(os.nrwd_Inds))]];
    % calculate reward roc curve
    glm = fitglm(spikeCounts, outcome, 'Distribution','binomial','Link','logit');
    score = glm.Fitted.Probability;
    [Xglm,Yglm,~,AUCglm] = perfcurve(outcome,score,1, 'NBoot',500, 'Alpha', 0.05);
    if (AUCglm(2)-0.5)*(AUCglm(3)-0.5)>0
        AUCglmSig = 1;
    else
        AUCglmSig = 0;
    end
    if glm.Coefficients.Estimate(2) < 0
        AUCglm = 1-AUCglm;
    end

    rwdAUC = [rwdAUC; AUCglm(1)];
    rwdAUCSig = [rwdAUCSig; AUCglmSig];
    % zscore focus
    spikeCountsZS = zscore(nansum(allTrial_spikeMatx_choice(:,p.Results.tb*1000+p.Results.focusWin(1):p.Results.tb*1000+p.Results.focusWin(2)),2));
    
    % count session number and unit numbers
    allSessions = [allSessions, session];
    allUnits = [allUnits, unit];



    %% calculate and combine pe
    % bin spikes into pe bins and calcualted bin 
    numBins = 10;
    currTarget = tbl.pe;
    edges = linspace(min(currTarget)-0.001, max(currTarget)+0.001, numBins+1);
    spikeMeans = zeros(1,numBins);
    targetMeans = zeros(1,numBins);
    for k = 1:numBins
        spikeNumsTemp = spikeCountsZS(currTarget >= edges(k) & currTarget < edges(k+1));
        targetMeans(k) = mean(currTarget(currTarget >= edges(k) & currTarget < edges(k+1)), 'omitnan');     
        spikeMeans(k) = mean(spikeNumsTemp, 'omitnan');
    end
    
    % collecting neuron data based on significance
    allPe = [allPe; targetMeans];
    focusSpikes = [focusSpikes; spikeMeans];
    
    targetVar = tbl.pe;
    if p.Results.sepOutcome
        % bin all trials into 2*3 groups by outcome and plotRegressor
        [~, tarInd_rwd] = sort(targetVar(os.rwd_Inds));
        [~, tarInd_nrwd] = sort(targetVar(os.nrwd_Inds));
        tarInd{1} = os.rwd_Inds(tarInd_rwd(1:floor(1/3*length(tarInd_rwd))));
        tarInd{2} = os.rwd_Inds(tarInd_rwd(floor(1/3*length(tarInd_rwd))+1:floor(2/3*length(tarInd_rwd))));
        tarInd{3} = os.rwd_Inds(tarInd_rwd(floor(2/3*length(tarInd_rwd))+1:end));
        tarInd{4} = os.nrwd_Inds(tarInd_nrwd(1:floor(1/3*length(tarInd_nrwd))));
        tarInd{5} = os.nrwd_Inds(tarInd_nrwd(floor(1/3*length(tarInd_nrwd))+1:floor(2/3*length(tarInd_nrwd))));
        tarInd{6} = os.nrwd_Inds(tarInd_nrwd(floor(2/3*length(tarInd_nrwd))+1:end)); 
    else
        % bin all trials into 3 groups by plotRegressor
        [~, tarInd_all] = sort(targetVar);
        tarInd{1} = tarInd_all(1:floor(1/3*length(targetVar)));
        tarInd{2} = tarInd_all(floor(1/3*length(targetVar))+1:floor(2/3*length(targetVar)));
        tarInd{3} = tarInd_all(floor(2/3*length(targetVar))+1:end);
    end
    
    % calculate mean and put into matrix
    spikes_target = zeros(length(tarInd), length(midPoints));
    for k = 1:length(tarInd)
        spikes_target(k,:) = mean(allTrial_spikeMatx_slide(tarInd{k},:), 1);
    end
    
    allSpikes = cat(3, allSpikes, spikes_target);   
    %% calculate and combine SW
    % bin trials into spike number bins and calcualted mean of bins
    numBins = 5;
    currTarget = spikeCountsZS;
    subject = tbl.svsNext;
    subject2 = tbl.choiceKernel;
    edges = linspace(min(currTarget)-0.001, max(currTarget)+0.001, numBins+1);
    subjectMeans = zeros(1,numBins);
    subjectMeans2 = zeros(1,numBins);
    targetMeans = zeros(1,numBins);
    for k = 1:numBins
        subjectTemp = subject(currTarget >= edges(k) & currTarget < edges(k+1));
        subjectTemp2 = subject2(currTarget >= edges(k) & currTarget < edges(k+1));
        targetMeans(k) = mean(currTarget(currTarget >= edges(k) & currTarget < edges(k+1)), 'omitnan');     
        subjectMeans(k) = mean(subjectTemp, 'omitnan');
        subjectMeans2(k) = mean(subjectTemp2, 'omitnan');
    end
    
    % collecting neuron data based on significance
    subjectMeans(~isnan(subjectMeans)) = zscore(subjectMeans(~isnan(subjectMeans)));
    allSW = [allSW; subjectMeans];
    subjectMeans2(~isnan(subjectMeans2)) = zscore(subjectMeans2(~isnan(subjectMeans2)));
    allKernel = [allKernel; subjectMeans2];
    focusSpikesBins = [focusSpikesBins; targetMeans];
       
    end 
end
%% classification
numCat = p.Results.numCat;
[coeff,score,latent, ~, explained, mu] = pca(zscore(waveformsSession, [], 1));
excitRatio = (excit(:,2)-excit(:,1))./excit(:,1);
matrix = [zscore(excit(:,1)) zscore(excit(:,2)) width' excitRatio score(:,1) score(:,2)];
names = {'baselineFreq', 'respFreq', 'width', 'excitRatio', 'PC1', 'PC2'};
matrix = zscore(matrix,0,1);
[coeffAll,scoreAll,latentAll, ~, explainedAll, muAll] = pca(matrix);
indAll = {};
dis = {};
for a = 1:10
    [indAll{a}, ~, dis{a}] = kmeans([score(:, 1:5), width'], numCat);
end
[~,optiInds] = min(cellfun(@mean, dis));
ind = indAll{optiInds};
%% plot everything
figure2;
colors = cool(numCat);
subplot(5,  size(matrix,2)+2, 1); hold on;
scatter3(scoreAll(:,1), scoreAll(:,2), scoreAll(:,3), 15, [0.7 0.7 0.7], 'filled');
for i = 1:numCat
    scatter3(scoreAll(ind==i,1), scoreAll(ind==i,2), scoreAll(ind==i,3), 12, colors(i,:), 'filled');
end
subplot(5,size(matrix,2)+2, 2); hold on;
for i = 1:numCat
    mat = waveformsSession(ind==i,:);
    plotFilled(1:size(mat,2), mat, colors(i,:));
end
plot([10 26], [-0.4 -0.4], 'k', 'LineWidth', 2)
text(12, -0.5, '0.5 ms', 'FontSize', 15)
set(gca, 'XColor', 'none');
set(gca, 'YColor', 'none');
for i = 1:size(matrix,2)
    subplot(5, size(matrix,2)+2, i+2);
    hold on
    edges = linspace(min(matrix(:,i))-0.001, max(matrix(:,i)+0.001),10);
    for j = 1:numCat
        histogram(matrix(ind==j,i), edges, 'FaceColor', colors(j,:), 'Normalization', 'probability');
    end
    title(names{i})
end
%% sp-pe
% subplot(5, 4, [5, 9]); hold on;
% for i = 1:numCat
%     currSpikes = focusSpikes(ind==i,:);
%     currTarget = allPe(ind==i, :);
%     spikeMeans = mean(currSpikes, 'omitnan');
%     spikeSems = sem(currSpikes);
%     targetMeans = mean(currTarget, 'omitnan');
% 
%     plot(targetMeans, spikeMeans, 'color', colors(i,:), 'linewidth', 2, 'Marker', 'none', 'LineStyle', '-');
%     fill([targetMeans fliplr(targetMeans)], [spikeMeans+spikeSems fliplr(spikeMeans-spikeSems)], colors(i,:), 'facealpha', 0.25, 'edgecolor', 'none')
% end
% xlabel('pe');
% ylabel('zscored(spikes/s)')


% colorsG = zeros(6,3);
% for c = 1:3
%     colorsG(c,:) = [1 1.2-c*2/5 1.2-c*2/5];
%     colorsG(c+3,:) = [(c-1)*2/5 (c-1)*2/5 1];
% end

% subplot(5, 4, [6]); hold on;
% currAllSpikes = allSpikes(:,:,ind==1);
% for k = 1:size(currAllSpikes,1)
%     spikesTemp = squeeze(currAllSpikes(k,:,:));% take all spikes from all neurons from first percental
%     plotFilled(slideTime, spikesTemp', colorsG(k,:));
% end
% title('1')
% subplot(5, 4, [10]); hold on;
% currAllSpikes = allSpikes(:,:,ind==2);
% for k = 1:size(currAllSpikes,1)
%     spikesTemp = squeeze(currAllSpikes(k,:,:));% take all spikes from all neurons from first percental
%     plotFilled(slideTime, spikesTemp', colorsG(k,:));
% end
% xlabel('time from lick');
% title('2')
%% SW
subplot(5, 4, [7, 11]); hold on;
for i = 1:numCat
    currSpikes = focusSpikesBins(ind==i,:);
    currTarget = allSW(ind==i, :);
    spikeMeans = mean(currSpikes, 'omitnan');
    targetSems = sem(currSpikes);
    targetMeans = mean(currTarget, 'omitnan');

    plot(spikeMeans, targetMeans, 'color', colors(i,:), 'linewidth', 2, 'Marker', 'none', 'LineStyle', '-');
    fill([spikeMeans fliplr(spikeMeans)], [targetMeans+targetSems fliplr(targetMeans-targetSems)], colors(i,:), 'facealpha', 0.25, 'edgecolor', 'none')
end
ylabel('p(sw)');
xlabel('zscored(spikes/s)')
% switchy
subplot(5, 4, [8, 12]); hold on;
for i = 1:numCat
    currSpikes = focusSpikesBins(ind==i,:);
    currTarget = allKernel(ind==i, :);
    spikeMeans = mean(currSpikes, 'omitnan');
    targetSems = sem(currSpikes);
    targetMeans = mean(currTarget, 'omitnan');

    plot(spikeMeans, targetMeans, 'color', colors(i,:), 'linewidth', 2, 'Marker', 'none', 'LineStyle', '-');
    fill([spikeMeans fliplr(spikeMeans)], [targetMeans+targetSems fliplr(targetMeans-targetSems)], colors(i,:), 'facealpha', 0.25, 'edgecolor', 'none')
end
ylabel('choiceEnt');
xlabel('zscored(spikes/s)')
%% postition
subplot(5, 4, [13, 17]); hold on;
edges = linspace(min(populationDepths)-0.000001, max(populationDepths)+0.000001, 15);
for i = 1:numCat
    histogram(populationDepths(ind==i), edges, 'FaceColor', colors(i,:), 'Normalization', 'probability');
end
xlabel('depth')

subplot(5, 4, [14, 18]); hold on;
edges = linspace(min(populationDepths)-0.000001, max(populationDepths)+0.000001, 15);
ratio = NaN(numCat, length(edges)-1);
for i = 1:numCat
    tempDepths = populationDepths(ind==i);
    for j = 1:length(edges)-1
        ratio(i, j) = sum(tempDepths>=edges(j)&tempDepths<edges(j+1))/sum(populationDepths>=edges(j)&populationDepths<edges(j+1));
    end
    plot(1:length(edges)-1, ratio(i, :), 'color', colors(i,:), 'linewidth', 2, 'Marker', 'none', 'LineStyle', '-');
end
xlabel('depth')
ylabel('ratio')
 
%% rwd vs nrwd
subplot(5, 4, [5, 9]); hold on;
top = max(rwdvsNrwd,[], 'all');
bottom = min(rwdvsNrwd,[], 'all');
scatter(rwdvsNrwd(:,1), rwdvsNrwd(:,2), 18, [0.7 0.7 0.7], 'filled');
for i = 1:numCat
    scatter(rwdvsNrwd(ind==i,1), rwdvsNrwd(ind==i,2), 15, colors(i,:), 'filled');
end
plot([-5 5], [-5 5], 'LineStyle', '--', 'Color', [0.7 0.7 0.7], 'LineWidth', 2);

xlim([bottom top])
ylim([bottom top])
xlabel('rwd')
ylabel('nRwd')

subplot(5, 4, [6]); hold on;
edges = 0:0.05:1;
for i = 1:numCat
    histogram(rwdAUC(ind==i), edges, 'FaceColor', colors(i,:), 'FaceAlpha', 0.4, 'Normalization', 'probability');
end
plot([0.5 0.5], [0 0.3], 'LineStyle', '--', 'Color', [0.2 0.2 0.2], 'LineWidth', 2);

xlabel('rwdAUC')
subplot(5, 4, [10]); hold on;
edges = 0:0.05:1;
% histogram(rwdAUC, edges, 'FaceColor', [0.7 0.7 0.7], 'Normalization', 'probability');
for i = 1:numCat
    histogram(rwdAUC(ind==i&rwdAUCSig==1), edges, 'FaceColor', colors(i,:), 'FaceAlpha', 0.4, 'Normalization', 'probability');
end
plot([0.5 0.5], [0 0.3], 'LineStyle', '--', 'Color', [0.2 0.2 0.2], 'LineWidth', 2);

xlabel('rwdAUC-sig')

%% evoke vs baseline
subplot(5, 4, [15, 19]); hold on;
scatter(excit(:,1), excit(:,2), 18, [0.7 0.7 0.7], 'filled');
for i = 1:numCat
    scatter(excit(ind==i,1), excit(ind==i,2), 15, colors(i,:), 'filled');
end

xlabel('baseline')
ylabel('evoke')

subplot(5, 4, [16, 20]); hold on;
edges = linspace(min((excit(:,2)-excit(:,1))./excit(:,1))-0.00001, max((excit(:,2)-excit(:,1))./excit(:,1))+0.00001, 15);
for i = 1:numCat
    histogram((excit(ind==i,2)-excit(ind==i,1))./excit(ind==i,1), edges, 'FaceColor', colors(i,:), 'FaceAlpha', 0.4, 'Normalization', 'probability');
end
plot([0 0], [0 0.3], 'LineStyle', '--', 'Color', [0.2 0.2 0.2], 'LineWidth', 2);

xlabel('evoke ratio')


screen = get(0,'Screensize');
screen(4) = screen(4) - 100;
set(gcf, 'Position', screen)














