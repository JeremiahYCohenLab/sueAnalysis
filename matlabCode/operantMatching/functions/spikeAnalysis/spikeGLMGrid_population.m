    function spikeGLMGrid_population(animalNames, category, varargin)
%task and model parameters
p = inputParser;
% default parameters if none given
p.addParameter('cellName', ['all']);
p.addParameter('plotFlag', 1);
p.addParameter('maxTrial', 1000);
% p.addParameter('modelName','7params_absPePeAN_scale_int_bias_ord')
p.addParameter('modelName','5params')
p.addParameter('regressors', '1+outcome+pe*rightSide')
p.addParameter('xRegressor', 'outcome')
p.addParameter('yRegressor', 'pe')
p.addParameter('plotRegressor', 'pe')
p.addParameter('binSize', 100)% in ms
p.addParameter('stepSize', 100)
p.addParameter('tb', 1)% in s
p.addParameter('tf', 3)% in s
p.addParameter('focusWin', [300 1300]); % in ms from respond;l;
p.addParameter('sepOutcome', 1)
p.addParameter('saveFigFlag', 1);
p.parse(varargin{:});
populationSig = []; % the matrix with 1 for positive beta, -1 for negative beta
populationTStats = []; % t statistics for each parameter
populationCoeffs =[]; % coeffs for each regressor
paramNames = getParamNames_dF(p.Results.modelName, 1);
maxTrial = p.Results.maxTrial;
% basic info
[root, sep] = currComputer();
% time window
time = -1000*p.Results.tb:1000*p.Results.tf;
midPoints = (0.5*p.Results.binSize + 1):p.Results.stepSize:(length(time)-0.5*p.Results.binSize);
slideTime = midPoints - p.Results.tb*1000;
% zscored spikeCount after rwd
focusSpikes = cell(3,3); % zscored within Session; vector 1*allTrials; 
allPe = cell(3,3); % zscored within session; vector 1*allTrials;
% zscored spikeCount with 100ms sliding window 50ms time step  
allSpikes = cell(3,3); % three curves for each neuron (curves are zscored before squeezed into 3) 3 * bins * neurons
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

    drift = unitsInfo(2:end,10);
    drift = contains(drift, 'drift')|contains(drift, 'Drift');
    quality = nums(:,1);
    sessionList = sessionList(quality<=0.05 & ~drift);
    unitList = unitList(quality<=0.05 & ~drift);
    optoUnitList = optoUnitList(quality<=0.05 & ~drift);
    subFolders = subFolders(quality<=0.05 & ~drift);
    
    % load good days
    [~ , goodDayList, ~] = xlsread([root xlFile], animalName);
    [~,col] = find(~cellfun(@isempty,strfind(goodDayList, category)) == 1);
    goodDayList = goodDayList(2:end,col);
    prevSession = [];
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
        if respLat > 15000 || met.isiV > 0.001 || met.distance>0.3 
            continue
        end
    end
    % load behavior and neurons
    if exist(neuralynxDataPath,'file')
        load(neuralynxDataPath)
    else
        sessionData = generateSessionData_nL_operantMatching(session);
    end
    if ~strcmp(session,prevSession) % avoiding recomputing model variables for different unit from same session
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
        [t,~,noSession] = getStanModelParams_samps(p.Results.modelName, [path sampFile '.mat'], 2000, 'sessionName', session);
        if noSession
            fprintf(['no good behavior in ' session '\n']);
            continue
        end
        
        % diff value
        Qdiff = abs(t.Q(:,2)-t.Q(:,1));
        % total value
        Qsum = sum(t.Q,2);
        % prepe
        prePe = [NaN; t.pe(1:end-1)];
        % pe
        pe = zscore(t.pe);
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
        preITI = os.timeBtwn';
        if contains(p.Results.modelName, '7params_absPePeAN_scale_int_bias_ord')
            aN = t.aN;
            peBar = t.peBar;
            pePe = t.pePe;
            scPe = pe.*(1-peBar);
            tbl = table(outcome, pe, prePe, preRwd, Qsum, Qdiff, choiceConf, biasSide, rightSide, timeInSession, lickLat, hmm, Qchosen, Qunchosen, QchosenUpdate, preITI, dawExp, scPe, aN, peBar, pePe);
        else
            tbl = table(outcome, pe, prePe, preRwd, Qsum, Qdiff, choiceConf, biasSide, rightSide, timeInSession, lickLat, hmm, Qchosen, Qunchosen, QchosenUpdate, preITI, dawExp);
        end
        names = tbl.Properties.VariableNames;
        % zscore all regressors
        for cols = 1:length(names)
            tmp = tbl.(names{cols});
            tmp(~isnan(tmp)) = zscore(tmp(~isnan(tmp)));
            tbl.(names{cols}) = tmp;
        end
        
        prevSession = session;
    end


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
    allTrial_lickMatx_slide = zeros(length(os.responseInds), length(midPoints));
    for w = 1:length(midPoints)
%         fprintf([num2str(w), '\n'])
        allTrial_lickMatx_slide(:,w) = ...
            sum(allTrial_lickMatx(:,midPoints(w)-0.5*p.Results.binSize:midPoints(w)+0.5*p.Results.binSize-1),2)*1000/p.Results.binSize;   
    end
    % spike matric for GLM
    allTrial_spikeMatx_choice = zeros(length(os.responseInds),length(time));         
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
    spikeCounts = zscore(nansum(allTrial_spikeMatx_choice(:,p.Results.tb*1000+p.Results.focusWin(1):p.Results.tb*1000+p.Results.focusWin(2)),2));

    % linear regression
    currTbl = addvars(tbl, spikeCounts);
    lm = fitlm(currTbl, ['spikeCounts~' p.Results.regressors]);
    sigs = zeros(1,length(lm.CoefficientNames)-1);
    for j = 1:(length(lm.CoefficientNames)-1)
        if lm.Coefficients.pValue(j+1)<0.05
           sigs(j) = sign(lm.Coefficients.Estimate(j+1));
        else
           sigs(j) = 0;
        end
    end
    
    xRegressorInd = find(strcmp(lm.CoefficientNames, p.Results.xRegressor))-1;
    yRegressorInd = find(strcmp(lm.CoefficientNames, p.Results.yRegressor))-1;
    
    % collecting neuron data based on significance
    allPe{2+sigs(yRegressorInd),2+sigs(xRegressorInd)} = [allPe{2+sigs(yRegressorInd),2+sigs(xRegressorInd)}; tbl.(p.Results.plotRegressor)];
    focusSpikes{2+sigs(yRegressorInd),2+sigs(xRegressorInd)} = [focusSpikes{2+sigs(yRegressorInd),2+sigs(xRegressorInd)}; spikeCounts];
    
    targetVar = tbl.(p.Results.plotRegressor);
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
    
    allSpikes{2+sigs(yRegressorInd),2+sigs(xRegressorInd)} = cat(3, allSpikes{2+sigs(yRegressorInd),2+sigs(xRegressorInd)}, spikes_target);
    end 
end

%% plot everything
figure;
numBins = 20;
colors = zeros(size(spikes_target,1),3);
if p.Results.sepOutcome
    for c = 1:0.5*size(spikes_target,1)
        colors(c,:) = [1 1.2-c*2/5 1.2-c*2/5];
        colors(c+3,:) = [(c-1)*2/5 (c-1)*2/5 1];
    end
else
    for c = 1:size(spikes_target,1)
        colors(c,:) = [1 1.2-c*0.4 1.2-c*0.4];
    end
end
for i = 1:3
    for j = 1:3
        currTarget = allPe{i,j};
        currSpikes = focusSpikes{i,j};
        currAllSpikes = allSpikes{i,j};
        edges = binEqualSize(currTarget, numBins);
        spikeMeans = zeros(numBins,1);
        spikeSems = zeros(numBins,1);
        targetMeans = zeros(numBins,1);
        for k = 1:numBins
            if k < numBins
                spikeNumsTemp = currSpikes(currTarget >= edges(k) & currTarget < edges(k+1));
                targetMeans(k) = mean(currTarget(currTarget >= edges(k) & currTarget < edges(k+1)));
            else
                spikeNumsTemp = currSpikes(currTarget >= edges(k) & currTarget <= edges(k+1));
                targetMeans(k) = mean(currTarget(currTarget >= edges(k) & currTarget <= edges(k+1)));
            end      
            spikeMeans(k) = mean(spikeNumsTemp);
            spikeSems(k) = sem(spikeNumsTemp);
        end
        subplot(3,6,(i-1)*6+2*j-1); hold on;
        plot(targetMeans, spikeMeans, 'color', [0 0 1], 'linewidth', 1, 'Marker', 'none', 'LineStyle', '-');
        fill([targetMeans' fliplr(targetMeans')], [spikeMeans'+spikeSems' fliplr(spikeMeans'-spikeSems')], [0 0 1], 'facealpha', 0.25, 'edgecolor', 'none')
        title(num2str(size(currAllSpikes,3)));
        subplot(3,6,(i-1)*6+2*j); hold on;
        if size(currAllSpikes,3)>1
            for k = 1:size(currAllSpikes,1)
                spikesTemp = squeeze(currAllSpikes(k,:,:));% take all spikes from all neurons from first percental
                plotFilled(slideTime, spikesTemp', colors(k,:));
            end
        else
            for k = 1:size(currAllSpikes,1)
                spikesTemp = squeeze(currAllSpikes(k,:,:));% take all spikes from all neurons from first percental
                plot(slideTime, spikesTemp, 'Color', colors(k,:));
            end
        end
        
        if i == 1 && j == 2
            ylabel('p.Results.yRegressor')
        end
        
        if i == 2 && j == 1
            xlabel('p.Results.xRegressor')
        end
    end
end














