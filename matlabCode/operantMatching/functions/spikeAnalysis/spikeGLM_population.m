    function spikeGLM_population(animalNames, category, varargin)
%task and model parameters
p = inputParser;
% default parameters if none given
p.addParameter('cellName', ['all']);
p.addParameter('plotFlag', 1);
p.addParameter('maxTrial', 1000);
% p.addParameter('modelName','7params_absPePeAN_scale_int_bias_ord')
p.addParameter('modelName','5params')
p.addParameter('regressors', '1 + outcome + Qchosen + rightSide')
p.addParameter('binSize', 1000)% in ms
p.addParameter('stepSize', 250)
p.addParameter('tb', 1.7)% in s 1.7
p.addParameter('tf', 3)% in s
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
allSessions = {};
allUnits = {};
allPe = {};
allChoices = {};
allSpikes = {};
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
    drift = contains(drift, 'drift')|contains(drift, 'Drift')|contains(drift, 'duplicate')|contains(drift, 'Duplicate');
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
        neuralynxDataPath = [pd.sortedFolder session '_sessionData_nL.mat'];
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
        if respLat > 15000 || met.isiV > 0.001 || met.distance>0.3 || isnan(met.distance) || isempty(met.distance)
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
        if length(os.behSessionData)~=length(sessionData)
            fprintf([session ' error \n'])
            return
        end
%         lickInds = os.lickInds
        choice = os.allChoices';
        choice(choice<0) = 0;
        outcome = abs(os.allRewards)';
        choice = choice(1:min(length(choice), maxTrial));
        outcome = outcome(1:length(choice));
        outcomeL = outcome;
        outcomeL(outcomeL==0) = -1;
        outcomeR = outcome;
        outcomeR(outcomeR==0) = -1;
        outcomeL(choice==1) = 0;
        outcomeR(choice==0) = 0;
        responseInds = os.responseInds(1:min(length(choice), maxTrial)); 
        preRwd = [NaN abs(os.allRewards(1:end-1))]';
        %% behavior
        % switch
        svs = zeros(length(os.responseInds),1);
        svs(os.changeChoice_Inds) = 1;
        svsNext = [svs(2:end); NaN];
        svsWhenNrwd = svsNext;
        svsWhenNrwd(os.rwd_Inds) = NaN;
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
        rightSide(os.allChoices<=0)=-1;
        preITI = os.timeBtwn';
        % consecutive no rewards
        conNrwds = zeros(size(pe));
        for j = 1:length(choice)
            if outcome(j) == 0
                k = 1;
                while j-k>0 
                    if outcome(j-k)==0
                        k = k+1;
                    else
                        break
                    end
                end                
                conNrwds(j) = k;
            end
        end
        
        if contains(p.Results.modelName, '7params_absPePeAN_scale_int_bias_ord')
            aN = t.aN;
            peBar = t.peBar;
            pePe = t.pePe;
            scPe = pe.*(1-peBar);
            tbl = table(outcome, outcomeL, outcomeR, pe, prePe, preRwd, Qsum, Qdiff, choiceConf, biasSide, rightSide, timeInSession, lickLat, hmm, Qchosen, Qunchosen, QchosenUpdate, preITI, dawExp, svs, svsNext, svsWhenNrwd, conNrwds, scPe, aN, peBar, pePe);
        else
            tbl = table(outcome, outcomeL, outcomeR, pe, prePe, preRwd, Qsum, Qdiff, choiceConf, biasSide, rightSide, timeInSession, lickLat, hmm, Qchosen, Qunchosen, QchosenUpdate, preITI, dawExp, svs, svsNext, svsWhenNrwd, conNrwds);
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
        
        prevSession = session;
    end


    %% create spike and lick cell
    allSessions = [allSessions; session];
    allUnits = [allUnits; unit];
    spikeFields = fields(sessionData);
    clust = find(contains(spikeFields,unit));
    allTrial_spike_choice = {};
    allTrial_lick = {};
    for k = 1:length(os.responseInds)
            if os.responseInds(k) == 1
                prevTrial_spike = [];
            else
                prevTrial_spikeInd = [sessionData(os.responseInds(k)-1).(spikeFields{clust})] > (sessionData(os.responseInds(k)).respondTime-p.Results.tb*1000);
%                 if contains(session,'mZS061d20210326')
%                     fprintf([num2str(k),'\n']);
%                 end
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

    sigs = [];
    tStats = [];
    coeffs = [];
    for k = 1:length(midPoints)
        currSpikes = zscore(allTrial_spikeMatx_slide(:,k));
        currTbl = addvars(tbl, currSpikes);
        lm = fitlm(currTbl, ['currSpikes~' p.Results.regressors]);
        sigTmp = zeros(1,length(lm.CoefficientNames)-1);
        tStatsTmp = zeros(1,length(lm.CoefficientNames)-1);
        coeffsTmp = zeros(1,length(lm.CoefficientNames)-1);
        for j = 1:(length(lm.CoefficientNames)-1)
            if lm.Coefficients.pValue(j+1)<0.05
               sigTmp(j) = sign(lm.Coefficients.Estimate(j+1));
            else
               sigTmp(j) = 0;
            end
            tStatsTmp(j) = lm.Coefficients.tStat(j+1);
            coeffsTmp(j) = lm.Coefficients.Estimate(j+1);
        end
        sigs = [sigs; sigTmp];
        tStats = [tStats; tStatsTmp];
        coeffs = [coeffs; coeffsTmp];
        
    end
    allPe = [allPe; pe];
    allChoices = [allChoices; choice];
    focusInd = find(slideTime>301,1);
    allSpikes = [allSpikes; allTrial_spikeMatx_slide(:,focusInd)];
    populationSig = cat(3, populationSig, sigs);
    populationTStats  = cat(3,populationTStats, tStats);
    populationCoeffs  = cat(3,populationCoeffs, coeffs);
    end 
end

%% plot everything
regressors = lm.CoefficientNames(2:end);
tFig = figure;
screen = get(0,'Screensize');
screen(4) = screen(4) - 100;slideTime
set(tFig, 'Position', screen)
suptitle('tStats distribution')
colors = cool(length(regressors));
subplot(length(regressors)+1,1,1); hold on;
allSig = abs(populationSig);
allSig = mean(allSig,3);
for i = 1:length(regressors)
    plot(slideTime,allSig(:,i),'Color', colors(i,:), 'LineStyle','-', 'Marker','none', 'linewidth', 2);
end
edges = [slideTime - 0.5*p.Results.stepSize slideTime(end)+0.5*p.Results.stepSize];
for i = 1:length(edges)
    line([edges(i) edges(i)], [0 1.2*max(allSig,[],'all')], 'color', [0.7 0.7 0.7], 'LineStyle','--')
end
line([300 300], [0 1.2*max(allSig,[],'all')], 'color', 'r', 'LineStyle','--');
legend(regressors)
ylim([0 1.2*max(allSig,[],'all')])
xlim(minmax(edges));
ylabel('ratio of sig untis')
xlabel('time from respond')


minT = min(populationTStats,[],'all');
maxT = max(populationTStats,[],'all');
bins = linspace(minT,maxT,40);
for k = 1:length(midPoints)
    for j = 1:length(regressors)
        subplot(length(regressors)+1, length(midPoints), length(midPoints)*j+k); hold on;
        tmpTStats = populationTStats(k,j,:);
        tmpSig = populationSig(k,j,:);
        nonSig = tmpTStats(tmpSig == 0);
        Sig = tmpTStats(tmpSig ~= 0);
        histogram(nonSig, bins, 'FaceColor', [0.5 0.5 0.5], 'Normalization', 'probability');
        histogram(Sig, bins, 'FaceColor', colors(j,:),  'Normalization', 'probability');
        if k == 1
            ylabel(regressors{j})
        end
    end
end

% for i = 1:length(midPoints)
%     titleStr = sprintf('From %d To %d', edges(i), edges(i+1)); 
%     figure;
%     scatterAll(squeeze(populationCoeffs(i,:,:))', regressors,7,'m');
%     suptitle(titleStr)
% end  












