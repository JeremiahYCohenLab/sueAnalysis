    function spikeGNG_populationCue(animalNames, category, varargin)
%task and model parameters
p = inputParser;
% default parameters if none given
p.addParameter('cellName', ['all']);
p.addParameter('plotFlag', 1);
p.addParameter('maxTrial', 1000);
p.addParameter('regressors', '1+pe+biasSide+pe*biasSide')
p.addParameter('binSize', 2000)% in ms
p.addParameter('stepSize', 500)
p.addParameter('tb', 4.5)% in s
p.addParameter('tf', 2.5)% in s
p.addParameter('saveFigFlag', 1);
p.parse(varargin{:});
populationSig = []; % the matrix with 1 for positive beta, -1 for negative beta
populationTStats = []; % t statistics for each parameter
populationCoeffs =[]; % coeffs for each regressor

maxTrial = p.Results.maxTrial;
% basic info
[root, sep] = currComputer();
% time window
time = -1000*p.Results.tb:1000*p.Results.tf;
midPoints = (0.5*p.Results.binSize + 1):p.Results.stepSize:(length(time)-0.5*p.Results.binSize);
slideTime = midPoints - p.Results.tb*1000;
prevSession = [];
%% animal loop 
for ani = 1:length(animalNames)
    % load model fitting results
    animalName = animalNames{ani};
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
        %% behavior
        % correct (hit CR) and wrong (miss, FA)
        correct = zeros(length(os.behSessionData),1);
        correct(setxor(find(os.CSplus>0), find(~isnan(os.lickSide)>0))) = 1;
       
        
        % respond
        respond = zeros(size(correct));
        respond(~isnan(os.lickSide)) = 1;
        
        % CS
        CS = zeros(size(correct));
        CS(os.CSplus) = 1;
        
        % chosen side
        rightSide = NaN(length(os.behSessionData),1);
        rightSide(os.lickSide==1) = 1;
        rightSide(os.lickSide==-1) = -1;
        
        tbl = table(respond, correct, rightSide, CS);
        names = tbl.Properties.VariableNames;
        
        prevSession = session;
    end


    %% create spike and lick cell
    spikeFields = fields(sessionData);
    clust = find(contains(spikeFields,unit));
    allTrial_spike_cue = {};
    for k = 1:length(os.behSessionData)
            if k == 1
                prevTrial_spike = [];
            else
                prevTrial_spikeInd = [sessionData(k-1).(spikeFields{clust})] > (sessionData(k).CSon-p.Results.tb*1000);
%                 if contains(session,'mZS061d20210326')
%                     fprintf([num2str(k),'\n']);
%                 end
                prevTrial_spike = sessionData(k-1).(spikeFields{clust})(prevTrial_spikeInd) - sessionData(k).CSon;
            end

            currTrial_spikeInd = sessionData(k).(spikeFields{clust}) < sessionData(k).CSon+p.Results.tf*1000 ... 
                & sessionData(k).(spikeFields{clust}) > sessionData(k).CSon-p.Results.tb*1000;
            currTrial_spike = sessionData(k).(spikeFields{clust})(currTrial_spikeInd) - sessionData(k).CSon;

            allTrial_spike_cue{k} = [prevTrial_spike currTrial_spike];

    end

    % sometimes no licks/spikes are considered 1x0 and sometimes they are []
    % plotSpikeRaster does not place nicely with [] so this converts all empty indices to 1x0
    allTrial_spike_cue(cellfun(@isempty,allTrial_spike_cue)) = {zeros(1,0)}; 
    %% initialize lick matrices
    for j = 1:length(os.behSessionData)
        trialDurDiff(j) = (sessionData(j).trialEnd - sessionData(j).CSon)- p.Results.tf*1000;
    end
    trialDurDiff(end) = 0; 

    % spike matric for GLM
    allTrial_spikeMatx_cue = zeros(length(os.behSessionData),length(time));         
   for j = 1:length(allTrial_spike_cue)
        tempSpike = allTrial_spike_cue{j};
        tempSpike = tempSpike + p.Results.tb*1000; % add this to pad time for SDF
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
    allTrial_spikeMatx_slide = zeros(length(os.behSessionData), length(midPoints));
    for w = 1:length(midPoints)
        allTrial_spikeMatx_slide(:,w) = ...
            nansum(allTrial_spikeMatx_cue(:,midPoints(w)-0.5*p.Results.binSize:midPoints(w)+0.5*p.Results.binSize-1),2)*1000/p.Results.binSize;
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

    populationSig = cat(3, populationSig, sigs);
    populationTStats  = cat(3,populationTStats, tStats);
    populationCoeffs  = cat(3,populationCoeffs, coeffs);
    end 
end

%% plot everything
regressors = lm.CoefficientNames(2:end);
tFig = figure;
screen = get(0,'Screensize');
screen(4) = screen(4) - 100;
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
xlabel('time from cue')


minT = min(populationTStats,[],'all'); 
maxT = max(populationTStats,[],'all');
bins = linspace(minT,maxT,30);
for k = 1:length(midPoints)
    for j = 1:length(regressors)
        subplot(length(regressors)+1, length(midPoints), length(midPoints)*j+k); hold on;
        tmpTStats = populationTStats(k,j,:);
        tmpSig = populationSig(k,j,:);
        nonSig = tmpTStats(tmpSig == 0);
        Sig = tmpTStats(tmpSig ~= 0);
        histogram(nonSig, bins, 'FaceColor', [0.5 0.5 0.5]);
        histogram(Sig, bins, 'FaceColor', colors(j,:));
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












