function compareLatVar_vSim_dF(xlFile, sheet, beh, varargin)

p = inputParser;
% default parameters if none given
p.addParameter('bernFlag', 1);
p.addParameter('plotFlag', 0);
p.addParameter('modelName', 'sixParam_absPePeAN_bi')
p.addParameter('varNames', [{'peBar'} {'pePe'} {'pe'}]);
p.addParameter('sessionParamFlag', 1);
p.addParameter('rwdProbs', [90 50 10]);
p.addParameter('randomSeed', 657);
p.addParameter('binWin', 10);
p.addParameter('numSims', 500);
p.parse(varargin{:});

[root, sep] = currComputer();

[numbers, sessionCellList, ~] = xlsread(xlFile, sheet);
revForFlag = numbers(:,1);
intanFlag = numbers(:,2);
cellList = sessionCellList(2:end, 1);
sessionList = sessionCellList(2:end, 2);

if size(numbers,2) > 2
    trialList = numbers(:,3:4);
else
    trialList = nan(length(cellList),2);
end 

tB = -1500;
tF = 5000;

intL = 5;
x = [1:intL];
intKern = -0.025*(x-1) + 0.2;

binWin = p.Results.binWin*1000;
numSims = p.Results.numSims;

numVars = length(p.Results.varNames);
sessionName = [];

numCells = length(sessionList);

sig = nan(numCells, numVars);
simT = nan(numCells, numVars, numSims);

rSeed = p.Results.randomSeed;

for currCell = 1:numCells
    
    fprintf('Analyzing cell %d of %d \n', currCell, length(sessionList));
    
    if strcmp(sessionName, sessionList{currCell}) == 0
        sessionName = sessionList{currCell};
        [animalName, date] = strtok(sessionName, 'd'); 
        animalName = animalName(2:end);

        %load spike data
        if isstrprop(sessionName(end), 'alpha')
            sortedFolderLocation = [root animalName sep sessionName(1:end-1) sep 'sorted' sep 'session ' sessionName(end) sep];
        else
            sortedFolderLocation = [root animalName sep sessionName sep 'sorted' sep 'session' sep];
        end
        sortedFolder = dir(sortedFolderLocation);

        if revForFlag(currCell) == 1
            if any(~cellfun(@isempty,strfind({sortedFolder.name},'_intan.mat'))) == 1
                sessionDataInd = ~cellfun(@isempty,strfind({sortedFolder.name},'_intan.mat'));
                load([sortedFolderLocation sortedFolder(sessionDataInd).name])
            else
                [sessionData] = generateSessionData_intan_operantMatching(sessionName);
            end
        else
            if any(~cellfun(@isempty,strfind({sortedFolder.name},'_nL.mat'))) == 1
                sessionDataInd = ~cellfun(@isempty,strfind({sortedFolder.name},'_nL.mat'));
                load([sortedFolderLocation sortedFolder(sessionDataInd).name])
            elseif any(~cellfun(@isempty,strfind({sortedFolder.name},'_intan.mat'))) == 1
                sessionDataInd = ~cellfun(@isempty,strfind({sortedFolder.name},'_intan.mat'));
                load([sortedFolderLocation sortedFolder(sessionDataInd).name])
            else
                if intanFlag(currCell) == 1
                    [sessionData] = generateSessionData_intan_operantMatching(sessionName);
                else
                    [sessionData] = generateSessionData_nL_operantMatching(sessionName);
                end
            end
        end
        [s] = behDataForStruct_opMD(sessionName, 'revForFlag', revForFlag(currCell));
        
        spikeFields = fields(sessionData);
        cellInds = find(~cellfun(@isempty,strfind(spikeFields, 'TT')) | ~cellfun(@isempty,strfind(spikeFields, 'C_')));
        sessionTime = (sessionData(end).CSon + tF) - (sessionData(1).CSon + tB);
        csTimes = [sessionData.CSon] -  (sessionData(1).CSon + tB);
        csTimes = [csTimes(s.responseInds) csTimes(end) + tF];
        
        %get model terms
        if p.Results.bernFlag
            modelPath = [root animalName sep animalName 'sorted' sep 'stan' sep 'bernoulli' sep p.Results.modelName sep animalName...
            beh '_' p.Results.modelName '.mat'];
        else
            modelPath = [root animalName sep animalName 'sorted' sep 'stan' sep p.Results.modelName sep animalName...
            beh '_' p.Results.modelName '.mat'];
        end
        t = generateStanModelTerms_opMD(p.Results.modelName, modelPath, sessionName, p.Results.sessionParamFlag, revForFlag(currCell));
        
        sumQ = [zscore(t.Q(:,1) + t.Q(:,2))]';
        diffQ = [zscore(t.Q(:,1) - t.Q(:,2))]';
        confQ = [zscore(abs(diffQ))]';
        pe = [zscore(abs(t.pe))]';

        allVars = [];
        for currVar = 1:length(p.Results.varNames)
            if isfield(t, p.Results.varNames{currVar})
                allVars(currVar, :) = zscore(eval(['t.' p.Results.varNames{currVar} ';']));
            else
                allVars(currVar, :) = eval([p.Results.varNames{currVar} ';']);
            end
        end

        allVarsSesh = zeros(size(allVars,1), sessionTime);
        for currVar = 1:numVars
            for currTrial = 1:length(csTimes)-1
                allVarsSesh(currVar, csTimes(currTrial):csTimes(currTrial+1)-1) = allVars(currVar, currTrial);
            end
        end
        
    end

    %% Sort all spikes into a raster-able matrix
    cellInd = find(~cellfun(@isempty,strfind(spikeFields, cellList{currCell})));
    cellNum = find(cellInds == cellInd);
    spikeTimes = sessionData(cellNum).allSpikes - (sessionData(1).CSon + tB);
    spikeTimes = spikeTimes(spikeTimes > 0 & spikeTimes <= sessionTime);
    sessionSpikes = zeros(1, sessionTime);
    sessionSpikes(spikeTimes) = 1;
    sessionSpikes = movsum(sessionSpikes, binWin);
    sessionSpikes = downsample(sessionSpikes, binWin, binWin-1);
    
    for currVar = 1:numVars
        sessionVar = movsum(allVarsSesh(currVar,:), binWin);
        sessionVar = downsample(sessionVar, binWin, binWin-1);
        mdl = fitlm(zscore(sessionVar)', zscore(sessionSpikes)'); 
        tStat(currCell, currVar) = mdl.Coefficients.tStat(2);
        
        varACF(currVar,:) = autocorr(sessionVar, 30);
        acfSum(currVar) = sum(varACF(currVar,:));
    end

    for currS = 1:numSims

        rSeed = rSeed + 1;
        switch p.Results.modelName
            case 'sixParam_absPePeAN_bi'
                [out, ~, ~, ~, ~] = qLearningModel_absPePeAN_bi_simNoPlot('params', t.params,...
                    'maxTrials', length(s.responseInds), 'randomSeed', rSeed, 'rwdProbs', p.Results.rwdProbs);
            case 'sevenParam_absPePeAN_bi_k'
                [out, ~, ~, ~, ~] = qLearningModel_absPePeAN_bi_k_simNoPlot('params', t.params,...
                    'maxTrials', length(s.responseInds), 'randomSeed', rSeed, 'rwdProbs', p.Results.rwdProbs);
        end

        s_sumQ = [zscore(out.Q(:,1) + out.Q(:,2))]';
        s_diffQ = [zscore(out.Q(:,1) - out.Q(:,2))]';
        s_confQ = [zscore(abs(s_diffQ))]';
        s_pe = [zscore(abs(out.pe))]';

        s_allVars = [];
        for currVar = 1:length(p.Results.varNames)
            if isfield(out, p.Results.varNames{currVar})
                s_allVars(currVar, :) = zscore(eval(['out.' p.Results.varNames{currVar} ';']));
            else
                s_allVars(currVar, :) = eval(['s_' p.Results.varNames{currVar} ';']);
            end
        end

        s_allVarsSesh = zeros(size(s_allVars,1), sessionTime);
        for currVar = 1:numVars
            for currTrial = 1:length(csTimes)-1
                s_allVarsSesh(currVar, csTimes(currTrial):csTimes(currTrial+1)-1) = s_allVars(currVar, currTrial);
            end
        end
        
        for currVar = 1:numVars
            s_sessionVar = movsum(s_allVarsSesh(currVar,:), binWin);
            s_sessionVar = downsample(s_sessionVar, binWin, binWin-1);

            mdl = fitlm(zscore(s_sessionVar)', zscore(sessionSpikes)');
            simT(currCell, currVar, currS) = mdl.Coefficients.tStat(2);

            simACF = autocorr(s_sessionVar, 30);
            diffTmp(currVar, currS) = sum(abs([varACF(currVar, :) - simACF]));
            sumTmp(currVar, currS) = sum(simACF);
        end
    end

    sig(currCell, currVar) = prctile(abs(simT(currCell, currVar, :)), 95);

    tmp = squeeze([abs(simT(currCell,currVar,:))]);
    mdl = fitlm(sumTmp', tmp);
    acfSumT(currCell, currVar) = mdl.Coefficients.tStat(2);

    mdl = fitlm(diffTmp', tmp);
    acfDiffT(currCell, currVar) = mdl.Coefficients.tStat(2);


    if p.Results.plotFlag == 1
        figure; hold on;
        histogram(abs(simT(currCell, currVar, :)), 20, 'normalization', 'probability', 'facecolor', 'c')
        yl = ylim;
        plot([sig(currCell, currVar) sig(currCell, currVar)], [0 yl(2)], ':k', 'linewidth', 3)
        plot([abs(tStat(currCell, currVar)) abs(tStat(currCell, currVar))], [0 yl(2)], '-m', 'linewidth', 3)
        xlabel('t-stat')
        set(gca, 'tickdir', 'out', 'box', 'off')
        set(gcf, 'renderer', 'painters')
        legend('sim', '5% sig boundary', 'actual')
    end

end

colors = cool(numVars);
figure;
for currVar = 1:numVars
    
    subplot(numVars,3,currVar*3-2)
    histogram(acfSumT(:,currVar), 20, 'normalization', 'probability', 'facecolor', colors(currVar,:))
    xlabel('t-stat')
    ylabel('probability')
    title('relationship between total autocorr and corr')
    set(gca, 'tickdir', 'out', 'box' , 'off')

    subplot(numVars,3,currVar*3-1)
    histogram(acfDiffT(:,currVar), 20, 'normalization', 'probability', 'facecolor', colors(currVar,:))
    xlabel('t-stat')
    ylabel('probability')
    title('relationship between acf similarity and corr')
    set(gca, 'tickdir', 'out', 'box' , 'off')
    set(gcf, 'renderer', 'painters')
    
    subplot(numVars,3,currVar*3)
    sigInds = abs(tStat(:,currVar)) > sig(:,currVar);
    h = pie([sum(sigInds) numCells-sum(sigInds)]);
    TextArr = findobj(h, 'Type','Text');   
    PatchArr = findobj(h, 'Type','patch');      
    PatchArr(2).FaceColor = [1 1 1];  
    PatchArr(1).FaceColor = colors(currVar,:);
    titleTxt = strrep(p.Results.varNames{currVar}, '_', ' ');
    title(titleTxt)
    
end

set(gcf, 'renderer', 'painters', 'position', [-1497 454 979 402])

end