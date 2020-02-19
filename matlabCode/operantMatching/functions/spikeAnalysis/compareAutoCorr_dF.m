function compareAutoCorr_dF(xlFile, sheet, beh, varargin)

p = inputParser;
% default parameters if none given
p.addParameter('bernFlag', 1);
p.addParameter('saveFigFlag', 0);
p.addParameter('modelName', 'sixParam_absPePeAN_bi')
p.addParameter('varNames', [{'peBar'} {'pePe'} {'pe'}]);
p.addParameter('sessionParamFlag', 1);
p.addParameter('binWin', 10);
p.addParameter('numLags', 30);
p.addParameter('numShifts', 60);
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


binWin = p.Results.binWin*1000;
numLags = p.Results.numLags;
numShifts = p.Results.numShifts;

numVars = length(p.Results.varNames);
sessionName = [];

numCells = length(sessionList);

spikesTau = nan(1,numCells);
varTau = nan(numCells, numVars);

rho = nan(numCells, numVars);
pVal = nan(numCells, numVars);
shiftP = nan(numCells, numVars, numShifts);

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
        
        csTimes = round(csTimes(1:end-1)/1000);
        iti = diff(csTimes);
    end

    %% Sort all spikes into a raster-able matrix
    cellInd = find(~cellfun(@isempty,strfind(spikeFields, cellList{currCell})));
    cellNum = find(cellInds == cellInd);
    spikeTimes = sessionData(cellNum).allSpikes - (sessionData(1).CSon + tB);
    spikeTimes = spikeTimes(spikeTimes > 0 & spikeTimes <= sessionTime);
    sessionSpikes = zeros(1, sessionTime);
    sessionSpikes(spikeTimes) = 1;
    sessionSpikesS = sessionSpikes;
    sessionSpikesS = movsum(sessionSpikesS, 1000);
    sessionSpikesS = downsample(sessionSpikesS, 1000, 999);
    sessionSpikes = movsum(sessionSpikes, binWin);
    sessionSpikes = downsample(sessionSpikes, binWin, binWin-1);
    
    spikesACF = autocorr(sessionSpikes, numLags);
    binSteps = [0:numLags] * (binWin)/1000;
    spikesACFfit = singleExpFitInt(spikesACF, binSteps);
    spikesTau(currCell) = 1/spikesACFfit.b;
    
    figure; plot(spikesACF)
    
    shiftSpikes = [sessionSpikes(2:end) sessionSpikes(1)];
    for currS = 1:numShifts-1
        shiftSpikes(currS+1,:) = [shiftSpikes(currS, 2:end) shiftSpikes(currS,1)];
    end
    
    %figure;
    for currVar = 1:numVars
        sessionVar = movsum(allVarsSesh(currVar,:), binWin);
        sessionVar = downsample(sessionVar, binWin, binWin-1);
        mdl = fitlm(zscore(sessionVar)', zscore(sessionSpikes)'); 
        pVal(currCell, currVar) = mdl.Coefficients.pValue(2);
        varACF = autocorr(sessionVar, numLags);
        varACFfit = singleExpFitInt(varACF, binSteps);
        varTau(currCell, currVar) = 1/varACFfit.b;
        
        for currS = 1:numShifts
            mdl = fitlm(zscore(sessionVar)', zscore(shiftSpikes(currS,:))');
            shiftP(currCell, currVar, currS) = mdl.Coefficients.pValue(2);
        end
        
%         sessionVar = movsum(allVarsSesh(currVar,:), 1000);
%         sessionVar = downsample(sessionVar, 1000, 999);
%         mdl = fitlm(zscore(sessionVar)', zscore(sessionSpikesS)');
%         subplot(2,numVars,(currVar-1)*2+1); hold on;
%         plot(mdl.Residuals.Raw, '-k')
%         scatter(csTimes(logical(s.allRewardsBinary)), mdl.Residuals.Raw(csTimes(logical(s.allRewardsBinary))), [], [0 0.7 1], 'filled')
%         scatter(csTimes(logical(~s.allRewardsBinary)), mdl.Residuals.Raw(csTimes(logical(~s.allRewardsBinary))), [], [1 0.3 0.3], 'filled')
%         xlim([0 length(sessionVar)])
%         set(gca, 'tickdir', 'out', 'box', 'off')
%         xlabel('time (s)')
%         ylabel('raw residuals')
        
%         itiRes = [];
%         for currT = 2:length(csTimes)
%             itiRes(currT-1) = mdl.Residuals.Raw(csTimes(currT)-1);
%         end
%         [binInds, binEdges] = discretize(iti, 10);
%         for currB = 1:10
%             tmpS = itiRes(binInds == currB);
%             meanRes(currB) = mean(tmpS);
%             semRes(currB) = std(tmpS) / sqrt(length(tmpS));
%         end
%         subplot(2,numVars,currVar*2)
%         x = binEdges(1:end-1) + diff(binEdges(1:2));
%         errorbar(x, meanRes, semRes, 'color', [1 0.4 0.2], 'linewidth', 2);
%         xlim([0 binEdges(end)])
%         xlabel('ITI length')
%         ylabel('raw residuals 1s pre CS')
%         set(gca, 'tickdir', 'out', 'box', 'off')
    end
%     titleTxt = [sessionName ' ' strrep(cellList{currCell}, '_', ' ')];
%     suptitle(titleTxt)
%     set(gcf, 'renderer', 'painters', 'position', [-1826 144 3457 701])
%     saveFigurePDF(gcf,['C:\Users\cooper\Desktop\analysis\autocorr\residualsTmp\' sessionName '_' cellList{currCell}])
%     close;
end

% dirTmp = dir('C:\Users\cooper\Desktop\analysis\autocorr\residualsTmp\');
% for currFig = 3:length(dirTmp)
%     append_pdfs(['C:\Users\cooper\Desktop\analysis\autocorr\residuals_' p.Results.modelName '_' sheet '.pdf'], ...
%         [dirTmp(currFig).folder '\' dirTmp(currFig).name]);
% end

figure;
for currVar = 1:numVars
    
    sigInds = find(pVal(:,currVar) < 0.05);
    nonSigInds = find(pVal(:,currVar) >= 0.05);
    
    subplot(numVars,2,(currVar-1)*2+1); hold on;
    scatter(varTau(nonSigInds), spikesTau(nonSigInds), [], 'markeredgecolor', [0 0 0])
    scatter(varTau(sigInds), spikesTau(sigInds), [], 'markeredgecolor', [0.75 0 1])
    set(gca, 'tickdir', 'out')
    xlabel([p.Results.varNames{currVar} ' \tau'])
    ylabel('firing rate \tau')
    xlim([0 300])
    ylim([0 300])
    
   subplot(numVars,2,currVar*2)
   numSig = length(sigInds)/numCells;
   for currS = 1:numShifts
       numSig = [numSig sum(shiftP(:,currVar,currS))/numCells];
   end
   plot([0:size(shiftP,3)], numSig, '-', 'color', [0.75 0 1],  'linewidth', 4)
   set(gca, 'tickdir', 'out', 'box', 'off')
   xlabel('number of spike bin shifts')
   ylabel('fraction of neurons with significant \beta coeff')
   ylim([0 max(numSig)])
   
end

set(gcf, 'renderer', 'painters', 'position', [-1497 454 979 402])



end