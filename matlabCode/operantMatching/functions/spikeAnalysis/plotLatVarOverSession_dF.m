function plotLatVarOverSession_dF(xlFile, sheet, beh, varargin)

p = inputParser;
% default parameters if none given
p.addParameter('biasFlag',0);
p.addParameter('bernFlag', 1);
p.addParameter('plotFlag', 0);
p.addParameter('saveFigFlag', 0);
p.addParameter('modelName', 'sixParam_absPePeAN_bi')
p.addParameter('varNames', [{'peBar'} {'pePe'} {'pe'}]);
p.addParameter('sessionParamFlag', 1);
p.addParameter('coefWindow', []);
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

smoothKernSize = 30000;
smoothKern = ones(1,smoothKernSize);
roughKernSize = 10000;
roughKern = ones(1,roughKernSize);
binWin = 10000;

numVars = length(p.Results.varNames);
sessionName = [];

numCells = length(sessionList);
rho = nan(numCells, numVars);
pVal = nan(numCells, numVars);
numBins = 6;
allMeanSpikes = nan(numCells, numVars, numBins);

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
    smoothSpikes = conv(sessionSpikes, smoothKern);
    smoothSpikes = smoothSpikes(1:(end-smoothKernSize+1));
    roughSpikes = conv(sessionSpikes, roughKern);
    roughSpikes = roughSpikes(1:(end-roughKernSize+1));
    
    binSpikes = movsum(sessionSpikes, binWin);
    binSpikes = downsample(binSpikes, binWin, binWin-1);
    
    zSpikes = zscore(smoothSpikes);
    zSpikes = downsample(zSpikes(smoothKernSize:end), 1000);
    smoothSpikes = smoothSpikes / (smoothKernSize/1000);
    smoothSpikes = downsample(smoothSpikes(smoothKernSize:end), 1000);
    roughSpikes = roughSpikes / (roughKernSize/1000);
    roughSpikes = downsample(roughSpikes(roughKernSize:end), 1000);
    
    
    colors = cool(length(p.Results.varNames));
    figure; hold on;
    subplot(numVars+1,6,[1:5]); hold on;
    plot(roughSpikes, '-', 'color', [0.6 0.6 0.6], 'linewidth', 1.3)
    plot(smoothSpikes, '-k', 'linewidth', 1.3)
    xticks([1:600:length(smoothSpikes)])    
    xticklabels([1:length([1:600:length(smoothSpikes)])] * 10 - 10)
    xlabel('Time (min)')
    ylabel('Spikes / s')
    xlim([0 length(smoothSpikes)])
    legend('10s kern', '30s kern')
    set(gca, 'tickdir', 'out')
    
    dVar = [];
    for currVar = 1:numVars
        subplot(numVars+1,6,[currVar*6+1:currVar*6+5]); hold on;
        dVar = downsample(allVarsSesh(currVar,smoothKernSize:end), 1000);
        plot(dVar, '-', 'color', colors(currVar,:), 'linewidth', 1.3)
        %[rho(currCell, currVar), pVal(currCell, currVar)] = corr(dVar', zSpikes', 'Type', 'Spearman');
        
        xticks([1:600:length(smoothSpikes)])
        xticklabels([1:length([1:600:length(smoothSpikes)])] * 10 - 10)
        xlabel('Time (min)')
        ylabel('z-scored values')
        xlim([0 length(smoothSpikes)])
        legend([p.Results.varNames{currVar}])
        set(gca, 'tickdir', 'out')
        xl = xlim; yl = ylim;
        
        binVar = movsum(allVarsSesh(currVar,:), binWin);
        binVar = downsample(binVar, binWin, binWin-1);
        [rho(currCell, currVar), pVal(currCell, currVar)] = corr(zscore(binVar'), zscore(binSpikes'));
        text('Position', [0 0.3 0.5], 'HorizontalAlignment', 'left', 'VerticalAlignment', 'top', 'Units', 'Normalized',...
            'string', sprintf(['rho = ' num2str(rho(currCell, currVar)) '\npVal = ' num2str(pVal(currCell, currVar))]));
        
        subplot(numVars+1,6,[currVar*6+6]); hold on
        [sortInds, edges] = discretize(dVar, numBins);
        for currBin = 1:numBins
            tmp = zSpikes(sortInds==currBin);
            meanSpikes(currBin) = mean(tmp);
            semSpikes(currBin) = std(tmp) / sqrt(length(tmp));
        end
        scatter(dVar, zSpikes, 'k', 'filled')
        x = edges(1:end-1) + diff(edges(1:2))/2;
        errorbar(x, meanSpikes, semSpikes, 'Color', colors(currVar,:),'linewidth',2);
        xlim([x(1)-0.2 x(end)+0.2]);
        xticks([])
        xlabel(p.Results.varNames{currVar})
        set(gca, 'tickdir', 'out')
        
        
        if pVal(currCell, currVar) < 0.05
            if meanSpikes(end) > meanSpikes(1)
                allMeanSpikes(currCell, currVar, :) = meanSpikes;
            else
                allMeanSpikes(currCell, currVar, :) = -meanSpikes;
            end
        end       
    end
    
    titleTxt = strrep([sessionName ' - ' cellList{currCell}], '_', ' ');
    suptitle(titleTxt)
    set(gcf, 'renderer', 'painters', 'position', [-1343 42 1012 954]);
    if p.Results.plotFlag == 0
        close;
    end
    
    if p.Results.saveFigFlag
        saveFigurePDF(gcf,['C:\Users\cooper\Desktop\analysis\population pdfs\latVarTmp\' sessionName '_' cellList{currCell}])
    end
    
    
end

%amend cell figures as one pdf
if p.Results.saveFigFlag
    dirTmp = dir('C:\Users\cooper\Desktop\analysis\population pdfs\latVarTmp\');
    for currFig = 3:length(dirTmp)
        append_pdfs(['C:\Users\cooper\Desktop\analysis\population pdfs\latVar_' p.Results.modelName '_' sheet '.pdf'], ...
            [dirTmp(currFig).folder '\' dirTmp(currFig).name]);
    end
end


%% plot averages

figure;
for currVar = 1:numVars
    subplot(numVars,3,(currVar-1)*3+1)
    title(p.Results.varNames{currVar})
    plotFilled(x, squeeze(allMeanSpikes(:,currVar,:)), colors(currVar,:));
    set(gca, 'tickdir', 'out', 'box', 'off')
    xlim([x(1) - 0.25 x(end) + 0.25])
    ylabel('z-scored firing rate')
    titleTxt = strrep(p.Results.varNames{currVar}, '_', ' ');
    xlabel(titleTxt)
    
    subplot(numVars,3,(currVar-1)*3+2)
    numSig = sum(pVal(:,currVar) < 0.05);
    h = pie([length(pVal(:,currVar))-numSig, numSig]);
    TextArr = findobj(h, 'Type','Text');   
    PatchArr = findobj(h, 'Type','patch');      
    PatchArr(1).FaceColor = [1 1 1];  
    PatchArr(2).FaceColor = colors(currVar,:);
    title(titleTxt)
    
    subplot(numVars,3,(currVar-1)*3+3)
    histogram(rho(:,currVar), 'FaceColor', colors(currVar,:))
    set(gca, 'tickdir', 'out', 'box', 'off')
    ylabel('count')
    xlabel('correlation coeff')

end
set(gcf, 'renderer', 'painters', 'position', [-1352 234 819 656])

end