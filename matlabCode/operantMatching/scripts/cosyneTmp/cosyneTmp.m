function cosyneTmp(xlFile, sheet, beh, varargin)

p = inputParser;
% default parameters if none given
p.addParameter('biasFlag',0);
p.addParameter('bernFlag', 1);
p.addParameter('saveFigFlag', 0);
p.addParameter('modelName', 'sixParam_absPePeAN_bi')
p.addParameter('varNames', [{'peBar'} {'pePe'} {'pe'}]);
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

% smoothKern = normpdf(0:30000, 0, 10000);
% smoothKern = smoothKern/sum(smoothKern);
smoothKernSize = 30000;
smoothKern = ones(1,smoothKernSize);
roughKernSize = 10000;
roughKern = ones(1,roughKernSize);

binWin = 30000;
respWin = 2000;

numVars = length(p.Results.varNames);
sessionName = [];

numCells = length(sessionList);
rho = nan(numCells, numVars);
pVal = nan(numCells, numVars);
numBins = 6;
allMeanSpikes = nan(numCells, numVars, numBins);

for currCell = 1:numCells
    
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
        cellInds = find(~cellfun(@isempty,strfind(spikeFields, 'TT')));
        sessionTime = (sessionData(end).CSon + tF) - (sessionData(1).CSon + tB);
        csTimes = [sessionData.CSon] -  (sessionData(1).CSon + tB);
        csTimes = [csTimes(s.responseInds) csTimes(end) + tF];
        
        oTime = nan(1, length(s.responseInds));
        for currT = 1:length(s.responseInds)
            if isnan(sessionData(s.responseInds(currT)).rewardR)
                if length(sessionData(s.responseInds(currT)).licksL) > 1 & ...
                        sessionData(s.responseInds(currT)).licksL(2) - sessionData(s.responseInds(currT)).CSon <  respWin;
                    oTime(currT) = sessionData(s.responseInds(currT)).licksL(2);
                else
                    oTime(currT) = sessionData(s.responseInds(currT)).licksL(1);
                end
            else
                if length(sessionData(s.responseInds(currT)).licksR) > 1 & ...
                        sessionData(s.responseInds(currT)).licksR(2) - sessionData(s.responseInds(currT)).CSon <  respWin;
                    oTime(currT) = sessionData(s.responseInds(currT)).licksR(2);
                else
                    oTime(currT) = sessionData(s.responseInds(currT)).licksR(1);
                end
            end
        end
        oTime = oTime - (sessionData(1).CSon + tB);

        %get model terms
        if p.Results.bernFlag
            modelPath = [root animalName sep animalName 'sorted' sep 'stan' sep 'bernoulli' sep p.Results.modelName sep animalName...
            beh '_' p.Results.modelName '.mat'];
        else
            modelPath = [root animalName sep animalName 'sorted' sep 'stan' sep p.Results.modelName sep animalName...
            beh '_' p.Results.modelName '.mat'];
        end
        t = generateStanModelTerms_opMD(p.Results.modelName, modelPath, sessionName, 1, revForFlag(currCell));
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
    
    zSpikes = zscore(smoothSpikes);
    zSpikes = downsample(zSpikes, 1000);
    smoothSpikes = smoothSpikes / (smoothKernSize/1000);
    smoothSpikes = downsample(smoothSpikes(smoothKernSize:end), 1000);
    roughSpikes = roughSpikes / (roughKernSize/1000);
    roughSpikes = downsample(roughSpikes(roughKernSize:end), 1000);
    
    
    itiSpikes = nan(1, length(csTimes)-1);
    for i = 1:length(csTimes)-1
        itiSpikes(i) = sum(sessionSpikes(csTimes(i)+2000 : csTimes(i+1))) / (csTimes(i+1) - (csTimes(i)+2000));
    end
    [itiRho(currCell), itiP(currCell)] = corr(allVars(2:end)', itiSpikes(1:end-1)');
       
    
    trialSpikes = nan(1, length(csTimes)-1);
    for i = 1:length(csTimes)-1
        trialSpikes(i) = sum(sessionSpikes(csTimes(i): csTimes(i)+2000)) / 2000;
    end
    [trialRho(currCell), trialP(currCell)] = corr(allVars(2:end)', trialSpikes(1:end-1)');
    
    
    csSpikes = nan(1, length(csTimes)-1);
    for i = 1:length(csTimes)-1
        csSpikes(i) = sum(sessionSpikes(csTimes(i): oTime(i))) / (oTime(i) - csTimes(i));
    end
    [csRho(currCell), csP(currCell)] = corr(allVars', csSpikes');
    
    
    oSpikes = nan(1, length(csTimes)-1);
    for i = 1:length(csTimes)-1
        oSpikes(i) = sum(sessionSpikes(oTime(i): csTimes(i)+2000)) / ((csTimes(i)+2000) - oTime(i));
    end
    [oRho(currCell), oP(currCell)] = corr(allVars(2:end)', oSpikes(1:end-1)');
    
    
%     binSpikes = [];
%     binVar = [];
%     timeBins = floor(sessionTime/binWin);
%     for i = 1:timeBins
%         binSpikes(i) = sum(sessionSpikes((i-1)*binWin+1 : i*binWin));
%         binVar(i) = sum(allVarsSesh((i-1)*binWin+1 : i*binWin));
%     end
%     binSpikes(i+1) = sum(sessionSpikes(i*binWin+1 : end));
%     binVar(i+1) = sum(allVarsSesh(i*binWin+1 : end));
%     [~, woof(currCell)] = corr(binVar', binSpikes');
    
    
    
end

itiSig = sum(itiP < 0.05) / length(itiP);
trialSig = sum(trialP < 0.05) / length(trialP);
csSig = sum(csP < 0.05) / length(csP);
oSig = sum(oP < 0.05) / length(oP);

colors = cool(4);
figure; 
subplot(1,2,1); hold on;
bar(1, itiSig, 'FaceColor', colors(1,:));
bar(2, trialSig, 'FaceColor', colors(2,:));
bar(3, csSig, 'FaceColor', colors(3,:));
bar(4, oSig, 'FaceColor', colors(4,:));
xticks([1:4])
xticklabels({'ITI', 'trial', 'CS', 'outcome'})
set(gca, 'tickdir', 'out')
ylabel('fraction of significant neurons')

edges = [-1:0.1:1];
subplot(1,2,2); hold on;
histogram(itiRho, edges, 'facecolor', colors(1,:))
histogram(trialRho, edges, 'facecolor', colors(2,:))
histogram(csRho, edges, 'facecolor', colors(3,:))
histogram(oRho, edges, 'facecolor', colors(4,:))
set(gca, 'tickdir', 'out')
ylabel('number of neurons')
xlabel('pearson \rho')

suptitle(p.Results.varNames)
set(gcf, 'renderer', 'painters', 'position', [-1636 260 1300 676])


figure
subplot(1,4,1); hold on;
mdl = fitlm(itiRho, oRho);
scatter(itiRho, oRho, 'k', 'filled')
xlim([-0.55 0.55])
ylim([-0.55 0.55])
x = [-0.55:0.01:0.55];
y = mdl.Coefficients.Estimate(2)*x + mdl.Coefficients.Estimate(1);
plot(x,y, 'k')
legTxt = [{'coeff: ' mdl.Coefficients.Estimate(2)}, {'pVal: ' num2str(mdl.Coefficients.pValue(2))}]
legend([{strcat('coeff: ', num2str(mdl.Coefficients.Estimate(2)))}, {strcat('pVal: ', num2str(mdl.Coefficients.pValue(2)))}])
xlabel('iti \rho')
ylabel('outcome \rho')
set(gca, 'tickdir', 'out')
xlim([-0.55 0.55])
ylim([-0.55 0.55])

subplot(1,4,2); hold on;
mdl = fitlm(itiRho, csRho);
scatter(itiRho, csRho, 'k', 'filled')
xlim([-0.55 0.55])
ylim([-0.55 0.55])
x = [-0.55:0.01:0.55];
y = mdl.Coefficients.Estimate(2)*x + mdl.Coefficients.Estimate(1);
plot(x,y, 'k')
legTxt = [{'coeff: ' mdl.Coefficients.Estimate(2)}, {'pVal: ' num2str(mdl.Coefficients.pValue(2))}]
legend([{strcat('coeff: ', num2str(mdl.Coefficients.Estimate(2)))}, {strcat('pVal: ', num2str(mdl.Coefficients.pValue(2)))}])
xlabel('iti \rho')
ylabel('cs \rho')
set(gca, 'tickdir', 'out')
xlim([-0.55 0.55])
ylim([-0.55 0.55])


subplot(1,4,3); hold on;
mdl = fitlm(itiRho, trialRho);
scatter(itiRho, trialRho, 'k', 'filled')
xlim([-0.55 0.55])
ylim([-0.55 0.55])
x = [-0.55:0.01:0.55];
y = mdl.Coefficients.Estimate(2)*x + mdl.Coefficients.Estimate(1);
plot(x,y, 'k')
legTxt = [{'coeff: ' mdl.Coefficients.Estimate(2)}, {'pVal: ' num2str(mdl.Coefficients.pValue(2))}]
legend([{strcat('coeff: ', num2str(mdl.Coefficients.Estimate(2)))}, {strcat('pVal: ', num2str(mdl.Coefficients.pValue(2)))}])
xlabel('iti \rho')
ylabel('trial \rho')
set(gca, 'tickdir', 'out')

subplot(1,4,4); hold on;
mdl = fitlm(csRho, oRho);
scatter(csRho, oRho, 'k', 'filled')
xlim([-0.55 0.55])
ylim([-0.55 0.55])
x = [-0.55:0.01:0.55];
y = mdl.Coefficients.Estimate(2)*x + mdl.Coefficients.Estimate(1);
plot(x,y, 'k')
legTxt = [{'coeff: ' mdl.Coefficients.Estimate(2)}, {'pVal: ' num2str(mdl.Coefficients.pValue(2))}]
legend([{strcat('coeff: ', num2str(mdl.Coefficients.Estimate(2)))}, {strcat('pVal: ', num2str(mdl.Coefficients.pValue(2)))}])
xlabel('cs \rho')
ylabel('outcome \rho')
set(gca, 'tickdir', 'out')
xlim([-0.55 0.55])
ylim([-0.55 0.55])


suptitle(p.Results.varNames)
set(gcf, 'renderer', 'painters', 'position', [-1636 260 1300 676])

end