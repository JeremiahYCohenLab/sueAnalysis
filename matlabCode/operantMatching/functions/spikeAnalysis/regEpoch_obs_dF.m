function regEpoch_obs_dF(xlFile, sheet)

[root, sep] = currComputer();

[numbers, sessionCellList, ~] = xlsread(xlFile, sheet);
revForFlagList = numbers(:,1);
intanFlagList = numbers(:,2);
cellList = sessionCellList(2:end, 1);
sessionList = sessionCellList(2:end, 2);

if size(numbers,2) > 2
    trialList = numbers(:,3:4);
else
    trialList = nan(length(cellList),2);
end 

%set epoch windows for spike analyses
tB = [-1000 0];
tF = [0 500];
tW = tF - tB;

if any(tB < 0)
    tL = min(tB(tB < 0));
else
    tL = 0;
end

smoothKern = normpdf(0:5000, 0, 250);
smoothKern = smoothKern/sum(smoothKern);

mdlStruct = [];
sessionName = [];
lastwarn('');

for currCell = 1:length(sessionList)
    
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

        if revForFlagList(currCell) == 1
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
                if intanFlagList(currCell) == 1
                    [sessionData] = generateSessionData_intan_operantMatching(sessionName);
                else
                    [sessionData] = generateSessionData_nL_operantMatching(sessionName);
                end
            end
        end
        [s] = behDataForStruct_opMD(sessionName, 'revForFlag', revForFlagList(currCell));
        O = s.allRewardsBinary;
        C = s.allChoices;
        X = O.*C;
        pO = [nan O(1:end-1)];
        
        spikeFields = fields(sessionData);
        cellInds = find(~cellfun(@isempty,strfind(spikeFields, 'TT')));
        sessionTime = (sessionData(end).CSon + max(tF)) - (sessionData(1).CSon + tL);
        csTimes = [sessionData.CSon] -  (sessionData(1).CSon + tL);
        csTimes = [csTimes(s.responseInds) csTimes(end) + max(tF)];
        tEndTimes = [0 csTimes(1:end-1)+2500];
    end

    %% Sort all spikes into a raster-able matrix
    cellInd = find(~cellfun(@isempty,strfind(spikeFields, cellList{currCell})));
    cellNum = find(cellInds == cellInd);
    spikeTimes = sessionData(cellNum).allSpikes - (sessionData(1).CSon + tL);
    spikeTimes = spikeTimes(spikeTimes > 0 & spikeTimes <= sessionTime);
    sessionSpikes = zeros(1, sessionTime);
    sessionSpikes(spikeTimes) = 1;
    smoothSpikes = conv(sessionSpikes, smoothKern);
    smoothSpikes = smoothSpikes(1:(end-length(smoothKern)+1));
    zSpikes = zscore(smoothSpikes);

    for eInd = 1:length(tB)
        spikeTmp = nan(1,length(s.responseInds));
        for trialInd = 1:length(s.responseInds)
            winInds = [(csTimes(trialInd) + tB(eInd) + 1) : (csTimes(trialInd) + tW(eInd))];
            winInds = winInds(winInds > tEndTimes(trialInd) & winInds <  csTimes(trialInd + 1));
            if ~isempty(winInds)
                spikeTmp(trialInd) = mean(zSpikes(winInds));
            end
        end
        
        glm = fitlm([O' C' X' pO'], spikeTmp);
        
        regName = ['reg_' num2str(tB(eInd))];
        regName = strrep(regName, '-', 'min');
        mdlStruct(currCell).(regName).coef = glm.Coefficients.Estimate(2:end);
        mdlStruct(currCell).(regName).pVal = glm.Coefficients.pValue(2:end);
        mdlStruct(currCell).(regName).tStat = glm.Coefficients.tStat(2:end);
        mdlStruct(currCell).(regName).warning = lastwarn;
        lastwarn('');
    end
end

figure; hold on;
O_sig = [];
C_sig = [];
X_sig = []; 
pO_sig = [];
warn_inds = logical([]);
for n = 1:length(sessionList)
    for t = tB
        regName = ['reg_' num2str(t)];
        regName = strrep(regName, '-', 'min');
        t_ind = find(tB == t);
        O_sig(t_ind, n) = mdlStruct(n).(regName).tStat(1);
        C_sig(t_ind, n) = mdlStruct(n).(regName).tStat(2);
        X_sig(t_ind, n) = mdlStruct(n).(regName).tStat(3);
        pO_sig(t_ind, n) = mdlStruct(n).(regName).tStat(4);
        warn_inds(t_ind, n) = ~isempty(mdlStruct(n).(regName).warning);
    end
end
O_sig(warn_inds) = NaN;
C_sig(warn_inds) = NaN;
X_sig(warn_inds) = NaN;
pO_sig(warn_inds) = NaN;

regNames = [{'O pre'} {'O cs'} {'C pre'} {'C cs'} {'X pre'} {'X cs'} {'pO pre'} {'pO cs'} ];

regMatx = [O_sig; C_sig; X_sig; pO_sig];
numRegs = size(regMatx,1);
combs = nchoosek(numRegs,2);
cols = round(combs/7);
comb = 1;
linetype = 'k';
tCut = 1.96;
for i = 1:numRegs-1
    for j = i+1:numRegs
        subplot(cols,7,comb); hold on;
        for k = 1:size(regMatx,2)
            if abs(regMatx(i,k)) > tCut &  abs(regMatx(j,k)) > tCut
                scatter(regMatx(i,k), regMatx(j,k), 'c', 'filled')
            else
                scatter(regMatx(i,k), regMatx(j,k), 'k', 'filled')
            end
        end 
        plot([-6 6], [tCut tCut], ':k')
        plot([-6 6], [-tCut -tCut], ':k')
        xlim([-6 6])
        ylim([-6 6])
        vline([-tCut tCut], linetype)
        xlabel(regNames{i})
        ylabel(regNames{j})
        set(gca, 'tickdir', 'out')
        comb = comb + 1;
    end
end
set(gcf, 'renderer', 'painters')
