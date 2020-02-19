function regOverTime_obs_dF(xlFile, sheet)

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

%set time window for spike analyses
tB = -1500;
tF = 5000;
tW = 500;
tS = 100;
timeToPlot = tB:tS:tF - tW;

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
        sessionTime = (sessionData(end).CSon + tF) - (sessionData(1).CSon + tB);
        csTimes = [sessionData.CSon] -  (sessionData(1).CSon + tB); 
        csTimes = [csTimes(s.responseInds) csTimes(end) + tF];  %add extra csTime for indexing purposes below
        tEndTimes = [0 csTimes(1:end-1)+1500];                  %end times of the previous trial to avoid taking spikes during previous trials
    end

    %% Sort all spikes into a raster-able matrix
    cellInd = find(~cellfun(@isempty,strfind(spikeFields, cellList{currCell})));
    cellNum = find(cellInds == cellInd);
    spikeTimes = sessionData(cellNum).allSpikes - (sessionData(1).CSon + tB);
    spikeTimes = spikeTimes(spikeTimes > 0 & spikeTimes <= sessionTime);
    sessionSpikes = zeros(1, sessionTime);
    sessionSpikes(spikeTimes) = 1;
    smoothSpikes = conv(sessionSpikes, smoothKern);
    smoothSpikes = smoothSpikes(1:(end-length(smoothKern)+1));
    zSpikes = zscore(smoothSpikes);
    
    for regTime = timeToPlot
        spikeTmp = zeros(1,length(s.responseInds));
        for trialInd = 1:length(s.responseInds)
            winInds = [(csTimes(trialInd) + regTime + 1) : (csTimes(trialInd) + regTime + tW)];
            winInds = winInds(winInds > tEndTimes(trialInd) & winInds <  csTimes(trialInd + 1));
            if ~isempty(winInds)
                spikeTmp(trialInd) = mean(zSpikes(winInds));
            end
        end
        
        glm = fitlm([O' C' X' pO'], spikeTmp);
        
        regName = ['reg_' num2str(regTime)];
        regName = strrep(regName, '-', 'min');
        mdlStruct(currCell).(regName).coef = glm.Coefficients.Estimate(2:end);
        mdlStruct(currCell).(regName).pVal = glm.Coefficients.pValue(2:end);
        mdlStruct(currCell).(regName).tStat = glm.Coefficients.tStat(2:end);
        mdlStruct(currCell).(regName).rSqr = glm.Rsquared(1).Ordinary(1);
        mdlStruct(currCell).(regName).warning = lastwarn;
        lastwarn('');
    end
end

O_sig = [];
C_sig = [];
X_sig = []; 
pO_sig = [];
O_coef = [];
C_coef = [];
X_coef = []; 
pO_coef = [];
rSqr = [];
warn_inds = logical([]);
for n = 1:length(sessionList)
    for t = timeToPlot
        regName = ['reg_' num2str(t)];
        regName = strrep(regName, '-', 'min');
        t_ind = find(timeToPlot == t);
        O_sig(t_ind, n) = mdlStruct(n).(regName).pVal(1);
        C_sig(t_ind, n) = mdlStruct(n).(regName).pVal(2);
        X_sig(t_ind, n) = mdlStruct(n).(regName).pVal(3);
        pO_sig(t_ind, n) = mdlStruct(n).(regName).pVal(4);
        O_coef(t_ind, n) = mdlStruct(n).(regName).coef(1);
        C_coef(t_ind, n) = mdlStruct(n).(regName).coef(2);
        X_coef(t_ind, n) = mdlStruct(n).(regName).coef(3);
        pO_coef(t_ind, n) = mdlStruct(n).(regName).coef(4);
        rSqr(t_ind, n) = mdlStruct(n).(regName).rSqr;
        warn_inds(t_ind, n) = ~isempty(mdlStruct(n).(regName).warning);
    end
end
O_sig(warn_inds) = NaN;
C_sig(warn_inds) = NaN;
X_sig(warn_inds) = NaN;
pO_sig(warn_inds) = NaN;
O_coef(warn_inds) = NaN;
C_coef(warn_inds) = NaN;
X_coef(warn_inds) = NaN;
pO_coef(warn_inds) = NaN;
rSqr(warn_inds) = NaN;


figure; 
subplot(1,2,1); hold on;
colors = cool(4);
pCut = 0.05;
plot(timeToPlot + tW, nanmean(O_sig < pCut, 2), 'linewidth', 2, 'Color', colors(1,:))
plot(timeToPlot + tW, nanmean(C_sig < pCut, 2), 'linewidth', 2, 'Color', colors(2,:))
plot(timeToPlot + tW, nanmean(X_sig < pCut, 2), 'linewidth', 2, 'Color', colors(3,:))
plot(timeToPlot + tW, nanmean(pO_sig < pCut, 2), 'linewidth', 2, 'Color', colors(4,:))
xlim([timeToPlot(1)+tW timeToPlot(end)+tW])
legend('Outcome','Choice','C x O','Previous Outcome')
xlabel('Time from cue (ms)')
ylabel('Fraction of significant neurons (P < 0.05)')
ylim([0 1])
set(gca, 'tickdir', 'out')

subplot(1,2,2); hold on;
colors = cool(size(rSqr,2));
for currCell = 1:size(rSqr,2)
    plot(timeToPlot + tW, rSqr(:,currCell), 'color', colors(currCell,:))
end
xlabel('Time from cue (ms)')
ylabel('R^2')
set(gca, 'tickdir', 'out')

suptitle(sheet)
set(gcf, 'renderer', 'painters')


figure; hold on; 
colors = cool(4);
subplot(1,4,1); hold on; 
sigInds = find(O_sig < 0.05);
histogram(O_coef, 20,  'FaceColor', 'k', 'FaceAlpha', 0.3)
histogram(O_coef(sigInds), 20,  'FaceColor', colors(1,:))
title('Outcome')
set(gca, 'tickdir', 'out')

subplot(1,4,2); hold on; 
sigInds = find(C_sig < 0.05);
histogram(C_coef, 20,  'FaceColor', 'k', 'FaceAlpha', 0.3)
histogram(C_coef(sigInds), 20,  'FaceColor', colors(2,:))
title('Choice')
set(gca, 'tickdir', 'out')

subplot(1,4,3); hold on; 
sigInds = find(X_sig < 0.05);
histogram(X_coef, 20,  'FaceColor', 'k', 'FaceAlpha', 0.3)
histogram(X_coef(sigInds), 20,  'FaceColor', colors(3,:))
title('Choice x Outcome')
set(gca, 'tickdir', 'out')

subplot(1,4,4); hold on; 
sigInds = find(pO_sig < 0.05);
histogram(pO_coef, 20,  'FaceColor', 'k', 'FaceAlpha', 0.3)
histogram(pO_coef(sigInds), 20,  'FaceColor', colors(4,:))
title('Previous Outcome')
set(gca, 'tickdir', 'out')

set(gcf, 'renderer', 'painters')
