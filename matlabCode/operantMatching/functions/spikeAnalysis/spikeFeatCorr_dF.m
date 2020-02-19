function [rhos, pVals] = spikeFeatCorr_dF(xlFile, sheet)

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
tb = 1.5;
tf = 5;
time = -1000*tb:1000*tf;
smoothWin = 250;
trialBeg = tb*1000;
CSoff = tb*1000 + 500;


rhos = struct;
pVals = struct;
rSqurs = struct;

sessionName = [];
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
        [s] = behAnalysisNoPlot_opMD(sessionName, 'revForFlag', revForFlagList(currCell));


    %     rwd_Inds = s.responseInds(s.allRewards==1);
    %     sessionTime = s.behSessionData(end).CSon + 3000 - s.behSessionData(1).CSon;       %pad time for reward on last trial
    %     sessionRwds = [s.behSessionData.rewardTime] - s.behSessionData(1).CSon;     %baseline to start time and convert to s from ms
    %     session_rwdsArray = zeros(1,sessionTime);
    %     sessionRwds = sessionRwds(rwd_Inds);
    %     session_rwdsArray(sessionRwds) = 1;
    % %     normKern = normpdf(-60000:0,0,30000);
    % %     normKern = normKern / sum(normKern);
    %     x = [-121000:0];
    %     expKern = 3.507*exp(0.00006308*x);
    %     sessionRwdsSmooth = conv(session_rwdsArray, expKern);
    %     CSonTimes = [s.behSessionData.CSon] - s.behSessionData(1).CSon;
    %     CSonTimes(1) = 1;
    %     trialRwdsSmooth = sessionRwdsSmooth(CSonTimes);
    %     trialRwdsSmooth = trialRwdsSmooth(s.responseInds);
        x = [-12:0];
        expKern = 2.861*exp(0.5821*x);
        allRewardsBin = s.allRewards;
        allRewardsBin(allRewardsBin==-1) = 1;
        trialRwdsSmooth = conv(allRewardsBin, expKern);
        trialRwdsSmooth = trialRwdsSmooth(1:end-(length(expKern)-1));
        trialRwdsSmooth = [0 trialRwdsSmooth(1:end-1)];

    %     if ~isempty(p.Results.trialList)
    %         s.responseInds = s.responseInds(p.Results.trialList);
    %         s.sessionRwds = s.sessionRwds(p.Results.trialList);
    %         rwdHx = rwdHx(p.Results.trialList);
    %     end
    end

    %% Sort all spikes into a raster-able matrix

    spikeFields = fields(sessionData);
    cellInd = find(~cellfun(@isempty,strfind(spikeFields, cellList{currCell})));

    allTrial_spike = {};
    for k = 1:length(sessionData)
        if k == 1
            prevTrial_spike = [];
            currTrial_lick = [];
        else
            prevTrial_spikeInd = sessionData(k-1).(spikeFields{cellInd}) > sessionData(k-1).trialEnd-tb*1000;
            prevTrial_spike = sessionData(k-1).(spikeFields{cellInd})(prevTrial_spikeInd) - sessionData(k).CSon;
        end

        currTrial_spikeInd = sessionData(k).(spikeFields{cellInd}) < sessionData(k).CSon+tf*1000;
        currTrial_spike = sessionData(k).(spikeFields{cellInd})(currTrial_spikeInd) - sessionData(k).CSon;

            allTrial_spike{k} = [prevTrial_spike currTrial_spike];
        if ~isnan(sessionData(k).licksL)
            currTrial_lickInd = sessionData(k).licksL < sessionData(k).CSon + tf*1000;
            currTrial_lick = sessionData(k).licksL(currTrial_lickInd) - sessionData(k).CSon;
        elseif ~isnan(sessionData(k).licksR)
            currTrial_lickInd = sessionData(k).licksR < sessionData(k).CSon + tf*1000;
            currTrial_lick = sessionData(k).licksR(currTrial_lickInd) - sessionData(k).CSon;  
        else
            currTrial_lick = 0;
        end
        allTrial_lick{k} = [currTrial_lick];
    end

    % sometimes no licks/spikes are considered 1x0 and sometimes they are []
    % plotSpikeRaster does not place nicely with [] so this converts all empty indices to 1x0
    allTrial_spike(cellfun(@isempty,allTrial_spike)) = {zeros(1,0)}; 


    %% set time window and smoothing parameters, run analysis for all cells

    for ind = 1:length(s.behSessionData)
        trialDurDiff(ind) = (s.behSessionData(ind).trialEnd - s.behSessionData(ind).CSon)- tf*1000;
    end

    allTrial_spikeMatx = NaN(length(sessionData),length(time));
    for j = 1:length(allTrial_spike)
        tempSpike = allTrial_spike{j};
        tempSpike = tempSpike + tb*1000; % add this to pad time for SDF
        allTrial_spikeMatx(j,tempSpike) = 1;
        if trialDurDiff(j) < 0
            allTrial_spikeMatx(j, isnan(allTrial_spikeMatx(j, 1:end+trialDurDiff(j)))) = 0;  %converts within trial duration NaNs to 0's
        else
            allTrial_spikeMatx(j, isnan(allTrial_spikeMatx(j,:))) = 0;
        end
        if sum(allTrial_spikeMatx(j,:)) == 0     %if there is no spike data for this trial, don't count it
            allTrial_spikeMatx(j,:) = NaN;
        end
    end

    %% find features of spike rate on each trial    
    preCSspikeCount = nan(1, length(allTrial_spike));
    postCSspikeCount = nan(1, length(allTrial_spike));
    maxFRcs = nan(1, length(allTrial_spike));
    minFRcs = nan(1, length(allTrial_spike));
    
    for j = 1:length(allTrial_spike)
        if ~isempty(allTrial_spikeMatx(j))
            preCSspikeCount(j) = sum(allTrial_spikeMatx(j, 1:trialBeg));              %find total spikes before CS on
            postCSspikeCount(j) = sum(allTrial_spikeMatx(j, trialBeg:CSoff));

            spikeTemp = fastsmooth(allTrial_spikeMatx(j,:)*1000, smoothWin, 3);         %smooth raw spikes to find features of spike rate
            maxFRcs(j) = max(spikeTemp(trialBeg:CSoff));
            minFRcs(j) = min(spikeTemp(trialBeg:CSoff));
            if ~isnan(maxFRcs(j))
                maxFRtime(j) = find(spikeTemp == max(spikeTemp(trialBeg:CSoff)), 1);
            else
                maxFRtime(j) = NaN;
            end

        else
            preCSspikeCount(j) = NaN;
            postCSspikeCount(j) = NaN;
            maxFRcs(j) = NaN;
            minFRcs(j) = NaN;
            maxFRtime(j) = NaN;
        end
    end
    
    preCSspikeCount = zscore(preCSspikeCount);
    postCSspikeCount = zscore(postCSspikeCount);
    maxFRcs = zscore(maxFRcs);
    minFRcs = zscore(minFRcs);
    
    
    [rho,pVal] = corr(maxFRcs(s.responseInds)', trialRwdsSmooth', 'Type', 'Spearman');
    rhos.maxFR_rwdHist(currCell) = rho;
    pVals.maxFR_rwdHist(currCell) = pVal;
    [rho,pVal] = corr(postCSspikeCount(s.responseInds)', trialRwdsSmooth', 'Type', 'Spearman');
    rhos.postCS_rwdHist(currCell) = rho;
    pVals.postCS_rwdHist(currCell) = pVal;
    [rho,pVal] = corr(preCSspikeCount(s.responseInds)', trialRwdsSmooth', 'Type', 'Spearman');
    rhos.preCS_rwdHist(currCell) = rho;
    pVals.preCS_rwdHist(currCell) = pVal;
    
    
    [rho,pVal] = corr(maxFRcs(s.responseInds(2:end))', s.allRewardsBinary(1:end-1)', 'Type', 'Spearman');
    rhos.maxFR_prevRwd(currCell) = rho;
    pVals.maxFR_prevRwd(currCell) = pVal;
    [rho,pVal] = corr(postCSspikeCount(s.responseInds(2:end))', s.allRewardsBinary(1:end-1)', 'Type', 'Spearman');
    rhos.postCS_prevRwd(currCell) = rho;
    pVals.postCS_prevRwd(currCell) = pVal;
    [rho,pVal] = corr(preCSspikeCount(s.responseInds(2:end))', s.allRewardsBinary(1:end-1)', 'Type', 'Spearman');
    rhos.preCS_prevRwd(currCell) = rho;
    pVals.preCS_prevRwd(currCell) = pVal;
    
    
    
end

sigInds_postCS = find(pVals.postCS_rwdHist < 0.05);
nonsigInds_postCS = find(pVals.postCS_rwdHist > 0.05);

sigInds_preCS = find(pVals.preCS_rwdHist < 0.05);
nonsigInds_preCS = find(pVals.preCS_rwdHist > 0.05);

bothSigInds = intersect(sigInds_preCS, sigInds_postCS);
sigInds_postCS(ismember(sigInds_postCS, bothSigInds)) = [];
sigInds_preCS(ismember(sigInds_preCS, bothSigInds)) = [];

bothNonSigInds = intersect(nonsigInds_preCS, nonsigInds_postCS);



figure;
subplot(2,2,1); hold on;
scatter(rhos.postCS_rwdHist(bothNonSigInds), rhos.preCS_rwdHist(bothNonSigInds), 'k', 'filled')
scatter(rhos.postCS_rwdHist(sigInds_postCS), rhos.preCS_rwdHist(sigInds_postCS), 'r', 'filled')
scatter(rhos.postCS_rwdHist(sigInds_preCS), rhos.preCS_rwdHist(sigInds_preCS), 'c', 'filled')
scatter(rhos.postCS_rwdHist(bothSigInds), rhos.preCS_rwdHist(bothSigInds), [], [0.7 0 1], 'filled')
title('Recency-weighted reward history')
set(gca, 'TickDir', 'out')

subplot(2,2,3)
pie([length(bothSigInds) length(sigInds_postCS), length(sigInds_preCS), length(bothNonSigInds)])
colormap([0.7 0 1; 1 0 0; 0 1 1; 0 0 0])
legend('Both significant', 'CS firing rate significant', 'Pre-CS firing rate significant', 'Neither significant')

sigInds_postCS = find(pVals.postCS_prevRwd < 0.05);
nonsigInds_postCS = find(pVals.postCS_prevRwd > 0.05);
sigInds_preCS = find(pVals.preCS_prevRwd < 0.05);
nonsigInds_preCS = find(pVals.preCS_prevRwd > 0.05);
bothSigInds = intersect(sigInds_preCS, sigInds_postCS);
sigInds_postCS(ismember(sigInds_postCS, bothSigInds)) = [];
sigInds_preCS(ismember(sigInds_preCS, bothSigInds)) = [];
bothNonSigInds = intersect(nonsigInds_preCS, nonsigInds_postCS);


subplot(2,2,2); hold on;
scatter(rhos.postCS_prevRwd(bothNonSigInds), rhos.preCS_prevRwd(bothNonSigInds), 'k', 'filled')
scatter(rhos.postCS_prevRwd(sigInds_postCS), rhos.preCS_prevRwd(sigInds_postCS), 'r', 'filled')
scatter(rhos.postCS_prevRwd(sigInds_preCS), rhos.preCS_prevRwd(sigInds_preCS), 'c', 'filled')
scatter(rhos.postCS_prevRwd(bothSigInds), rhos.preCS_prevRwd(bothSigInds), [], [0.7 0 1], 'filled')
title('Previous reward')

subplot(2,2,4)
pie([length(bothSigInds) length(sigInds_postCS), length(sigInds_preCS), length(bothNonSigInds)])
colormap([0.7 0 1; 1 0 0; 0 1 1; 0 0 0])


for i=1:2
    subplot(2,2,i)
    set(gca, 'TickDir', 'out')
    xlabel('Spearman''s \rho CS firing rate')
    ylabel('Spearman''s \rho pre-CS firing rate')
end

end