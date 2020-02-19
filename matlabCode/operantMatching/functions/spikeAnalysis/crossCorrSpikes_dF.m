function [corrCoeffs] = crossCorrSpikes_dF(xlFile, sheet, varargin)

p = inputParser;
% default parameters if none given
p.addParameter('trialFlag', 1);
p.addParameter('figFlag', 1)
p.parse(varargin{:});

[root, sep] = currComputer();

[revForFlagList, sessionCellList, ~] = xlsread(xlFile, sheet);
revForFlagList = revForFlagList(:,1);
cellList = sessionCellList(2:end, 1);
sessionList = sessionCellList(2:end, 2);

timeMax = 181000;
binSize = 30000;
timeBinEdges = [1000:binSize:timeMax];  %no trials shorter than 1s between outcome and CS on
tMax = length(timeBinEdges) - 1;
corrCoeffs = nan(1,length(cellList));

for currCell = 1:length(sessionList)
    
    fprintf('Analyzing cell %d of %d \n', currCell, length(sessionList));
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
        else
            [sessionData] = generateSessionData_nL_operantMatching(sessionName);
        end
    end


    %create arrays for choices and rewards
    responseInds = find(~isnan([sessionData.rewardTime])); % find CS+ trials with a response in the lick window
    allReward_R = [sessionData(responseInds).rewardR]; 
    allReward_L = [sessionData(responseInds).rewardL]; 
    allChoices = NaN(1,length(sessionData(responseInds)));
    allChoices(~isnan(allReward_R)) = 1;
    allChoices(~isnan(allReward_L)) = -1;
    allChoice_R = double(allChoices == 1);
    allChoice_L = double(allChoices == -1);

    allReward_R(isnan(allReward_R)) = 0;
    allReward_L(isnan(allReward_L)) = 0;
    allRewards = zeros(1,length(allChoices));
    allRewards(logical(allReward_R)) = 1;
    allRewards(logical(allReward_L)) = 1;
    
    rwd_Inds = responseInds(allRewards==1);
    sessionTime = sessionData(end).CSon + 3000 - sessionData(1).CSon;       %pad time for reward on last trial
    sessionRwds = [sessionData.rewardTime] - sessionData(1).CSon;     %baseline to start time and convert to s from ms
    session_rwdsArray = zeros(1,sessionTime);
    sessionRwds = sessionRwds(rwd_Inds);
    session_rwdsArray(sessionRwds) = 1;
    
    binEnd = ceil(sessionTime / 15000);
    binInds = reshape([1:binEnd*15000], [15000, binEnd]);
    sessionRwdsCounts = zeros(size(binInds));
    sessionRwdsCounts(binInds <= sessionTime) = session_rwdsArray(binInds <= sessionTime);
    sessionRwdsCounts = sum(sessionRwdsCounts);
    
    
    spikeFields = fields(sessionData);
    clust = find(~cellfun(@isempty,strfind(spikeFields,'C_')) | ~cellfun(@isempty,strfind(spikeFields,'TT')));
    cellInd = find(~cellfun(@isempty,strfind(spikeFields(clust), cellList{currCell})));
    
    %spikes
    sessionSpikes = sessionData(cellInd).allSpikes - sessionData(1).CSon;
    sessionSpikeInd = sessionSpikes(sessionSpikes > 0);
    session_spikeArray = zeros(1,sessionTime);      
    session_spikeArray(sessionSpikeInd(find(sessionSpikeInd<sessionTime))) = 1;
    
    sessionSpikeCounts = zeros(size(binInds));
    sessionSpikeCounts(binInds <= sessionTime) = session_spikeArray(binInds <= sessionTime);
    sessionSpikeCounts = sum(sessionSpikeCounts);
    
    [rho,~] = corr(sessionSpikeCounts', sessionRwdsCounts', 'Type', 'Spearman');
    corrCoeffs(currCell) = rho;
    
end

 figure; histogram(corrCoeffs, 10, 'FaceColor', 'k', 'FaceAlpha', 1)

