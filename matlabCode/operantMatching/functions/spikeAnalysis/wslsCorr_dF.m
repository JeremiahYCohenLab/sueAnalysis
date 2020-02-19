function [rhos, pVals] = wslsCorr_dF(xlFile, sheet)

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
trialEnd = tb*1000 + 1500;


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
        
        changeChoice = [false abs(diff(s.allChoices)) > 0];
        changeChoice = changeChoice(2:end);
        rwdStay = changeChoice == 0 & s.allRewardsBinary(1:end-1)==1;
        rwdSwitch = changeChoice == 1 & s.allRewardsBinary(1:end-1)==1;
        noRwdStay = changeChoice == 0 & s.allRewardsBinary(1:end-1)==0;
        noRwdSwitch = changeChoice == 1 & s.allRewardsBinary(1:end-1)==0;

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
            postCSspikeCount(j) = sum(allTrial_spikeMatx(j, CSoff:trialEnd));

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
   
   
    spikeTmp = postCSspikeCount(s.responseInds);
    spikeTmp = spikeTmp(1:end-1);
    mdl = fitlm([rwdStay' rwdSwitch' noRwdStay' noRwdSwitch'], spikeTmp);

    coefVals = mdl.Coefficients.Estimate(2:5);
    CIbands = coefCI(mdl);
    errorL = abs(coefVals - CIbands(2:5,1));
    errorU = abs(coefVals - CIbands(2:5,2));
    figure; hold on
    colors = cool(4);
    for currCoeff = 1:4
        errorbar(currCoeff,coefVals(currCoeff),errorL(currCoeff),errorU(currCoeff),'-o', 'Color', colors(currCoeff,:),'linewidth',2)
        if mdl.Coefficients.pValue(currCoeff+1) < 0.05
            plot(currCoeff, [mdl.Coefficients.Estimate(currCoeff+1) + errorL(currCoeff) + 0.5], '*', 'Color', 'k')
        end
    end
    set(gca, 'tickdir', 'out')
    xlim([0.75 4.25])
    xticks([1:4])
    xticklabels([{'rwd stay'}, {'rwd switch'}, {'no rwd stay'}, {'no rwd switch'}])
    titleTxt = strrep([sessionName ' - ' cellList{currCell}], '_', ' ');
    suptitle(titleTxt)
   
   
    
%    spikeTmp = maxFRcs(s.responseInds);
%    figure; 
%    subplot(1,2,1)
%    errorbar([1 2], [mean(spikeTmp(rwdStay)) mean(spikeTmp(rwdSwitch))],...
%        [sem(spikeTmp(rwdStay)) sem(spikeTmp(rwdSwitch))], 'Color', [0 0.6 1], 'linewidth', 2);
%    set(gca, 'tickdir', 'out')
%    xticks([1 2])
%    xticklabels([{'rwd-stay'} {'rwd-shift'}])
%     xlim([0.75 2.25])
%    
%    subplot(1,2,2)
%    errorbar([1 2], [mean(spikeTmp(noRwdStay)) mean(spikeTmp(noRwdSwitch))],...
%        [sem(spikeTmp(noRwdStay)) sem(spikeTmp(noRwdSwitch))], 'Color', [1 0.1 0], 'linewidth', 2);
%    set(gca, 'tickdir', 'out')
%    xticks([1 2])
%    xticklabels([{'no-rwd-stay'} {'no-rwd-shift'}])
%    xlim([0.75 2.25])
    
    
    
end

end