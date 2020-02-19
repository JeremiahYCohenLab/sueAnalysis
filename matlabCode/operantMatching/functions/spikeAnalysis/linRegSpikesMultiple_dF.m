function [mdlStruct] = spikeAnalysisPrevRwd_dF(xlFile, sheet)


% path info
[root,sep] = currComputer();

% load data from excel sheet
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

%set time window for spike analyses
tb = 1.5;
tf = 5;
time = -1000*tb:1000*tf;
trialBeg = tb*1000;
CSoff = tb*1000 + 500;
smoothWin = 250;


sessionName = [];

for currCell = 1:length(cellList)
    
    fprintf('Analyzing cell %d of %d \n', currCell, length(sessionList));
    
    if strcmp(sessionName, sessionList{currCell}) == 0
        sessionName = sessionList{currCell};
        [animalName, date] = strtok(sessionName, 'd'); 
        animalName = animalName(2:end);

        [animalName] = strtok(sessionName, 'd');
        animalName = animalName(2:end);

        if isstrprop(sessionName(end), 'alpha')
            sortedFolderLocation = [root animalName sep sessionName(1:end-1) sep 'sorted' sep 'session ' sessionName(end) sep];
        else
            sortedFolderLocation = [root animalName sep sessionName sep 'sorted' sep 'session' sep];
        end
        sortedFolder = dir(sortedFolderLocation);


        if intanFlag(currCell)
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
        [s] = behAnalysisNoPlot_opMD(sessionName, 'revForFlag', revForFlag(currCell));
    end


    %% Sort all spikes into a raster-able matrix

    spikeFields = fields(sessionData);
    cellInd = find(~cellfun(@isempty,strfind(spikeFields,cellList(currCell))));

    allTrial_spike = {};
    for k = 1:length(sessionData)
        if k == 1
            prevTrial_spike = [];
        else
            prevTrial_spikeInd = sessionData(k-1).(spikeFields{cellInd}) > sessionData(k-1).trialEnd-tb*1000;
            prevTrial_spike = sessionData(k-1).(spikeFields{cellInd})(prevTrial_spikeInd) - sessionData(k).CSon;
        end

        currTrial_spikeInd = sessionData(k).(spikeFields{cellInd}) < sessionData(k).CSon+tf*1000;
        currTrial_spike = sessionData(k).(spikeFields{cellInd})(currTrial_spikeInd) - sessionData(k).CSon;

        allTrial_spike{k} = [prevTrial_spike currTrial_spike];
    end

    % sometimes no licks/spikes are considered 1x0 and sometimes they are []
    % plotSpikeRaster does not place nicely with [] so this converts all empty indices to 1x0
    allTrial_spike(cellfun(@isempty,allTrial_spike)) = {zeros(1,0)}; 


    %% set time window and smoothing parameters, run analysis for all cells


    for i = 1:length(s.behSessionData)
        trialDurDiff(i) = (s.behSessionData(i).trialEnd - s.behSessionData(i).CSon)- tf*1000;
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
    for j = 1:length(allTrial_spike)
        if ~isempty(allTrial_spikeMatx(i,j))
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

%% create models and put them into the structure
    mdl_maxCS = fitlm(s.rwdTimeMatxBin', maxFRcs(s.responseInds));
    mdl_preCS = fitlm(s.rwdTimeMatxBin', preCSspikeCount(s.responseInds));
    
    rSqr_maxCS(currCell,:) = [mdl_maxCS.Rsquared.Ordinary mdl_maxCS.Rsquared.Adjusted];
    rSqr_preCS(currCell,:) = [mdl_preCS.Rsquared.Ordinary mdl_preCS.Rsquared.Adjusted];
    prevRwd_maxCS_z(currCell,:) = mdl_maxCS.Coefficients.tStat(2);
    prevRwd_preCS_z(currCell,:) = mdl_preCS.Coefficients.tStat(2);
    
    xVals = [1:s.timeBinSize/1000:s.timeMax/1000];
    xVals = xVals(1:end-1) + diff(xVals)/2;
    if mdl_maxCS.Coefficients.Estimate(1) < 0
        expFitTime_maxCS = singleExpFit(-mdl_maxCS.Coefficients.Estimate(2:end), xVals');
    else
        expFitTime_maxCS = singleExpFit(mdl_maxCS.Coefficients.Estimate(2:end), xVals');
    end
    tau_maxCS(currCell) = expFitTime_maxCS.b;
    
    if mdl_maxCS.Coefficients.Estimate(1) < 0
        expFitTime_preCS = singleExpFit(-mdl_preCS.Coefficients.Estimate(2:end), xVals');
    else
        expFitTime_preCS = singleExpFit(mdl_preCS.Coefficients.Estimate(2:end), xVals');
    end
    tau_preCS(currCell) = expFitTime_maxCS.b;
end


figure;
subplot(2,2,1); hold on
tmpInds_ns = find(abs(prevRwd_maxCS_z) < 1.96 & abs(prevRwd_preCS_z) < 1.96);
scatter(prevRwd_maxCS_z(tmpInds_ns), prevRwd_preCS_z(tmpInds_ns), [], 'k', 'filled')
tmpInds_s = [1:length(prevRwd_maxCS_z)];
tmpInds_s(tmpInds_ns) = [];
scatter(prevRwd_maxCS_z(tmpInds_s), prevRwd_preCS_z(tmpInds_s), [], 'c', 'filled')
xlabel('Max FR Z Value')
ylabel('Pre-CS Z Value')

subplot(2,2,2)
scatter(rSqr_maxCS(:,1), rSqr_preCS(:,1), 'm', 'filled')
xlabel('Max FR R^2 Value')
ylabel('Pre-CS R^2 Value')

subplot(2,2,3)
tmpInds_s = find(abs(prevRwd_maxCS_z) > 1.96);
histogram(tau_maxCS(tau_maxCS(tmpInds_s)< 100), 30, 'FaceColor', 'k', 'FaceAlpha', 1, 'Normalization', 'probability')
ylabel('Probability')
xlabel('\tau of exponential fit to \beta coefficients')
title('Max FR during CS')


subplot(2,2,4)
tmpInds_s = find(abs(prevRwd_preCS_z) > 1.96);
histogram(tau_preCS(tau_preCS(tmpInds_s) < 100), 30, 'FaceColor', 'k',  'FaceAlpha', 1,'Normalization', 'probability')
ylabel('Probability')
xlabel('\tau of exponential fit to \beta coefficients')
title('FR Pre-CS')

for i = 1:4
    subplot(2,2,i)
    set(gca, 'TickDir', 'out')
    set(gca, 'Box', 'off')
end

suptitle('Linear regression: Rewards in time on firing rates')

    