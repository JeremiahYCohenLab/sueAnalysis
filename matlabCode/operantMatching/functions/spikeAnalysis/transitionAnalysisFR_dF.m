function [transChoiceMatxLow, transChoiceMatxHigh, mdlFitLow, mdlFitHigh] = transitionAnalysisFR_dF(xlFile, sheet, spikeFeatName)

if nargin < 3
    spikeFeatName = ['postCSspikeCount'];
end

[root, sep] = currComputer();

%initialize choice and FR matrices, set range
transChoiceMatxLow = []; 
transChoiceMatxHigh = [];
transFRmatxHigh = [];
transFRmatxLow = [];
range = 20;

%set time window for spike analyses
tb = 1.5;
tf = 5;
time = -1000*tb:1000*tf;
trialBeg = tb*1000;
CSoff = tb*1000 + 500;
smoothWin = 250;

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
    
    spikeFeat = eval(spikeFeatName);
    
    %% analyze choice behavior and spike feature at transitions  

    rwdProb_R = [s.behSessionData(s.responseInds).rewardProbR]; 
    rwdProb_L = [s.behSessionData(s.responseInds).rewardProbL]; 

%     rwdProbs = sort(unique([rwdProb_R rwdProb_L]));
%     pHigh = rwdProbs(3);
%     pLow = rwdProbs(2);
%     pTrans = rwdProbs(3);
    if max(([rwdProb_R rwdProb_L])) == 70 || max(([rwdProb_R rwdProb_L])) == 40
        pHigh = 70;
        pLow = 40;
        pTrans = 70;
    else
        pHigh = 90;
        pLow = 50;
        pTrans = 90;
    end
    
    transFRmatxHighTmp = [];
    transFRmatxLowTmp = [];
    for j = 2:(length(s.blockSwitch) - 1)
        tmpInd = s.blockSwitch(j);
        if rwdProb_R(tmpInd) == pHigh & rwdProb_L(tmpInd) == 10 & rwdProb_R(tmpInd+1) == 10 & rwdProb_L(tmpInd+1) == pTrans
                if (tmpInd - range - 1) > 0 & length(s.allChoices) > (tmpInd + range)
                    transChoiceMatxHigh = [transChoiceMatxHigh; s.allChoices((tmpInd-range+1):(tmpInd+range))];
                    transFRmatxHighTmp = [transFRmatxHighTmp; spikeFeat((tmpInd-range+1):(tmpInd+range))];
                end
        elseif rwdProb_R(tmpInd) == pLow & rwdProb_L(tmpInd) == 10 & rwdProb_R(tmpInd+1) == 10 & rwdProb_L(tmpInd+1) == pTrans
                if (tmpInd - range - 1) > 0 & length(s.allChoices) > (tmpInd + range)
                    transChoiceMatxLow = [transChoiceMatxLow; s.allChoices((tmpInd-range+1):(tmpInd+range))];
                    transFRmatxLowTmp = [transFRmatxLowTmp; spikeFeat((tmpInd-range+1):(tmpInd+range))];
                end
        elseif rwdProb_L(tmpInd) == pHigh & rwdProb_R(tmpInd) == 10 & rwdProb_L(tmpInd+1) == 10 & rwdProb_R(tmpInd+1) == pTrans
                if (tmpInd - range - 1) > 0 & length(s.allChoices) > (tmpInd + range)
                    transChoiceMatxHigh = [transChoiceMatxHigh; (s.allChoices((tmpInd-range+1):(tmpInd+range))*-1)];
                    transFRmatxHighTmp = [transFRmatxHighTmp; spikeFeat((tmpInd-range+1):(tmpInd+range))];
                end 
        elseif rwdProb_L(tmpInd) == pLow & rwdProb_R(tmpInd) == 10 & rwdProb_L(tmpInd+1) == 10 & rwdProb_R(tmpInd+1) == pTrans
                if (tmpInd - range - 1) > 0 & length(s.allChoices) > (tmpInd + range)
                    transChoiceMatxLow = [transChoiceMatxLow; (s.allChoices((tmpInd-range+1):(tmpInd+range))*-1)];
                    transFRmatxLowTmp = [transFRmatxLowTmp; spikeFeat((tmpInd-range+1):(tmpInd+range))];
                end
        end   
    end
    
    if ~isempty(transFRmatxHighTmp)
        if size(transFRmatxHighTmp,1) > 1
            transFRmatxHighTmp = zscore(mean(transFRmatxHighTmp,1));
        else
            transFRmatxHighTmp = zscore(transFRmatxHighTmp);
        end
        if mean(transFRmatxHighTmp(range:end)) > mean(transFRmatxHighTmp(1:range-1))
            transFRmatxHighTmp = -transFRmatxHighTmp;
        end
        transFRmatxHigh = [transFRmatxHigh; transFRmatxHighTmp];
    end
    if ~isempty(transFRmatxLowTmp)
        if size(transFRmatxLowTmp,1) > 1
            transFRmatxLowTmp = zscore(mean(transFRmatxLowTmp,1));
        else
            transFRmatxLowTmp = zscore(transFRmatxLowTmp);
        end
        if mean(transFRmatxLowTmp(range:end)) > mean(transFRmatxLowTmp(1:range-1))
            transFRmatxLowTmp = -transFRmatxLowTmp;
        end
        transFRmatxLow = [transFRmatxLow; transFRmatxLowTmp];
    end
    
end   


lowAvg = mean(transChoiceMatxLow,1);
highAvg = mean(transChoiceMatxHigh,1);

figure; 
set(gcf, 'renderer', 'painters')
suptitle([ 'Choice at block transitions'])
subplot(1,3,1); hold on
x = [-range+1:range];
plot(x,lowAvg, 'r')
plot(x,highAvg, 'b')
ylabel('Choice average')
xlabel('Trials from switch')
legend('medium -> low', 'high -> low')
plot([0 0], [-1 1], '--k')

xx = [1:range+1];
mdlFitLow = singleExpFitInt(xx,lowAvg(range:end));
mdlFitHigh = singleExpFitInt(xx,highAvg(range:end));
expConvLow = mdlFitLow.a*exp(-(mdlFitLow.b)*(1:range+1))+mdlFitLow.c;
expConvHigh = mdlFitHigh.a*exp(-(mdlFitHigh.b)*(1:range+1))+mdlFitHigh.c;
plot([0:range], expConvLow, '--r'); plot([0:range], expConvHigh,'--b')
ylim([-1 1])
set(gca, 'tickdir', 'out')

subplot(1,3,2); hold on
plot(x,[highAvg - lowAvg], 'k')
ylabel('Choice average difference')
xlabel('Trials from switch')
y = get(gca, 'ylim');
plot([0 0], y, '--k');
legend('high - medium');
ylim([y]);
set(gca, 'tickdir', 'out')

subplot(1,3,3); hold on
plotFilled(x, transFRmatxLow, [1 0 0])
plotFilled(x, transFRmatxHigh, [0 0 1])
ylabel('Z-scored FR')
xlabel('Trials from switch')
title(spikeFeatName)
y = get(gca, 'ylim');
plot([0 0], y, '--k');
ylim([y]);
set(gca, 'tickdir', 'out')


    