function [cellChoice, matChoice, matChoiceSlide, slideTime] = getUnitMatChoice(session, unit, tb, tf, stepSize, binSize)
[root, sep] = currComputer();
cellChoice = [];
matChoice = [];
matChoiceSlide = [];
slideTime = [];
time = -1000*tb:1000*tf;
midPoints = (0.5*binSize + 1):stepSize:(length(time)-0.5*binSize);
slideTime = midPoints - tb*1000;
% paths
pd = parseSessionString_df(session, root, sep);
neuralynxDataPath = [pd.sortedFolder session '_sessionData_nL.mat'];
load(neuralynxDataPath)
os = behAnalysisNoPlot_opMD(session, 'simpleFlag', 1);
spikeFields = fields(sessionData);
clust = find(contains(spikeFields,unit));
if length(os.behSessionData) ~= length(sessionData)
    fprintf([session ' realign \n', ]);
    return
end
% cell
allTrial_spike_choice = {};
for k = 1:length(os.responseInds)
    if  os.responseInds(k) == 1
        prevTrial_spike = [];
    else
        prevTrial_spikeInd = [sessionData(os.responseInds(k)-1).(spikeFields{clust})] > (sessionData(os.responseInds(k)).respondTime-tb*1000);
        prevTrial_spike = sessionData(os.responseInds(k)-1).(spikeFields{clust})(prevTrial_spikeInd) - sessionData(os.responseInds(k)).respondTime;
    end

    currTrial_spikeInd = sessionData(os.responseInds(k)).(spikeFields{clust}) < sessionData(os.responseInds(k)).respondTime+tf*1000 ... 
        & sessionData(os.responseInds(k)).(spikeFields{clust}) > sessionData(os.responseInds(k)).respondTime-tb*1000;
    currTrial_spike = sessionData(os.responseInds(k)).(spikeFields{clust})(currTrial_spikeInd) - sessionData(os.responseInds(k)).respondTime;

    allTrial_spike_choice{k} = [prevTrial_spike currTrial_spike];
end
allTrial_spike_choice(cellfun(@isempty,allTrial_spike_choice)) = {zeros(1,0)};
cellChoice = allTrial_spike_choice;
% mat
for j = 1:length(sessionData)
    trialDurDiff(j) = (sessionData(j).trialEnd - (sessionData(j).rewardTime - os.rwdDelay))- tf*1000;
end
trialDurDiff(end) = 0; 
allTrial_spikeMatx_choice = zeros(length(os.responseInds),length(time));         
for j = 1:length(allTrial_spike_choice)
    tempSpike = allTrial_spike_choice{j};
    tempSpike = tempSpike + tb*1000; % add this to pad time for SDF
    if any(tempSpike == 0)
        tempSpike(tempSpike == 0) = 1;
    end
    allTrial_spikeMatx_choice(j,tempSpike) = 1;
    if trialDurDiff(j) < 0
        allTrial_spikeMatx_choice(j, isnan(allTrial_spikeMatx_choice(j, 1:end+trialDurDiff(j)))) = 0;  %converts within trial duration NaNs to 0's
    else
        allTrial_spikeMatx_choice(j, isnan(allTrial_spikeMatx_choice(j,:))) = 0;
    end
end
matChoice = allTrial_spikeMatx_choice;
% slide window
allTrial_spikeMatx_slide = zeros(length(os.responseInds), length(midPoints));
for w = 1:length(midPoints)
    allTrial_spikeMatx_slide(:,w) = ...
        nansum(allTrial_spikeMatx_choice(:,midPoints(w)-0.5*binSize:midPoints(w)+0.5*binSize-1),2)*1000/binSize;
end
matChoiceSlide = allTrial_spikeMatx_slide;
end