function [cellCue, matCue, matCueSlide, slideTime] = getUnitMatCue(session, unit, tb, tf, stepSize, binSize)
[root, sep] = currComputer();
time = -1000*tb:1000*tf;
midPoints = (0.5*binSize + 1):stepSize:(length(time)-0.5*binSize);
slideTime = midPoints - tb*1000;
% paths
pd = parseSessionString_df(session, root, sep);
neuralynxDataPath = [pd.sortedFolder session '_sessionData_nL.mat'];
load(neuralynxDataPath)
spikeFields = fields(sessionData);
clust = find(contains(spikeFields,unit));
% cell
allTrial_spike_choice = {};
for k = 1:length(sessionData)
    if k == 1
        prevTrial_spike = [];
    else
        prevTrial_spikeInd = [sessionData(k-1).(spikeFields{clust})] > (sessionData(k).CSon-tb*1000);
        prevTrial_spike = sessionData(k-1).(spikeFields{clust})(prevTrial_spikeInd) - sessionData(k).CSon;
    end

    currTrial_spikeInd = sessionData(k).(spikeFields{clust}) < sessionData(k).CSon+tf*1000 ... 
        & sessionData(k).(spikeFields{clust}) > sessionData(k).CSon-tb*1000;
    currTrial_spike = sessionData(k).(spikeFields{clust})(currTrial_spikeInd) - sessionData(k).CSon;

    allTrial_spike_choice{k} = [prevTrial_spike currTrial_spike];
end
allTrial_spike_choice(cellfun(@isempty,allTrial_spike_choice)) = {zeros(1,0)};
cellCue = allTrial_spike_choice;
% mat
for j = 1:length(sessionData)
    trialDurDiff(j) = (sessionData(j).trialEnd - sessionData(j).CSon)- tf*1000;
end
trialDurDiff(end) = 0; 
allTrial_spikeMatx_choice = zeros(length(sessionData),length(time));         
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
matCue = allTrial_spikeMatx_choice;
% slide window
allTrial_spikeMatx_slide = zeros(length(sessionData), length(midPoints));
for w = 1:length(midPoints)
    allTrial_spikeMatx_slide(:,w) = ...
        nansum(allTrial_spikeMatx_choice(:,midPoints(w)-0.5*binSize:midPoints(w)+0.5*binSize-1),2)*1000/binSize;
end
matCueSlide = allTrial_spikeMatx_slide;
end