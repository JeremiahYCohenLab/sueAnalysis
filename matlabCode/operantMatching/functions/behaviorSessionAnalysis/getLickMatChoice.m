function [cellChoice, matChoice, matChoiceSlide, slideTime] = getLickMatChoice(session, tb, tf, stepSize, binSize)
[root, sep] = currComputer();
cellChoice = [];
matChoice = [];
matChoiceSlide = [];
slideTime = [];
time = -1000*tb:1000*tf;
midPoints = (0.5*binSize + 1):stepSize:(length(time)-0.5*binSize);
slideTime = midPoints - tb*1000;
% paths
os = behAnalysisNoPlot_opMD(session, 'simpleFlag', 1);
sessionData = os.behSessionData;

% cell
allTrial_lick_choice = {};
for k = 1:length(os.responseInds)
    if  os.responseInds(k) == 1
        prevTrial_lick = [];
    else
        if os.allChoices(k)>0
            prevLicksTemp = sessionData(os.responseInds(k)-1).licksR;
        else
            prevLicksTemp = sessionData(os.responseInds(k)-1).licksL;
        end
        prevTrial_lickInd = prevLicksTemp > (sessionData(os.responseInds(k)).respondTime-tb*1000);
        prevTrial_lick = prevLicksTemp(prevTrial_lickInd) - sessionData(os.responseInds(k)).respondTime;
    end
    if os.allChoices(k)>0
        currLicksTemp = sessionData(os.responseInds(k)).licksR;
    else
        currLicksTemp = sessionData(os.responseInds(k)).licksL;
    end
    currTrial_spikeInd = currLicksTemp < sessionData(os.responseInds(k)).respondTime+tf*1000 ... 
        & currLicksTemp > sessionData(os.responseInds(k)).respondTime-tb*1000;
    currTrial_lick = currLicksTemp(currTrial_spikeInd) - sessionData(os.responseInds(k)).respondTime;

    allTrial_lick_choice{k} = [prevTrial_lick currTrial_lick];
end
allTrial_lick_choice(cellfun(@isempty,allTrial_lick_choice)) = {zeros(1,0)};
cellChoice = allTrial_lick_choice;
% mat
for j = 1:length(sessionData)
    trialDurDiff(j) = (sessionData(j).trialEnd - (sessionData(j).rewardTime - os.rwdDelay))- tf*1000;
end
trialDurDiff(end) = 0; 
allTrial_lickMatx_choice = zeros(length(os.responseInds),length(time));         
for j = 1:length(allTrial_lick_choice)
    tempSpike = allTrial_lick_choice{j};
    tempSpike = tempSpike + tb*1000; % add this to pad time for SDF
    if any(tempSpike == 0)
        tempSpike(tempSpike == 0) = 1;
    end
    allTrial_lickMatx_choice(j,tempSpike) = 1;
    if trialDurDiff(j) < 0
        allTrial_lickMatx_choice(j, isnan(allTrial_lickMatx_choice(j, 1:end+trialDurDiff(j)))) = 0;  %converts within trial duration NaNs to 0's
    else
        allTrial_lickMatx_choice(j, isnan(allTrial_lickMatx_choice(j,:))) = 0;
    end
end
matChoice = allTrial_lickMatx_choice;
% slide window
allTrial_lickMatx_slide = zeros(length(os.responseInds), length(midPoints));
for w = 1:length(midPoints)
    allTrial_lickMatx_slide(:,w) = ...
        nansum(allTrial_lickMatx_choice(:,midPoints(w)-0.5*binSize:midPoints(w)+0.5*binSize-1),2)*1000/binSize;
end
matChoiceSlide = allTrial_lickMatx_slide;
end