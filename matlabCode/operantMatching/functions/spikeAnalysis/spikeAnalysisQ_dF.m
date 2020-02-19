function [mdlStruct] = spikeAnalysisQ_dF(sessionName, varargin)


p = inputParser;
% default parameters if none given
p.addParameter('intanFlag', 0);
p.addParameter('revForFlag', 0);
p.addParameter('biasFlag',0);
p.addParameter('cellName', ['all']);
p.parse(varargin{:});

cellName = p.Results.cellName;

% Path
[root,sep] = currComputer();

[animalName] = strtok(sessionName, 'd');
animalName = animalName(2:end);

if isstrprop(sessionName(end), 'alpha')
    sortedFolderLocation = [root animalName sep sessionName(1:end-1) sep 'sorted' sep 'session ' sessionName(end) sep];
else
    sortedFolderLocation = [root animalName sep sessionName sep 'sorted' sep 'session' sep];
end
sortedFolder = dir(sortedFolderLocation);


if p.Results.intanFlag
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
[s] = behAnalysisNoPlot_opMD(sessionName, 'revForFlag', p.Results.revForFlag);


Qmdls = qLearning_fitOpponency([sessionName '.asc'], 'revForFlag', p.Results.revForFlag);
rBar = Qmdls.fiveParams_opponency.rBar;
if p.Results.biasFlag
    sumQ = Qmdls.fiveParams_twoLearnRates_alphaForget_bias.Q(:,1) +  Qmdls.fiveParams_twoLearnRates_alphaForget_bias.Q(:,2);
    diffQ = Qmdls.fiveParams_twoLearnRates_alphaForget_bias.Q(:,1) -  Qmdls.fiveParams_twoLearnRates_alphaForget_bias.Q(:,2);
    pe = Qmdls.fiveParams_twoLearnRates_alphaForget_bias.pe;
else
    sumQ = Qmdls.fourParams_twoLearnRates_alphaForget.Q(:,1) +  Qmdls.fourParams_twoLearnRates_alphaForget.Q(:,2);
    diffQ = Qmdls.fourParams_twoLearnRates_alphaForget.Q(:,1) -  Qmdls.fourParams_twoLearnRates_alphaForget.Q(:,2);
    pe = Qmdls.fourParams_twoLearnRates_alphaForget.pe;
end
confQ = abs(diffQ);
choiceQ = NaN(length(s.allChoices), 1);
choiceQ(s.allChoices == 1) = Qmdls.fourParams_twoLearnRates_alphaForget.Q((s.allChoices == 1),1);
choiceQ(s.allChoices == -1) = Qmdls.fourParams_twoLearnRates_alphaForget.Q((s.allChoices == -1),2);



%% Sort all spikes into a raster-able matrix

%set time window for spike analyses
tb = 1.5;
tf = 5;
time = -1000*tb:1000*tf;

spikeFields = fields(sessionData);
if iscell(cellName)
    for i = 1:length(cellName)
        cellInd(i) = find(~cellfun(@isempty,strfind(spikeFields,cellName{i})));
    end
elseif regexp(cellName, 'all')
    cellInd = find(~cellfun(@isempty,strfind(spikeFields,'C_')) | ~cellfun(@isempty,strfind(spikeFields,'TT')));
else
    cellInd = find(~cellfun(@isempty,strfind(spikeFields,cellName)));
end
   
allTrial_spike = {};
for k = 1:length(sessionData)
    for i = 1:length(cellInd)
        if k == 1
            prevTrial_spike = [];
        else
            prevTrial_spikeInd = sessionData(k-1).(spikeFields{cellInd(i)}) > sessionData(k-1).trialEnd-tb*1000;
            prevTrial_spike = sessionData(k-1).(spikeFields{cellInd(i)})(prevTrial_spikeInd) - sessionData(k).CSon;
        end
        
        currTrial_spikeInd = sessionData(k).(spikeFields{cellInd(i)}) < sessionData(k).CSon+tf*1000;
        currTrial_spike = sessionData(k).(spikeFields{cellInd(i)})(currTrial_spikeInd) - sessionData(k).CSon;
        
        allTrial_spike{i,k} = [prevTrial_spike currTrial_spike];

    end
end

% sometimes no licks/spikes are considered 1x0 and sometimes they are []
% plotSpikeRaster does not place nicely with [] so this converts all empty indices to 1x0
allTrial_spike(cellfun(@isempty,allTrial_spike)) = {zeros(1,0)}; 


%% set time window and smoothing parameters, run analysis for all cells

smoothWin = 250;
trialBeg = tb*1000;
CSoff = tb*1000 + 500;
for i = 1:length(s.behSessionData)
    trialDurDiff(i) = (s.behSessionData(i).trialEnd - s.behSessionData(i).CSon)- tf*1000;
end

mdlStruct = struct;

for i = 1:length(cellInd)
    allTrial_spikeMatx = NaN(length(sessionData),length(time));
    for j = 1:length(allTrial_spike)
        tempSpike = allTrial_spike{i,j};
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
            preCSspikeCount(i,j) = sum(allTrial_spikeMatx(j, 1:trialBeg));              %find total spikes before CS on
            postCSspikeCount(i,j) = sum(allTrial_spikeMatx(j, trialBeg:CSoff));

            spikeTemp = fastsmooth(allTrial_spikeMatx(j,:)*1000, smoothWin, 3);         %smooth raw spikes to find features of spike rate
            maxFRcs(i,j) = max(spikeTemp(trialBeg:CSoff));
            minFRcs(i,j) = min(spikeTemp(trialBeg:CSoff));
            if ~isnan(maxFRcs(i,j))
                maxFRtime(i,j) = find(spikeTemp == max(spikeTemp(trialBeg:CSoff)), 1);
            else
                maxFRtime(i,j) = NaN;
            end

        else
            preCSspikeCount(i,j) = NaN;
            postCSspikeCount(i,j) = NaN;
            maxFRcs(i,j) = NaN;
            minFRcs(i,j) = NaN;
            maxFRtime(i,j) = NaN;
        end
    end
       
%% create models and put them into the structure
    
    if length(cellInd) > 1
        cellNameTemp = cellName{i};
    else
        cellNameTemp = cellName;
    end
    mdlStruct.(cellNameTemp).maxFRtrialQ = fitlm([sumQ diffQ confQ choiceQ rBar], maxFRcs(i,s.responseInds));
    mdlStruct.(cellNameTemp).preCSspikeCountQ = fitlm([sumQ diffQ confQ choiceQ rBar], preCSspikeCount(i,s.responseInds));
    

    
end

save([sortedFolderLocation sessionName '_spikeMdls.mat'], 'mdlStruct');


    