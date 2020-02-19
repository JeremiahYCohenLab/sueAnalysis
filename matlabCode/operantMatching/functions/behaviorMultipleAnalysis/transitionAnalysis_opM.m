function [transChoiceMatx, mdlFit] = transitionAnalysis_opM(xlFile, sheet, category, probs)

if nargin < 4
    pHigh = 90;
else
    pHigh = probs(1);
end

[root, sep] = currComputer();
transChoiceMatx = [];
range = 20;

[~, dayList, ~] = xlsread(xlFile, sheet);
[~,col] = find(~cellfun(@isempty,strfind(dayList, category)) == 1);
dayList = dayList(2:end,col);
endInd = find(cellfun(@isempty,dayList),1);
if ~isempty(endInd)
    dayList = dayList(1:endInd-1,:);
end

for i = 1: length(dayList)
    sessionName = dayList{i};
    [animalName, date] = strtok(sessionName, 'd'); 
    animalName = animalName(2:end);
    date = date(1:9);
    sessionFolder = ['m' animalName date];

    if isstrprop(sessionName(end), 'alpha')
        sessionDataPath = [root animalName sep sessionFolder sep 'sorted' sep 'session ' sessionName(end) sep sessionName '_sessionData_behav.mat'];
    else
        sessionDataPath = [root animalName sep sessionFolder sep 'sorted' sep 'session' sep sessionName '_sessionData_behav.mat'];
    end

    if exist(sessionDataPath,'file')
        delete(sessionDataPath)
    end
    
    if exist(sessionDataPath,'file')
        load(sessionDataPath)
    else
        [sessionData, blockSwitch, blockProbs] = generateSessionData_behav_operantMatching(sessionName);
    end

    responseInds = find(~isnan([sessionData.rewardTime])); % find CS+ trials with a response in the lick window
    omitInds = isnan([sessionData.rewardTime]); 

    tempBlockSwitch = blockSwitch;
    for i = 2:length(blockSwitch)
        subVal = sum(omitInds(tempBlockSwitch(i-1):tempBlockSwitch(i)));
        blockSwitch(i:end) = blockSwitch(i:end) - subVal;
    end

    allReward_R = [sessionData(responseInds).rewardR]; 
    allReward_L = [sessionData(responseInds).rewardL]; 
    allChoices = NaN(1,length(sessionData(responseInds)));
    allChoices(~isnan(allReward_R)) = 1;
    allChoices(~isnan(allReward_L)) = -1;

    for j = 2:(length(blockSwitch) - 1)
        tmpInd = blockSwitch(j);
        if str2double(blockProbs{j}(1:2)) == pHigh
                if (tmpInd - range - 1) > 0 & length(allChoices) > (tmpInd + range)
                    transChoiceMatx = [transChoiceMatx; allChoices((tmpInd-range+1):(tmpInd+range))];
                end
        elseif str2double(blockProbs{j}(4:5)) == pHigh
                if (tmpInd - range - 1) > 0 & length(allChoices) > (tmpInd + range)
                    transChoiceMatx = [transChoiceMatx; (allChoices((tmpInd-range+1):(tmpInd+range)))*-1];
                end
        end   
    end
end   


choiceAvg = mean(transChoiceMatx,1);

figure; hold on
x = [-range+1:range];
plot(x,choiceAvg, 'r')
title(['Choice at block transitions'])
ylabel('Choice average')
xlabel('Trials from switch')

xx = [1:range+1];
mdlFit = singleExpFitInt(xx,choiceAvg(range:end));
expConv = mdlFit.a*exp(-(mdlFit.b)*(1:range+1))+mdlFit.c;
plot([0:range], expConv, '--r')
linetype = 'k';
vline(0, linetype)
ylim([-1 1])
legend('actual', 'exp fit')
    