load D:\computerDataBackUp\tmpData\catWithOutcome.mat
ind = cats;
[root, sep] = currComputer();
%% 
stepSize = 100;
binSize = 200;
tb = 2;
tf = 5;
modelName = '5params';
category = 'good';

prevSession = [];
for i = 1:length(allUnits)
    session = allSessions{i};
    unit = allUnits{i};
    % choice aligned
    currStruct = struct;
    [cellChoice, matChoice, matChoiceSlide, slideTime] = getUnitMatChoice(session, unit, tb, tf, stepSize, binSize);
     currStruct.matChoice = matChoice;
%     % cue aligned
%     [cellCue, matCue, matCueSlide] = getUnitMatCue(session, unit, tb, tf, stepSize, binSize);
%     allData.([session unit]).matCue = matCue;
    % beh
    pd = parseSessionString_df(session, root, sep);
    neuralynxDataPath = [pd.sortedFolder session '_sessionData_nL.mat'];
    load(neuralynxDataPath)
    currStruct.sessionData = sessionData;
    % model
    if ~strcmp(session, prevSession)
        sampFile = [pd.animalName category '_', modelName];
        path = [root pd.animalName sep pd.animalName 'sorted' sep 'stan' sep 'bernoulli' sep modelName sep category sep];
        [t,~,noSession] = getStanModelParams_samps(modelName, [path sampFile '.mat'], 4000, 'sessionName', session);
        prevSession = session;
        currStruct.model = t;
    else
        currStruct.model = t;
    end
    eval([session unit '= currStruct';]);
end
%%

%%