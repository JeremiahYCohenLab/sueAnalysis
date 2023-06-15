dayList = getDayList('allDBh-cre', 'all-DBh', 'allWithManip');
% remove duplicate
dayListNew = {};
for i = 1:length(dayList)
    if sum(cellfun(@(x) strcmp(x, dayList{i}), dayListNew))==0
        dayListNew = [dayListNew; dayList{i}];
    end
end
dayList = dayListNew;
%% mean sessions per animal
aniList = cellfun(@(x) x(2:6), dayList, 'UniformOutput', false);
aniList = unique(aniList);

sessionNum = zeros(size(aniList));

for i = 1:length(aniList)
    expression = ['\w*' aniList{i} '\w*'];
    aniSessInd = cellfun(@(x) regexp(x,expression), dayList, 'UniformOutput', false);
    sessionNum(i) = sum(cellfun(@(x) ~isempty(x), aniSessInd));
end

figure2;
histogram(sessionNum, 30, 'FaceColor', [0.4 0.4 0.4], 'EdgeColor', 'none');
xlabel('No. sessions/animal')

%%
sessionLen = zeros(size(dayList));
sessionTime = zeros(size(dayList));
itiCombine = [];
lickLatZStay = NaN(size(dayList));
lickLatZSwitch = NaN(size(dayList));
for i = 1:length(dayList)
    s = behAnalysisNoPlot_opMD(dayList{i}, 'simpleFlag', 1);
    timeBtwn = [s.behSessionData(2:end).CSon] - [s.behSessionData(1:end-1).CSon] - 2800 - s.rwdDelay;
    itiCombine = [itiCombine timeBtwn];
    
    sessionLen(i) = length(s.allChoices);
    sessionTime(i) = 1/1000*1/60*(s.behSessionData(end).CSon - s.behSessionData(1).CSon);
    
    lickLatZStay(i) = mean(s.lickLat(s.stayChoice_Inds));
    lickLatZSwitch(i) = mean(s.lickLat(s.changeChoice_Inds));
end
%%
figure2;
histogram(sessionLen, 20, 'FaceColor', [0.4 0.4 0.4], 'EdgeColor', 'none');
xlabel('No. trials/session')

figure2;
histogram(sessionTime, 20, 'FaceColor', [0.4 0.4 0.4], 'EdgeColor', 'none');
xlabel('Time/session')

figure2; hold on;
scatter(lickLatZStay, lickLatZSwitch, 15, 'k');
plot([50 550], [50 550], 'LineStyle', '--', 'Color', 'r');
xlabel('stay (ms)');
ylabel('switch (ms)')
set(gca, 'TickDir', 'out');
%% parameter distribution
[root, sep] = currComputer();
[~, dayList, ~] = xlsread([root 'aniModel.xlsx'], 'all');
allAnis = dayList(2:end, 1);
allCols = dayList(2:end, 2);
allFile = dayList(2:end, 3);
allSheet = dayList(2:end, 4);
% check model fitting
modelName = '5params';
for i = 1:length(allAnis)
    animalName = allAnis{i};
    category = allCols{i};
    sampFile = [animalName category '_', modelName];
    path = [root animalName sep animalName 'sorted' sep 'stan' sep 'bernoulli' sep modelName sep category sep];
    modelPath = [path sampFile '.mat'];

    if ~exist(modelPath, 'file')
        fprintf([animalName ' ' category, ' not fitted yet with ' modelName '\n'])
        stan_qLearningFit(allFile{i}, allAnis{i}, allCols{i}, 'modelName', modelName, 'iter', 10000);
    end
end
%% calculate MAP
paramNames = getParamNames_dF(modelName, 0);
aniParams = NaN(length(allAnis), length(paramNames));
aniMeanLL = NaN(length(allAnis), 1);
numBin = 50;
for i = 1:length(allAnis)
    animalName = allAnis{i};
    category = allCols{i};
    sampFile = [animalName category '_', modelName];
    path = [root animalName sep animalName 'sorted' sep 'stan' sep 'bernoulli' sep modelName sep category sep];
    modelPath = [path sampFile '.mat'];
    load(modelPath);
    eval(['samps = ' sampFile]);
    eval(['clear ' sampFile]);
    % get map of all params
    for j = 1:length(paramNames)
        tmp = samps.(['mu_' paramNames{j}]);
        noDInd = samps.divergent__~= 1;
        tmp = tmp(noDInd);
        [counts, edges] = histcounts(tmp, 50);
        [~, maxInd] = max(counts);
        if maxInd<length(counts)
            aniParams(i, j) = mean(tmp(tmp>=edges(maxInd)&tmp<edges(maxInd+1)));
        else
            aniParams(i, j) = mean(tmp(tmp>=edges(maxInd)&tmp<=edges(maxInd+1)));
        end
    end
end
%%
colors = cool(length(paramNames));
figure2;
for j = 1:length(paramNames)
    subplot(1,length(paramNames), j); hold on;
    histogram(aniParams(:,j), 10, 'FaceColor', colors(j,:), 'EdgeColor', 'none');
    title(paramNames{j})
end
%% session by session trialLL
[root, sep] = currComputer();
[~, dayList, ~] = xlsread([root 'aniModel.xlsx'], 'all');
allAnis = dayList(2:end, 1);
allCols = dayList(2:end, 2);
allFile = dayList(2:end, 3);
allSheet = dayList(2:end, 4);
dayList = getDayList('allDBh-cre', 'all-DBh', 'allWithManip');
%%
trialLL = NaN(size(dayList));
trialLength = NaN(size(dayList));
probChoice = cell(size(dayList));
clear allFit
[root, sep] = currComputer();
modelName = '5params';
numSamps = 2000;
for i = 1:length(dayList)
    pd = parseSessionString_df(dayList{i}, root, sep);
    aniInd = strcmp(allAnis, pd.animalName);
    col = allCols{aniInd};
    params = getStanModelParams_sampsOnly(pd.animalName, col, modelName, numSamps, 'sessionName', dayList{i}, 'sessionParamsFlag', 1);
    if isnan(params)
        continue
    end
    t = inferModelVar(dayList{i}, params, modelName, 'perturb', 0);
    s = behAnalysisNoPlot_opMD(dayList{i}, 'simpleFlag', 1);
    trialLL(i) = t.LH;
    trialLength(i) = length(s.allChoices);
    probChoice{i} = t.probChoice;
    allFit(i) = t; 
end
%%
figure2;
histogram(trialLL./trialLength)
%% animal level
% dayList = getDayList('allDBh-cre', 'all-DBh', 'allWithManip');
correct = cellfun(@(x) x>=0.5, probChoice, 'UniformOutput', false);
correct = cellfun(@sum, correct);
aniTrialLL = NaN(size(allAnis));
aniCorrect = NaN(size(allAnis));
for i = 1:length(allAnis)
    currInd = contains(dayList, allAnis{i});
    aniTrialLL(i) = sum(trialLL(currInd))/sum(trialLength(currInd));
    aniCorrect(i) = sum(correct(currInd))/sum(trialLength(currInd));
end
%%
figure2;
histogram(aniTrialLL, 10, 'FaceColor', [0.4 0.4 0.4], 'EdgeColor', 'none')
title('aniTrialLL')
set(gca, 'Box', 'off');
set(gca, 'TickDir', 'out');
%% animal level correct
figure2;
histogram(aniCorrect, 10, 'FaceColor', [0.4 0.4 0.4], 'EdgeColor', 'none')
title('aniCorrect')
set(gca, 'Box', 'off');
set(gca, 'TickDir', 'out');
%% model fitting
[root, sep] = currComputer();
[~, dayList, ~] = xlsread([root 'aniModel.xlsx'], 'all');
allAnis = dayList(2:end, 1);
allCols = dayList(2:end, 2);
allFile = dayList(2:end, 3);
allSheet = dayList(2:end, 4);
% check model fitting
modelName = '5params_k_bias';
for i = 1:length(allAnis)
    animalName = allAnis{i};
    category = allCols{i};
    sampFile = [animalName category '_', modelName];
    path = [root animalName sep animalName 'sorted' sep 'stan' sep 'bernoulli' sep modelName sep category sep];
    modelPath = [path sampFile '.mat'];

%     if ~exist(modelPath, 'file')
        fprintf([animalName ' ' category, ' not fitted yet with ' modelName '\n'])
        stan_qLearningFit(allFile{i}, allAnis{i}, allCols{i}, 'modelName', modelName, 'iter', 10000);
%     end
end
%% model fitting
allAnis = {'ZS059', 'ZS060', 'ZS061', 'ZS062'};
col = 'good';
modelName = '5params_k_bias';
[root, sep] = currComputer();
for i = 1:length(aniNames)
    animalName = allAnis{i};
    category = col;
    sampFile = [animalName category '_', modelName];
    path = [root animalName sep animalName 'sorted' sep 'stan' sep 'bernoulli' sep modelName sep category sep];
    modelPath = [path sampFile '.mat'];

    fprintf([animalName ' ' category, ' not fitted yet with ' modelName '\n'])
    stan_qLearningFit(allAnis{i}, allAnis{i}, col, 'modelName', modelName, 'iter', 10000);   
end
%%