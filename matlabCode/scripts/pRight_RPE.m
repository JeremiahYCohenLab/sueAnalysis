%% get all animal lists and dayLists
[root, sep] = currComputer();
[~, dayList, ~] = xlsread([root 'aniModel.xlsx'], 'all');
allAnis = dayList(2:end, 1);
allCols = dayList(2:end, 2);
allFile = dayList(2:end, 3);
allSheet = dayList(2:end, 4);
dayList = getDayList('allDBh-cre', 'all-DBh', 'allWithManip');
% remove duplicate
dayListNew = {};
for i = 1:length(dayList)
    if sum(cellfun(@(x) strcmp(x, dayList{i}), dayListNew))==0
        dayListNew = [dayListNew; dayList{i}];
    end
end
dayList = dayListNew;

%%
modelName = '5params';
sampNum = 2000;
allPe = {};
allChoices = {};
allPredR = {};

peCombined = [];
choiceCombined = [];
predRCombined = [];

for i = 1:length(dayList)
    session = dayList{i};
    pd = parseSessionString_df(session, root, sep);
    colInd = strcmp(allAnis, pd.aniName);
    category = allCols{colInd};
    params = getStanModelParams_sampsOnly(pd.animalName, category, modelName, sampNum, 'sessionName', session, 'biasFlag', 1, 'sessionParamsFlag', 1);
    t = inferModelVar(session, params, modelName, 'perturb', 0);
    s = behAnalysisNoPlot_opMD(session, 'simpleFlag', 1);
    pRight = t.probChoice;
    pRight(s.allChoices<0) = 1 - pRight(s.allChoices<0);
    allPe{i} = t.pe;
    allChoices{i} = 0.5*(s.allChoices + 1)';
    allPredR{i} = pRight;
    allT(i) = t;
    
    peCombined = [peCombined; t.pe];
    choiceCombined = [choiceCombined; 0.5*(s.allChoices + 1)'];
    predRCombined = [predRCombined; pRight];
 
end
%% get all RPEs
load F:\tmpData\allSessionModel.mat
%% plot predicted P(right)-rpe
rightRpe = peCombined;
rightRpe(choiceCombined==0) = -rightRpe(choiceCombined==0);

pRightChange = predRCombined(2:end) - predRCombined(1:end-1);
pRightChange = [pRightChange; NaN];

figure2;
scatter(rightRpe,pRightChange, 3, 'k', 'filled');
title('combined')
figure2Wide;
subplot(1, 2, 1)
scatter(peCombined(choiceCombined==0), pRightChange(choiceCombined==0), 3, 'k', 'filled');
title('L')
subplot(1, 2, 2)
scatter(peCombined(choiceCombined==1), pRightChange(choiceCombined==1), 3, 'k', 'filled');
title('R')
%% load data
% load F:\tmpData\allSessionRPE.mat
%% plot P(right) - rpe
numBins = 6;
Lmeans = NaN(length(allPe), numBins);
Rmeans = NaN(length(allPe), numBins);
allMeans = NaN(length(allPe), numBins);
meanPe = NaN(length(allPe), numBins);

for i = 1:length(allChoices)
    edges = linspace(-1-0.0001, 1+0.0001, numBins+1);
    rpeR = allPe{i};
    choices = allChoices{i};
    pRPred = allPredR{i};
    rpeR(choices==0) = - rpeR(choices==0);
    rpeR(abs(rpeR)<0.07) = NaN;
    
    for j = 1:numBins
        currInd = find(rpeR>=edges(j) & rpeR<edges(j+1));
        pRcurr = mean(pRPred(currInd), 'omitnan');
        pRnext = mean(choices(intersect(currInd+1, 1:length(choices))), 'omitnan');
        allMeans(i,j) = pRnext-pRcurr;
        meanPe(i,j) = mean(rpeR(currInd), 'omitnan');
    end
end
% plot
figure2;
errorbar(mean(meanPe, 'omitnan'), mean(allMeans, 'omitnan'), sem(allMeans), 'LineWidth', 2, 'Color', [0.4 0.4 0.4]);
xlabel('rpe');
ylabel('deltaP(right)')
set(gca,'Box', 'off');
set(gca, 'TickDir', 'out')
%% plot pStay - rpe
numBins = 7;
allMeans = NaN(length(allPe), numBins);
meanPe = NaN(length(allPe), numBins);

for i = 1:length(allChoices)
    choices = allChoices{i};
    rpe = allPe{i};
    edges = quantile(rpe, linspace(0, 1, numBins+1));
    edges(1) = edges(1) - 0.001;
    edges(end) = edges(end) + 0.001;
    edges = linspace(min(rpe)-0.001, max(rpe)+0.001, numBins+1);
    stay = ones(length(choices)-1,1);
    stay(choices(1:end-1)~=choices(2:end)) = 0;
    stay = [NaN; stay];
    pStay = allPredR{i};
    pStay(find(choices(1:end-1)==0)+1) = 1 - pStay(find(choices(1:end-1)==0)+1);
    pStay(1) = NaN;
    pChoicePred = allPredR{i};
    pChoicePred(choices==0) = 1 - pChoicePred(choices==0);
    meanChange = mean(stay(2:end) - pChoicePred(1:end-1));
    for j = 1:numBins
        currInd = find(rpe>=edges(j) & rpe<edges(j+1));
        pRcurr = mean(pChoicePred(currInd), 'omitnan');
        pRnext = mean(stay(intersect(currInd+1, 1:length(choices))), 'omitnan');
        allMeans(i,j) = pRnext-pRcurr;
        meanPe(i,j) = mean(rpe(currInd), 'omitnan');
    end
%     allMeans(i,:) = allMeans(i,:) - meanChange;
end
% plot
figure2;
errorbar(mean(meanPe, 'omitnan'), mean(allMeans, 'omitnan'), sem(allMeans), 'LineWidth', 2, 'Color', [0.4 0.4 0.4]);
xlabel('rpe');
ylabel('deltaP(currChoice)')
set(gca,'Box', 'off');
set(gca, 'TickDir', 'out')
%% spearman correlation 
combineAllMeans = reshape(allMeans, [], 1);
combineMeanPe = reshape(meanPe, [], 1);
valInds = ~isnan(combineAllMeans)& ~isnan(combineMeanPe);
[rho,pval] = corr(combineAllMeans(valInds), combineMeanPe(valInds), 'type', 'Spearman');
%% plot predicted P(right)-rpe
% rightRpe = peCombined;
% rightRpe(choiceCombined==0) = -rightRpe(choiceCombined==0);
% 
% pRightChange = predRCombined(2:end) - predRCombined(1:end-1);
% pRightChange = [pRightChange; NaN];
% 
% numBins = 6;
% edges = linspace(-1-0.001, 1+0.001, numBins+1);
% pRightPredMean = NaN(1, numBins);
% pRightPredSem = NaN(1, numBins);
% pRightMean = NaN(1, numBins);
% pRightSem = NaN(1, numBins);
% peMean = NaN(1, numBins);
% 
% rightRpe(abs(rightRpe)<0.05) = NaN;
% for j = 1:numBins
%     rpeInd = find(rightRpe>=edges(j) & rightRpe<edges(j+1));
%     peMean(j) = mean(rightRpe(rpeInd), 'omitnan');
%     pRightPredMean(j) =  mean(pRightChange(rpeInd), 'omitnan');
%     pRightPredSem(j) = sem(pRightChange(rpeInd));
%     pRightMean(j) = mean(choiceCombined(intersect(rpeInd+1, 1:length(choiceCombined)))) - mean(choiceCombined(rpeInd));
%     pRightPost(j) = mean(choiceCombined(intersect(rpeInd+1, 1:length(choiceCombined))));
%     pRightPre(j) = mean(predRCombined(rpeInd));
%     pRightSem(j) = sqrt(sem_bern(choiceCombined(intersect(rpeInd+1, 1:length(choiceCombined))))^2 + sem(predRCombined(rpeInd))^2);
% end
%% use ones with neuron data
load('F:\tmpData\allUnitAUC.mat');
% get all RPEs
modelName = '5params';
sampNum = 2000;
allPe = {};
allChoices = {};
allPredR = {};

peCombined = [];
choiceCombined = [];
predRCombined = [];

for i = 1:length(allSessions)
    session = allSessions{i};
    pd = parseSessionString_df(session, root, sep);
    colInd = strcmp(allAnis, pd.aniName);
    category = allCols{colInd};
    params = getStanModelParams_sampsOnly(pd.animalName, category, modelName, sampNum, 'sessionName', session, 'biasFlag', 1, 'sessionParamsFlag', 1);
    t = inferModelVar(session, params, modelName, 'perturb', 0);
    s = behAnalysisNoPlot_opMD(session, 'simpleFlag', 1);
    pRight = t.probChoice;
    pRight(s.allChoices<0) = 1 - pRight(s.allChoices<0);
    allPe{i} = t.pe;
    allChoices{i} = 0.5*(s.allChoices + 1)';
    allPredR{i} = pRight;
 
end
%% deltaRight-RPE
numBins = 6;
allMeans = NaN(length(allSessions), numBins);
meanPe = NaN(length(allSessions), numBins);

for i = 1:length(allSessions)
    spikeRPE = spikeFocus{i};
    choices = allChoices{i};
    pRPred = allPredR{i};
    spikeRPE(~isnan(spikeRPE)) = zscore(spikeRPE(~isnan(spikeRPE)));
    spikeRPE(choices==0) = -spikeRPE(choices==0);
    edges = linspace(min(spikeRPE)-0.0001, max(spikeRPE)+0.0001, numBins+1);
    
    for j = 1:numBins
        currInd = find(spikeRPE>=edges(j) & spikeRPE<edges(j+1));
        pRcurr = mean(pRPred(currInd), 'omitnan');
        pRnext = mean(pRPred(intersect(currInd+1, 1:length(choices))), 'omitnan');
        allMeans(i,j) = pRnext-pRcurr;
        meanPe(i,j) = mean(spikeRPE(currInd), 'omitnan');
    end
    
end
%%
coeffsMax = [lmMaxAll.coeffs]';
tStatsMax = [lmMaxAll.tStats]';
sigMax = [lmMaxAll.ps]'<0.05;

figure2Wide;
subplot(1,2,1);
plotFilled(mean(meanPe(sigMax(:,1)==1&coeffsMax(:,1)>0&ind==2,:), 'omitnan'), allMeans(sigMax(:, 1)==1&coeffsMax(:,1)>0&ind==2,:), 'm');
title('deltaRight-spike, typeI')
subplot(1,2,2);
plotFilled(mean(meanPe(sigMax(:,1)==1&coeffsMax(:,1)<0&ind==1,:), 'omitnan'), allMeans(sigMax(:, 1)==1&coeffsMax(:,1)<0&ind==1,:), 'c');
title('deltaRight-spike, typeII')
%% deltaStay-spikes

numBins = 4;
allMeans = NaN(length(allSessions), numBins);
meanPe = NaN(length(allSessions), numBins);
allPeI = [];
allSpikesI = [];
allPeII = [];
allSpikesII = [];
for i = 1:length(allSessions)
    choices = allChoices{i};
    spikeRPE = spikeFocus{i};
    spikeRPE(~isnan(spikeRPE)) = zscore(spikeRPE(~isnan(spikeRPE)));
    spikeRPE(choices==0) = spikeRPE(choices==0) - mean(spikeRPE(choices==0), 'omitnan');
    spikeRPE(choices==1) = spikeRPE(choices==1) - mean(spikeRPE(choices==1), 'omitnan');
    edges = linspace(min(spikeRPE)-0.0001, max(spikeRPE)+0.0001, numBins+1);
    edges = quantile(spikeRPE, linspace(0, 1, numBins+1));
    edges(1) = edges(1) - 0.001;
    edges(end) = edges(end) + 0.001;
    stay = ones(length(choices)-1,1);
    stay(choices(1:end-1)~=choices(2:end)) = 0;
    stay = [NaN; stay];
    pStay = allPredR{i};
    pStay(find(choices(1:end-1)==0)+1) = 1 - pStay(find(choices(1:end-1)==0)+1);
    pStay(1) = NaN;
    pChoicePred = allPredR{i};
    pChoicePred(choices==0) = 1 - pChoicePred(choices==0);
    meanChange = mean(pStay(2:end) - pChoicePred(1:end-1));
    for j = 1:numBins
        currInd = find(spikeRPE>=edges(j) & spikeRPE<edges(j+1));
        pRcurr = mean(pChoicePred(currInd), 'omitnan');
        pRnext = mean(pStay(intersect(currInd+1, 1:length(choices))), 'omitnan');
        allMeans(i,j) = pRnext-pRcurr;
        meanPe(i,j) = mean(spikeRPE(currInd), 'omitnan');
    end
%     temp =  allMeans(i,:);
%     temp(~isnan(temp)) = zscore(temp(~isnan(temp)));
%     allMeans(i,:) = temp;
    allMeans(i,:) = allMeans(i,:) - meanChange;
end
%%
coeffsMax = [lmMaxAll.coeffs]';
tStatsMax = [lmMaxAll.tStats]';
sigMax = [lmMaxAll.ps]'<0.05;

color1 = [0 0 0];
color2 = [0.5 0.5 0.5];

figure2Wide;
subplot(1,2,1);
errorbar(mean(meanPe(sigMax(:,1)==1&coeffsMax(:,1)>0&ind==2,:), 'omitnan'), mean(allMeans(sigMax(:, 1)==1&coeffsMax(:,1)>0&ind==2,:), 'omitnan'), ...
    sem(allMeans(sigMax(:, 1)==1&coeffsMax(:,1)>0&ind==2,:)), 'LineWidth', 2, 'Color', color2);
title('deltaP(curr)-spike, typeI')
set(gca,'Box', 'off');
set(gca, 'TickDir', 'out')
xlabel('spike/s (zscored)')
ylabel('deltaP(currC)')
ylim([-0.05 0.05])
subplot(1,2,2);
errorbar(mean(meanPe(sigMax(:,1)==1&coeffsMax(:,1)<0&ind==1,:), 'omitnan'), mean(allMeans(sigMax(:, 1)==1&coeffsMax(:,1)<0&ind==1,:), 'omitnan'), ...
    sem(allMeans(sigMax(:, 1)==1&coeffsMax(:,1)<0&ind==1,:)), 'LineWidth', 2, 'Color', color1);
title('deltaP(currC)-spike, typeII')
set(gca,'Box', 'off');
set(gca, 'TickDir', 'out')
xlabel('spike/s (zscored)')
ylim([-0.05 0.05])
%% calculate coeff
% type I 
spike = reshape(meanPe(sigMax(:,1)==1&coeffsMax(:,1)>0&ind==2,:), [], 1);
change = reshape(allMeans(sigMax(:,1)==1&coeffsMax(:,1)>0&ind==2,:), [], 1);
lm = fitlm(spike, change);
% type II
spike = reshape(meanPe(sigMax(:,1)==1&coeffsMax(:,1)<0&ind==1,:), [], 1);
change = reshape(allMeans(sigMax(:,1)==1&coeffsMax(:,1)<0&ind==1,:), [], 1);
lm = fitlm(spike, change);
%% no selection
figure2Wide;
subplot(1,2,1);
errorbar(mean(meanPe(ind==2,:), 'omitnan'), mean(allMeans(ind==2,:), 'omitnan'), ...
    sem(allMeans(ind==2,:)), 'LineWidth', 2, 'Color', color2);
title('deltaP(curr)-spike, typeI')
set(gca,'Box', 'off');
set(gca, 'TickDir', 'out')
xlabel('spike/s (zscored)')
ylabel('deltaP(currC)')
ylim([-0.04 0.03])
subplot(1,2,2);
errorbar(mean(meanPe(ind==1,:), 'omitnan'), mean(allMeans(ind==1,:), 'omitnan'), ...
    sem(allMeans(ind==1,:)), 'LineWidth', 2, 'Color', color1);
title('deltaP(currC)-spike, typeII')
set(gca,'Box', 'off');
set(gca, 'TickDir', 'out')
xlabel('spike/s (zscored)')
ylim([-0.04 0.03])
%%









