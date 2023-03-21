%% get all animal lists and dayLists
[root, sep] = currComputer();
xlFile = 'activationAll';
sheet = 'all';
col = 'activation';
dayList = getDayList(xlFile, sheet, col);
%% get all RPEs;
modelName = '5params_k_bias_LaserNegRPE';
modelName = '5params';
sampNum = 1000;
allPe = {};
allChoices = {};
allPredR = {};

peCombined = [];
choiceCombined = [];
predRCombined = [];

for i = 1:length(dayList)
    session = dayList{i};
    pd = parseSessionString_df(session, root, sep);
    params = getStanModelParams_sampsOnly(pd.animalName, col, modelName, sampNum, 'sessionName', session, 'biasFlag', 1, 'sessionParamsFlag', 1);
    t = inferModelVar(session, params, modelName, 'perturb', 1);
    s = behAnalysisNoPlot_opMD(session, 'simpleFlag', 1);
    pRight = t.probChoice;
    pRight(s.allChoices<0) = 1 - pRight(s.allChoices<0);
    allPe{i} = t.pe;
    allChoices{i} = 0.5*(s.allChoices + 1)';
    allPredR{i} = pRight;
    
    
    peCombined = [peCombined; t.pe];
    choiceCombined = [choiceCombined; 0.5*(s.allChoices + 1)'];
    predRCombined = [predRCombined; pRight];
 
end
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
errorbar(mean(meanPe, 'omitnan'), mean(allMeans, 'omitnan'), sem(allMeans), 'LineWidth', 2, 'Color', [0.3 0.3 0.3]);
xlabel('rpe');
ylabel('deltaP(right)')
set(gca,'Box', 'off');
set(gca, 'TickDir', 'out')
%% plot pStay - rpe
numBins = 5;
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
    allMeans(i,:) = allMeans(i,:) - meanChange;
end
% plot
figure2;
errorbar(mean(meanPe, 'omitnan'), mean(allMeans, 'omitnan'), sem(allMeans), 'LineWidth', 2, 'Color', [0.4 0.4 0.4]);
xlabel('rpe');
ylabel('deltaP(currChoice)')
set(gca,'Box', 'off');
set(gca, 'TickDir', 'out')
%% separate laser vs no laser, plot pStay - rpe
numBins = 3;
allMeansN = NaN(length(allPe), numBins);
meanPeN = NaN(length(allPe), numBins);
allMeansL = NaN(length(allPe), numBins);
meanPeL = NaN(length(allPe), numBins);

for i = 1:length(allChoices)
    choicesAll = allChoices{i};
    rpeAll = allPe{i};
    s = behAnalysisNoPlot_opMD(dayList{i}, 'simpleFlag', 1);
    stay = ones(length(choicesAll)-1,1);
    stay(choicesAll(1:end-1)~=choicesAll(2:end)) = 0;
    stay = [NaN; stay];
    stayNext = [stay(2:end); NaN];
    pStay = allPredR{i};
    pStay(find(choicesAll(1:end-1)==0)+1) = 1 - pStay(find(choicesAll(1:end-1)==0)+1);
    pStay(1) = NaN;
    pStayNext = [pStay(2:end); NaN];
    pChoicePred = allPredR{i};
    pChoicePred(choicesAll==0) = 1 - pChoicePred(choicesAll==0);
    meanChange = mean(stay(2:end) - pChoicePred(1:end-1));
    edges = quantile(rpeAll, linspace(0, 1, numBins+1));
    edges(1) = edges(1) - 0.001;
    edges(end) = edges(end) + 0.001;
    edges = linspace(min(rpeAll)-0.001, max(rpeAll)+0.001, numBins+1);
    % no laser
    rpe = rpeAll(s.laser==0);
    pChoicePredN = pChoicePred(s.laser==0);
    stayNextN = stayNext(s.laser==0);
    pStayNextN = pStayNext(s.laser==0);
%     edges = quantile(rpe, linspace(0, 1, numBins+1));
%     edges(1) = edges(1) - 0.001;
%     edges(end) = edges(end) + 0.001;
%     
%     edges = linspace(min(rpe)-0.001, max(rpe)+0.001, numBins+1);
    for j = 1:numBins
        currInd = find(rpe>=edges(j) & rpe<edges(j+1));
        pRcurr = mean(pChoicePredN(currInd), 'omitnan');
        pRnext = mean(stayNextN(currInd), 'omitnan');
        allMeansN(i,j) = pRnext-pRcurr;
        meanPeN(i,j) = mean(rpe(currInd), 'omitnan');
    end
    allMeansN(i,:) = allMeansN(i,:) - meanChange;
    % laser
    rpe = rpeAll(s.laser==1);
    pChoicePredL = pChoicePred(s.laser==1);
    stayNextL = stayNext(s.laser==1);
    pStayNextL = pStayNext(s.laser==1);
%     edges = quantile(rpe, linspace(0, 1, numBins+1));
%     edges(1) = edges(1) - 0.001;
%     edges(end) = edges(end) + 0.001;
%     
    edges = linspace(min(rpe)-0.001, max(rpe)+0.001, numBins+1);
    for j = 1:numBins
        currInd = find(rpe>=edges(j) & rpe<edges(j+1));
        pRcurr = mean(pChoicePredL(currInd), 'omitnan');
        pRnext = mean(stayNextL(currInd), 'omitnan');
        allMeansL(i,j) = pRnext-pRcurr;
        meanPeL(i,j) = mean(rpe(currInd), 'omitnan');
    end
    allMeansL(i,:) = allMeansL(i,:) - meanChange;    
end

%% plot
figure2; hold on;
errorbar(mean(meanPeN, 'omitnan'), mean(allMeansN, 'omitnan'), sem(allMeansN), 'LineWidth', 2, 'Color', [0.4 0.4 0.4]);
errorbar(mean(meanPeL, 'omitnan'), mean(allMeansL, 'omitnan'), sem(allMeansL), 'LineWidth', 2, 'Color', [0.2 0.2 1]);
legend('no laser', 'laser')
xlabel('rpe');
ylabel('deltaP(currChoice)')
set(gca,'Box', 'off');
set(gca, 'TickDir', 'out')
%% paired test and correction
p = NaN(1, numBins);
tstats = NaN(1, numBins);
for i = 1:numBins
    valInds = ~isnan(allMeansN(:,i))& ~isnan(allMeansL(:,i));
%     [h,p(i),ci,stats] = ttest(allMeansN(valInds,i), allMeansL(valInds,i)); 
%     tstats(i) = stats.tstat;
    [p(i), h] = signrank(allMeansN(valInds,i), allMeansL(valInds,i));
end
%%


%%



