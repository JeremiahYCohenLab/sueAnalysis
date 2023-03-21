load('F:\tmpData\allSessionUnit.mat');
load('F:\tmpData\waveformFiltered.mat');
tb = 1;
tf = 2;
stepSize = 200;
binSize = 1000;

transIdx = zeros(length(waveformsSession), 1);
for i = 1:length(transIdx)
    x = cellfun(@(x)strcmp(allSessionWF{i}, x), allSessions);
    y = cellfun(@(x)strcmp(allUnitWF{i}, x), allUnits);
    if ~isempty(find(x+y>=2,1))
        transIdx(i) = 1;
    end
end
waveformsSession = waveformsSession(logical(transIdx),:);
color1 = [0.2 0.2 0.2];
color2 = [0.8 0.8 0.8];
%% clustering
numCat = 2;
[coeff,score,latent, ~, explained, mu] = pca(zscore(waveformsSession, [], 1));
indAll = {};
dis = {};
for a = 1:10
    [indAll{a}, ~, dis{a}] = kmeans([score(:, 1:5)], numCat);
end
[~,optiInds] = min(cellfun(@mean, dis));
ind = indAll{optiInds};
figure2;
hold on;
plotFilled(1:size(waveformsSession,2), waveformsSession(ind==1,:), color1);
plotFilled(1:size(waveformsSession,2), waveformsSession(ind==2,:), color2);

plot([5 21], [-0.3 -0.3], 'LineWidth', 3,'Color', 'k');
text(10, -0.4, '0.5 ms', 'FontSize', 14)

set(gca, 'Box', 'off')
set(gca,'tickdir', 'out')
xlabel('expected value tstats', 'FontSize', 18)
set(gca, 'XTick', [-10:5:5], 'FontSize', 14)
set(gca, 'YTick', [0 0.1 0.2], 'FontSize', 14)
set(gca, 'XColor', 'none')
set(gca, 'YColor', 'none')
%% beh analysis
uniqSessions = unique(allSessions);
FAs = [];
MISSs = [];
csPlus = [];
csMinus = [];
for i = 1:length(uniqSessions)
    session = uniqSessions{i};
    s = behAnalysisNoPlot_opMD(session, 'simpleFlag', 1);
    cue = 2*(double(s.CSplus)'-0.5);
    resp = zeros(size(cue));
    resp(~isnan(s.lickSide)) = 1;
    resp(isnan(s.lickSide)) = -1;
    corre = resp.*cue;
    FAs(i) = sum(cue==-1 & resp==1); 
    MISSs(i) = sum(cue==1 & resp==-1); 
    len(i) = length(s.behSessionData);
    csPlus(i) = sum(s.CSplus);
    csMinus(i) = sum(s.CSminus);
end
%% histogram
% miss
figure2;
histogram(MISSs./csPlus,  10, 'FaceColor', [0.6 0.6 0.6], 'EdgeColor', 'none', 'Normalization', 'probability');
title('miss rate')
set(gca, 'TickDir', 'out');
set(gca, 'Box', 'off');
% FA
figure2;
histogram(FAs./csMinus,  10, 'FaceColor', [0.6 0.6 0.6], 'EdgeColor', 'none', 'Normalization', 'probability');
title('FA rate')
set(gca, 'TickDir', 'out');
set(gca, 'Box', 'off');
% hit vs FA
figure2; hold on;
scatter(FAs./csMinus, 1-(MISSs./csPlus), 25, 'k', 'MarkerEdgeAlpha', 0.6, 'LineWidth', 1.5);
plot([0 1], [0 1], 'LineStyle', '--', 'LineWidth', 2, 'Color', 'r');
xlabel('P(resp|csMinus)')
ylabel('P(resp|csPlus)')
%% neuron correlation
tb = 2; % in s
tf = 0.3; % in s
fomula = 'resp ~ 1 + go + baseline';
blStats = [];
corrPre = [];
corrPost = [];
pPre = [];
pPost = [];
allTstats = [];
allSigs = [];
allCoeffs = [];
allRsq = [];
allLL = [];
allBIC = [];
indCorrect = [];
FAs = [];
MISSs = [];
csPlus = [];
csMinus = [];
for i = 1:length(allSessions)
    session = allSessions{i};
    unit = allUnits{i};
    [~, matCue] = getUnitMatCue(session, unit, tb, 0, 100, 100);
    baseline = zscore(sum(matCue, 2)/tb);
    [~, matCue] = getUnitMatCue(session, unit, 0, tf, 100, 100);
    go = zscore(sum(matCue, 2)/tf);
    s = behAnalysisNoPlot_opMD(session, 'simpleFlag', 1);
    cue = 2*(double(s.CSplus)'-0.5);
    resp = zeros(size(cue));
    resp(~isnan(s.lickSide)) = 1;
    resp(isnan(s.lickSide)) = 0;
    corre = resp.*cue;
    FAs(i) = sum(cue~=1 & resp==1); 
    MISSs(i) = sum(cue==1 & resp~=1); 
    len(i) = length(s.behSessionData);
    csPlus(i) = sum(s.CSplus);
    csMinus(i) = sum(s.CSminus);
    if (sum(cue==1 & resp~=1))>=0.01*length(s.CSplus)
        baseline = baseline(s.CSplus);
        go = go(s.CSplus);
        lmTmp = fitlm(baseline, go);
        blStats = [blStats, lmTmp.Coefficients.tStat(2)];
        [R, P] = corrcoef(baseline, go);
        corrPre = [corrPre; R(1,2)];   
        pPre = [pPre; P(1,2)];
        goR = lmTmp.Residuals.Raw;
        [R, P] = corrcoef(baseline, goR);
        corrPost = [corrPost; R(1,2)];
        pPost = [pPost; P(1,2)];
        resp = resp(s.CSplus);
        indCorrect = [indCorrect; ind(i)];
        tbl = table(baseline, go, resp, goR);
        lm = fitglm(tbl, fomula, 'Distribution','binomial');
        allTstats = [allTstats; lm.Coefficients.tStat(2:end)'];
        allSigs = [allSigs; lm.Coefficients.pValue(2:end)'];
        allCoeffs = [allCoeffs; lm.Coefficients.Estimate(2:end)'];
        allRsq = [allRsq; lm.Rsquared.Adjusted'];  
        allLL = [allLL; lm.LogLikelihood/length(s.CSplus)];
        allBIC = [allBIC; lm.ModelCriterion.BIC];
    end
end
%%
regressors = lm.CoefficientNames(2:end);
figure2;
for i = 1:length(regressors)
    subplot(length(regressors), 1, i); hold on;
    edges = linspace(min(allTstats(:, i))-0.01, max(allTstats(:, i))+0.01, 10);
    histogram(allTstats(indCorrect==1, i), edges, 'FaceColor', color1, 'EdgeColor', 'none', 'Normalization', 'probability');
    histogram(allTstats(indCorrect==2, i), edges, 'FaceColor', color2, 'EdgeColor', 'none', 'Normalization', 'probability');
    [h, p, ~, stats] = ttest2(allTstats(indCorrect==1, i),allTstats(indCorrect==2, i));
    title([regressors{i} 'p=' num2str(p) 'tstat=' num2str(stats.tstat), 'df=' num2str(stats.df)])
    legend({'Type I', 'Type II'})
end
% figure2;
% histogram(allRsq, 'FaceColor', [0.6 0.6 0.6], 'EdgeColor', 'none', 'Normalization', 'probability');
% title('Rsq')
%%
regressors = lm.CoefficientNames(2:end);
figure2;
for i = 1:length(regressors)
    subplot(length(regressors), 1, i); hold on;
    edges = linspace(min(allTstats(:, i))-0.01, max(allTstats(:, i))+0.01, 10);
    histogram(allTstats(indCorrect==1&allSigs(:, i)<0.05, i), edges, 'FaceColor', color1, 'EdgeColor', 'none', 'Normalization', 'probability');
    histogram(allTstats(indCorrect==2&allSigs(:, i)<0.05, i), edges, 'FaceColor', color2, 'EdgeColor', 'none', 'Normalization', 'probability');
    title(regressors{i})
    legend({'Type I', 'Type II'})
end
%%
figure2; hold on;
scatter(allTstats(indCorrect==1, 1), allTstats(indCorrect==1, 2), 25, color1, 'LineWidth', 2, 'MarkerEdgeAlpha', 0.8);
scatter(allTstats(indCorrect==2, 1), allTstats(indCorrect==2, 2), 25, color2, 'LineWidth', 2, 'MarkerEdgeAlpha', 0.8);
xlabel('baseline')
ylabel('goCueResidual')
legend({'TypeII', 'TypeI'})
%%
figure2; hold on;
edges = linspace(min(corrPre)-0.001, max(corrPre)+0.001, 10);
histogram(corrPre, edges, 'FaceColor', color1, 'Normalization', 'probability');
xlabel('corrCoeff(blGo)')
[h, p, ~, stats] = ttest(corrPre);
title(['bl-go' 'p=' num2str(p) 'tstat=' num2str(stats.tstat), 'df=' num2str(stats.df), 'sig=' num2str(sum(pPre<0.05)/length(pPre))])

%%
figure2; hold on;
histogram(corrPost, 'FaceColor', color1, 'Normalization', 'probability');
xlabel('corrCoeff(blGoR)')
[h, p, ~, stats] = ttest(corrPost);
title(['bl-go' 'p=' num2str(p) 'tstat=' num2str(stats.tstat), 'df=' num2str(stats.df), 'sig=' num2str(sum(pPost<0.05)/length(pPost))])
%% neuron baseline
tb = 2; % in s
tf = 0.3; % in s
stepSize = 100;
binSize = 500;
blHit = [];
blMiss = [];
goHit = [];
goMiss = [];
indCorrect = [];

for i = 1:length(allSessions)
    session = allSessions{i};
    unit = allUnits{i};
    s = behAnalysisNoPlot_opMD(session, 'simpleFlag', 1);
    cue = 2*(double(s.CSplus)'-0.5);
    resp = zeros(size(cue));
    resp(~isnan(s.lickSide)) = 1;
    resp(isnan(s.lickSide)) = 0;
    corre = resp.*cue;
    
    [~, matCue, matCueSlide, slidetime] = getUnitMatCue(session, unit, tb, -0.1, stepSize, binSize);
    baseline = matCueSlide;
    [~, matCue, matCueSlide, slidetime] = getUnitMatCue(session, unit, 0, tf, stepSize, binSize);
    goCue = matCue;
    if (sum(cue==1 & resp~=1))>=0.01*length(s.CSplus)
        tempHit = mean(baseline(cue==1 & resp==1, :));
        tempMiss = mean(baseline(cue==1 & resp==0, :));
        blHit = [blHit; mean(tempHit)];
        blMiss = [blMiss; mean(tempMiss)];
        tempHit = mean(goCue(cue==1 & resp==1, :));
        tempMiss = mean(goCue(cue==1 & resp==0, :));    
        goHit = [goHit; sum(tempHit)/tf];
        goMiss = [goMiss; sum(tempMiss)/tf];
        indCorrect = [indCorrect; ind(i)];
    end
end
%% all neurons baseline
% figure2;
% plotFilled(slidetime, blHit,'k')
% hold on
% plotFilled(slidetime, blMiss,'r')
% xlabel('time from go cue')
% ylabel('spikes/s zscored')
% legend('Hit', '', 'Miss')

figure2;
subplot(1,2,1); hold on; scatter(blHit(indCorrect==1), blMiss(indCorrect==1), 15, color1, 'filled');
hold on; plot([0 10], [0 10]);
ylabel('Type II')
[h,p,ci,stats] = ttest(blHit(indCorrect==1), blMiss(indCorrect==1));
title(['p=' num2str(p), 'tState=' num2str(stats.tstat), 'df=' num2str(stats.df)]);
xlabel('baseline Hit')
ylabel('Type II: baseline Miss')

subplot(1,2,2); hold on; scatter(blHit(indCorrect==2), blMiss(indCorrect==2), 15, color2, 'filled');
hold on; plot([0 10], [0 10]);
ylabel('Type I')
[h,p,ci,stats] = ttest(blHit(indCorrect==2), blMiss(indCorrect==2));
title(['p=' num2str(p), 'tState=' num2str(stats.tstat), 'df=' num2str(stats.df)]);
xlabel('baseline Hit')
ylabel('Type I: baseline Miss')
%% all neurons go cue
figure2;
subplot(1,2,1); hold on; scatter(goHit(indCorrect==1), goMiss(indCorrect==1), 15, color1, 'filled');
hold on; plot([0 10], [0 10]);
ylabel('Type II')
[h,p,ci,stats] = ttest(goHit(indCorrect==1), goMiss(indCorrect==1));
title(['p=' num2str(p), 'tState=' num2str(stats.tstat), 'df=' num2str(stats.df)]);
xlabel('go Hit')
ylabel('Type II: go Miss')

subplot(1,2,2); hold on; scatter(goHit(indCorrect==2), goMiss(indCorrect==2), 15, color2, 'filled');
hold on; plot([0 10], [0 10]);
ylabel('Type I')
[h,p,ci,stats] = ttest(goHit(indCorrect==2), goMiss(indCorrect==2));
title(['p=' num2str(p), 'tState=' num2str(stats.tstat), 'df=' num2str(stats.df)]);
xlabel('go Hit')
ylabel('Type I: go Miss')
%% sep by type
figure2Wide;
subplot(1,2,1); hold on;
plotFilled(slidetime, blHit(indCorrect==1,:),'k')
plotFilled(slidetime, blMiss(indCorrect==1,:),'r')
xlabel('time from go cue')
ylabel('spikes/s zscored')                                                                             
legend('Hit', '', 'Miss')
title('type II')
subplot(1,2,2); hold on;
plotFilled(slidetime, blHit(indCorrect==2,:),'k')
plotFilled(slidetime, blMiss(indCorrect==2,:),'r')
xlabel('time from go cue')
ylabel('spikes/s zscored')
legend('Hit', '', 'Miss')
title('Type I')
%%