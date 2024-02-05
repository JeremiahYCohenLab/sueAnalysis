dayList = getDayList('allDBh-cre', 'all-DBh', 'noRetractionDelay');
%%
allDDs = [];
allQc = [];
allQalt = [];
allQsum = [];
allQdiff = [];
allBias = [];
allSwitch = [];
allSwitchNext = [];
modelName = '5params';
category = 'good';
sampNum = 2000;
[root, sep] = currComputer();

for i = 1:length(dayList)
    session = dayList{i};
    pd = parseSessionString_df(session, root, sep);
    s = behAnalysisNoPlot_opMD(session,'simpleFlag', 1);
    params = getStanModelParams_sampsOnly(pd.animalName, category, modelName, sampNum, 'sessionName', session, 'biasFlag', 1, 'sessionParamsFlag', 1);
    t = inferModelVar(session, params, modelName);
    Qsum = sum(t.Q, 2);
    Qchosen = zeros(length(s.allChoices),1);
    Qchosen(s.allChoices==1) = t.Q(s.allChoices==1, 2);
    Qchosen(s.allChoices==-1) = t.Q(s.allChoices==-1, 1);

    Qalt = zeros(length(s.allChoices),1);
    Qalt(s.allChoices==1) = t.Q(s.allChoices==1, 1);
    Qalt(s.allChoices==-1) = t.Q(s.allChoices==-1, 2);

    Qdiff = Qchosen - Qalt;
    
    biasSide = zeros(size(s.allChoices))';
    if mean(params(:, end)) > 0
        biasSide(s.allChoices==1) = 1;  
    else
        biasSide(s.allChoices==-1) = 1;
    end
    svs = zeros(size(s.allChoices))';
    svs(s.changeChoice_Inds) = 1;
    svsNext = [svs(2:end); NaN];
    allDDs = [allDDs; abs(s.dd)];
    allQc = [allQc; Qchosen];
    allQalt = [allQalt; Qalt];
    allQsum = [allQsum; Qsum];
    allQdiff = [allQdiff; Qdiff];
    allBias = [allBias; biasSide];
    allSwitch = [allSwitch; svs];
    allSwitchNext = [allSwitchNext; svsNext];
end
%%
fit = fitglm([allQdiff, allQsum, allBias, allSwitch], allDDs);
%%
coeffs = fit.Coefficients.Estimate(2:end);
CI = coefCI(fit);
CI = CI(2:end, :);
upper = CI(:, 2) - coeffs;
lower = CI(:, 1) - coeffs;
figure; hold on;
bar(1:length(coeffs), coeffs, 'FaceColor', [0.7 0.7 0.7])
errorbar(1:length(coeffs), coeffs, lower, upper, 'LineStyle', 'none', 'LineWidth',2, 'Color', [0.4 0.4 0.4])
%%  
fitNext = fitglm([allQdiff, allQsum, allSwitch, allDDs], allSwitchNext);
%%
coeffs = fitNext.Coefficients.Estimate(2:end);
CI = coefCI(fitNext);
CI = CI(2:end, :);
upper = CI(:, 2) - coeffs;
lower = CI(:, 1) - coeffs;
figure; hold on;
bar(1:length(coeffs), coeffs, 'FaceColor', [0.7 0.7 0.7])
errorbar(1:length(coeffs), coeffs, lower, upper, 'LineStyle', 'none', 'LineWidth',2, 'Color', [0.4 0.4 0.4])
%%