% load data and params
load('F:\tmpData\allUnitAUC.mat');
modelName = '5params';
[root, sep] = currComputer();
stepSize = 100;
binSize = 100;
numSamps = 2000;
%%
clear allLM
for i = 1:length(allSessions)
    session = allSessions{i};
    unit = allUnits{i};
    s = behAnalysisNoPlot_opMD(session, 'simpleFlag', 1);
    rightSide = zeros(1, length(s.allChoices))';
    rightSide(s.allChoices>0) = 1;
    tf = 0.1*(s.rwdDelay + max(50, min(s.lickLat)));
    if tf>=0.5
        tb = 0.5;
        tf = 0;
    else
        tb = 0.5 - tf;
    end
    pd = parseSessionString_df(session, root, sep);
    params = getStanModelParams_sampsOnly(pd.animalName, 'good', modelName, numSamps);
    t = inferModelVar(session, params, modelName);
    
    [spikeCell] = getUnitMatCue(session, unit, tf, tb, stepSize, binSize);
    spikeNumPre = cellfun(@length, spikeCell);
    spikeNumPre = spikeNumPre(s.responseInds)';
    lickDiff = s.lickDiff';
    lickLat = s.lickLat';
    Qchosen = t.Q(:,2);
    Qchosen(s.allChoices<0) = t.Q(s.allChoices<0,1);
    svs = NaN(length(s.allChoices),1);
    svs(s.changeChoice_Inds) = 1;
    svs(s.stayChoice_Inds) = 0;
    tbl = table(spikeNumPre, rightSide, lickDiff, lickLat, Qchosen, svs);
    lm = fitlm(tbl, 'spikeNumPre ~ 1 + Qchosen + rightSide + svs');
    allLM(i).coeffs = lm.Coefficients.Estimate(2:end);
    allLM(i).tStats = lm.Coefficients.tStat(2:end);
    allLM(i).pSig = lm.Coefficients.pValue(2:end);
end
%%
preTstats = [allLM.tStats];
preSig = [allLM.pSig];
preSig = preSig<=0.05;
regressor = lm.CoefficientNames(2:end);
colors = cool(length(regressor));


figure2;

for j = 1:length(regressor)
    edges = linspace(min(preTstats(j,:)), max(preTstats(j,:)), 20);
    subplot(length(regressor), 1, j); hold on;
    histogram(preTstats(j,:), edges, 'FaceColor', [0.7 0.7 0.7], 'EdgeColor', 'none');
    histogram(preTstats(j,preSig(j,:)), edges, 'FaceColor', colors(j,:), 'EdgeColor', 'none');
    title(regressor{j})
end
%%


