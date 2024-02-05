col = 'cueOnGood';
modelName = '5params';
sampNum = 2000;
[root,sep] = currComputer();
dayList = getDayList('inhibitionGt', 'allGt', col);
%% velocity
numBins = 4;

varT = {'QdiffChosen', 'Qsum'};
varK = {'maxS', 'speedFixInter', 'vXFixInterC', 'vYFixInter', 'speedPrePeak', 'xMaxC'};

matAllK = NaN(length(varT), length(varK), numBins, length(dayList));
matAllT = NaN(length(varT), length(varK), numBins, length(dayList));

for sess = 1:length(dayList)
    session = dayList{sess};
    kine = lickPlotSessionNoPlot(session);
    
    pd = parseSessionString_df(session, root, sep);
    params = getStanModelParams_sampsOnly(pd.aniName, col, modelName, sampNum, 'sessionName', session, 'biasFlag',1);
    t = inferModelVar(session, params, modelName, 'perturb', 0);
    s = behAnalysisNoPlot_opMD(session, 'simpleFlag', 1);
    Qdiff = t.Q(:,2) - t.Q(:,1);
    Qsum = t.Q(:,2) + t.Q(:,1);
    Qsum(s.allChoices == 1) = zscore(Qsum(s.allChoices == 1));
    Qsum(s.allChoices == -1) = zscore(Qsum(s.allChoices == -1));
    QdiffChosen = Qdiff;
    QdiffChosen(s.allChoices==-1) = - QdiffChosen(s.allChoices==-1);
    pRight = t.probChoice;
    pRight(s.allChoices==1) = 1 - t.probChoice(s.allChoices==1);
    conf = t.probChoice;

    for i = 1:length(varT)
        eval(['currT = ' varT{i} ';']);
        for j = 1:length(varK)
            eval(['currK = kine.' varK{j} ';']);
            edges = quantile(currT, linspace(0, 1, numBins+1));
            edges(1) = edges(1) - 0.0001;
            edges(end) = edges(end) + 0.0001;
            for b = 1:numBins
                currMT = mean(currT(currT > edges(b) & currT <= edges(b+1)), "all", 'omitnan');
                currMK = mean(currK(currT > edges(b) & currT <= edges(b+1)), "all", 'omitnan');
                matAllT(i, j, b, sess) = currMT;
                matAllK(i, j, b, sess) = currMK;
            end
        end
    end

end
%%
figure2;
for i = 1:length(varT)
    for j = 1:length(varK)
        currMatT = squeeze(matAllT(i, j, :, :));
        currMatK = squeeze(matAllK(i, j, :, :));

        currMatTM = mean(currMatT, 2, 'omitmissing')';
        currMatKM = mean(currMatK, 2, 'omitmissing')';
        currMatKSem = sem(currMatK');

        subplot(length(varK), length(varT), length(varT)*(j-1)+i);
        hold on;
        plot(currMatTM, currMatKM, 'Color', 'k', 'LineWidth', 2);
        patch([currMatTM, flip(currMatTM)], [currMatKM-currMatKSem, flip(currMatKM+currMatKSem)], [0.5 0.5 0.5], 'FaceAlpha', 0.4,  'EdgeColor','none')
        xlabel(varT{i});
        ylabel(varK{j});
    end
end
%% variance
numBins = 5;

varT = {'QdiffChosen', 'Qsum'};
varK = {'xMaxC'};

matAllK = NaN(length(varT), length(varK), numBins, length(dayList));
matAllT = NaN(length(varT), length(varK), numBins, length(dayList));

for sess = 1:length(dayList)
    session = dayList{sess};
    kine = lickPlotSessionNoPlot(session);
    
    pd = parseSessionString_df(session, root, sep);
    params = getStanModelParams_sampsOnly(pd.aniName, col, modelName, sampNum, 'sessionName', session, 'biasFlag',1);
    t = inferModelVar(session, params, modelName, 'perturb', 0);
    s = behAnalysisNoPlot_opMD(session, 'simpleFlag', 1);
    Qdiff = t.Q(:,2) - t.Q(:,1);
    Qsum = t.Q(:,2) + t.Q(:,1);
    Qsum(s.allChoices == 1) = zscore(Qsum(s.allChoices == 1));
    Qsum(s.allChoices == -1) = zscore(Qsum(s.allChoices == -1));
    QdiffChosen = Qdiff;
    QdiffChosen(s.allChoices==-1) = - QdiffChosen(s.allChoices==-1);
    pRight = t.probChoice;
    pRight(s.allChoices==1) = 1 - t.probChoice(s.allChoices==1);
    conf = t.probChoice;

    for i = 1:length(varT)
        eval(['currT = ' varT{i} ';']);
        for j = 1:length(varK)
            currX = kine.xMaxC;
            currY = kine.yMax;

            % C = cov(currX(~isnan(currX)), currY(~isnan(currX)));

            eval(['currK = kine.' varK{j} ';']);
            edges = quantile(currT, linspace(0, 1, numBins+1));
            edges(1) = edges(1) - 0.0001;
            edges(end) = edges(end) + 0.0001;
            % edges = linspace(min(currT)-0.001, max(currT)+0.001, numBins+1);
            for b = 1:numBins
                currMT = mean(currT(currT > edges(b) & currT <= edges(b+1)), "all", 'omitmissing');
                currMK = std(currK(currT > edges(b) & currT <= edges(b+1)), 'omitmissing');
                matAllT(i, j, b, sess) = currMT;
                matAllK(i, j, b, sess) = currMK;        
                % mah
                tempX = currX(currT > edges(b) & currT <= edges(b+1));
                tempY = currY(currT > edges(b) & currT <= edges(b+1));
                disM = mahDis([tempX', tempY']);

            end
        end
    end

end
%%
figure2;
for i = 1:length(varT)
    for j = 1:length(varK)
        currMatT = squeeze(matAllT(i, j, :, :));
        currMatK = squeeze(matAllK(i, j, :, :));

        currMatTM = mean(currMatT, 2, 'omitmissing')';
        currMatKM = mean(currMatK, 2, 'omitmissing')';
        currMatKSem = sem(currMatK');

        subplot(length(varK), length(varT), length(varT)*(j-1)+i);
        hold on;
        plot(currMatTM, currMatKM, 'Color', 'k', 'LineWidth', 2);
        patch([currMatTM, flip(currMatTM)], [currMatKM-currMatKSem, flip(currMatKM+currMatKSem)], [0.5 0.5 0.5], 'FaceAlpha', 0.4,  'EdgeColor','none')
        xlabel(varT{i});
        ylabel(varK{j});
    end
end
%% malDis
numBins = 4;

varT = {'QdiffChosen', 'Qsum'};
matAllK = NaN(length(varT), 1, numBins, length(dayList));
matAllT = NaN(length(varT), 1, numBins, length(dayList));

for sess = 1:length(dayList)
    session = dayList{sess};
    kine = lickPlotSessionNoPlot(session);
    
    pd = parseSessionString_df(session, root, sep);
    params = getStanModelParams_sampsOnly(pd.aniName, col, modelName, sampNum, 'sessionName', session, 'biasFlag',1);
    t = inferModelVar(session, params, modelName, 'perturb', 0);
    s = behAnalysisNoPlot_opMD(session, 'simpleFlag', 1);
    Qdiff = t.Q(:,2) - t.Q(:,1);
    Qsum = t.Q(:,2) + t.Q(:,1);
    Qsum(s.allChoices == 1) = zscore(Qsum(s.allChoices == 1));
    Qsum(s.allChoices == -1) = zscore(Qsum(s.allChoices == -1));
    mean(Qsum)
    figure;
    histogram(Qsum)
    QdiffChosen = Qdiff;
    QdiffChosen(s.allChoices==-1) = - QdiffChosen(s.allChoices==-1);
    pRight = t.probChoice;
    pRight(s.allChoices==1) = 1 - t.probChoice(s.allChoices==1);
    conf = t.probChoice;

    for i = 1:length(varT)
        eval(['currT = ' varT{i} ';']);
        currX = kine.xMaxC;
        currY = kine.yMaxZ;
        currMat = [currX', currY'];
        currMat = currMat(~isnan(sum(currMat, 2)), :);
        currCov = cov(currMat);
        edges = quantile(currT, linspace(0, 1, numBins+1));
        edges(1) = edges(1) - 0.0001;
        edges(end) = edges(end) + 0.0001;
        % edges = linspace(min(currT)-0.001, max(currT)+0.001, numBins+1);
        for b = 1:numBins
            currMT = mean(currT(currT > edges(b) & currT <= edges(b+1)), "all", 'omitmissing');
            matAllT(i, 1, b, sess) = currMT;
                   
            % mah
            tempX = currX(currT > edges(b) & currT <= edges(b+1));
            tempY = currY(currT > edges(b) & currT <= edges(b+1));
            disM = mahDis([tempX', tempY'], 'cov', currCov);
            matAllK(i, 1, b, sess) = mean(disM); 

        end
    end

end
%%
figure;
for i = 1:length(varT)
    currMatT = squeeze(matAllT(i, 1, :, :));
    currMatK = squeeze(matAllK(i, 1, :, :));

    currMatTM = mean(currMatT, 2, 'omitmissing')';
    currMatKM = mean(currMatK, 2, 'omitmissing')';
    currMatKSem = sem(currMatK');

    subplot(1, length(varT), i);
    hold on;
    plot(currMatTM, currMatKM, 'Color', 'k', 'LineWidth', 2);
    patch([currMatTM, flip(currMatTM)], [currMatKM-currMatKSem, flip(currMatKM+currMatKSem)], [0.5 0.5 0.5], 'FaceAlpha', 0.4,  'EdgeColor','none')
    xlabel(varT{i});
    ylabel('disM');
end
%% compare lasered and not lasered trials
numBins = 4;

varT = {'QdiffChosen'};
varK = {'maxS', 'vXFixInterC', 'xMaxC'};

matAllK = NaN(length(varT), length(varK), numBins, length(dayList));
matAllKL = NaN(length(varT), length(varK), numBins, length(dayList));
matAllT = NaN(length(varT), length(varK), numBins, length(dayList));

for sess = 1:length(dayList)
    session = dayList{sess};
    kine = lickPlotSessionNoPlot(session);
    
    pd = parseSessionString_df(session, root, sep);
    params = getStanModelParams_sampsOnly(pd.aniName, col, modelName, sampNum, 'sessionName', session, 'biasFlag',1);
    t = inferModelVar(session, params, modelName, 'perturb', 0);
    s = behAnalysisNoPlot_opMD(session, 'simpleFlag', 1);
    Qdiff = t.Q(:,2) - t.Q(:,1);
    Qsum = t.Q(:,2) + t.Q(:,1);
    Qsum(s.allChoices == 1) = zscore(Qsum(s.allChoices == 1));
    Qsum(s.allChoices == -1) = zscore(Qsum(s.allChoices == -1));
    QdiffChosen = Qdiff;
    QdiffChosen(s.allChoices==-1) = - QdiffChosen(s.allChoices==-1);
    pRight = t.probChoice;
    pRight(s.allChoices==1) = 1 - t.probChoice(s.allChoices==1);
    conf = t.probChoice;

    for i = 1:length(varT)
        eval(['currT = ' varT{i} ';']);
        for j = 1:length(varK)
            eval(['currK = kine.' varK{j} ';']);
            edges = quantile(currT, linspace(0, 1, numBins+1));
            edges(1) = edges(1) - 0.0001;
            edges(end) = edges(end) + 0.0001;
            for b = 1:numBins
                currMT = mean(currT(currT > edges(b) & currT <= edges(b+1)), "all", 'omitnan');
                currMK = mean(currK(currT > edges(b) & currT <= edges(b+1) & s.laser' == 0), "all", 'omitnan');
                currMKL = mean(currK(currT > edges(b) & currT <= edges(b+1) & s.laser' == 1), "all", 'omitnan');
                
                matAllT(i, j, b, sess) = currMT;
                matAllK(i, j, b, sess) = currMK;
                matAllKL(i, j, b, sess) = currMKL;
            end
        end
    end

end
%%
figure2;
for i = 1:length(varT)
    for j = 1:length(varK)
        currMatT = squeeze(matAllT(i, j, :, :));
        currMatK = squeeze(matAllK(i, j, :, :));

        currMatTM = mean(currMatT, 2, 'omitmissing')';
        currMatKM = mean(currMatK, 2, 'omitmissing')';
        currMatKSem = sem(currMatK');

        currMatK = squeeze(matAllKL(i, j, :, :));
        currMatKML = mean(currMatK, 2, 'omitmissing')';
        currMatKLSem = sem(currMatK');

        subplot(length(varK), length(varT), length(varT)*(j-1)+i);
        hold on;
        plot(currMatTM, currMatKM, 'Color', 'k', 'LineWidth', 2);
        patch([currMatTM, flip(currMatTM)], [currMatKM-currMatKSem, flip(currMatKM+currMatKSem)], [0.5 0.5 0.5], 'FaceAlpha', 0.4,  'EdgeColor','none')
        plot(currMatTM, currMatKML, 'Color', 'r', 'LineWidth', 2);
        patch([currMatTM, flip(currMatTM)], [currMatKML-currMatKLSem, flip(currMatKML+currMatKLSem)], [1 0.5 0.5], 'FaceAlpha', 0.4,  'EdgeColor','none')

        xlabel(varT{i});
        ylabel(varK{j});
    end
end
%%
