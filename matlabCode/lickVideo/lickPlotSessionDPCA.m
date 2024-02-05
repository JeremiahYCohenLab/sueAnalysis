function lickPlotSessionDPCA(session)
    [root, sep] = currComputer();
    pd = parseSessionString_df(session, root, sep);
    savepath = [pd.sortedFolder session '_tongue.mat'];
    if exist(savepath, "file")
        load(savepath);
    else
        fprintf([session ' no tongue. \n'])
        return
    end
    
    % load([pd.sortedFolder 'lickSession.mat'], 'lickSession');
    modelName = '5params'; 
    col = 'cueOnGood';
    sampNum = 2000;
    [root,sep] = currComputer();
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
    
    % load kinematics
    kine = lickPlotSessionNoPlot(session);
    
    %%
    minPre = 10;
    minPost = 10;
    
    allLen = cellfun(@length, kine.allY);
    preLen = kine.maxInd-1;
    postLen = allLen - kine.maxInd';
    valTrajInd = preLen'>=minPre & postLen>=minPost;
    
    peakTrajRaw = NaN(length(kine.allY), 2*(minPre+minPost+1));
    peakTrajZ = NaN(length(kine.allY), 2*(minPre+minPost+1));
    for i = 1:length(kine.allY)
        if valTrajInd(i)
            if s.allChoices(i) == 1
                currX = -kine.allX{i};
            else
                currX = kine.allX{i};
            end
            currY = kine.allY{i};
            currTraj = [currX(kine.maxInd(i)-minPre:kine.maxInd(i)+minPost); currY(kine.maxInd(i)-minPre:kine.maxInd(i)+minPost)]';
            % currTraj = [currX(kine.maxInd(i)-minPre:kine.maxInd(i)+minPost) - currX(kine.maxInd(i)); currY(kine.maxInd(i)-minPre:kine.maxInd(i)+minPost)- currY(kine.maxInd(i))]';

            peakTrajRaw(i, :) = currTraj;
        end
    end
    peakTraj = peakTrajRaw;
    peakTraj(valTrajInd & [s.allChoices]' == 1, 1:minPre+minPost+1) = zscore(peakTraj(valTrajInd & [s.allChoices]' == 1, 1:minPre+minPost+1), [], 'all');
    peakTraj(valTrajInd & [s.allChoices]' == -1, 1:minPre+minPost+1) = zscore(peakTraj(valTrajInd & [s.allChoices]' == -1, 1:minPre+minPost+1), [], 'all');
    peakTraj(valTrajInd & [s.allChoices]' == 1, minPre+minPost+2:2*(minPre+minPost+1)) = zscore(peakTraj(valTrajInd & [s.allChoices]' == 1, minPre+minPost+2:2*(minPre+minPost+1)), [], 'all');
    peakTraj(valTrajInd & [s.allChoices]' == -1, minPre+minPost+2:2*(minPre+minPost+1)) = zscore(peakTraj(valTrajInd & [s.allChoices]' == -1, minPre+minPost+2:2*(minPre+minPost+1)), [], 'all');
    peakTrajZ(valTrajInd,:) = zscore(peakTraj(valTrajInd, :), 0, 1);
    
    [coeff,scores,latent,tsquared,explained,mu] = pca(peakTraj);
    scores(:, 7:end) = 0;
    peakTrajPCA = scores* inv(coeff) + mu;
    %%  linear model fit 
    laser = s.laser';
    conf = t.probChoice;
    vars = {'QdiffChosen', 'conf', 'Qsum'};
    fig = figure;
    for i = 1:length(vars)
        eval(['currVar = ' vars{i} ';'])
        lm = fitlm(scores(:,1:6), currVar);
        currCoeff = lm.Coefficients.Estimate(2:end);
        varCoeff = coeff(:, 1:6) * currCoeff;
        subplot(4, length(vars), [i, i + length(vars)]);
        imagesc(varCoeff); 
        maxPos = max(varCoeff, [], 'all');
        maxNeg = -(min(varCoeff, [], 'all'));
        
        myMap = [[linspace(0, 1, 100)', linspace(0, 1, 100)', linspace(1, 1, 100)']; ...
                [linspace(1, 1, 100*(maxPos/maxNeg))', linspace(1, 0, 100*(maxPos/maxNeg))', linspace(1, 0, 100*(maxPos/maxNeg))']];

        colorbar
        colormap(myMap)
        title(vars{i})  

        
        numBins = 4;
        colors = cool(numBins);

        eval(['target = ' vars{i} ';'])
        edges = quantile(target, linspace(0, 1, numBins+1));
        edges(1) = edges(1) + 0.001;
        edges(end) = edges(end) + 0.001;

        subplot(4, length(vars), [i + 2*length(vars)]); hold on;
        colors = [linspace(0, 1, numBins); linspace(0, 0, numBins); linspace(0, 0, numBins)]';
        for j = 1:numBins
            currInd = target>edges(j) & target<=edges(j+1); 
            meanTraj = mean(peakTrajPCA(currInd, :), 'omitmissing');
            % semTraj = sem(peakTraj(currInd, :));
            plot(meanTraj(1:minPost+minPre+1), meanTraj(minPost+minPre+2:2*(minPost+minPre+1)), 'Color', colors(j, :), 'LineWidth', 3);
            scatter(meanTraj(1), meanTraj(0.5*length(meanTraj)+1), 10, 'r', 'filled');
        end

        subplot(4, length(vars), [i + 3*length(vars)]); hold on;
        yPred =   lm.Fitted;
        edges = linspace(min(yPred), max(yPred), 10);
        histogram(yPred(s.laser==0), edges, 'Normalization', 'probability', 'FaceColor', 'k', 'FaceAlpha', 0.25);
        histogram(yPred(s.laser==1), edges, 'Normalization', 'probability', 'FaceColor', 'b', 'FaceAlpha', 0.25);
         
        legend({'none', 'laser'})
    end
    screen = get(0,'Screensize');
    screen(4) = screen(4) - 100;
    set(fig, 'Position', screen);
    sgtitle(session);
    saveFigurePDF(fig, [pd.saveFigFolder sep session '_lickTrackVarProj.pdf'])
end