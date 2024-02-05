function lickPlotSessionPCA(session)
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
    [coeff,scores,latent,tsquared,explained,mu] = pca(peakTrajZ);

    %% plot feature distribution
    pFig = figure;
    subplot(4,6,[1,2]); hold on;
    plot(cumsum(explained));
    scatter(1:length(mu), cumsum(explained), 10, 'filled');
    plot([0 length(mu)], [95 95]);
    subplot(4,6,[3, 4, 5, 9, 10, 11]); hold on;
    imagesc(coeff);
    maxPos = max(coeff, [], 'all');
    maxNeg = abs(min(coeff, [], 'all'));
    
    myMap = [[linspace(0, 1, 100)', linspace(0, 1, 100)', linspace(1, 1, 100)']; ...
            [linspace(1, 1, 100*(maxPos/maxNeg))', linspace(1, 0, 100*(maxPos/maxNeg))', linspace(1, 0, 100*(maxPos/maxNeg))']];
    xlim([0, length(mu)]);
    ylim([0, length(mu)]);
    colorbar
    colormap(myMap)
    %% plot score distribution 
    figure(pFig);
    subplot(4, 6, 7); hold on;
    scatter(scores(s.allChoices == 1,1), scores(s.allChoices == 1,2), 14, 'filled');
    scatter(scores(s.allChoices == -1,1), scores(s.allChoices == -1,2), 14, 'filled');
    xlabel('1')
    ylabel('2')
    legend({'R', 'L'})
    subplot(4, 6, 8); hold on;
    scatter(scores(s.allChoices == 1,3), scores(s.allChoices == 1,4), 14, 'filled');
    scatter(scores(s.allChoices == -1,3), scores(s.allChoices == -1,4), 14, 'filled');
    xlabel('3')
    ylabel('4')

    %% plot PCA features;
    numBins = 6;
    numFs = 6;
    colors = cool(numBins);
    for i = 1:numFs
        target = scores(:,i);
        edges = quantile(target, linspace(0, 1, numBins+1));
        edges(1) = edges(1) + 0.001;
        edges(end) = edges(end) + 0.001;
        subplot(4, numFs, 2*(numFs)+i); hold on;
        for j = 1:numBins
            currInd = target>edges(j) & target<=edges(j+1); 
            meanTraj = mean(peakTraj(currInd, :));
            % semTraj = sem(peakTraj(currInd, :));
            plot(meanTraj(1:minPost+minPre+1), meanTraj(minPost+minPre+2:2*(minPost+minPre+1)), 'Color', colors(j, :), 'LineWidth', 3);
            scatter(meanTraj(1), meanTraj(0.5*length(meanTraj)+1), 10, 'r', 'filled');
        end
        
    end
    sgtitle(session)
    screen = get(0,'Screensize');
    screen(4) = screen(4) - 100;
    set(pFig, 'Position', screen);

    %% plot PCA with behavior values
    laser = s.laser';
    tFig = figure;
    vars = {'QdiffChosen', 't.probChoice', 'Qsum', 'laser'};

    numBins = 6;
    for i = 1:numFs
        target = scores(:,i);
        edges = quantile(target, linspace(0, 1, numBins+1));
        edges(1) = edges(1) + 0.001;
        edges(end) = edges(end) + 0.001;
        subplot(length(vars)+1, numFs, i); hold on;
        for j = 1:numBins
            hold on
            currInd = target>edges(j) & target<=edges(j+1); 
            meanTraj = mean(peakTraj(currInd, :));
            % semTraj = sem(peakTraj(currInd, :));
            plot(meanTraj(1:minPost+minPre+1), meanTraj(minPost+minPre+2:2*(minPost+minPre+1)), 'Color', colors(j, :), 'LineWidth', 3);
            scatter(meanTraj(1), meanTraj(0.5*length(meanTraj)+1), 10, 'r', 'filled');    
        end
    
    end

    target = cell(length(vars), 1);
    for i = 1:length(vars)
        eval(['target{' num2str(i) '}=' vars{i} ';']);
    end

    numBins = 2;
    for i = 1:length(vars)
        currVar = target{i};
        edges = quantile(currVar, linspace(0, 1, numBins+1));
        edges(1) = edges(1) - 0.0001;
        edges(end) = edges(end) + 0.0001;
        for j = 1:numFs
            subplot(length(vars)+1, numFs, numFs*(i) + j);
            currM = NaN(1, numBins);
            currSem = NaN(1, numBins);
            currTM = NaN(1, numBins);
            for b = 1:numBins
                currInd = currVar>edges(b) & currVar<=edges(b+1);
                currM(b) = mean(scores(currInd, j), "omitmissing");
                currSem(b) = sem(scores(currInd, j));
                currTM(b) = mean(currVar(currInd), "omitmissing");
            end
            plot(currTM, currM, 'Color', 'k', 'LineWidth', 2); hold on;
            patch([currTM, flip(currTM)], [currM-currSem, flip(currM + currSem)], 'k', 'FaceAlpha', 0.25, 'EdgeColor', 'none'); 
            xlabel(vars{i});
            ylabel(['PC' num2str(j)]);

            tmpInd = ~isnan(currVar) & ~isnan(scores(:, j));
            currTarget = currVar(tmpInd);
            currKine = scores(tmpInd, j);
            [R, P] = corrcoef(currTarget, currKine);
            title(sprintf('R=%3.2f P=%3.2f', R(1,2), P(1,2)));
            
        end
    end
    sgtitle(session)
    screen = get(0,'Screensize');
    screen(4) = screen(4) - 100;
    set(tFig, 'Position', screen);

    saveFigurePDF(pFig, [pd.saveFigFolder sep session '_lickTrackPCA_laser.pdf'])
    saveFigurePDF(tFig, [pd.saveFigFolder sep session '_lickTuningPCA_laser.pdf'])
end