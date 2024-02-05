function lickPlotSession(session)
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

%%
ind  = 1:length(allLicks);
meanX = NaN(1, length(ind));
meanY = NaN(1, length(ind));
meanXPre = NaN(1, length(ind));
meanXMax = NaN(1, length(ind));
meanYMax = NaN(1, length(ind));
maxS = NaN(1, length(ind));
speedAtPort = NaN(1, length(ind));
vXAtPort = NaN(1, length(ind));
vYAtPort = NaN(1, length(ind));
vXPrePeak = NaN(1, length(ind));
vYPrePeak = NaN(1, length(ind));
vXFixInter = NaN(1, length(ind));
vYFixInter = NaN(1, length(ind));
speedFixInter = NaN(1, length(ind));
speedPrePeak = NaN(1, length(ind));
crossInd = NaN(1, length(ind));
maxIndAll = NaN(1, length(ind));
crossX = NaN(1, length(ind));
crossY = NaN(1, length(ind));
allX = cell(length(ind), 1);
allY = cell(length(ind), 1);
thresh = 0;
preBins = 6;
colors = cool(length(s.allChoices));
lickInd = ~cellfun(@isempty, {allLicks.decisionID});
for j = 1:length(ind)
    if allLicks(ind(j)).decisionID
        decisionID = allLicks(ind(j)).decisionID;
        interval = allLicks(ind(j)).time(2) - allLicks(ind(j)).time(1);
        disRL = norm(allLicks(ind(j)).portL - allLicks(ind(j)).portR);
        temp = allLicks(ind(j)).allLicksFilt{decisionID};
        tempX = (temp(:,1) - allLicks(ind(j)).midPoint(1)) * 4/disRL;
        tempY = -(temp(:,2) - allLicks(ind(j)).midPoint(2)) * 4/disRL;
        allX{j} = tempX;
        allY{j} = tempY;
        [meanYMax(ind(j)), maxInd] = max(tempY);
        maxIndAll(ind(j)) = maxInd;
        % fprintf([num2str(maxInd - length(tempY)) '\n']);
        meanXMax(ind(j)) = tempX(maxInd);
        meanXPre(ind(j)) = mean(tempX(1:maxInd));
        
        meanX(ind(j)) = mean(tempX);
        meanY(ind(j)) = mean(tempY);
        time = allLicks(ind(j)).time(allLicks(ind(j)).windows(decisionID,1): allLicks(ind(j)).windows(decisionID,2));

        lickDisTmpX = diff(tempX);
        lickDisTmpY = diff(tempY);
        lickSpeed = sqrt(lickDisTmpX.^2 + lickDisTmpY.^2)/(interval);
        maxS(ind(j)) = max(lickSpeed);
        lickVX = lickDisTmpX/(allLicks(ind(j)).time(2) - allLicks(ind(j)).time(1));
        lickVY = lickDisTmpY/(allLicks(ind(j)).time(2) - allLicks(ind(j)).time(1));
       
        if maxInd >= preBins+1 && maxInd < length(tempX)-2
            vXPrePeak(ind(j)) = (tempX(maxInd) - tempX(maxInd-preBins))/(preBins*interval);
            vYPrePeak(ind(j)) = (tempY(maxInd) - tempY(maxInd-preBins))/(preBins*interval);
            speedPrePeak(ind(j)) = mean(lickSpeed(maxInd-1:maxInd));
        end
        % find threshold crossing
        cross = tempY >thresh;
        cross = find(cross(1:end-1)==0 & cross(2:end)==1, 1);
        if ~isempty(cross)
            crossInd(ind(j)) = cross;
            speedAtPort(ind(j)) = lickSpeed(cross);
            vXAtPort(ind(j)) = lickVX(cross);
            vYAtPort(ind(j)) = lickVY(cross);
            crossX(ind(j)) = mean(tempX(cross));
            crossY(ind(j)) = mean(tempY(cross));
        end
        
        % use fixed inter
        if ~isempty(cross) && maxInd >= preBins+1 && maxInd < length(tempX)-2
            vXFixInter(ind(j)) = (tempX(maxInd) - tempX(cross))/((maxInd-cross)*interval);
            vYFixInter(ind(j)) = (tempY(maxInd) - tempY(cross))/((maxInd-cross)*interval);
            speedFixInter(ind(j)) = norm([vXFixInter(ind(j)), vYFixInter(ind(j))]);
        end

    end
end
%% regularization
vXFixInter(s.allChoices==1 & ~isnan(vXFixInter)) = zscore(vXFixInter(s.allChoices==1 & ~isnan(vXFixInter)));
vXFixInter(s.allChoices==-1 & ~isnan(vXFixInter)) = zscore(vXFixInter(s.allChoices==-1 & ~isnan(vXFixInter)));
vXFixInterAbs = vXFixInter;
vXFixInterAbs(s.allChoices==1) = -vXFixInterAbs(s.allChoices==1);

meanXMax(~isnan(meanXMax) & s.allChoices==1) = meanXMax(~isnan(meanXMax) & s.allChoices==1)/abs(mean(meanXMax(~isnan(meanXMax) & s.allChoices==1), 'all', 'omitnan'));
meanXMax(~isnan(meanXMax) & s.allChoices==-1) = meanXMax(~isnan(meanXMax) & s.allChoices==-1)/abs(mean(meanXMax(~isnan(meanXMax) & s.allChoices==-1), 'all', 'omitnan'));


meanXMaxZ(~isnan(meanXMax) & s.allChoices==1) = zscore(meanXMax(~isnan(meanXMax) & s.allChoices==1));
meanXMaxZ(~isnan(meanXMax) & s.allChoices==-1) = zscore(meanXMax(~isnan(meanXMax) & s.allChoices==-1));

meanXMaxC = meanXMaxZ;
meanXMaxC(s.allChoices==1) = - meanXMaxC(s.allChoices==1);

meanYMax(~isnan(meanYMax) & s.allChoices==1) = meanYMax(~isnan(meanYMax) & s.allChoices==1)/abs(mean(meanYMax(~isnan(meanYMax) & s.allChoices==1), 'all', 'omitnan'));
meanYMax(~isnan(meanYMax) & s.allChoices==-1) = meanYMax(~isnan(meanYMax) & s.allChoices==-1)/abs(mean(meanYMax(~isnan(meanYMax) & s.allChoices==-1), 'all', 'omitnan'));

vYFixInter(s.allChoices==1 & ~isnan(vYFixInter)) = zscore(vYFixInter(s.allChoices==1 & ~isnan(vYFixInter)));
vYFixInter(s.allChoices==-1 & ~isnan(vYFixInter)) = zscore(vYFixInter(s.allChoices==-1 & ~isnan(vYFixInter)));

speedFixInter(s.allChoices==1 & ~isnan(speedFixInter)) = zscore(speedFixInter(s.allChoices==1 & ~isnan(speedFixInter)));
speedFixInter(s.allChoices==-1 & ~isnan(speedFixInter)) = zscore(speedFixInter(s.allChoices==-1 & ~isnan(speedFixInter)));

maxS(s.allChoices==1 & ~isnan(maxS)) = zscore(maxS(s.allChoices==1 & ~isnan(maxS)));
maxS(s.allChoices==-1 & ~isnan(maxS)) = zscore(maxS(s.allChoices==-1 & ~isnan(maxS)));

meanXMax = -meanXMax;


%% plot trajectories
fig = figure;
[~, indQ] = sort(Qsum);
[~, indConf] = sort(t.probChoice);
[~, indX] = sort(vXFixInterAbs(~isnan(vXFixInterAbs)));
[~, indY] = sort(vYFixInter(~isnan(vYFixInter)));
[~, indS] = sort(maxS(~isnan(maxS)));
colors = cool(length(indQ));
colorsX = cool(length(indX));
colorsY = cool(length(indY));
colorsS = cool(length(indS));

for j = 1:length(s.allChoices)
    currX = -allX{j};
    currY = allY{j};

    subplot(4, 4, 1); hold on;
    newInd = find(indQ == j);
    scatter(currX, currY, 5, colors(newInd, :), "filled", 'MarkerFaceAlpha', 0.25, 'MarkerEdgeColor','none');
    title('Qsum')

    subplot(4, 4, 2); hold on;
    newInd = find(indConf == j);
    scatter(currX, currY, 5, colors(newInd, :), 'filled', 'MarkerFaceAlpha', 0.25, 'MarkerEdgeColor','none')
    title('Conf')

    if ~isnan(vXFixInterAbs(j))
        subplot(4, 4, 3); hold on;
        title('vXC')
        newInd = sum(~isnan(vXFixInterAbs(1:j)));
        newInd = indX == newInd;
        scatter(currX, currY, 5, colorsX(newInd, :), 'filled', 'MarkerFaceAlpha', 0.25, 'MarkerEdgeColor','none')
    end

    if ~isnan(vYFixInter(j))
        subplot(4, 4, 4); hold on;
        title('vY')
        newInd = sum(~isnan(vYFixInter(1:j)));
        newInd = indY == newInd;
        scatter(currX, currY, 5, colorsY(newInd, :), 'filled', 'MarkerFaceAlpha', 0.25, 'MarkerEdgeColor','none')
    end


end

for j = 1:length(allLicks)
    currX = -allX{j};
    currY = allY{j};
    if ~isnan(maxS(j))
        subplot(4, 4, [9:12]); hold on;
        title('speed')
        newInd = sum(~isnan(maxS(1:j)));
        newInd = indS == newInd;
        qInd = find(indQ == j);
        if s.allChoices(j) > 0
            plot(0.3 * qInd + currX, currY, 'color', colorsS(newInd, :));
        else
            plot(0.3 * qInd -currX, currY, 'color', colorsS(newInd, :))
        end
    end
end

%% mean trajectory
numBins = 5;
minPre = 10;
minPost = 10;

meanXTraj = NaN(numBins, minPre + minPost + 1);
meanYTraj = NaN(numBins, minPre + minPost + 1);
meanMaxS = NaN(numBins, 1);
meanQ = NaN(numBins, 1);

allLen = cellfun(@length, allY);
preLen = maxIndAll-1;
postLen = allLen - maxIndAll';
valTrajInd = preLen'>=minPre & postLen>=minPost;
edges = quantile(Qsum, linspace(0, 1, numBins+1));
edges(1) = edges(1) - 0.0001;
edges(end) = edges(end) + 0.0001;
colorsQS = cool(numBins);
figure(fig); hold on;
for i = 1:numBins
    currInd = Qsum > edges(i) & Qsum <= edges(i+1);
    currMaxInd = maxIndAll(currInd & valTrajInd);
    currX = allX(currInd & valTrajInd);
    currY = allY(currInd & valTrajInd);
    currChoice = s.allChoices(currInd & valTrajInd);
    currMaxS = maxS(currInd & valTrajInd);
    currQsum = Qsum(currInd & valTrajInd);
    lickMatX = NaN(length(currX), minPre+minPost+1);
    lickMatY = NaN(length(currX), minPre+minPost+1);
    for j = 1:length(currX)
        if currChoice(j) == 1
            tempX = -currX{j};
        else
            tempX = currX{j};
        end
        tempY = currY{j};
        lickMatX(j, :) = tempX(currMaxInd(j)-minPre:currMaxInd(j)+minPost);
        lickMatY(j, :) = tempY(currMaxInd(j)-minPre:currMaxInd(j)+minPost);
    end
    meanXTraj(i,:) = mean(lickMatX);
    meanYTraj(i,:) = mean(lickMatY);
    meanMaxS(i) = mean(currMaxS);
    meanQ(i) = mean(currQsum);
end
[~, speedInd] = sort(meanMaxS);
[~, QInd] = sort(meanQ);
for j = 1:length(speedInd)
    i = speedInd(j);
    Qrank = find(QInd==i);
    subplot(4, 4, 13:16); hold on;
    plot(meanXTraj(i, :)+0.2*Qrank, meanYTraj(i, :), 'color', colorsQS(j,:), 'LineWidth', 2);
end
%% feature extraction of tongue traj
%% mean trajectory
minPre = 10;
minPost = 10;

allLen = cellfun(@length, allY);
preLen = maxIndAll-1;
postLen = allLen - maxIndAll';
valTrajInd = preLen'>=minPre & postLen>=minPost;

peakTraj = NaN(length(allY), 2*(minPre+minPost+1));

for i = 1:length(allY)
    if valTrajInd(i)
        if s.allChoices(i) == 1
            currX = -allX{i};
        else
            currX = allX{i};
        end
        currY = allY{i};
        currTraj = [currX(maxIndAll(i)-minPre:maxIndAll(i)+minPost); currY(maxIndAll(i)-minPre:maxIndAll(i)+minPost)]';
        peakTraj(i, :) = currTraj;
    end
end
peakTraj(valTrajInd & [s.allChoices]' == 1, 1:minPre+minPost+1) = zscore(peakTraj(valTrajInd & [s.allChoices]' == 1, 1:minPre+minPost+1), [], 'all');
peakTraj(valTrajInd & [s.allChoices]' == -1, 1:minPre+minPost+1) = zscore(peakTraj(valTrajInd & [s.allChoices]' == -1, 1:minPre+minPost+1), [], 'all');
peakTraj(valTrajInd & [s.allChoices]' == 1, minPre+minPost+2:2*(minPre+minPost+1)) = zscore(peakTraj(valTrajInd & [s.allChoices]' == 1, minPre+minPost+2:2*(minPre+minPost+1)), [], 'all');
peakTraj(valTrajInd & [s.allChoices]' == -1, minPre+minPost+2:2*(minPre+minPost+1)) = zscore(peakTraj(valTrajInd & [s.allChoices]' == -1, minPre+minPost+2:2*(minPre+minPost+1)), [], 'all');
figure2;
peakTraj(valTrajInd,:) = zscore(peakTraj(valTrajInd, :), 0, 1);
[coeff,scores,latent,tsquared,explained,mu] = pca(peakTraj); 
%%
numBins = 10;
  
%% colorCodeQdiffC on maxPosition
figure(fig);

subplot(4, 4, 5);
[~, ind] = sort(Qsum);
scatter(meanXMax(ind), meanYMax(ind), 15, colors,  'filled', 'MarkerFaceAlpha', 0.75);
title('QdiffChosen')

subplot(4, 4, 6);
scatter(meanXMax(indConf), meanYMax(indConf), 15, colors,  'filled', 'MarkerFaceAlpha', 0.75);
title('conf')

subplot(4, 4, 7);
tmpX = meanXMax(~isnan(vXFixInterAbs));
tmpY = meanYMax(~isnan(vXFixInterAbs));
scatter(tmpX(indX), tmpY(indX), 15, colorsX,  'filled', 'MarkerFaceAlpha', 0.75);
title('vXC')

subplot(4, 4, 8);
tmpX = meanXMax(~isnan(vYFixInter));
tmpY = meanYMax(~isnan(vYFixInter));
scatter(tmpX(indY), tmpY(indY), 15, colorsY,  'filled', 'MarkerFaceAlpha', 0.75);
title('vY')

sgtitle(session);
%% tuning with location and velocity
figT = figure;
numBins = 3;
vars = {'Qdiff', 'QdiffChosen', 't.probChoice', 'Qsum'};
for i = 1:length(vars)
    eval(['target{' num2str(i) '}=' vars{i} ';']);
end
varsK = {'meanXMax', 'meanXMaxC', 'vXFixInterAbs', 'vYFixInter', 'speedFixInter', 'maxS'};
for i = 1:length(varsK)
    eval(['kines{' num2str(i) '}=' varsK{i} ';']);
end

for i = 1:length(target)
    for j = 1:length(kines)
        subplot(length(target), length(kines), length(kines)*(i-1)+j); hold on; 

        currTarget = target{i};
        currKine = kines{j};

        edges = linspace(min(currTarget)-0.001, max(currTarget)+0.001, numBins+1);
        edges = quantile(currTarget, linspace(0, 1, numBins+1));
        edges(1) = edges(1) - 0.0001;   
        edges(end) = edges(end) + 0.0001;
        meanTarget = NaN(1, numBins);
        meanKine = NaN(1, numBins);
        semKine = NaN(1, numBins);
        
        for k = 1:numBins
            meanTarget(k) = mean(currTarget(currTarget>=edges(k) & currTarget<edges(k+1)), 'omitnan');
            meanKine(k) = mean(currKine(currTarget>=edges(k) & currTarget<edges(k+1)), 'omitnan');
            semKine(k) = sem(currKine(currTarget>=edges(k) & currTarget<edges(k+1)));
        end
        
        plot(meanTarget, meanKine, 'Color', 'k', 'LineWidth', 2);
        patch([meanTarget, flip(meanTarget)], [meanKine-semKine, flip(meanKine+semKine)], [0.5 0.5 0.5], 'FaceAlpha', 0.4,  'EdgeColor','none')
        
        tmpInd = ~isnan(currTarget') & ~isnan(currKine);
        currTarget = currTarget(tmpInd);
        currKine = currKine(tmpInd);
        [R, P] = corrcoef(currTarget, currKine);
        title(sprintf('R=%3.2f P=%3.2f', R(1,2), P(1,2)));
        xlabel(vars{i})
        ylabel(varsK{j})
        % lm = fitlm(currTarget, currKine);
        % R = lm.Coefficients.Estimate(2);
        % P = lm.Coefficients.pValue(2);
        % title(sprintf('R=%3.2f P=%3.2f', R, P));
    end
end
sgtitle(session)
%%
figC = figure;
matC = [];
for j = 1:length(kines)
    matC(j,:) = kines{j};
end
scatterAll(matC', varsK, 10, 'k');
sgtitle(session);
%%
screen = get(0,'Screensize');
screen(4) = screen(4) - 100;
set(fig, 'Position', screen);
set(figT, 'Position', screen);
set(figC, 'Position', screen);

saveFigurePDF(fig, [pd.saveFigFolder sep session '_lickTrack.pdf'])
saveFigurePDF(figT, [pd.saveFigFolder sep session '_lickTuning.pdf'])
saveFigurePDF(figC, [pd.saveFigFolder sep session '_lickCorr.pdf'])
end 