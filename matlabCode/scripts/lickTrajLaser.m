file = 'inhibitionGt';
sheet = 'allGt';
col = 'cueOnGood';
dayList = getDayList(file, sheet, col);
[root, sep] = currComputer();
minPre = 10;
minPost = 10;
modelName = '5params'; 
sampNum = 2000;
vars = {'QdiffChosen', 'conf', 'Qsum'}; 
%%
allMeansControl = NaN(length(dayList), length(vars));
allMeansLaser = NaN(length(dayList), length(vars));
allSemsControl = NaN(length(dayList), length(vars));
allSemsLaser = NaN(length(dayList), length(vars));
allSigs = NaN(length(dayList), length(vars));
aucCurve = NaN(length(dayList), 3);
aucSig = NaN(length(dayList), 3);
for sess = 1:length(dayList)
    session = dayList{sess};
    pd = parseSessionString_df(session, root, sep);
    savepath = [pd.sortedFolder session '_tongue.mat'];
    if exist(savepath, "file")
        load(savepath);
    else
        fprintf([session ' no tongue. \n'])
        return
    end
    
    % load([pd.sortedFolder 'lickSession.mat'], 'lickSession');
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
    for i = 1:length(vars)
        eval(['currVar = ' vars{i} ';'])
        lm = fitlm(scores(:,1:6), currVar);
        currCoeff = lm.Coefficients.Estimate(2:end);
        varCoeff = coeff(:, 1:6) * currCoeff;
      
  
        yPred = lm.Fitted;
        allMeansLaser(sess, i) = mean(yPred(s.laser == 1), 'omitmissing');
        allMeansControl(sess, i)  = mean(yPred(s.laser == 0), 'omitmissing');
        allSemsLaser(sess, i) = sem(yPred(s.laser == 1));
        allSemsControl(sess, i)  = sem(yPred(s.laser == 0));
        [h] = ttest2(yPred(s.laser == 1), yPred(s.laser == 0));
        allSigs(sess, i) = h;

    end
    glm = fitglm(scores(:,1:6), s.laser, 'Distribution', 'binomial', 'Link','logit');
    score = glm.Fitted.Probability;
    [Xglm,Yglm,~,AUCglm] = perfcurve(s.laser, score, 1, 'NBoot',500, 'Alpha', 0.05);
    aucCurve(sess, :) = AUCglm;
    aucSig(sess) = (AUCglm(2)-0.5) * (AUCglm(3)-0.5) > 0;
end 

%% plot everything

for i = 1:length(vars)
    figure2;
    hold on;
    errorbar(allMeansControl(:, i), allMeansLaser(:,i), allSemsLaser(:,i), allSemsLaser(:,i), allSemsControl(:,i), allSemsControl(:,i), 'LineStyle', 'none', 'Marker','o', 'MarkerSize', 10, 'MarkerEdgeColor', 'k')
    scatter(allMeansControl(allSigs(:,i)==1, i), allMeansLaser(allSigs(:,i)==1, i), 15, 'red', 'filled', 'o')
    title(vars{i})
    plot(minmax(allMeansControl(:, i)), minmax(allMeansControl(:, i)), 'LineStyle','--', 'Color', [0.5 0.5 0.5])
    xlabel('Control')
    ylabel('Laser')
end
%%

[~, indAUC] = sort(aucCurve(:,1));
sortedAUC = aucCurve(indAUC, :);  
sortedSig = aucSig(indAUC);
figure2; hold on;
errorbar(1:length(dayList), sortedAUC(:,1), sortedAUC(:,1) - sortedAUC(:,2), sortedAUC(:,3) - sortedAUC(:,1), 'LineStyle', 'none', 'Marker','o')
plot([1, length(dayList)], [0.5 0.5], 'LineWidth', 2, 'LineStyle', '--')
scatter(find(sortedSig == 1), sortedAUC(sortedSig==1, 1), 10, 'red', 'filled');
xlabel('session')
ylabel('AUC')
%%
