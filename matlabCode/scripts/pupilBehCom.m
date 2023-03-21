%% load animal list
[root, sep] = currComputer();
[~, dayList, ~] = xlsread([root 'aniModel.xlsx'], 'all');
allAnis = dayList(2:11, 5);
allCols = dayList(2:11, 6);
allFile = dayList(2:11, 7);
allSheet = dayList(2:11, 8);


%% parameter distribution
modelName = '5params';
%% calculate MAP
paramNames = getParamNames_dF(modelName, 0);
aniParams = NaN(length(allAnis), length(paramNames));
aniMeanLL = NaN(length(allAnis), 1);
sessionParams = cell(length(allAnis), 1);
numBin = 50;
for i = 1:length(allAnis)
    animalName = allAnis{i};
    category = allCols{i};
    sampFile = [animalName category '_', modelName];
    path = [root animalName sep animalName 'sorted' sep 'stan' sep 'bernoulli' sep modelName sep category sep];
    modelPath = [path sampFile '.mat'];
    load(modelPath);
    eval(['samps = ' sampFile]);
    eval(['clear ' sampFile]);
    % get map of all params
    for j = 1:length(paramNames)
        tmp = samps.(['mu_' paramNames{j}]);
        noDInd = samps.divergent__~= 1;
        tmp = tmp(noDInd);
        [counts, edges] = histcounts(tmp, 50);
        [~, maxInd] = max(counts);
        if maxInd<length(counts)
            aniParams(i, j) = mean(tmp(tmp>=edges(maxInd)&tmp<edges(maxInd+1)));
        else
            aniParams(i, j) = mean(tmp(tmp>=edges(maxInd)&tmp<=edges(maxInd+1)));
        end
    end
    paramMatx = plotStanSessionParams({animalName}, 'modelName', '5params', 'beh', category, 'plotFlag', 0);
    sessionParams{i} = paramMatx;
end

%% lickRate
numSamps = 2000;
coeffsAni = NaN(1, length(allAnis));
tstatsAni = NaN(1, length(allAnis));
pValsAni = NaN(1, length(allAnis));
sessionCoeffs = cell(length(allAnis), 1);
sessionTstats = cell(length(allAnis), 1);
sessionpVals = cell(length(allAnis), 1);
for i = 1:length(allAnis)
    animalName = allAnis{i};
    category = allCols{i};
    sampFile = [animalName category '_', modelName];
    path = [root animalName sep animalName 'sorted' sep 'stan' sep 'bernoulli' sep modelName sep category sep];
    modelPath = [path sampFile '.mat'];
    load(modelPath, 'dayList');
    coeffsTemp = NaN(length(dayList), 1);
    tStatsTemp = NaN(length(dayList), 1);
    pValsTemp = NaN(length(dayList), 1);
    lickRateCmb = [];
    negRPECmb = [];
    for j = 1:length(dayList)
        s = behAnalysisNoPlot_opMD(dayList{j}, 'simpleFlag', 1);
        lickRateCurr = zscore(s.lickNumRwd(s.allRewards==0));
        params = getStanModelParams_sampsOnly(animalName, category, modelName, numSamps, 'sessionName', dayList{j}, 'sessionParamsFlag', 1);
        t = inferModelVar(dayList{j}, params, modelName, 'perturb', 0);
        negRPECurr = zscore(t.pe(s.allRewards==0));
        lmTmp = fitlm(negRPECurr, lickRateCurr);
        coeffsTemp(j,:) = lmTmp.Coefficients.Estimate(2);
        tStatsTemp(j,:) = lmTmp.Coefficients.tStat(2);
        pValsTemp(j,:) = lmTmp.Coefficients.pValue(2);
        lickRateCmb = [lickRateCmb; lickRateCurr'];
        negRPECmb = [negRPECmb; negRPECurr];
    end
    lm = fitlm(negRPECmb, lickRateCmb);
    coeffsAni(i) = lm.Coefficients.Estimate(2);
    tstatsAni(i) = lm.Coefficients.tStat(2);
    pValsAni(i) = lm.Coefficients.pValue(2);
    
    sessionCoeffs{i} = coeffsTemp;
    sessionTstats{i} = tStatsTemp;
    sessionpVals{i} = pValsTemp;
end
%% Pupil
maxSteps = 30;
numSamps = 2000;
coeffsPupilAni = NaN(maxSteps, 2, length(allAnis));
tstatsPupilAni = NaN(maxSteps, 2, length(allAnis));
pValsPupilAni = NaN(maxSteps, 2, length(allAnis));
sessionPupilCoeffs = cell(length(allAnis), 1);
sessionPupilTstats = cell(length(allAnis), 1);
sessionPupilpVals = cell(length(allAnis), 1);
%% 
for i = 1:length(allAnis)
    animalName = allAnis{i};
    category = allCols{i};
    sampFile = [animalName category '_', modelName];
    path = [root animalName sep animalName 'sorted' sep 'stan' sep 'bernoulli' sep modelName sep category sep];
    modelPath = [path sampFile '.mat'];
    load(modelPath, 'dayList');
    coeffsTemp = NaN(maxSteps, 2, length(dayList));
    tStatsTemp = NaN(maxSteps, 2, length(dayList));
    pValsTemp = NaN(maxSteps, 2, length(dayList));
    outcomeCmb = [];
    QchosenCmb = [];
    pupilMaxCmb = [];
    for j = 1:length(dayList)
        s = behAnalysisNoPlot_opMD(dayList{j}, 'simpleFlag', 1);
        outcomeCurr = abs(s.allRewards);
        % load pupil
        pd = parseSessionString_df(dayList{j}, root, sep);
        pupilPath = [pd.sortedFolder dayList{j} '_pupil.mat'];
        if exist(pupilPath, 'file')
            load(pupilPath);
            if errorProp < 0.7 && sum(qualInd)>=10
                % behavior
                params = getStanModelParams_sampsOnly(animalName, category, modelName, numSamps, 'sessionName', dayList{j}, 'sessionParamsFlag', 1);
                t = inferModelVar(dayList{j}, params, modelName, 'perturb', 0);
                QchosenCurr = t.Q(:,2);
                QchosenCurr(s.allChoices<=0) = t.Q(s.allChoices<=0,1);
                % pupil
                stepSize = 4; % in frames
                binSize = 10;  % in frames
                midpoints = round(0.5*binSize+1):stepSize:size(sessionPupilCue,2)-0.5*binSize+1;
                midpoints = midpoints(1:maxSteps);
                pupilSlide = zeros(size(sessionPupilCue,1), length(midpoints));

                for a = 1:length(midpoints)
                    pupilSlide(:,a) = mean(sessionPupilCue(:,midpoints(a)-0.5*binSize:midpoints(a)+0.5*binSize-1),2,'omitnan');
                end
                coeff = NaN(length(midpoints), 2);
                tstat = NaN(length(midpoints), 2);
                sigs = NaN(length(midpoints), 2);
                combineMat = [outcomeCurr', QchosenCurr];
                for a = 1:length(midpoints)
                    lm = fitlm(combineMat(ismember(s.responseInds, find(qualInd>0)),:),pupilSlide(s.responseInds(ismember(s.responseInds, find(qualInd>0))),a));
                    coeff(a,:) = lm.Coefficients.Estimate(2:end);
                    tstat(a,:) = lm.Coefficients.tStat(2:end);
                    sigs(a,:) = double(lm.Coefficients.pValue(2:end)<0.05);
                end
                

                coeffsTemp(:,:,j) = coeff;
                tStatsTemp(:,:,j) = tstat;
                pValsTemp(:,:,j) = sigs;   
                
                QchosenCmb = [QchosenCmb; combineMat(ismember(s.responseInds, find(qualInd>0)), 2)];
                outcomeCmb = [outcomeCmb; combineMat(ismember(s.responseInds, find(qualInd>0)), 1)];
                pupilMaxCmb = [pupilMaxCmb; pupilSlide(s.responseInds(ismember(s.responseInds, find(qualInd>0))),:)];
            else
                fprintf([dayList{j} ' not well aligned \n'])
            end
           
        else
            fprintf([dayList{j} ' no pupil \n'])           
        end
    end
    
    for a = 1:size(pupilMaxCmb, 2)
        lm = fitlm([QchosenCmb, outcomeCmb], pupilMaxCmb(:,a));
        coeffsPupilAni(a,:,i) = lm.Coefficients.Estimate(2:end);
        tstatsPupilAni(a,:,i) = lm.Coefficients.tStat(2:end);
        pValsPupilAni(a,:,i) = double(lm.Coefficients.pValue(2:end)<0.05);
    end
    sessionPupilCoeffs{i} = coeffsTemp;
    sessionPupilTstats{i} = tStatsTemp;
    sessionPupilpVals{i} = pValsTemp;
end
%% plot animal by animal analysis
for i = 1:length(allAnis)
    tstatsCurr = sessionPupilTstats{i};
    time = linspace(-2, 10, 60);
    time = time(1:maxSteps);
    figure2;
    % outcome
    for t = 1:maxSteps
        subplot(5, maxSteps, t); hold on;
        histogram(tstatsCurr(t, 1, :), 10, 'FaceColor', [0.6 0.6 0.6], 'EdgeColor', 'none');
        plot([0 0], [0 5], 'r')
        if t==1
            ylabel('outcome')
        end
    end
    % Qchosen
    for t = 1:maxSteps
        subplot(5, maxSteps, t+maxSteps); hold on;
        histogram(tstatsCurr(t, 2, :), 10, 'FaceColor', [0.6 0.6 0.6], 'EdgeColor', 'none');
        plot([0 0], [0 5], 'r')
        if t==1
            ylabel('Qchosen')
        end
    end    
    
    % whole animal
    subplot(5, 1, 3); hold on;
    yyaxis left
    plot(time, tstatsPupilAni(:, 2, i), 'LineWidth', 2, 'Color', 'b');
    scatter(time(pValsPupilAni(:,2,i)==1), tstatsPupilAni(pValsPupilAni(:,2,i)==1, 2, i), 10, 'r')
    ylabel('outcome')
    yyaxis right
    plot(time, tstatsPupilAni(:, 1, i),'LineWidth', 2, 'Color', 'r');
    scatter(time(pValsPupilAni(:,1,i)==1), tstatsPupilAni(pValsPupilAni(:,1,i)==1, 1, i), 10, 'r')
    ylabel('Qchosen')    
    
    % corr with lickCoeff
    
    % find maxTstat
    [~, maxInd] = max(abs(tstatsPupilAni(11:end, 2, i)));
    maxInd = maxInd + 10;
    subplot(5, 1, 3); hold on;
    plot([time(maxInd), time(maxInd)], [-10 10], 'm', 'LineWidth', 3);
    statsMaxSession = tstatsCurr(maxInd, 1, :);
    subplot(5, 3, [10 13]); hold on;
    scatter(statsMaxSession, sessionTstats{i}, 15, 'filled');
    
    x = statsMaxSession;
    y = sessionTstats{i};
    y = y(~isnan(x));
    x = x(~isnan(x));
    [R, P] = corrcoef(x, y);
    
    title(['p=' num2str(P(2,1)) 'Coeff=' num2str(R(2,1)), 'n=' num2str(length(x))]);
    ylabel('lick to rpe')
    
    subplot(5, 3, 11); hold on;
    % find params
    paramsCurr = sessionParams{i};
    scatter(statsMaxSession, paramsCurr(:, 1)./(paramsCurr(:, 1) + paramsCurr(:, 2)), 15, 'filled');
    
    x = statsMaxSession;
    y = paramsCurr(:, 1)./(paramsCurr(:, 1) + paramsCurr(:, 2));
    y = y(~isnan(x));
    x = x(~isnan(x));
    [R, P] = corrcoef(x, y);
    
    title(['p=' num2str(P(2,1)) 'Coeff=' num2str(R(2,1)), 'n=' num2str(length(x))]);
    ylabel('aN/(aN + aP)')
    
    subplot(5, 3, 12); hold on;
    scatter(statsMaxSession, paramsCurr(:, 2)./(paramsCurr(:, 1) + paramsCurr(:, 2)), 15, 'filled');
    
    x = statsMaxSession;
    y = paramsCurr(:, 2)./(paramsCurr(:, 1) + paramsCurr(:, 2));
    y = y(~isnan(x));
    x = x(~isnan(x));
    [R, P] = corrcoef(x, y);
    
    title(['p=' num2str(P(2,1)) 'Coeff=' num2str(R(2,1)), 'n=' num2str(length(x))]);
    ylabel('aP/(aN + aP)')
    
    subplot(5, 3, 14); hold on;
    scatter(statsMaxSession, paramsCurr(:, 4).*(paramsCurr(:, 2)+paramsCurr(:,1)), 15, 'filled');
    
    x = statsMaxSession;
    y = paramsCurr(:, 4).*(paramsCurr(:, 2)+paramsCurr(:,1));
    y = y(~isnan(x));
    x = x(~isnan(x));
    [R, P] = corrcoef(x, y);
    
    title(['p=' num2str(P(2,1)) 'Coeff=' num2str(R(2,1)), 'n=' num2str(length(x))]);
    ylabel('beta/(aN + aP)')
    
    subplot(5, 3, 15); hold on;
    scatter(statsMaxSession, paramsCurr(:, 4), 15, 'filled');
    
    x = statsMaxSession;
    y = paramsCurr(:, 4);
    y = y(~isnan(x));
    x = x(~isnan(x));
    [R, P] = corrcoef(x, y);
    
    title(['p=' num2str(P(2,1)) 'Coeff=' num2str(R(2,1)), 'n=' num2str(length(x))]);
    ylabel('beta')
    
    sgtitle(allAnis{i});
    
    screen = get(0,'Screensize');
    screen(4) = screen(4) - 100;
    set(gcf, 'Position', screen);
end
%% all animals togther
% tstats
negRPE = tstatsAni;
[~, maxInd] = max(abs(squeeze(tstatsPupilAni(11:end, 2, :))),[],1);
maxInd = maxInd + 10;
pupilT = NaN(size(allAnis));
for i = 1:length(allAnis)
    pupilT(i) = tstatsPupilAni(maxInd(i), 2, i);
end
[R, P] = corrcoef(pupilT, negRPE);
figure2;
scatter(pupilT, negRPE, 20, 'filled');
title(['p = ' num2str(P(1,2)) 'R=' num2str(R(1,2))]);
xlabel('pupil-outcome tstats')
ylabel('lick-rpe tstats')
set(gca, 'tickdir', 'out')
figure2;
subplot(2,1,1)
histogram(pupilT, 5)
title('pupilTstats')
subplot(2,1,2)
histogram(negRPE, 5);
title('lickTstats')
%% coeff
negRPE = coeffsAni;
[~, maxInd] = max(abs(squeeze(coeffsPupilAni(11:end, 2, :))),[],1);
maxInd = maxInd + 10;
pupilT = NaN(size(allAnis));
for i = 1:length(allAnis)
    pupilT(i) = coeffsPupilAni(maxInd(i), 2, i);
end
figure2;
scatter(pupilT, negRPE, 20, 'filled');
[R, P] = corrcoef(pupilT, negRPE);
title(['p = ' num2str(P(1,2)) 'R=' num2str(R(1,2))]);
xlabel('pupil-outcome Coeff')
ylabel('lick-rpe Coeff')
set(gca, 'tickdir', 'out')
%%