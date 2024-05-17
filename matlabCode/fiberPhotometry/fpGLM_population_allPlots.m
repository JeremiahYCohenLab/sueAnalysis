function fpGLM_population_allPlots(xlFile, sheet, category, region, varargin)
%task and model parameters
p = inputParser;
% default parameters if none given
p.addParameter('cellName', ['all']);
p.addParameter('plotFlag', 1);
p.addParameter('maxTrial', 1000);
p.addParameter('modelName','5params')
p.addParameter('regressors', '1+outcome+Qchosen+rightSide')
p.addParameter('binSize', 1500)% in ms
p.addParameter('stepSize', 500)
p.addParameter('focusWin', [0 1500])% in ms, from time of reward
p.addParameter('saveFigFlag', 1);
p.parse(varargin{:});
dayList = getDayList(xlFile, sheet, category);
tb = 2;
tf = 5;
numBins = 6;
focusBin = [0 1500]; % from reward delivery
% focusBin = [0 2000]; % from reward delivery
paramNames = getParamNames_dF(p.Results.modelName, 1);
maxTrial = p.Results.maxTrial;
% basic info
[root, sep] = currComputer();

% all regressors
sigs = [];
tStats = [];
coeffs = [];
uppers = [];
lowers = [];

meanPEs = NaN(length(dayList), numBins);
meanSignal = NaN(length(dayList), numBins);
psthMat = [];
%% session loop 
for sess = 1:length(dayList)
    % load model fitting results
    session = dayList{sess};
    pd = parseSessionString_df(session, root, sep);
    sampFile = [pd.animalName 'good' '_', p.Results.modelName];
    path = [root pd.animalName sep pd.animalName 'sorted' sep 'stan' sep 'bernoulli' sep p.Results.modelName sep 'good' sep];
    fprintf([session '\n']);
    % paths
    fpDataPath = [pd.sortedFolder session '_photometryCombinewithKH.mat'];
    % load behavior and neurons
    if exist(fpDataPath,'file')
        load(fpDataPath)
    else
        fprintf([session 'no fp file' '\n']);
        continue
    end
        midPoints = 1000*midPoints;

        %% behavior preparation 
        % parse behavior
        os = behAnalysisNoPlot_opMD(session, 'simpleFlag', 1);
        choice = os.allChoices';
        choice(choice<0) = 0;
        outcome = abs(os.allRewards)';
        choice = choice(1:min(length(choice), maxTrial));
        outcome = outcome(1:length(choice));
        outcomeL = outcome;
        outcomeL(outcomeL==0) = -1;
        outcomeR = outcome;
        outcomeR(outcomeR==0) = -1;
        outcomeL(choice==1) = 0;
        outcomeR(choice==0) = 0;
        responseInds = os.responseInds(1:min(length(choice), maxTrial)); 
        preRwd = [NaN abs(os.allRewards(1:end-1))]';
        %% behavior
        % switch
        svs = zeros(length(os.responseInds),1);
        svs(os.changeChoice_Inds) = 1;
        svsNext = [svs(2:end); NaN];
        svsWhenNrwd = svsNext;
        svsWhenNrwd(os.rwd_Inds) = NaN;
        params = getStanModelParams_sampsOnly(pd.animalName, 'good', p.Results.modelName, 2000, 'sessionName', session);
        t = inferModelVar(session, params, p.Results.modelName);
        
        % diff value
        Qdiff = abs(t.Q(:,2)-t.Q(:,1));
        % total value
        Qsum = sum(t.Q,2);
        % prepe
        prePe = [NaN; t.pe(1:end-1)];
        % pe
        pe = t.pe;
        % dawExp
        dawExp = double(t.probChoice <= 0.5);
        % confidence
        choiceConf = 2.*t.probChoice - 1;
        % chosen valie
        Qchosen  = zeros(length(choice),1);
        Qunchosen  = zeros(length(choice),1);
        QchosenUpdate = NaN(length(choice),1);
        for j = 1:length(choice)
            if j < length(choice)
                if choice(j)>0
                    Qchosen(j) = t.Q(j,2);
                    Qunchosen(j) = t.Q(j,1);
                    QchosenUpdate(j) = t.Q(j+2);
                else
                    Qchosen(j) = t.Q(j,1);
                    Qunchosen(j) = t.Q(j,2);
                    QchosenUpdate(j) = t.Q(j+1);
                end
            else                
                if choice(j)>0
                    Qchosen(j) = t.Q(j,2);
                    Qunchosen(j) = t.Q(j,1);
                else
                    Qchosen(j) = t.Q(j,1);
                    Qunchosen(j) = t.Q(j,2);
                end
                
            end
        end
        QchosenUpdate(~isnan(QchosenUpdate)) = zscore(QchosenUpdate(~isnan(QchosenUpdate)));
        % bias side
        biasSide = zeros(size(responseInds))';
        biasInd = contains(paramNames, 'bias');
        if mean(t.params(:,biasInd))>0
            biasSide(os.lickR_Inds)=1;
        else
            biasSide(os.lickL_Inds)=1;
        end
        % hmm = double(os.hmmStates==1)';
        lickLat = os.lickLatLogZ';
        rightSide = zeros(size(pe));
        rightSide(os.allChoices>0)=1;
        rightSide(os.allChoices<=0)=-1;
        preITI = os.timeBtwn';
        preITI(~isnan(preITI)) = zscore(preITI(~isnan(preITI)));
        % consecutive no rewards
        conNrwds = zeros(size(pe));
        for j = 1:length(choice)
            if outcome(j) == 0
                k = 1;
                while j-k>0 
                    if outcome(j-k)==0
                        k = k+1;
                    else
                        break
                    end
                end                
                conNrwds(j) = k;
            end
        end
       

    % get signal
    if strcmp(region, 'LC')
        signalMat = LCmatG(responseInds, :);
    else
        if strcmp(region, 'mPFC')
            signalMat = mPFCmatChoiceG(:, :);
        else 
            if strcmp(region, 'LCN')
                signalMat = LCNmatChoiceG(:, :);
            end
        end
    end

    
    signalMat = zscore(signalMat, [], 'all');
    bl = mean(signalMat(:,midPoints<-500 & midPoints>-1500), 2);
    % signalMat = signalMat - bl;
    % 
    % fprintf([session ' removed baseline \n'])
    

    % bin signal by PE
    focusInd = (midPoints >= focusBin(1) + os.rwdDelay & midPoints <= focusBin(2) + os.rwdDelay);
    signalFocus = zscore(mean(signalMat(:,focusInd),2));
    edges = [linspace(-1-0.001, 0, 0.5*numBins+1) linspace(0, 1+0.001, 0.5*numBins+1)];
    edges = [edges(1:0.5*numBins), edges(0.5*numBins+2:end)];
    for k = 1:numBins
        meanPEs(sess, k) = mean(pe(pe >= edges(k) & pe < edges(k+1)));
        meanSignal(sess, k) = mean(signalFocus(pe >= edges(k) & pe < edges(k+1)));
    end
    
    psthTemp = zeros(numBins, length(midPoints));
    target = t.pe;
    edges = quantile(target, linspace(0, 1, numBins+1));
    edges = linspace(min(target)-0.001, max(target)+0.001, numBins+1);
    edges = [linspace(-1-0.001, 0, 0.5*numBins+1) linspace(0, 1+0.001, 0.5*numBins+1)];
    edges = [edges(1:0.5*numBins), edges(0.5*numBins+2:end)];
    for k = 1:numBins
        tmpInd = (target >= edges(k) & target < edges(k+1));
        psthTemp(k,:) = mean(signalMat(tmpInd,:), 'omitmissing');
    end   
    psthMat = cat(3, psthMat, psthTemp);
    
    % linear model
    
    tbl = table(outcome, Qchosen, Qunchosen, Qsum, Qdiff, svs, rightSide);
    currSignal = signalFocus;
    currTbl = addvars(tbl, currSignal);
    lm = fitlm(currTbl, ['currSignal~' p.Results.regressors]);
    sigs(sess,:) = (lm.Coefficients.pValue(2:end)<0.05)';
    tStats(sess,:) = lm.Coefficients.tStat(2:end)';
    coeffs(sess,:) = lm.Coefficients.Estimate(2:end)';
    
end
%% plot everything
regressors = lm.CoefficientNames(2:end);
% tuning
figure2;
subplot(2,2,1)
plotFilled(mean(meanPEs), meanSignal, 'r');
xlabel('RPE')
ylabel('zscored signal')
title([sheet region ' photometry signal'])

colorsPSTH = zeros(numBins, 3);
for i = 1:0.5*numBins
    colorsPSTH(i, :) = [1, 0.3*(i-1), 0.3*(i-1)]; 
    colorsPSTH(numBins-(i-1), :) = [0.3*(i-1) , 0.3*(i-1), 1];
end
% psth

subplot(2,2,2)
hold on
for i = 1:numBins
    temp = squeeze(psthMat(i,:,:))';
    plotFilled(midPoints, temp, colorsPSTH(i,:));
end

xlim([-1000 4000])
set(gca, 'Box', 'off' )
set(gca, 'TickDir', 'Out')
xlabel('time from respond','FontSize', 15)
ylabel('zscored (dF/F - baseline)', 'FontSize', 15)
title([sheet ' ' region], 'FontSize', 15 )

% polar histogram
subplot(2,2,3); hold on;
tmpCoeffX = squeeze(coeffs(:,1));
tmpCoeffY = squeeze(coeffs(:,2));
tmpTStatsX = squeeze(tStats(:,1));
tmpTStatsY = squeeze(tStats(:,2));
scatter(tmpTStatsX, tmpTStatsY, 'm', 'filled');
plot([0 0], minmax(tmpTStatsY'), "Color", 'k', 'LineStyle', '--', 'LineWidth', 2);
plot(minmax(tmpTStatsX'), [0 0], "Color", 'k', 'LineStyle', '--', 'LineWidth', 2);
xlabel(regressors{1});
ylabel(regressors{2});


allVec = [tmpCoeffX, tmpCoeffY];
% allVec = [coeffsMax(:, outcomeInd), coeffsMax(:,qInd)];
[theta, rho] = cart2pol(allVec(:,1), allVec(:,2));
subplot(2,2,4) 
edges = linspace(-pi, pi, 20);
polarhistogram(theta,edges, 'FaceColor', [0.1 0.1 0.1], 'FaceAlpha',.7, 'EdgeColor', 'none', 'Normalization', 'Probability');
sgtitle(region);
%%


% for i = 1:length(midPoints)
%     titleStr = sprintf('From %d To %d', edges(i), edges(i+1)); 
%     figure;
%     scatterAll(squeeze(populationCoeffs(i,:,:))', regressors,7,'m');
%     suptitle(titleStr)
% end  












