%task and model parameters
clear all
xlFile = 'photometry'; 
sheet = 'all'; 
category = 'goodFP';
region = 'mPFC';
regressors = 'outcome + Qchosen + (Qchosen-1|group)';
modelName = '5params';
binSize = 1500; % in ms
stepSize = 500;
maxTrial = 1000;

dayList = getDayList(xlFile, sheet, category);
tb = 2;
tf = 5;
numBins = 6;
focusBin = [0 1500]; % from reward delivery
populationSig = []; % the matrix with 1 for positive beta, -1 for negative beta
populationTStats = []; % t statistics for each parameter
populationCoeffs =[]; % coeffs for each regressor
paramNames = getParamNames_dF(modelName, 1);
% basic info
[root, sep] = currComputer();
% time window
midPointsGLM = (-tb*1000 + 0.5 * binSize):stepSize:(1000*tf - 0.5 * binSize);
slideTime = midPointsGLM;
allPe = {};
allOutcome = {};
allRightSide = {};
allQchosen = {};
allSlides = {};
allSignalFocus = {};
focusCombined = [];
% all regressors
outcomeCombined = [];
QchosenCombined = []; 
peCombined = [];
rightSideCombined = [];
allTrial_Matx_slideCombined = [];
group = [];
%% session loop 
for sess = 1:length(dayList)
    % load model fitting results
    session = dayList{sess};
    pd = parseSessionString_df(session, root, sep);
    sampFile = [pd.animalName 'good' '_', modelName];
    path = [root pd.animalName sep pd.animalName 'sorted' sep 'stan' sep 'bernoulli' sep modelName sep 'good' sep];
    fprintf([session '\n']);
    % paths
    fpDataPath = [pd.sortedFolder session '_photometry.mat'];
    % load behavior and neurons
    if exist(fpDataPath,'file')
        load(fpDataPath)
    else
        fprintf([session 'no fp file' '\n']);
    end
        midPoints = 1000*midPoints;

        %% behavior preparation 
        % parse behavior
        os = behAnalysisNoPlot_opMD(session, 'simpleFlag', 1);
        choice = os.allChoices;
        choice(choice==-1) = 0;
        outcome = abs(os.allRewards)';
        responseInds = os.responseInds(1:min(length(choice), maxTrial)); 
        %% behavior
        [t,~,noSession] = getStanModelParams_samps(modelName, [path sampFile '.mat'], 2000, 'sessionName', session, 'sessionParamsFlag', 1);
        if noSession
            fprintf(['no good behavior in ' session '\n']);
            continue
        end
        pe = t.pe;
        % chosen valie
        Qchosen  = zeros(length(choice),1);
        for j = 1:length(choice)
            if choice(j)>0
                Qchosen(j) = t.Q(j,2);
            else
                Qchosen(j) = t.Q(j,1);
            end
        end
        rightSide = zeros(size(pe));
        rightSide(os.allChoices>0)=1;
        rightSide(os.allChoices<=0)=-1;
        
    % slide window
    allTrial_Matx_slide = zeros(length(os.responseInds), length(midPointsGLM));
    if strcmp(region, 'LC')
        signalMat = LCmatChoice;
    else
        if strcmp(region, 'mPFC')
            signalMat = mPFCmatChoice;
        end
    end
    for w = 1:length(midPointsGLM)
        timeInd = (midPoints >= midPointsGLM(w)-0.5*binSize) & (midPoints <= (midPointsGLM(w)+0.5*binSize -1));
        allTrial_Matx_slide(:,w) = nanmean(signalMat(:,timeInd),2);
    end
    
    % bin signal by PE
    focusInd = (midPoints >= focusBin(1) + os.rwdDelay & midPoints <= focusBin(2) + os.rwdDelay);
    signalFocus = nanmean(signalMat(:,focusInd),2);
    
    outcomeCombined = [outcomeCombined; outcome];
    peCombined = [peCombined; pe];
    QchosenCombined = [QchosenCombined; Qchosen];
    rightSideCombined = [rightSideCombined; rightSide];
    group = [group; repmat(session, length(outcome), 1)];
    allTrial_Matx_slideCombined = [allTrial_Matx_slideCombined; allTrial_Matx_slide];
    focusCombined = [focusCombined; signalFocus];
    
    allOutcome{sess} = outcome;
    allQchosen{sess} = Qchosen;
    allPe{sess} = pe;
    allRightSide{sess} = rightSide;
    allSlides{sess} = allTrial_Matx_slide;
    allSignalFocus{sess} = signalFocus;
    
end
% %% regression through time
%     outcome = outcomeCombined;
%     Qchosen = QchosenCombined;
%     
%     tbl = table(group, outcome, Qchosen);
% 
%     sigs = [];
%     tStats = [];
%     coeffs = [];
%     uppers = [];
%     lowers = [];
%     for k = 1:length(midPointsGLM)
%         currSignal = allTrial_Matx_slideCombined(:,k);
%         currTbl = addvars(tbl, currSignal);
%         lm = fitlme(currTbl, ['currSignal~' p.Results.regressors]);
%         sigs(k,:) = (lm.Coefficients.pValue(2:end)<0.05)';
%         tStats(k,:) = lm.Coefficients.tStat(2:end)';
%         coeffs(k,:) = lm.Coefficients.Estimate(2:end)';
%         lowers(k,:) = lm.Coefficients.Lower(2:end)';
%         uppers(k,:) = lm.Coefficients.Upper(2:end)';        
%     end
%% mix model on focus window
outcome = outcomeCombined;
Qchosen = QchosenCombined;

tbl = table(group, outcome, Qchosen);

currTbl = addvars(tbl, focusCombined);
lm = fitlme(currTbl, 'focusCombined ~ outcome + Qchosen + (outcome-1|group) + (Qchosen-1|group) + (1|group)');
[B,Bnames,stats] = randomEffects(lm);
outcomeInds = strcmp(Bnames.Name, 'outcome');
qChosenInds = strcmp(Bnames.Name, 'Qchosen');
outcomeRand = B(outcomeInds);
qChosenRand = B(qChosenInds);
allCoeff = lm.Coefficients.Estimate(2:end);
outcomeCoeff = allCoeff(1) + outcomeRand;
QchosenCoeff = allCoeff(2) + qChosenRand;
figure2;
scatter(outcomeCoeff, QchosenCoeff)
%% fit all models separately
allT = zeros(length(allOutcome),2);
allSig = zeros(length(allOutcome),2);
allCoeff = zeros(length(allOutcome),2);
for i = 1:length(allOutcome)
    mat = [allOutcome{i}, allQchosen{i}];
    lm = fitlm(mat, allSignalFocus{i});
    allT(i,:) = lm.Coefficients.tStat(2:end);
    allSig(i,:) = lm.Coefficients.pValue(2:end)<0.05;
    allCoeff(i,:) = lm.Coefficients.Estimate(2:end);
end
%%
figure2;
hold on; 
scatter(allT(:,1), allT(:,2))
plot([-2 8], [0 0], 'LineWidth', 1.5, 'LineStyle', '--', 'Color', [0.6 0.6 0.6]);
plot([0 0], [-4 1], 'LineWidth', 1.5, 'LineStyle', '--', 'Color', [0.6 0.6 0.6]);
title('tstats')
%%
figure2;
hold on; 
scatter(allCoeff(:,1), allCoeff(:,2))
plot([-0.4 1], [0 0], 'LineWidth', 1.5, 'LineStyle', '--', 'Color', [0.6 0.6 0.6]);
plot([0 0], [-0.6 0.2], 'LineWidth', 1.5, 'LineStyle', '--', 'Color', [0.6 0.6 0.6]);
title('coeff')
%%
edges = linspace(-pi, pi, 20);
% max tStats
outcomeInd = 1;
qInd = 2;

allVec = allT;
[theta, rho] = cart2pol(allVec(:,1), allVec(:,2));
figure2; 
polarhistogram(theta,edges, 'FaceColor', [0.1 0.1 0.1], 'FaceAlpha',.7, 'EdgeColor', 'none', 'Normalization', 'Probability');

title('maxWin, all')
%% plot everything
regressors = lm.CoefficientNames(2:end);
tFig = figure;
screen = get(0,'Screensize');
screen(4) = screen(4) - 100;
set(tFig, 'Position', screen)
colors = cool(length(regressors));
subplot(2,1,1); hold on;

for i = 1:length(regressors)
    plot(midPointsGLM,coeffs(:,i),'Color', colors(i,:), 'LineStyle','-', 'Marker','none', 'linewidth', 2);
end
for i = 1:length(regressors)
    fill([midPointsGLM, flip(midPointsGLM)], [lowers(:,i); flip(uppers(:,i))], colors(i,:), 'FaceAlpha', 0.2, 'EdgeColor', 'none');
    tmpCoeff = coeffs(:,i);
    tmpCoeff(sigs(:,i)==0) = NaN;
    scatter(midPointsGLM, tmpCoeff, 20, 'r');    
end
edges = [midPointsGLM - 0.5*p.Results.stepSize midPointsGLM(end)+0.5*p.Results.stepSize];
for i = 1:length(edges)
    line([edges(i) edges(i)], [0 1.2*max(uppers,[],'all')], 'color', [0.7 0.7 0.7], 'LineStyle','--')
end
line([os.rwdDelay os.rwdDelay], [0 1.2*max(uppers,[],'all')], 'color', 'r', 'LineStyle','--');
line([midPointsGLM(1), midPointsGLM(end)], [0 0], 'color', [0.7 0.7 0.7], 'LineStyle','--');
legend(regressors)
ylim([1.2*min(lowers,[],'all') 1.2*max(uppers,[],'all')])
xlim(minmax(edges));
title('coeffs')


subplot(2,1,2); hold on;

for i = 1:length(regressors)
    plot(midPointsGLM,tStats(:,i),'Color', colors(i,:), 'LineStyle','-', 'Marker','none', 'linewidth', 2);
end
for i = 1:length(regressors)
    tmptStats = tStats(:,i);
    tmptStats(sigs(:,i)==0) = NaN;
    scatter(midPointsGLM, tmptStats, 20, 'r');    
end
edges = [midPointsGLM - 0.5*p.Results.stepSize midPointsGLM(end)+0.5*p.Results.stepSize];
for i = 1:length(edges)
    line([edges(i) edges(i)], [0 1.2*max(tStats,[],'all')], 'color', [0.7 0.7 0.7], 'LineStyle','--')
end
line([os.rwdDelay os.rwdDelay], [0 1.2*max(tStats,[],'all')], 'color', 'r', 'LineStyle','--');
line([midPointsGLM(1), midPointsGLM(end)], [0 0], 'color', [0.7 0.7 0.7], 'LineStyle','--');
legend(regressors)
ylim([1.2*min(tStats,[],'all') 1.2*max(tStats,[],'all')])
xlim(minmax(edges));
xlabel('time from respond')
title('tStats')

sgtitle([sheet]);

figure2;
plotFilled(mean(meanPEs), meanSignal, 'r');
xlabel('RPE')
ylabel('zscored signal')
title([sheet region ' photometry signal'])

colorsPSTH = zeros(numBins, 3);
for i = 1:0.5*numBins
    colorsPSTH(i, :) = [0.3*(i-1) , 0.3*(i-1), 1];
    colorsPSTH(numBins-(i-1), :) = [1, 0.3*(i-1), 0.3*(i-1)];
end

figure2;
for i = 1:numBins
    temp = squeeze(psthMat(i,:,:))';
    plotFilled(midPoints, temp, colorsPSTH(i,:));
end

xlim([-1000 2500])
ylabel('dF/F')
set(gca, 'Box', 'off' )
set(gca, 'TickDir', 'Out')
xlabel('time from respond','FontSize', 15)
ylabel('dF/F', 'FontSize', 15)
title([sheet ' ' region], 'FontSize', 15 )

% for i = 1:length(midPoints)
%     titleStr = sprintf('From %d To %d', edges(i), edges(i+1)); 
%     figure;
%     scatterAll(squeeze(populationCoeffs(i,:,:))', regressors,7,'m');
%     suptitle(titleStr)
% end  












