function fpGLMmix_population(xlFile, sheet, category, region, varargin)
%task and model parameters
p = inputParser;
% default parameters if none given
p.addParameter('cellName', ['all']);
p.addParameter('plotFlag', 1);
p.addParameter('maxTrial', 1000);
% p.addParameter('modelName','7params_absPePeAN_scale_int_bias_ord')
p.addParameter('modelName','5params')
p.addParameter('regressors', '1+outcome+Qchosen')
p.addParameter('binSize', 1500)% in ms
p.addParameter('stepSize', 500)
p.addParameter('saveFigFlag', 1);
p.parse(varargin{:});
dayList = getDayList(xlFile, sheet, category);
tb = 2;
tf = 5;
numBins = 6;
focusBin = [0 1500]; % from reward delivery
populationSig = []; % the matrix with 1 for positive beta, -1 for negative beta
populationTStats = []; % t statistics for each parameter
populationCoeffs =[]; % coeffs for each regressor
paramNames = getParamNames_dF(p.Results.modelName, 1);
maxTrial = p.Results.maxTrial;
% basic info
[root, sep] = currComputer();
% time window
midPointsGLM = (-tb*1000 + 0.5 * p.Results.binSize):p.Results.stepSize:(1000*tf - 0.5 * p.Results.binSize);
slideTime = midPointsGLM;
allPe = {};
allChoices = {};
allSignal = {};
% all regressors
outcomeCombined = [];
peCombined = [];
prePeCombined = [];
preRwdCombined = []; 
QsumCombined = []; 
QdiffCombined = []; 
choiceConfCombined = []; 
biasSideCombined = []; 
rightSideCombined = []; 
lickLatCombined = []; 
hmmCombined = []; 
QchosenCombined = []; 
QunchosenCombined = []; 
QchosenUpdateCombined = []; 
preITICombined = []; 
dawExpCombined = []; 
svsCombined = []; 
svsNextCombined = []; 
svsWhenNrwdCombined = []; 
conNrwdsCombined = [];
allTrial_Matx_slideCombined = [];
group = [];
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
    fpDataPath = [pd.sortedFolder session '_photometry.mat'];
    % load behavior and neurons
    if exist(fpDataPath,'file')
        load(fpDataPath)
    else
        fprintf([session 'no fp file' '\n']);;
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
        [t,~,noSession] = getStanModelParams_samps(p.Results.modelName, [path sampFile '.mat'], 2000, 'sessionName', session);
        if noSession
            fprintf(['no good behavior in ' session '\n']);
            continue
        end
        
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
        hmm = double(os.hmmStates==1)';
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
       

    % slide window
    allTrial_Matx_slide = zeros(length(os.responseInds), length(midPointsGLM));
    if strcmp(region, 'LC')
        signalMat = LCmatChoice;
    else
        if strcmp(region, 'mPFC')
            signalMat = mPFCmatChoice;
        else 
            if strcmp(region, 'LCN')
                signalMat = LCNmatChoice;
            end
        end
    end
    for w = 1:length(midPointsGLM)
        timeInd = (midPoints >= midPointsGLM(w)-0.5*p.Results.binSize) & (midPoints <= (midPointsGLM(w)+0.5*p.Results.binSize -1));
        allTrial_Matx_slide(:,w) = nanmean(signalMat(:,timeInd),2);
    end
    
    outcomeCombined = [outcomeCombined; outcome];
    peCombined = [peCombined; pe];
    prePeCombined = [prePeCombined; prePe];
    preRwdCombined = [preRwdCombined; preRwd]; 
    QsumCombined = [QsumCombined; zscore(Qsum)]; 
    QdiffCombined = [QdiffCombined; zscore(Qdiff)]; 
    choiceConfCombined = [choiceConfCombined; zscore(choiceConf)]; 
    biasSideCombined = [biasSideCombined; biasSide]; 
    rightSideCombined = [rightSideCombined; rightSide]; 
    lickLatCombined = [lickLatCombined; lickLat]; 
    hmmCombined = [hmmCombined; hmm]; 
    QchosenCombined = [QchosenCombined; Qchosen]; 
    QunchosenCombined = [QunchosenCombined; zscore(Qunchosen)]; 
    QchosenUpdateCombined = [QchosenUpdateCombined; QchosenUpdate]; 
    preITICombined = [preITICombined; preITI]; 
    svsCombined = [svsCombined; svs]; 
    svsNextCombined = [svsNextCombined; svsNext]; 
    svsWhenNrwdCombined = [svsWhenNrwdCombined; svsWhenNrwd]; 
    conNrwdsCombined = [conNrwdsCombined; conNrwds]; 
    group = [group; repmat(session, length(outcome), 1)];
    allTrial_Matx_slideCombined = [allTrial_Matx_slideCombined; allTrial_Matx_slide];
    
    % bin signal by PE
    focusInd = (midPoints >= focusBin(1) + os.rwdDelay & midPoints <= focusBin(2) + os.rwdDelay);
    signalFocus = zscore(nanmean(signalMat(:,focusInd),2));
    edges = [linspace(-1-0.001, 0, 0.5*numBins+1) linspace(0, 1+0.001, 0.5*numBins+1)];
    edges = [edges(1:0.5*numBins), edges(0.5*numBins+2:end)];
    for k = 1:numBins
        meanPEs(sess, k) = mean(pe(pe >= edges(k) & pe < edges(k+1)));
        meanSignal(sess, k) = mean(signalFocus(pe >= edges(k) & pe < edges(k+1)));
    end
    
    psthTemp = zeros(numBins, length(midPoints));
    for k = 1:numBins
        tmpInd = (pe >= edges(k) & pe < edges(k+1));
        psthTemp(k,:) = mean(signalMat(tmpInd,:));
    end   
    psthMat = cat(3, psthMat, psthTemp);
    
end
    % change name for convenience
    outcome = outcomeCombined;
    pe = peCombined;
    prePe = prePeCombined;
    preRwd = preRwdCombined; 
    Qsum = QsumCombined; 
    Qdiff = QdiffCombined;
    choiceConf = choiceConfCombined;
    biasSide = biasSideCombined;
    rightSide = rightSideCombined; 
    lickLat = lickLatCombined; 
    hmm = hmmCombined;
    Qchosen = QchosenCombined;
    Qunchosen = QunchosenCombined;
    QchosenUpdate = QchosenUpdateCombined;
    preITI = preITICombined;
    svs = svsCombined; 
    svsNext = svsNextCombined;
    svsWhenNrwd = svsWhenNrwdCombined;
    conNrwds = conNrwdsCombined;
    
    tbl = table(group, outcome, pe, prePe, preRwd, Qsum, Qdiff, choiceConf, biasSide, rightSide, lickLat, hmm, Qchosen, Qunchosen, QchosenUpdate, preITI, svs, svsNext, svsWhenNrwd, conNrwds);

    sigs = [];
    tStats = [];
    coeffs = [];
    uppers = [];
    lowers = [];
    for k = 1:length(midPointsGLM)
        currSignal = allTrial_Matx_slideCombined(:,k);
        currTbl = addvars(tbl, currSignal);
        lm = fitlme(currTbl, ['currSignal~' p.Results.regressors]);
        sigs(k,:) = (lm.Coefficients.pValue(2:end)<0.05)';
        tStats(k,:) = lm.Coefficients.tStat(2:end)';
        coeffs(k,:) = lm.Coefficients.Estimate(2:end)';
        lowers(k,:) = lm.Coefficients.Lower(2:end)';
        uppers(k,:) = lm.Coefficients.Upper(2:end)';        
    end
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












