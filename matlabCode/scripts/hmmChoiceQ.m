%% hidden markov model fitting
trans = [0.95,0.05;
      0.10,0.90];
emis = [1/6, 1/6, 1/6, 1/6, 1/6, 1/6;
   1/10, 1/10, 1/10, 1/10, 1/10, 1/2];

seq1 = hmmgenerate(100,trans,emis);
seq2 = hmmgenerate(200,trans,emis);
seqs = {seq1,seq2};
[estTR,estE] = hmmtrain(seqs,trans,emis);
state1 = hmmviterbi(seq1,estTR,estE);
state2 = hmmviterbi(seq2,estTR,estE);
%% focus on recording days
% prior = 10, 10
xlFile = {'ZS059', 'ZS060', 'ZS061', 'ZS062'};
sheet = {'ZS059', 'ZS060', 'ZS061', 'ZS062'};
col = 'good';
[root, sep] = currComputer();
modelName = '5params';
%% all animals
[root, sep] = currComputer();
[~, dayList, ~] = xlsread([root 'aniModel.xlsx'], 'all');
allAnis = dayList(2:end, 1);
allCols = dayList(2:end, 2);
allFile = dayList(2:end, 3);
allSheet = dayList(2:end, 4);

modelName = '5params';
xlFile = allFile;
sheet = allSheet;
cols = allCols;
%%
allEmis = {};
allQ = {};
allPe = {};
allStates = {};
allSession = {};
allRewards = {};
allSvSNext = {};
transGuess = [0.4 0.3 0.3;
              0.3 0.7 0;
              0.3 0 0.7];
emisGuess = [0.5 0.5;
             1 0;
             0 1];
for ani = 1:length(sheet)
    col = cols{ani};
    dayList = getDayList(xlFile{ani}, sheet{ani}, col);
    Emis = cell(length(dayList),1);
    Q = cell(length(dayList),1);
    Pe = cell(length(dayList),1);
    Rewards = cell(length(dayList),1);
    sN = cell(length(dayList),1);
    for i = 1:length(dayList)
        session = dayList{i};
        s = behAnalysisNoPlot_opMD(session, 'simpleFlag', 1);
        allChoices = 0.5*(s.allChoices + 3);
        Emis{i} = allChoices;
        svsNext = NaN(size(s.allChoices));
        svsNext(s.stayChoice_Inds-1) = 0;
        svsNext(s.changeChoice_Inds-1) = 1;
    
        pd = parseSessionString_df(session, root, sep);
        sampFile = [pd.animalName col '_' modelName];
        path = [root pd.animalName sep pd.animalName 'sorted' sep 'stan' sep 'bernoulli' sep modelName sep col sep];  
        [t,~,noSession] = getStanModelParams_samps(modelName, [path sampFile '.mat'], 4000, 'sessionName', session);

        Q{i} = t.Q;
        Pe{i} = t.pe;
        Rewards{i} = s.allRewards;
        sN{i} = svsNext;
        
%         eval([dayList{i} '.sessionData = s.behSessionData';])
    end
    % fitting
      
    % Estimation
    States = cell(length(dayList),1);
    [estTR,estE] = hmmtrainPrior(Emis,transGuess,emisGuess);
    for i = 1:length(dayList)   
%         [estTR,estE] = hmmtrainPrior(Emis{i},transGuess,emisGuess);
        States{i} = hmmviterbi(Emis{i}, estTR,estE);
%         eval([dayList{i} '.hmm = States{i}';])
    end
    
    allEmis = [allEmis; Emis];
    allQ = [allQ; Q];
    allPe = [allPe; Pe];
    allStates = [allStates; States];
    allSession = [allSession; dayList];
    allRewards = [allRewards; Rewards];
    allSvSNext = [allSvSNext; sN];
  
end
%% plotting and checking
for i = 1:length(allEmis)
    choices = allEmis{i};
    rewards = allRewards{i};
    states = allStates{i};
    % removing transition points
    transEInd = find(states(2:end-1)==1 & states(1:end-2)~=states(3:end) & states(1:end-2)~= 1 & states(3:end)~= 1)+1;
    states(transEInd) = NaN;
    figure;
    subplot(3,1,1);
    plot([1:length(choices); 1:length(choices)], [zeros(1,length(choices)); choices - 1.5 + 0.5*rewards], 'Color', 'k');
    hold on;
    EInd = find(states==1);
    patch([EInd-0.5;EInd+0.5;EInd+0.5;EInd-0.5], [-ones(1,length(EInd));-ones(1,length(EInd));ones(1,length(EInd));ones(1,length(EInd))], 'r', 'FaceAlpha', 0.4, 'EdgeColor', 'none');
    sgtitle(allSession{i})
    screen = get(0,'Screensize');
    screen(4) = screen(4) - 100;
    set(gcf, 'Position', screen);
    
end
%% plot p(R|Qdiff)
numBins = 6;
Qexploit = zeros(length(allEmis), numBins);
Qexplore = zeros(length(allEmis), numBins);
PRexploit = zeros(length(allEmis), numBins);
PRexplore = zeros(length(allEmis), numBins);
allQt = [];
allQe = [];
allChoicesT = [];
allChoicesE = [];
for i = 1:length(allSession)
    states = allStates{i};
    % remove those single trials in trainsition states
    transEInd = find(states(2:end-1)==1 & states(1:end-2)~=states(3:end) & states(1:end-2)~= 1 & states(3:end)~= 1)+1;
    states(transEInd) = NaN;
    if isempty(find(states==1,1))
        Qexploit(i, :) = NaN;
        Qexplore(i, :) = NaN;
        PRexploit(i, :) = NaN;
        PRexplore(i, :) = NaN;        
        continue
    end
    Q = allQ{i};
    choices = allEmis{i} - 1;
    Qdiff = Q(:,2) - Q(:,1);
    Qt = Qdiff(states==2 | states==3);
    Qe = Qdiff(states==1);
    choiceT = choices(states==2 | states==3);
    choiceE = choices(states==1);
    edgesT = linspace(min(Qt)-0.01, max(Qt)+0.01, numBins+1);
    edgesE = linspace(min(Qe)-0.01, max(Qe)+0.01, numBins+1);

    for j = 1:numBins
        
        Qexploit(i, j) = mean(Qt(Qt>=edgesT(j) & Qt<edgesT(j+1)), 'omitnan');
        Qexplore(i, j) = mean(Qe(Qe>=edgesE(j) & Qe<edgesE(j+1)), 'omitnan');
        PRexploit(i, j) = mean(choiceT(Qt>=edgesT(j) & Qt<edgesT(j+1)), 'omitnan');
        PRexplore(i, j) = mean(choiceE(Qe>=edgesE(j) & Qe<edgesE(j+1)), 'omitnan');
    end
    
    allQt = [allQt; Qt];
    allQe = [allQe; Qe];
    allChoicesT = [allChoicesT choiceT];
    allChoicesE = [allChoicesE choiceE];
end
%% plot mean of all sessions
semE = sem(PRexplore);
semT = sem(PRexploit);
meanCe = mean(PRexplore, 'omitnan');
meanCt = mean(PRexploit, 'omitnan');
meanQe = mean(Qexplore, 'omitnan');
meanQt = mean(Qexploit, 'omitnan');

figure2;
hold on;
plot(meanQe, meanCe, 'Color', 'm', 'LineWidth', 2);
patch([meanQe, flip(meanQe)], [meanCe-semE flip(meanCe+semE)], 'm', 'EdgeColor', 'none', 'FaceAlpha', 0.4);
plot(meanQt, meanCt, 'Color', 'b', 'LineWidth', 2);
patch([meanQt, flip(meanQt)], [meanCt-semT flip(meanCt+semT)], 'b', 'EdgeColor', 'none', 'FaceAlpha', 0.4);
%% plot p(R) after pulling all trials
edgesT = linspace(min(allQt)-0.01, max(allQt)+0.01, numBins+1);
edgesE = linspace(min(allQe)-0.01, max(allQe)+0.01, numBins+1);

meanQt = zeros(numBins,1);
meanQe = zeros(numBins,1);
meanCt = zeros(numBins,1);
meanCe = zeros(numBins,1);
semT = zeros(numBins,1);
semE = zeros(numBins,1);
for j = 1:numBins
    meanQt(j) = mean(allQt(allQt>=edgesT(j) & allQt<edgesT(j+1)), 'omitnan');
    meanQe(j) = mean(allQe(allQe>=edgesE(j) & allQe<edgesE(j+1)), 'omitnan');
    meanCt(j) = mean(allChoicesT(allQt>=edgesT(j) & allQt<edgesT(j+1)), 'omitnan');
    meanCe(j) = mean(allChoicesE(allQe>=edgesE(j) & allQe<edgesE(j+1)), 'omitnan');
    semT(j) = sem(allChoicesT(allQt>=edgesT(j) & allQt<edgesT(j+1)));
    semE(j) = sem(allChoicesE(allQe>=edgesE(j) & allQe<edgesE(j+1)));
end

figure2;
hold on;
plot(meanQe, meanCe, 'Color', 'm', 'LineWidth', 2);
patch([meanQe; flip(meanQe)]', [meanCe-semE; flip(meanCe+semE)]', 'm', 'EdgeColor', 'none', 'FaceAlpha', 0.4);
plot(meanQt, meanCt, 'Color', 'b', 'LineWidth', 2);
patch([meanQt; flip(meanQt)]', [meanCt-semT; flip(meanCt+semT)]', 'b', 'EdgeColor', 'none', 'FaceAlpha', 0.4);
%% plot rewards before start and before end of the exploration
trialS = 6;
trialE = 5;
rwdStart = [];
rwdEnd = [];
svsEnd = [];
lickStart = [];
lickEnd = [];
lickT = [];
lickE = [];
lickTStay = [];
lickEStay = [];
lickTSwitch = [];
lickESwitch = [];
for i = 1:length(allSession)
    states = allStates{i};
    % remove those single trials in trainsition states
    transEInd = find(states(2:end-1)==1 & states(1:end-2)~=states(3:end) & states(1:end-2)~= 1 & states(3:end)~= 1)+1;
    states(transEInd) = NaN;
    states(states==2|states==3) = 0;
    s = behAnalysisNoPlot_opMD(allSession{i}, 'simpleFlag', 1);
    outcomes = abs(s.allRewards);
    % before exploration start
    currInds = find(states(1:end-1) == 0 & states(2:end) == 1) + 1;
    % check animal did not change state while in trialS before exploration
    sumS = conv(states, ones(1, trialS));
    sumS = sumS(1:length(states));
    currInds = currInds(sumS(currInds)==1);
    % gather all starts together
    outcomeTemp = [NaN(1, trialS-1), outcomes];
    lickTemp = [NaN(1, trialS-1), s.lickLatZ];
    currInds = currInds + trialS - 1; % added trialS -1 to take care of appended trials
    if ~isempty(currInds)
        rwdMatTemp = NaN(length(currInds), trialS);
        lickMatTemp = NaN(length(currInds), trialS);
        for j = 1:trialS
            rwdMatTemp(:,j) = outcomeTemp(currInds - (trialS - j))';
            lickMatTemp(:,j) = lickTemp(currInds - (trialS - j))';
        end
        rwdStart = [rwdStart; rwdMatTemp];
        lickStart = [lickStart; lickMatTemp];
    end
    % before exploration end
    currInds = find(states(1:end-1) == 1 & states(2:end) == 0) + 1;
    % check animal did not change state while in trialS before exploration
    sumE = conv(states, ones(1, trialE));
    sumE = sumE(1:length(states));
    currInds = currInds(sumE(currInds)==trialE-1);
    % gather all starts together
    outcomeTemp = [NaN(1, trialE-1), outcomes];
    lickTemp = [NaN(1, trialE-1), s.lickLatZ];
    svsTemp = NaN(1, trialE-1+length(states));
    svsTemp(s.stayChoice_Inds+trialE-1) = 0;
    svsTemp(s.changeChoice_Inds+trialE-1) = 1;
    currInds = currInds + trialE - 1; % added trialS -1 to take care of appended trials
    if ~isempty(currInds)
        rwdMatTemp = NaN(length(currInds), trialE);
        lickMatTemp = NaN(length(currInds), trialE);
        svsMatTemp = NaN(length(currInds), trialE);
        for j = 1:trialE
            rwdMatTemp(:,j) = outcomeTemp(currInds - (trialE - j))';
            lickMatTemp(:,j) = lickTemp(currInds - (trialE - j))';
            svsMatTemp(:,j) = svsTemp(currInds - (trialE - j))';
        end
        rwdEnd = [rwdEnd; rwdMatTemp];
        lickEnd = [lickEnd; lickMatTemp];
        svsEnd = [svsEnd; svsMatTemp];
    end
    % get licks
    lickT(i) = mean(s.lickLat(states==0));
    lickE(i) = mean(s.lickLat(states==1));
    % stay licks
    lickTStay(i) = mean(s.lickLat(intersect(find(states==0), s.stayChoice_Inds)));
    lickEStay(i) = mean(s.lickLat(intersect(find(states==1), s.stayChoice_Inds)));
    % switch licks
    lickTSwitch(i) = mean(s.lickLat(intersect(find(states==0), s.changeChoice_Inds)));
    lickESwitch(i) = mean(s.lickLat(intersect(find(states==1), s.changeChoice_Inds)));
end
%% plot results
% rwd
figure2;
meanS = mean(rwdStart, 'omitnan');
semS = sem_bern(rwdStart);
errorbar([1:trialS]-trialS, meanS, semS, 'r', 'LineWidth', 2);
xlabel('from start of exploration')
ylabel('P(rwd)')
figure2;
meanE = mean(rwdEnd, 'omitnan');
semE = sem_bern(rwdEnd);
errorbar([1:trialE]-trialE, meanE, semE, 'b', 'LineWidth', 2);
xlabel('from start of exploitation')
ylabel('P(rwd)')
% svs
figure2;
meanE = mean(svsEnd, 'omitnan');
semE = sem_bern(svsEnd);
errorbar([1:trialE]-trialE, meanE, semE, 'b', 'LineWidth', 2);
xlabel('from start of exploitation')
ylabel('P(switch)')
% lick
figure2;
meanS = mean(lickStart, 'omitnan');
semS = sem(lickStart);
errorbar([1:trialS]-trialS, meanS, semS, 'r', 'LineWidth', 2);
xlabel('from start of exploration')
ylabel('lickLat (zscored)')
figure2;
meanE = mean(lickEnd, 'omitnan');
semE = sem(lickEnd);
errorbar([1:trialE]-trialE, meanE, semE, 'b', 'LineWidth', 2);
xlabel('from start of exploitation')
ylabel('lickLat (zscored)')
%
figure2; hold on;
scatter(lickT(lickT<=400), lickE(lickT<=400), 15, [0.5 0.5 0.5], 'MarkerEdgeAlpha', 0.6);
xlabel('exploitation');
ylabel('exploration')
plot([0 500], [0 500], 'LineWidth', 2, 'LineStyle', '--');
[h, p, ~, stats] = ttest(lickT(lickT<=400), lickE(lickT<=400));
title(['E vs T' 'p=' num2str(p) 'tstat=' num2str(stats.tstat), 'df=' num2str(stats.df)])
%
figure2; hold on;
scatter(lickTStay(lickT<=400), lickEStay(lickT<=400), 15, [0.5 0.5 0.5], 'MarkerEdgeAlpha', 0.6);
xlabel('exploitation');
ylabel('exploration')
plot([0 500], [0 500], 'LineWidth', 2, 'LineStyle', '--');
[h, p, ~, stats] = ttest(lickTStay(lickT<=400), lickEStay(lickT<=400));
title(['Stay: E vs T' 'p=' num2str(p) 'tstat=' num2str(stats.tstat), 'df=' num2str(stats.df)])

%
figure2; hold on;
scatter(lickTSwitch(lickT<=400), lickESwitch(lickT<=400), 15, [0.5 0.5 0.5], 'MarkerEdgeAlpha', 0.6);
xlabel('exploitation');
ylabel('exploration')
plot([0 500], [0 500], 'LineWidth', 2, 'LineStyle', '--');
[h, p, ~, stats] = ttest(lickTSwitch(lickT<=400), lickESwitch(lickT<=400));
title(['Switch: E vs T' 'p=' num2str(p) 'tstat=' num2str(stats.tstat), 'df=' num2str(stats.df)])
%% hmm-dependent laser effect;
numBins = 2;
QtMeansL = zeros(length(allSession), numBins);
QtMeansN = zeros(length(allSession), numBins);
QeMeansL = zeros(length(allSession), numBins);
QeMeansN = zeros(length(allSession), numBins);
QMeansL = zeros(length(allSession), numBins);
QMeansN = zeros(length(allSession), numBins);

pSwTL = zeros(length(allSession), numBins);
pSwTN = zeros(length(allSession), numBins);
pSwEL = zeros(length(allSession), numBins);
pSwEN = zeros(length(allSession), numBins);
pSwL = zeros(length(allSession), numBins);
pSwN = zeros(length(allSession), numBins);

combinedChoices = [];
combinedRewards = [];
combinedSvsNext = [];
combinedLaser = [];
combinedStates = [];
combinedQchosen = [];
combinedChoiceNext = [];
combinedQdiff = [];
combinedLickLat = [];
combinedLickLatPre = [];

for i = 1:length(allSession)
    states = allStates{i};
    session = allSession{i};
    % remove those single trials in trainsition states
    transEInd = find(states(2:end-1)==1 & states(1:end-2)~=states(3:end) & states(1:end-2)~= 1 & states(3:end)~= 1)+1;
    states(transEInd) = NaN;
    Q = allQ{i};
    s = behAnalysisNoPlot_opMD(session, 'simpleFlag', 1);
 
    Qdiff = Q(:,2) - Q(:,1);
    Qchosen = zeros(size(s.allChoices));
    choice = s.allChoices';
    choice(choice<0) = 0;
    for j = 1:length(choice)
        if choice(j)>0
            Qchosen(j) = Q(j,2);
        else
            Qchosen(j) = Q(j,1);
        end
    end
%     Qchosen(s.allChoices == 1) = Q(s.allChoices == 1,2);
%     Qchosen(s.allChoices == -1) = Q(s.allChoices == -1,1);
    Qchosen = zscore(Qchosen);
    Qt = Qchosen(states==2 | states==3);
    Qe = Qchosen(states == 1);
    
%     edgesT = linspace(min(Qt)-0.01, max(Qt)+0.01, numBins+1);
%     edgesE = linspace(min(Qe)-0.01, max(Qe)+0.01, numBins+1);
    edgesT = quantile(Qt, linspace(0, 1, numBins+1));
    edgesE = quantile(Qe, linspace(0, 1, numBins+1));
    edges = linspace(min(Qchosen)-0.01, max(Qchosen)+0.01, numBins+1);
    
    edgesT(end) = edgesT(end)+0.001;
    edgesE(end) = edgesE(end)+0.001;
    
    svsNext = NaN(size(s.allChoices));
    svsNext(s.stayChoice_Inds-1) = 0;
    svsNext(s.changeChoice_Inds-1) = 1;
    
    svsNextT = svsNext(states==2 | states==3);
    svsNextE = svsNext(states == 1);
    
    laser = s.laser;
    laserT = s.laser(states==2 | states==3);
    laserE = s.laser(states == 1);
    
    rwd = abs(s.allRewards);
    rwdT = abs(s.allRewards(states==2 | states == 3));
    rwdE = abs(s.allRewards(states==1));
    
    selectRwd = 0;
    
    for j = 1:numBins
        QtMeansL(i,j) = mean(Qt(Qt>=edgesT(j)&Qt<edgesT(j+1)&laserT==1&rwdT==selectRwd));
        QtMeansN(i,j) = mean(Qt(Qt>=edgesT(j)&Qt<edgesT(j+1)&laserT==0&rwdT==selectRwd));
        QeMeansL(i,j) = mean(Qe(Qe>=edgesE(j)&Qe<edgesE(j+1)&laserE==1&rwdE==selectRwd));
        QeMeansN(i,j) = mean(Qe(Qe>=edgesE(j)&Qe<edgesE(j+1)&laserE==0&rwdE==selectRwd));
        QMeansL(i,j) = mean(Qchosen(Qchosen>=edges(j)&Qchosen<edges(j+1)&laser==1&rwd==selectRwd));
        QMeansN(i,j) = mean(Qchosen(Qchosen>=edges(j)&Qchosen<edges(j+1)&laser==0&rwd==selectRwd));
        
        pSwTL(i,j) = mean(svsNextT(Qt>=edgesT(j)&Qt<edgesT(j+1)&laserT==1&rwdT==selectRwd), 'omitnan');
        pSwTN(i,j) = mean(svsNextT(Qt>=edgesT(j)&Qt<edgesT(j+1)&laserT==0&rwdT==selectRwd), 'omitnan');
        pSwEL(i,j) = mean(svsNextE(Qe>=edgesE(j)&Qe<edgesE(j+1)&laserE==1&rwdE==selectRwd), 'omitnan');
        pSwEN(i,j) = mean(svsNextE(Qe>=edgesE(j)&Qe<edgesE(j+1)&laserE==0&rwdE==selectRwd), 'omitnan');
        pSwL(i,j) = mean(svsNext(Qchosen>=edges(j)&Qchosen<edges(j+1)&laser==1&rwd==selectRwd), 'omitnan');
        pSwN(i,j) = mean(svsNext(Qchosen>=edges(j)&Qchosen<edges(j+1)&laser==0&rwd==selectRwd), 'omitnan');
    end
    choiceNext = [0.5*(s.allChoices(2:end)+1), NaN];
    states(states==2|states==3) = 0;
    combinedStates = [combinedStates; states'];
    combinedQchosen = [combinedQchosen; Qchosen'];
    combinedQdiff = [combinedQdiff; Qdiff];
    combinedRewards = [combinedRewards; abs(s.allRewards)'];
    combinedChoices = [combinedChoices; 0.5*(s.allChoices+1)'];
    combinedChoiceNext = [combinedChoiceNext; choiceNext'];
    combinedSvsNext = [combinedSvsNext; svsNext'];
    combinedLaser = [combinedLaser; s.laser'];
end
%% glm & lm test p(switch)
tbl = table(combinedStates, combinedQchosen, -combinedRewards+1, combinedChoices, combinedChoiceNext, combinedLaser, combinedSvsNext,...
    'VariableNames', {'states', 'Qchosen', 'rwd', 'choice', 'choiceNext', 'laser', 'svs'});
glm = fitglm(tbl, 'svs~1 + rwd*states + Qchosen + states', 'Distribution', 'binomial', 'Link', 'probit');
lm = fitlm(tbl, 'svs~1 + laser*rwd*states + Qchosen + states');
%% glm & lm test p(choice=R)
tbl = table(combinedStates, combinedQdiff, combinedRewards-0.5, combinedChoices, combinedChoiceNext, combinedLaser-0.5, combinedSvsNext,...
    'VariableNames', {'states', 'Qdiff', 'rwd', 'choice', 'choiceNext', 'laser', 'svs'});
glmChoice = fitglm(tbl, 'choiceNext~1 + choice + Qdiff + states*choice*rwd - rwd:states - rwd - states + states:laser:rwd:choice', 'Distribution', 'binomial', 'Link', 'probit');
%%
colorL = [0.2 0.2 1];
colorN = [0.5 0.5 0.5];
QtL = mean(QtMeansL, 'omitnan');
QtN = mean(QtMeansN, 'omitnan');
QeL = mean(QeMeansL, 'omitnan');
QeN = mean(QeMeansN, 'omitnan');
QL = mean(QMeansL, 'omitnan');
QN = mean(QMeansN, 'omitnan');
meanSwTL = mean(pSwTL, 'omitnan');
meanSwTN = mean(pSwTN, 'omitnan');
meanSwEL = mean(pSwEL, 'omitnan');
meanSwEN = mean(pSwEN, 'omitnan');
meanSwL = mean(pSwL, 'omitnan');
meanSwN = mean(pSwN, 'omitnan');
semSwTL = sem(pSwTL);
semSwTN = sem(pSwTN);
semSwEL = sem(pSwEL);
semSwEN = sem(pSwEN);
semSwL = sem(pSwL);
semSwN = sem(pSwN);

figure2;
subplot(1,3,1); hold on;
plot(QeL, meanSwEL, 'Color', colorL, 'LineWidth', 2);
patch([QeL, flip(QeL)], [meanSwEL-semSwEL flip(meanSwEL+semSwEL)], colorL, 'edgeColor', 'none', 'FaceAlpha', 0.4);
plot(QeN, meanSwEN, 'Color', colorN, 'LineWidth', 2);
patch([QeN, flip(QeN)], [meanSwEN-semSwEN flip(meanSwEN+semSwEN)], colorN, 'edgeColor', 'none', 'FaceAlpha', 0.4);

subplot(1,3,2); hold on;
plot(QtL, meanSwTL, 'Color', colorL, 'LineWidth', 2);
patch([QtL, flip(QtL)], [meanSwTL-semSwTL flip(meanSwTL+semSwTL)], colorL, 'edgeColor', 'none', 'FaceAlpha', 0.4);
plot(QtN, meanSwTN, 'Color', colorN, 'LineWidth', 2);
patch([QtN, flip(QtN)], [meanSwTN-semSwTN flip(meanSwTN+semSwTN)], colorN, 'edgeColor', 'none', 'FaceAlpha', 0.4);

subplot(1,3,3); hold on;
plot(QL, meanSwL, 'Color', colorL, 'LineWidth', 2);
patch([QL, flip(QL)], [meanSwL-semSwL flip(meanSwL+semSwL)], colorL, 'edgeColor', 'none', 'FaceAlpha', 0.4);
plot(QN, meanSwN, 'Color', colorN, 'LineWidth', 2);
patch([QN, flip(QN)], [meanSwN-semSwN flip(meanSwN+semSwN)], colorN, 'edgeColor', 'none', 'FaceAlpha', 0.4);
%% without binning 
% hmm-dependent laser effect;
% switch after no rwd
pSwNrwdTL = zeros(length(allSession), 1);
pSwNrwdTN = zeros(length(allSession), 1);
pSwNrwdEL = zeros(length(allSession), 1);
pSwNrwdEN = zeros(length(allSession), 1);
pSwNrwdL = zeros(length(allSession), 1);
pSwNrwdN = zeros(length(allSession), 1);
% switch after rwd
pSwRwdTL = zeros(length(allSession), 1);
pSwRwdTN = zeros(length(allSession), 1);
pSwRwdEL = zeros(length(allSession), 1);
pSwRwdEN = zeros(length(allSession), 1);
pSwRwdL = zeros(length(allSession), 1);
pSwRwdN = zeros(length(allSession), 1);

for i = 1:length(allSession)
    states = allStates{i};
    session = allSession{i};
    % remove those single trials in trainsition states
    transEInd = find(states(2:end-1)==1 & states(1:end-2)~=states(3:end) & states(1:end-2)~= 1 & states(3:end)~= 1)+1;
    states(transEInd) = NaN;
    Q = allQ{i};
    s = behAnalysisNoPlot_opMD(session, 'simpleFlag', 1);
 
    Qdiff = Q(:,2) - Q(:,1);
    Qchosen = zeros(size(s.allChoices));
    choice = s.allChoices';
    choice(choice<0) = 0;
    for j = 1:length(choice)
        if choice(j)>0
            Qchosen(j) = Q(j,2);
        else
            Qchosen(j) = Q(j,1);
        end
    end
%     Qchosen(s.allChoices == 1) = Q(s.allChoices == 1,2);
%     Qchosen(s.allChoices == -1) = Q(s.allChoices == -1,1);
    Qchosen = zscore(Qchosen);
    Qt = Qchosen(states==2 | states==3);
    Qe = Qchosen(states == 1);
    
    svsNext = NaN(size(s.allChoices));
    svsNext(s.stayChoice_Inds-1) = 0;
    svsNext(s.changeChoice_Inds-1) = 1;
    
    svsNextT = svsNext(states==2 | states==3);
    svsNextE = svsNext(states == 1);
    
    laser = s.laser;
    laserT = s.laser(states==2 | states==3);
    laserE = s.laser(states == 1);
    
    rwd = abs(s.allRewards);
    rwdT = abs(s.allRewards(states==2 | states == 3));
    rwdE = abs(s.allRewards(states==1));
    
    selectRwd = 0;
           
    pSwNrwdTL(i) = mean(svsNextT(laserT==1&rwdT==selectRwd), 'omitnan');
    pSwNrwdTN(i) = mean(svsNextT(laserT==0&rwdT==selectRwd), 'omitnan');
    pSwNrwdEL(i) = mean(svsNextE(laserE==1&rwdE==selectRwd), 'omitnan');
    pSwNrwdEN(i) = mean(svsNextE(laserE==0&rwdE==selectRwd), 'omitnan');
    pSwNrwdL(i) = mean(svsNext(laser==1&rwd==selectRwd), 'omitnan');
    pSwNrwdN(i) = mean(svsNext(laser==0&rwd==selectRwd), 'omitnan');

    selectRwd = 1;
           
    pSwRwdTL(i) = mean(svsNextT(laserT==1&rwdT==selectRwd), 'omitnan');
    pSwRwdTN(i) = mean(svsNextT(laserT==0&rwdT==selectRwd), 'omitnan');
    pSwRwdEL(i) = mean(svsNextE(laserE==1&rwdE==selectRwd), 'omitnan');
    pSwRwdEN(i) = mean(svsNextE(laserE==0&rwdE==selectRwd), 'omitnan');
    pSwRwdL(i) = mean(svsNext(laser==1&rwd==selectRwd), 'omitnan');
    pSwRwdN(i) = mean(svsNext(laser==0&rwd==selectRwd), 'omitnan');
end


%% 
meanpSwNrwdTL = mean(pSwNrwdTL, 'omitnan');
semSwNrwdTL = sem(pSwNrwdTL);
meanpSwNrwdTN = mean(pSwNrwdTN, 'omitnan');
semSwNrwdTN = sem(pSwNrwdTN);
meanpSwNrwdEL = mean(pSwNrwdEL, 'omitnan');
semSwNrwdEL = sem(pSwNrwdEL);
meanpSwNrwdEN = mean(pSwNrwdEN, 'omitnan');
semSwNrwdEN = sem(pSwNrwdEN);

meanpSwRwdTL = mean(pSwRwdTL, 'omitnan');
semSwRwdTL = sem(pSwRwdTL);
meanpSwRwdTN = mean(pSwRwdTN, 'omitnan');
semSwRwdTN = sem(pSwRwdTN);
meanpSwRwdEL = mean(pSwRwdEL, 'omitnan');
semSwRwdEL = sem(pSwRwdEL);
meanpSwRwdEN = mean(pSwRwdEN, 'omitnan');
semSwRwdEN = sem(pSwRwdEN);

figure2;
subplot(1,2,1);
hold on;
errorbar([meanpSwNrwdTN meanpSwNrwdTL], [semSwNrwdTN semSwNrwdTL], 'Color', colorN, 'LineWidth', 2);
errorbar([meanpSwNrwdEN meanpSwNrwdEL], [semSwNrwdEN semSwNrwdEL], 'Color', colorL, 'LineWidth', 2);
pT = signrank(pSwNrwdTL, pSwNrwdTN);
pE = signrank(pSwNrwdEL, pSwNrwdEN);
text(1.5, 0.3, num2str(pT));
text(1.5, 0.8, num2str(pE));
xlim([0.75 2.25])
ylim([-0.1 0.9])
title('P(switch|nrwd)', 'FontSize', 18)
xlabel('control    laser', 'FontSize', 18)
set(gca,'tickdir', 'out')
set(gca, 'XTick', [])
set(gca, 'YTick', [0:0.2:0.8], 'FontSize',14)


subplot(1,2,2);
hold on;
errorbar([meanpSwRwdTN meanpSwRwdTL], [semSwRwdTN semSwRwdTL], 'Color', colorN, 'LineWidth', 2);
errorbar([meanpSwRwdEN meanpSwRwdEL], [semSwRwdEN semSwRwdEL], 'Color', colorL, 'LineWidth', 2);
pT = signrank(pSwRwdTL, pSwRwdTN);
pE = signrank(pSwRwdEL, pSwRwdEN);
text(1.5, 0.35, num2str(pT));
text(1.5, 0.05, num2str(pE));
xlim([0.75 2.25])
ylim([-0.1 0.9])
title('P(switch|rwd)', 'FontSize', 18)
xlabel('control    laser', 'FontSize', 18)
set(gca,'tickdir', 'out')
set(gca, 'XTick', [])
set(gca, 'YTick', [0:0.2:0.8], 'FontSize',14)
legend({'exploit', 'explore'})
%% load data
load('F:\tmpData\hmmDataAllAni.mat')
%% glm rwd history on p(R)
combinedRwdMat = [];
combinedNrwdMat = [];
combinedStatesMat = [];
combinedLatMat = [];

combinedChoices = [];
combinedStates = [];
combinedSvS = [];
combinedLat = [];

numTrial = 7;

for i = 1:length(allSession)
    session = allSession{i};
    states = allStates{i};
    % remove those single trials in trainsition states
    transEInd = find(states(2:end-1)==1 & states(1:end-2)~=states(3:end) & states(1:end-2)~= 1 & states(3:end)~= 1)+1;
    states(transEInd) = NaN;
    states(states==2 | states==3) = 0;
    s = behAnalysisNoPlot_opMD(session, 'simpleFlag', 1);
    svs = NaN(1, length(s.allChoices));
    svs(s.changeChoice_Inds) = 1;
    svs(s.stayChoice_Inds) = 0;
    choicesTmp = s.allChoices;
    choicesTmp(choicesTmp==-1) = 0;
    combinedChoices = [combinedChoices choicesTmp];
    combinedStates = [combinedStates states];
    combinedLat = [combinedLat s.lickLatZ];
%     combinedLatPre = [combinedLatPre [NaN s.lickLatLogZ(1:end-1)]];
    combinedSvS = [combinedSvS svs];
    
    rwdMat = NaN(numTrial, length(choicesTmp));
    nRwdMat = NaN(numTrial, length(choicesTmp));
    stateMat = NaN(numTrial, length(choicesTmp));
    lickMat = NaN(numTrial, length(choicesTmp));
    for j = 1:numTrial
        rwdMat(j, j+1:end) = s.allRewards(1:end-j);
        nRwdMat(j, j+1:end) = s.allNoRewards(1:end-j);
        stateMat(j, j+1:end) = states(1:end-j);
        lickMat(j, j+1:end) = s.lickLatZ(1:end-j);
    end
    combinedRwdMat = [combinedRwdMat rwdMat];
    combinedNrwdMat = [combinedNrwdMat nRwdMat];
    combinedStatesMat = [combinedStatesMat stateMat];
    combinedLatMat = [combinedLatMat lickMat];
end
%%
glm = fitglm([combinedRwdMat(:,combinedStates==0)', combinedNrwdMat(:,combinedStates==0)'], combinedChoices(combinedStates==0)');

%%
glm = fitglm([combinedRwdMat', combinedStatesMat'.*combinedRwdMat'], combinedChoices');
%%
combinedRwdMatE = combinedRwdMat;
combinedRwdMatE(combinedStatesMat==0) = 0;
combinedRwdMatE(isnan(combinedStatesMat)) = NaN;
combinedRwdMatT = combinedRwdMat;
combinedRwdMatT(combinedStatesMat==1) = 0;
combinedRwdMatT(isnan(combinedStatesMat)) = NaN;

combinednRwdMatE = combinedNrwdMat;
combinednRwdMatE(combinedStatesMat==0) = 0;
combinednRwdMatE(isnan(combinedStatesMat)) = NaN;
combinednRwdMatT = combinedNrwdMat;
combinednRwdMatT(combinedStatesMat==1) = 0;
combinednRwdMatT(isnan(combinedStatesMat)) = NaN;

combinedChoicesMat = combinedRwdMat + combinednRwdMatT;
combinedChoicesMatE = combinedChoicesMat;
combinedChoicesMatE(combinedStatesMat==0) = 0;
combinedChoicesMatE(isnan(combinedStatesMat)) = NaN;
combinedChoicesMatT = combinedChoicesMat;
combinedChoicesMatT(combinedStatesMat==1) = 0;
combinedChoicesMatT(isnan(combinedStatesMat)) = NaN;


mat = [combinednRwdMatT' combinednRwdMatE' combinedChoicesMatT' combinedChoicesMatE'];

glm = fitglm(mat, combinedChoices');
%% for first numTrials
relevInds = 2:numTrial+1;
coefVals = glm.Coefficients.Estimate(relevInds);
CIbands = coefCI(glm);
errorL = abs(coefVals - CIbands(relevInds,1));
errorU = abs(coefVals - CIbands(relevInds,2));
figure;
errorbar((1:numTrial),coefVals,errorL,errorU,'Color', [0.7 0 1],'linewidth',2)
%% for two numTrials
relevInds = 2:numTrial+1;
coefVals = glm.Coefficients.Estimate(relevInds);
CIbands = coefCI(glm);
errorL = abs(coefVals - CIbands(relevInds,1));
errorU = abs(coefVals - CIbands(relevInds,2));
figure; hold on;
errorbar((1:numTrial),coefVals,errorL,errorU,'Color', [1 0.3 0.3],'linewidth',2)

relevInds = (2:numTrial+1) + numTrial;
coefVals = glm.Coefficients.Estimate(relevInds);
CIbands = coefCI(glm);
errorL = abs(coefVals - CIbands(relevInds,1));
errorU = abs(coefVals - CIbands(relevInds,2));
errorbar((1:numTrial),coefVals,errorL,errorU,'Color', [0 0 0],'linewidth',2)

%% for four numTrials
% first two trials
relevInds = 2:numTrial+1;
coefVals = glm.Coefficients.Estimate(relevInds);
CIbands = coefCI(glm);
errorL = abs(coefVals - CIbands(relevInds,1));
errorU = abs(coefVals - CIbands(relevInds,2));
figure; hold on;
errorbar((1:numTrial),coefVals,errorL,errorU,'Color', [0 0 0],'linewidth',2)

relevInds = relevInds + numTrial;
coefVals = glm.Coefficients.Estimate(relevInds);
CIbands = coefCI(glm);
errorL = abs(coefVals - CIbands(relevInds,1));
errorU = abs(coefVals - CIbands(relevInds,2));
errorbar((1:numTrial),coefVals,errorL,errorU,'Color', [1 0 0],'linewidth',2)
legend({'1', '2'})

% last two trials
relevInds = relevInds + numTrial;
coefVals = glm.Coefficients.Estimate(relevInds);
CIbands = coefCI(glm);
errorL = abs(coefVals - CIbands(relevInds,1));
errorU = abs(coefVals - CIbands(relevInds,2));
figure; hold on;
errorbar((1:numTrial),coefVals,errorL,errorU,'Color', [0 0 0],'linewidth',2)

relevInds = relevInds + numTrial;
coefVals = glm.Coefficients.Estimate(relevInds);
CIbands = coefCI(glm);
errorL = abs(coefVals - CIbands(relevInds,1));
errorU = abs(coefVals - CIbands(relevInds,2));
errorbar((1:numTrial),coefVals,errorL,errorU,'Color', [1 0 0],'linewidth',2)
legend({'3', '4'})
%% for interactions
% terms
relevInds = 2:numTrial+1;
coefVals = glm.Coefficients.Estimate(relevInds);
CIbands = coefCI(glm);
errorL = abs(coefVals - CIbands(relevInds,1));
errorU = abs(coefVals - CIbands(relevInds,2));
figure; 
subplot(1,2,1);hold on;
errorbar((1:numTrial),coefVals,errorL,errorU,'Color', [0.7 0 1],'linewidth',2)
legend({'1'})

relevInds = relevInds + numTrial;
coefVals = glm.Coefficients.Estimate(relevInds);
CIbands = coefCI(glm);
errorL = abs(coefVals - CIbands(relevInds,1));
errorU = abs(coefVals - CIbands(relevInds,2));
subplot(1,2,2);hold on;
errorbar((1:numTrial),coefVals,errorL,errorU,'Color', [0 0 1],'linewidth',2)
legend({'2'})

% interaction terms
relevInds = relevInds + numTrial;
sigs = glm.Coefficients.pValue(relevInds);
sigs(sigs>=0.05) = NaN;
sigs(sigs<0.05) = -0.1;
subplot(1,2,1); hold on;
scatter(1:numTrial, sigs, 15, "magenta", "filled")

relevInds = relevInds + numTrial;
sigs = glm.Coefficients.pValue(relevInds);
sigs(sigs>=0.05) = NaN;
sigs(sigs<0.05) = 0.05;
subplot(1,2,2); hold on;
scatter(1:numTrial, sigs, 15, 'k', '*')
%% lickLat linear model
combinedLatMatE = combinedLatMat;
combinedLatMatE(combinedStatesMat==0) = 0;
combinedLatMatE(isnan(combinedStatesMat)) = NaN;
combinedLatMatT = combinedLatMat;
combinedLatMatT(combinedStatesMat==1) = 0;
combinedLatMatT(isnan(combinedStatesMat)) = NaN;

combinedRwdMatE = abs(combinedRwdMat);
combinedRwdMatE(combinedStatesMat==0) = 0;
combinedRwdMatE(isnan(combinedStatesMat)) = NaN;
combinedRwdMatT = abs(combinedRwdMat);
combinedRwdMatT(combinedStatesMat==1) = 0;
combinedRwdMatT(isnan(combinedStatesMat)) = NaN;

lm = fitlm([combinedRwdMatE', combinedRwdMatT', combinedLatMat', combinedStates',  combinedSvS', combinedStates'.*combinedSvS'], combinedLat');
%% lick interaction
    
lm2 = fitlm([abs(combinedRwdMat)', abs(combinedRwdMat)'.*combinedStatesMat', combinedLatMat', combinedStatesMat', combinedSvS', combinedStates'.*combinedSvS'], combinedLat');

%% 
figure2;
relevInds = (size(lm.CoefficientNames,2)-2):(size(lm.CoefficientNames,2));
coefVals = glm.Coefficients.Estimate(relevInds);
CIbands = coefCI(glm);
errorL = abs(coefVals - CIbands(relevInds,1));
errorU = abs(coefVals - CIbands(relevInds,2));
errorbar((0.2*1:length(relevInds)),coefVals,errorL,errorU,'Color', [0 0 0],'linewidth',2, 'LineStyle', 'none')
ylabel('Coeff')
xlabel(['States       ', 'Svs       ', 'Inter      '])
set(gca, 'Box', 'off')
set(gca, 'TickDir', 'out')
%%
lm = fitlm([abs(combinedRwdMat)', abs(combinedRwdMat)'.*combinedStatesMat', combinedStatesMat', combinedLatMat', combinedSvS', combinedStates'.*combinedSvS'], combinedLat');
%% holm-bonferroni correction
% interaction terms
mat = [combinedNrwdMat', combinedChoicesMat', combinedNrwdMat'.*combinedStatesMat', combinedChoicesMat'.*combinedStatesMat', combinedStatesMat'];
glm = fitglm(mat, combinedChoices');
rwdInds =  (2:(numTrial+1)) + 2*numTrial;
choiceInds = (2:(numTrial+1)) + 3*numTrial;
%%
% rwd correction
trialInd = 1:numTrial;
alpha = 0.05./((numTrial+1) - (1:numTrial));
pVRwd = glm.Coefficients.pValue(rwdInds);
[pVRwd, sortInd] = sort(pVRwd);
trialInd = trialInd(sortInd);
sigs = pVRwd' < alpha;
stopInd = find(sigs == 0, 1);
if stopInd
    sigs(stopInd:end) = 0;
end
trialInd = trialInd(sigs);

scatter(trialInd, 0.05*ones(size(trialInd)), 15, 'k', 'Marker', '*');
%%
% choice correction
trialInd = 1:numTrial;
alpha = 0.05./((numTrial+1) - (1:numTrial));
pVCh = glm.Coefficients.pValue(choiceInds);
[pVCh, sortInd] = sort(pVCh);
trialInd = trialInd(sortInd);
sigs = pVCh' < alpha;
stopInd = find(sigs == 0, 1);
if stopInd
    sigs(stopInd:end) = 0;
end
trialInd = trialInd(sigs);
scatter(trialInd, 0.05*ones(size(trialInd)), 15, 'k', 'Marker', '*');
%% calculate beta 
% load data
clear all;
load('F:\tmpData\hmmDataAllAni.mat', 'allAnis', 'cols');
modelName = '5paramsHmm';
[root,sep] = currComputer();
%%
betaE = NaN(length(allAnis), 3);
betaT = NaN(length(allAnis), 3);
betaDiff = NaN(length(allAnis), 3);
betaSum = NaN(length(allAnis), 3);

for i = 1:length(allAnis)
    currAni = allAnis{i};
    currCol = cols{i};
    fileName = [currAni currCol '_' modelName];
    modelPath = [root currAni sep currAni 'sorted' sep 'stan' sep 'bernoulli' sep modelName sep currCol sep fileName '.mat'];
    load(modelPath);
    eval(['fit =' fileName ';']);
    eval(['clear ' fileName ';']);
    mu_E = fit.mu_betaE(fit.divergent__==0);
    mu_T = fit.mu_betaT(fit.divergent__==0);
    betaE(i,:) = quantile(mu_E, [0.05 0.5 0.95]);
    betaT(i,:) = quantile(mu_T, [0.05 0.5 0.95]);
    betaDiff(i,:) = quantile(0.5*(mu_E-mu_T), [0.05 0.5 0.95]);
    betaSum(i,:) = quantile(0.5*(mu_E+mu_T), [0.05 0.5 0.95]);
end
%% plotting
yNeg = betaE(:,2) - betaE(:,1);
yPos = betaE(:,3) - betaE(:,2);

xNeg = betaT(:,2) - betaT(:,1);
xPos = betaT(:,3) - betaT(:,2);
figure2; hold on;
errorbar(betaT(:,2), betaE(:,2), yNeg, yPos, xNeg, xPos, 'LineStyle', 'none', 'Color', [0.5 0.5 0.5]);
scatter(betaT(:,2), betaE(:,2), 18, 'MarkerEdgeColor', 'k', 'MarkerEdgeAlpha', 0.7, 'LineWidth', 1.5, 'MarkerFaceColor', 'w');
plot([min(betaE, [], 'all'), max(betaT, [], 'all')], [min(betaE, [], 'all'), max(betaT, [], 'all')], 'LineStyle', '--', 'Color', '0.7 0.7 0.7')
xlabel('betaT')
ylabel('betaE')
%%
%% plotting
yNeg = betaSum(:,2) - betaSum(:,1);
yPos = betaSum(:,3) - betaSum(:,2);

xNeg = betaDiff(:,2) - betaDiff(:,1);
xPos = betaDiff(:,3) - betaDiff(:,2);
figure2; hold on;
errorbar(betaDiff(:,2), betaSum(:,2), yNeg, yPos, xNeg, xPos, 'LineStyle', 'none', 'Color', [0.5 0.5 0.5]);
scatter(betaDiff(:,2), betaSum(:,2), 18, 'MarkerEdgeColor', 'k', 'MarkerEdgeAlpha', 0.7, 'LineWidth', 1.5, 'MarkerFaceColor', 'w');

xlabel('betaDiff')
ylabel('betaSum')
% plot([min(betaDiff, [], 'all'), max(betaSum, [], 'all')], [min(betaDiff, [], 'all'), max(betaSum, [], 'all')], 'LineStyle', '--', 'Color', '0.7 0.7 0.7')
%%















