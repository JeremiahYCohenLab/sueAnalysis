load('F:\tmpData\allUnitAUC.mat');
load('F:\tmpData\hmmNeuronData.mat');
%% get all states
% prior = 10, 10
xlFile = {'ZS059', 'ZS060', 'ZS061', 'ZS062'};
sheet = {'ZS059', 'ZS060', 'ZS061', 'ZS062'};
col = 'good';
[root, sep] = currComputer();
modelName = '5params';

allEmis = {};
allQ = {};
allPe = {};
allStates = {};
allSessionHmm = {};
allRewards = {};
allSvSNext = {};
transGuess = [0.4 0.3 0.3;
              0.3 0.7 0;
              0.3 0 0.7];
emisGuess = [0.5 0.5;
             1 0;
             0 1];
for ani = 1:length(sheet)
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
        spikeCount{i} = spikeFocus;
        
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
    allSessionHmm = [allSessionHmm; dayList];
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
%     states(transEInd) = NaN;
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
spikeExploit = zeros(length(allSessions), numBins);
spikeExplore = zeros(length(allSessions), numBins);
spike = zeros(length(allSessions), numBins);
RPE = zeros(length(allSessions), numBins);
RPEExploit = zeros(length(allSessions), numBins);
RPEExplore = zeros(length(allSessions), numBins);


for i = 1:length(allSessions)
    session = allSessions{i};
    hmmInd = cellfun(@(x)strcmp(session, x), allSessionHmm);
    states = allStates{hmmInd};
    % remove those single trials in trainsition states
    transEInd = find(states(2:end-1)==1 & states(1:end-2)~=states(3:end) & states(1:end-2)~= 1 & states(3:end)~= 1)+1;
    states(transEInd) = NaN;
    pe = allPe{hmmInd};
    spikeCounts = spikeFocus{i};
    spikeCounts(~isnan(spikeCounts)) = zscore(spikeCounts(~isnan(spikeCounts)));
    spikeT = spikeCounts(states==2 | states==3);
    spikeE = spikeCounts(states==1);
    peT = pe(states==2 | states==3);
    peE = pe(states==1);
    target = peT;
    edgesT = [linspace(min(target)- 0.01, 0, 0.5*numBins+1) linspace(0, max(target)+ 0.01,0.5*numBins+1)];
    edgesT = [edgesT(1:0.5*numBins+1), edgesT(0.5*numBins+3:end)];
    target = peE;
    edgesE = [linspace(min(target)- 0.01, 0, 0.5*numBins+1) linspace(0, max(target)+ 0.01,0.5*numBins+1)];
    edgesE = [edgesE(1:0.5*numBins+1), edgesE(0.5*numBins+3:end)];
    target = pe;
    edges = [linspace(min(target)- 0.01, 0, 0.5*numBins+1) linspace(0, max(target)+ 0.01,0.5*numBins+1)];
    edges = [edges(1:0.5*numBins+1), edges(0.5*numBins+3:end)];
    for j = 1:numBins
        spikeExploit(i, j) = mean(spikeT(peT>=edgesT(j) & peT<edgesT(j+1)), 'omitnan');
        spikeExplore(i, j) = mean(spikeE(peE>=edgesE(j) & peE<edgesE(j+1)), 'omitnan');
        RPEExploit(i, j) = mean(peT(peT>=edgesT(j) & peT<edgesT(j+1)), 'omitnan');
        RPEExplore(i, j) = mean(peE(peE>=edgesE(j) & peE<edgesE(j+1)), 'omitnan');
        spike(i, j) = mean(spikeCounts(pe>=edges(j) & pe<edges(j+1)), 'omitnan');
        RPE(i,j) =  mean(pe(pe>=edges(j) & pe<edges(j+1)), 'omitnan');
    end
%     allQt = [allQt; spikeT];
%     allQe = [allQe; spikeE];
%     allChoicesT = [allChoicesT peT];
%     allChoicesE = [allChoicesE peE];
end
%% plot mean of all sessions, both states
coeffsMax = [lmMaxAll.coeffs]';
tStatsMax = [lmMaxAll.tStats]';
sigMax = [lmMaxAll.ps]'<0.05;
color1 = [0.2 0.2 0.2];
color2 = [0.8 0.8 0.8];
figure2;
hold on;
plotFilled(mean(RPE(sigMax(:,1) & tStatsMax(:,1)>0 & ind==2,:), 'omitnan'), spike(sigMax(:,1) & tStatsMax(:,1)>0 & ind==2,:), color2);
plotFilled(mean(RPE(sigMax(:,1) & tStatsMax(:,1)<0 & ind==1,:), 'omitnan'), spike(sigMax(:,1) & tStatsMax(:,1)<0 & ind==1,:), color1);

set(gca, 'Box', 'off')
set(gca,'tickdir', 'out')
set(gca, 'XTick', [-1 0 1])
set(gca, 'YTick', [-0.5 0 0.5])
xlabel('rpe', 'FontSize', 18)
ylabel('spikes/s (zscored)', 'FontSize', 18)
%% plot mean of all sessions, separate states
figure2;
% exploitation
subplot(1,2,1);
hold on;
plotFilled(mean(RPEExploit(sigMax(:,1) & tStatsMax(:,1)>0 & ind==2,:), 'omitnan'), spikeExploit(sigMax(:,1) & tStatsMax(:,1)>0 & ind==2,:), color2);
plotFilled(mean(RPEExploit(sigMax(:,1) & tStatsMax(:,1)<0 & ind==1,:), 'omitnan'), spikeExploit(sigMax(:,1) & tStatsMax(:,1)<0 & ind==1,:), color1);

set(gca, 'Box', 'off')
set(gca,'tickdir', 'out')
set(gca, 'XTick', [-1 0 1])
set(gca, 'YTick', [-0.5 0 0.5])
xlabel('rpe', 'FontSize', 18)
ylabel('spikes/s (zscored)', 'FontSize', 18)
title('exploit', 'Color', 'k', 'FontSize', 18)

subplot(1,2,2);
hold on;
plotFilled(mean(RPEExplore(sigMax(:,1) & tStatsMax(:,1)>0 & ind==2,:), 'omitnan'), spikeExplore(sigMax(:,1) & tStatsMax(:,1)>0 & ind==2,:), color2);
plotFilled(mean(RPEExplore(sigMax(:,1) & tStatsMax(:,1)<0 & ind==1,:), 'omitnan'), spikeExplore(sigMax(:,1) & tStatsMax(:,1)<0 & ind==1,:), color1);

set(gca, 'Box', 'off')
set(gca,'tickdir', 'out')
set(gca, 'XTick', [-1 0 1])
set(gca, 'YTick', [-0.5 0 0.5])
xlabel('rpe', 'FontSize', 18)
ylabel('spikes/s (zscored)', 'FontSize', 18)
title('explore', 'Color', 'r', 'FontSize', 18)

legend({'TypeI','', 'TypeII',''})
%% Sep by two plots
figure2;
colorE = [1 0 0];
colorT = [0 0 0];
% Type I 
subplot(1,2,1);
hold on;
plotFilled(mean(RPEExploit(sigMax(:,1) & tStatsMax(:,1)>0 & ind==2,:), 'omitnan'), spikeExploit(sigMax(:,1) & tStatsMax(:,1)>0 & ind==2,:), colorT);
plotFilled(mean(RPEExplore(sigMax(:,1) & tStatsMax(:,1)>0 & ind==2,:), 'omitnan'), spikeExplore(sigMax(:,1) & tStatsMax(:,1)>0 & ind==2,:), colorE);

set(gca, 'Box', 'off')
set(gca,'tickdir', 'out')
set(gca, 'XTick', [-1 0 1])
set(gca, 'YTick', [-0.5 0 0.5])
xlabel('rpe', 'FontSize', 18)
ylabel('spikes/s (zscored)', 'FontSize', 18)
title('Type I', 'FontSize', 18, 'Color', color1)

% Type II
subplot(1,2,2);
hold on;
plotFilled(mean(RPEExploit(sigMax(:,1) & tStatsMax(:,1)<0 & ind==1,:), 'omitnan'), spikeExploit(sigMax(:,1) & tStatsMax(:,1)<0 & ind==1,:), colorT);
plotFilled(mean(RPEExplore(sigMax(:,1) & tStatsMax(:,1)<0 & ind==1,:), 'omitnan'), spikeExplore(sigMax(:,1) & tStatsMax(:,1)<0 & ind==1,:), colorE);

set(gca, 'Box', 'off')
set(gca,'tickdir', 'out')
set(gca, 'XTick', [-1 0 1])
set(gca, 'YTick', [-0.5 0 0.5])
xlabel('rpe', 'FontSize', 18)
ylabel('spikes/s (zscored)', 'FontSize', 18)
title('Type II', 'FontSize', 18, 'Color', color2)
legend({'exploit','', 'explore',''})
%% baseline analysis
tb = 1;
tf = 0;
stepSize = 50;
binSize = 50;
numBins = 6;

spikeExploit = zeros(length(allSessions), numBins);
spikeExplore = zeros(length(allSessions), numBins);
spike = zeros(length(allSessions), numBins);

QExploit = zeros(length(allSessions), numBins);
QExplore = zeros(length(allSessions), numBins);
QchosenAll = zeros(length(allSessions), numBins);
auc = zeros(length(allSessions),1);
for i = 1:length(allSessions)
    session = allSessions{i};
    unit = allUnits{i};
    s = behAnalysisNoPlot_opMD(session, 'simpleFlag', 1);
    hmmInd = cellfun(@(x)strcmp(session, x), allSessionHmm);
    Q = allQ{hmmInd};
    Qchosen = NaN(length(s.allChoices),1);
    Qchosen(s.allChoices==1) = Q(s.allChoices==1, 2);
    Qchosen(s.allChoices==-1) = Q(s.allChoices==-1, 1);
    [~, matCue, ~] = getUnitMatCue(session, unit, tb, tf, stepSize, binSize);
    spikeCounts = sum(matCue, 2);

    states = allStates{hmmInd};
    % remove those single trials in trainsition states
    transEInd = find(states(2:end-1)==1 & states(1:end-2)~=states(3:end) & states(1:end-2)~= 1 & states(3:end)~= 1)+1;
    states(transEInd) = NaN;
    spikeCounts = zscore(spikeCounts);
    spikeT = spikeCounts(states==2 | states==3);
    spikeE = spikeCounts(states==1);
    Qt = Qchosen(states==2 | states==3);
    Qe = Qchosen(states==1);
    target = Qt;
    edgesT = [linspace(min(target)- 0.01, max(target)+ 0.01, numBins+1)];
    target = Qe;
    edgesE = [linspace(min(target)- 0.01, max(target)+ 0.01, numBins+1)];
    target = Qchosen;
    edges = [linspace(min(target)- 0.01, max(target)+ 0.01, numBins+1)];
    auc(i) = auROCZS(spikeT, spikeE);
    for j = 1:numBins
        spikeExploit(i, j) = mean(spikeT(Qt>=edgesT(j) & Qt<edgesT(j+1)), 'omitnan');
        spikeExplore(i, j) = mean(spikeE(Qe>=edgesE(j) & Qe<edgesE(j+1)), 'omitnan');
        QExploit(i, j) = mean(Qt(Qt>=edgesT(j) & Qt<edgesT(j+1)), 'omitnan');
        QExplore(i, j) = mean(Qe(Qe>=edgesE(j) & Qe<edgesE(j+1)), 'omitnan');
        spike(i, j) = mean(spikeCounts(Qchosen>=edges(j) & Qchosen<edges(j+1)), 'omitnan');
        QchosenAll(i,j) =  mean(Qchosen(Qchosen>=edges(j) & Qchosen<edges(j+1)), 'omitnan');
    end
end

%%
%% plot mean of all sessions, both states
coeffsMax = [lmMaxAll.coeffs]';
tStatsMax = [lmMaxAll.tStats]';
sigMax = [lmMaxAll.ps]'<0.05;
color1 = [0.2 0.2 0.2];
color2 = [0.8 0.8 0.8];
figure2;
hold on;
plotFilled(mean(QchosenAll(ind==2,:), 'omitnan'), spike(sigMax(:,1) & tStatsMax(:,1)>0 & ind==2,:), color2);
plotFilled(mean(QchosenAll(ind==1,:), 'omitnan'), spike(sigMax(:,1) & tStatsMax(:,1)<0 & ind==1,:), color1);

set(gca, 'Box', 'off')
set(gca,'tickdir', 'out')
set(gca, 'XTick', [ 0 1])
set(gca, 'YTick', [-0.5 0 0.5])
xlabel('Qchosen', 'FontSize', 18)
ylabel('spikes/s (zscored)', 'FontSize', 18)
%% plot mean of all sessions, separate states
figure2;
% exploitation
subplot(1,2,1);
hold on;
plotFilled(mean(QExploit(sigMax(:,1) & tStatsMax(:,1)>0 & ind==2,:), 'omitnan'), spikeExploit(sigMax(:,1) & tStatsMax(:,1)>0 & ind==2,:), color2);
plotFilled(mean(QExploit(sigMax(:,1) & tStatsMax(:,1)<0 & ind==1,:), 'omitnan'), spikeExploit(sigMax(:,1) & tStatsMax(:,1)<0 & ind==1,:), color1);

set(gca, 'Box', 'off')
set(gca,'tickdir', 'out')
set(gca, 'XTick', [-1 0 1])
set(gca, 'YTick', [0:1:5])
xlabel('Qchosen', 'FontSize', 18)
ylabel('spikes/s (zscored)', 'FontSize', 18)
title('exploit', 'Color', 'k', 'FontSize', 18)

subplot(1,2,2);
hold on;
plotFilled(mean(QExplore(sigMax(:,1) & tStatsMax(:,1)>0 & ind==2,:), 'omitnan'), spikeExplore(sigMax(:,1) & tStatsMax(:,1)>0 & ind==2,:), color2);
plotFilled(mean(QExplore(sigMax(:,1) & tStatsMax(:,1)<0 & ind==1,:), 'omitnan'), spikeExplore(sigMax(:,1) & tStatsMax(:,1)<0 & ind==1,:), color1);

set(gca, 'Box', 'off')
set(gca,'tickdir', 'out')
set(gca, 'XTick', [-1 0 1])
set(gca, 'YTick', [0:1:5])
xlabel('Qchosen', 'FontSize', 18)
ylabel('spikes/s (zscored)', 'FontSize', 18)
title('explore', 'Color', 'r', 'FontSize', 18)

legend({'TypeI','', 'TypeII',''})
%% Sep by two plots
figure2;
colorE = [1 0 0];
colorT = [0 0 0];
% Type I 
subplot(1,2,1);
hold on;
plotFilled(mean(RPEExploit(sigMax(:,1) & tStatsMax(:,1)>0 & ind==2,:), 'omitnan'), spikeExploit(sigMax(:,1) & tStatsMax(:,1)>0 & ind==2,:), colorT);
plotFilled(mean(RPEExplore(sigMax(:,1) & tStatsMax(:,1)>0 & ind==2,:), 'omitnan'), spikeExplore(sigMax(:,1) & tStatsMax(:,1)>0 & ind==2,:), colorE);

set(gca, 'Box', 'off')
set(gca,'tickdir', 'out')
set(gca, 'XTick', [-1 0 1])
set(gca, 'YTick', [-0.5 0 0.5])
xlabel('rpe', 'FontSize', 18)
ylabel('spikes/s (zscored)', 'FontSize', 18)
title('Type I', 'FontSize', 18, 'Color', color1)

% Type II
subplot(1,2,2);
hold on;
plotFilled(mean(RPEExploit(sigMax(:,1) & tStatsMax(:,1)<0 & ind==1,:), 'omitnan'), spikeExploit(sigMax(:,1) & tStatsMax(:,1)<0 & ind==1,:), colorT);
plotFilled(mean(RPEExplore(sigMax(:,1) & tStatsMax(:,1)<0 & ind==1,:), 'omitnan'), spikeExplore(sigMax(:,1) & tStatsMax(:,1)<0 & ind==1,:), colorE);

set(gca, 'Box', 'off')
set(gca,'tickdir', 'out')
set(gca, 'XTick', [-1 0 1])
set(gca, 'YTick', [-0.5 0 0.5])
xlabel('rpe', 'FontSize', 18)
ylabel('spikes/s (zscored)', 'FontSize', 18)
title('Type II', 'FontSize', 18, 'Color', color2)
legend({'exploit','', 'explore',''})
%% colormap
myMap = [[linspace(0, 1, 100)', linspace(0, 1, 100)', linspace(1, 1, 100)'];
         [linspace(1, 1, 100)', linspace(1, 0, 100)', linspace(1, 0, 100)']];
figure2;
subplot(1,2,1)
imagesc(auc(ind==2));
ylim([0 max([sum(ind==1) sum(ind==2)])])
title('TypeI', 'Color', color1)
colormap(myMap)
caxis([0 1])

subplot(1,2,2)
imagesc(auc(ind==1));
ylim([0 max([sum(ind==1) sum(ind==2)])])
title('TypeII', 'Color', color2)
colormap(myMap)
caxis([0 1])
%% baseline analysis
% baseline E vs T
%% baseline analysis
tb = 1;
tf = 1;
stepSize = 100;
binSize = 400;

meanSpikeE = NaN(length(allSessions), 1);
meanSpikeT = NaN(length(allSessions), 1);

preTrialTrendE = [];
preTrialTrendT = [];

auc = zeros(length(allSessions),1);
for i = 1:length(allSessions)
    session = allSessions{i};
    unit = allUnits{i};
    s = behAnalysisNoPlot_opMD(session, 'simpleFlag', 1);
    hmmInd = cellfun(@(x)strcmp(session, x), allSessionHmm);
    Q = allQ{hmmInd};
    Qchosen = NaN(length(s.allChoices),1);
    Qchosen(s.allChoices==1) = Q(s.allChoices==1, 2);
    Qchosen(s.allChoices==-1) = Q(s.allChoices==-1, 1);
    [~, matCue, matCueSlide, timeSlide] = getUnitMatCue(session, unit, tb, 0, stepSize, binSize);
    spikeCounts = sum(matCue(:, 1:(tb*1000-400)), 2);
    
    states = allStates{hmmInd};
    % remove those single trials in trainsition states
    transEInd = find(states(2:end-1)==1 & states(1:end-2)~=states(3:end) & states(1:end-2)~= 1 & states(3:end)~= 1)+1;
    states(transEInd) = NaN;
    spikeCounts = spikeCounts;
    spikeT = spikeCounts(states==2 | states==3);
    spikeE = spikeCounts(states==1);
    % diff between E vs T
    auc(i) = auROCZS(spikeT, spikeE);
    meanSpikeE(i) = mean(spikeE);
    meanSpikeT(i) = mean(spikeT);
    % append trend before trial
    preTrialTrendE(i,:) = zscore(mean(matCueSlide(states==1,:), 'omitnan'));
    preTrialTrendT(i,:) = zscore(mean(matCueSlide(states==2|states==3,:), 'omitnan'));
end
%%
% mean baseline
figure2Wide;
% subplot(1,3,1);
% scatter(meanSpikeT, meanSpikeE)
% [h, p, ~, stats] = ttest(meanSpikeT, meanSpikeE);
% title(['p =', num2str(p, 2)]);
% hold on; plot([0 8], [0 8])
subplot(1,2,1);
scatter(meanSpikeT(ind==1), meanSpikeE(ind==1), 20, 'k');
[h, p, ~, stats] = ttest(meanSpikeT(ind==1), meanSpikeE(ind==1));
title(['p=', num2str(p, 2), 't=', num2str(stats.tstat), 'df=', num2str(stats.df)]);
hold on; plot([0 6], [0 6], 'Color', [0.6 0.6 0.6], 'LineStyle', '--');
set(gca, 'XTick', [0:2:6])
set(gca, 'tickdir', 'out')
set(gca, 'FontSize', 12)
set(gca, 'YTick', [0:2:6])
subplot(1,2,2);
scatter(meanSpikeT(ind==2), meanSpikeE(ind==2), 20, 'k')
[h, p, ~, stats] = ttest(meanSpikeT(ind==2), meanSpikeE(ind==2));
title(['p=', num2str(p, 2), 't=', num2str(stats.tstat), 'df=', num2str(stats.df)]);
hold on; plot([0 6], [0 6], 'Color', [0.6 0.6 0.6], 'LineStyle', '--')
set(gca, 'XTick', [0:2:6])
set(gca, 'tickdir', 'out')
set(gca, 'FontSize', 12)
set(gca, 'YTick', [0:2:6])
%% trend before cue
figure2; 
subplot(1,3,1); hold on; 
plotFilled(timeSlide, preTrialTrendE, 'r')
plotFilled(timeSlide, preTrialTrendT, 'k')
title('all')
subplot(1,3,2); hold on; 
plotFilled(timeSlide, preTrialTrendE(ind==1, :), 'r')
plotFilled(timeSlide, preTrialTrendT(ind==1, :), 'k')
title('typeII')
subplot(1,3,3); hold on; 
plotFilled(timeSlide, preTrialTrendE(ind==2, :), 'r')
plotFilled(timeSlide, preTrialTrendT(ind==2, :), 'k')
title('typeI')
%% calculate interactions
time = zscore(repmat(timeSlide(timeSlide+0.5*binSize<0)', 2*length(allSessions), 1));
spikeBl = [preTrialTrendE(:, timeSlide+0.5*binSize<0), preTrialTrendT(:, timeSlide+0.5*binSize<0)];
spikeBl = reshape(spikeBl', [], 1);
EvT = repmat([ones(sum(timeSlide+0.5*binSize<0), 1); zeros(sum(timeSlide+0.5*binSize<0), 1)], length(allSessions),1);
type = repmat(ind, 1, 2*sum(timeSlide+0.5*binSize<0));
type = reshape(type', [], 1);
tbl = table(time, spikeBl, EvT, type);
%%
lm = stepwiselm(tbl, 'spikeBl ~ 1+time*type*EvT');

%% plot rewards before start and before end of the exploration
trialS = 6;
trialE = 5;
spikeStart = [];
spikeEnd = [];
startType = [];
endType = [];

tb = 2;
tf = -0.4;

for i = 1:length(allSessions)
    session = allSessions{i};
    unit = allUnits{i};
    s = behAnalysisNoPlot_opMD(session, 'simpleFlag', 1);
    hmmInd = cellfun(@(x)strcmp(session, x), allSessionHmm);
    Q = allQ{hmmInd};
    Qchosen = NaN(length(s.allChoices),1);
    Qchosen(s.allChoices==1) = Q(s.allChoices==1, 2);
    Qchosen(s.allChoices==-1) = Q(s.allChoices==-1, 1);
    [~, matCue, matCueSlide, timeSlide] = getUnitMatCue(session, unit, tb, tf, stepSize, binSize);
    spikeCounts = zscore(sum(matCue, 2));
    
    states = allStates{hmmInd};
    % remove those single trials in trainsition states
    transEInd = find(states(2:end-1)==1 & states(1:end-2)~=states(3:end) & states(1:end-2)~= 1 & states(3:end)~= 1)+1;
    states(transEInd) = NaN;
    states(states==2|states==3) = 0;
    % before exploration start
    currInds = find(states(1:end-1) == 0 & states(2:end) == 1) + 1;
    % check animal did not change state while in trialS before exploration
    sumS = conv(states, ones(1, trialS));
    sumS = sumS(1:length(states));
    currInds = currInds(sumS(currInds)==1);
    % gather all starts together
    spikesTemp = [NaN(1, trialS-1), spikeCounts'];
    currInds = currInds + trialS - 1; % added trialS -1 to take care of appended trials
    if ~isempty(currInds)
        spikeMatTemp = NaN(length(currInds), trialS);
        for j = 1:trialS
            spikeMatTemp(:,j) = spikesTemp(currInds - (trialS - j))';
        end
        spikeStart = [spikeStart; spikeMatTemp];
        startType = [startType; ind(i)];
    end
    % before exploration end
    currInds = find(states(1:end-1) == 1 & states(2:end) == 0) + 1;
    % check animal did not change state while in trialS before exploration
    sumE = conv(states, ones(1, trialE));
    sumE = sumE(1:length(states));
    currInds = currInds(sumE(currInds)==trialE-1);
    % gather all ends together
    spikesTemp = [NaN(1, trialE-1), spikeCounts'];
    currInds = currInds + trialE - 1; % added trialS -1 to take care of appended trials
    if ~isempty(currInds)
        spikeMatTemp = NaN(length(currInds), trialE);
        for j = 1:trialE
            spikeMatTemp(:,j) = spikesTemp(currInds - (trialE - j))';
        end
        spikeEnd = [spikeEnd; spikeMatTemp];
        endType = [endType; ind(i)];
    end
end
%% plot
figure2Wide;
subplot(1,2,1); hold on;
meanS = mean(spikeStart(startType==1,:), 'omitnan');
semS = sem(spikeStart(startType==1,:));
errorbar([1:trialS]-trialS, meanS, semS, 'Color', color1, 'LineWidth', 2);

meanS = mean(spikeStart(startType==2,:), 'omitnan');
semS = sem(spikeStart(startType==2,:));
errorbar([1:trialS]-trialS, meanS, semS, 'Color', color2, 'LineWidth', 2);

legend({'II', 'I'})

xlabel('from start of exploration')
ylabel('spikes (zscored)')

subplot(1,2,2); hold on;
meanS = mean(spikeEnd(endType==1,:), 'omitnan');
semS = sem(spikeEnd(endType==1,:));
errorbar([1:trialE]-trialE, meanS, semS, 'Color', color1, 'LineWidth', 2);

meanS = mean(spikeEnd(endType==2,:), 'omitnan');
semS = sem(spikeEnd(endType==2,:));
errorbar([1:trialE]-trialE, meanS, semS, 'Color', color2, 'LineWidth', 2);

legend({'II', 'I'})

xlabel('from start of exploitation')
ylabel('spikes (zscored)')
%%