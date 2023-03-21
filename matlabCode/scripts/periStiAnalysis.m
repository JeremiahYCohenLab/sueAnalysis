% pari-light analysis
load('F:\tmpData\allSessionHmmInhibition.mat');
ani = 'allGt';
sheet = 'inhibitionGt';
col = 'cueOnShamLate';
dayList = getDayList(sheet, ani, col);
[root, sep] = currComputer();
trialB = 3;
trialF = 3;
%% prepare figures
choiceFig = figure2;
swFig = figure2;
lcFig = figure2;
wcFig = figure2;


%%
% pre-location all params
periStiChoices = [];% choices aligned at laser, previous choice as 1, opposite as zeros
periStiChoicesNrwd = [];% choices aligned at laser, previous choice as 1, opposite as zeros, laser is on nrwd trials
periStiSwitches = [];% switches aligned at laser
periStiSwitchesNrwd = [];% switches aligned at laser, laser is on nrwd trials
periStiSwitchesNext = [];% switches aligned at laser
periStiSwitchesNextNrwd = [];% switches aligned at laser, laser is on nrwd trials
periStiRwds = [];% switches aligned at laser
periStiRwdsNrwd = [];% switches aligned at laser, laser is on nrwd trials
periStiHmm = [];% hmm states around light stimulations 
periStiSti = [];
periStiLatNrwd = [];

colorL = 'm';
% find all lasers
for i = 1:length(dayList)
    % getChoices, Switches and Hmm
    session = dayList{i};
    s = behAnalysisNoPlot_opMD(session, 'simpleFlag', 1);
%     hmmInd = cellfun(@(x) strcmp(session, x), allSession);
%     states = allStates{hmmInd};
%     transEInd = find(states(2:end-1)==1 & states(1:end-2)~=states(3:end) & states(1:end-2)~= 1 & states(3:end)~= 1)+1;
%     states(transEInd) = NaN;
    choices = s.allChoices;
    choices(choices==-1) = 0;
    svs = NaN(size(choices));
    svs(s.changeChoice_Inds) = 1;
    svs(s.stayChoice_Inds) = 0;
    svsNext = NaN(size(choices));
    svsNext(s.changeChoice_Inds-1) = 1;
    svsNext(s.stayChoice_Inds-1) = 0;
    % find trial Inds with laser
    laserInd = find(s.laser==1);
    noLaserInd = find(s.laser==0);
    
    % focus
    focus = noLaserInd;
    currColor = 'k';
    currTitle = 'laser, nrwd, next change';
    currInd = mintersect(focus, s.nrwd_Inds, s.changeChoice_Inds-1);
    %
    currInd = currInd(currInd<=length(choices));
    currInd = currInd + trialB;
    targetTemp = [NaN(1, trialB), choices, NaN(1, trialF)];
    laserTemp = [NaN(1, trialB), s.laser, NaN(1, trialF)];
    svsTemp = [NaN(1, trialB), svs, NaN(1, trialF)];
    svsNextTemp = [NaN(1, trialB), svsNext, NaN(1, trialF)];
    rwdTemp = [NaN(1, trialB), abs(s.allRewards), NaN(1, trialF)];
    lickLattemp = [NaN(1, trialB), s.lickLatZ, NaN(1, trialF)];
    for j = 1:length(currInd)
        % choice
        currSeq = targetTemp(currInd(j)-trialB:currInd(j)+trialF);
        if targetTemp(currInd(j))==0
            currSeq = 1-currSeq;
        end
        periStiChoicesNrwd = [periStiChoicesNrwd; currSeq];
        % sti
        currSeq = laserTemp(currInd(j)-trialB:currInd(j)+trialF);
        periStiSti = [periStiSti; currSeq];
        % Switch
        currSeq = svsTemp(currInd(j)-trialB:currInd(j)+trialF);
        periStiSwitchesNrwd = [periStiSwitchesNrwd; currSeq];
        % rwd
        currSeq = rwdTemp(currInd(j)-trialB:currInd(j)+trialF);
        periStiRwdsNrwd = [periStiRwdsNrwd; currSeq];
        % swtich next 
        currSeq = svsNextTemp(currInd(j)-trialB:currInd(j)+trialF);
        periStiSwitchesNextNrwd = [periStiSwitchesNextNrwd; currSeq];
        % lick lat
        currSeq = lickLattemp(currInd(j)-trialB:currInd(j)+trialF);
        periStiLatNrwd = [periStiLatNrwd; currSeq];
    end
    
end
%% plot summary plots
% choices
figure(choiceFig);
hold on; plotFilledBern(-trialB:1:trialF, periStiChoicesNrwd, currColor);
ylabel('P(choice=C(0))')
title(currTitle)
%% switches
figure(swFig);
hold on; plotFilledBern(-trialB:1:trialF, periStiSwitchesNrwd, currColor);
ylabel('P(switch)')
title(currTitle)
%% P(switch|nrwd)
nrwdInd = periStiRwdsNrwd==0;
nrwdSwNextMat = periStiSwitchesNextNrwd;
nrwdSwNextMat(~nrwdInd) = 0;
nrwdNum = sum(nrwdInd, 1, 'omitnan');
nrwdSwNum = sum(nrwdSwNextMat, 1, 'omitnan');

meanLC = nrwdSwNum./nrwdNum;
semLC = sem_bernoulli(nrwdSwNum, nrwdNum);
%%
figure(lcFig); hold on;
errorbar(-trialB:1:trialF,  meanLC, semLC, 'Color', currColor, 'LineWidth', 2)

% plot(-trialB:1:trialF, meanLC, 'Color', currColor);
% patch([-trialB:1:trialF, trialF:-1:-trialB], [meanLC-semLC, flip(meanLC+semLC)], currColor, 'edgeColor', 'none', 'FaceAlpha',0.3);
% ylabel('P(switch|nrwd)')
% title(currTitle)
%% 
% 
% plot(-trialB:1:trialF, meanLC, 'Color', 'k');
% patch([-trialB:1:trialF, trialF:-1:-trialB], [meanLC-semLC, flip(meanLC+semLC)], 'k', 'edgeColor', 'none', 'FaceAlpha',0.3);
%% P(switch|rwd)
rwdInd = periStiRwdsNrwd==1;
rwdSwNextMat = periStiSwitchesNextNrwd;
rwdSwNextMat(~rwdInd) = 0;
rwdNum = sum(rwdInd, 1, 'omitnan');
rwdSwNum = sum(rwdSwNextMat, 1, 'omitnan');

meanWC = rwdSwNum./rwdNum;
semWC = sem_bernoulli(rwdSwNum, rwdNum);
%%
figure(wcFig); hold on;
errorbar(-trialB:1:trialF,  meanWC, semWC, 'Color', currColor, 'LineWidth', 2)
% plot(-trialB:1:trialF, meanWC, 'Color', colorL);
% patch([-trialB:1:trialF, trialF:-1:-trialB], [meanWC-semWC, flip(meanWC+semWC)], colorL, 'edgeColor', 'none', 'FaceAlpha',0.3);
ylabel('P(switch|rwd)')
title(currTitle)
%% 
% errorbar(-trialB:1:trialF,  meanWC, semWC, 'Color', 'k', 'LineWidth', 2)
%%
% plot(-trialB:1:trialF, meanWC, 'Color', 'k');
% patch([-trialB:1:trialF, trialF:-1:-trialB], [meanWC-semWC, flip(meanWC+semWC)], 'k', 'edgeColor', 'none', 'FaceAlpha',0.3);
%%






%% peri-laser lick analysis

load('F:\tmpData\allSessionHmmInhibition.mat');
ani = 'allGt';
sheet = 'inhibitionGt';
col = 'cueOnGood';
dayList = getDayList(sheet, ani, col);
[root, sep] = currComputer();
trialB = 3;
trialF = 3;
%%
% pre-location all params
periStiSwitchesNrwd = [];
periStiLatNrwd = [];
periStiRate = [];
periStiRwds = [];
colorL = 'm';
% find all lasers
for i = 1:length(dayList)
    % getChoices, Switches and Hmm
    session = dayList{i};
    s = behAnalysisNoPlot_opMD(session, 'simpleFlag', 1);
%     hmmInd = cellfun(@(x) strcmp(session, x), allSession);
%     states = allStates{hmmInd};
%     transEInd = find(states(2:end-1)==1 & states(1:end-2)~=states(3:end) & states(1:end-2)~= 1 & states(3:end)~= 1)+1;
%     states(transEInd) = NaN;
    choices = s.allChoices;
    choices(choices==-1) = 0;
    svs = NaN(size(choices));
    svs(s.changeChoice_Inds) = 1;
    svs(s.stayChoice_Inds) = 0;
    svsNext = NaN(size(choices));
    svsNext(s.changeChoice_Inds-1) = 1;
    svsNext(s.stayChoice_Inds-1) = 0;
    % find trial Inds with laser
    laserInd = find(s.laser==1);
    noLaserInd = find(s.laser==0);
    
    % focus
    focus = noLaserInd;
    currColor = 'k';
%     focus = randperm(length(choices), floor(0.3*length(choices)));
    % find peri-laser choices & lasers.
    % focus on lasered + nrwd trials
    currInd = mintersect(focus);
    currTitle = 'laser, nrwd';
    currInd = currInd(currInd<=length(choices));
    currInd = currInd + trialB;
    targetTemp = [NaN(1, trialB), choices, NaN(1, trialF)];
    laserTemp = [NaN(1, trialB), s.laser, NaN(1, trialF)];
    svsTemp = [NaN(1, trialB), svs, NaN(1, trialF)];
    svsNextTemp = [NaN(1, trialB), svsNext, NaN(1, trialF)];
    rwdTemp = [NaN(1, trialB), abs(s.allRewards), NaN(1, trialF)];
    lickLattemp = [NaN(1, trialB), s.lickLatZ, NaN(1, trialF)];
    lickRatetemp = [NaN(1, trialB), zscore(s.lickNumRwd), NaN(1,trialF)];
    for j = 1:length(currInd)
        % choice
        currSeq = rwdTemp(currInd(j)-trialB:currInd(j)+trialF);
        periStiRwds = [periStiRwds; currSeq];

        % Switch
        currSeq = svsTemp(currInd(j)-trialB:currInd(j)+trialF);
        periStiSwitchesNrwd = [periStiSwitchesNrwd; currSeq];

        % lick lat
        currSeq = lickLattemp(currInd(j)-trialB:currInd(j)+trialF);
        periStiLatNrwd = [periStiLatNrwd; currSeq];
        
        % lick Rate
        currSeq = lickRatetemp(currInd(j)-trialB:currInd(j)+trialF);
        periStiRate = [periStiRate; currSeq];
    end
    
end
%% prepare figures
rawLatFig = figure2;
stayFig = figure2;
switchFig = figure2;
%% all lick Lats
figure(rawLatFig);hold on; 
meanLat = mean(periStiLatNrwd, 'omitnan');
semLat = sem(periStiLatNrwd);
errorbar(-trialB:1:trialF, meanLat, semLat, 'Color', currColor, 'LineWidth', 2);
title('allTrials')

%% all stay lats/all switch lats

for i = 1:(trialF+trialB+1)
    meanLatSw(i) = mean(periStiLatNrwd(periStiSwitchesNrwd(:,i)==1,i), 'omitnan');
    meanLatStay(i) = mean(periStiLatNrwd(periStiSwitchesNrwd(:,i)==0,i), 'omitnan');
    semLatSw(i) = sem(periStiLatNrwd(periStiSwitchesNrwd(:,i)==1,i));
    semLatStay(i) = sem(periStiLatNrwd(periStiSwitchesNrwd(:,i)==0,i));
end
%% 
figure(stayFig); hold on;
errorbar(-trialB:1:trialF, meanLatStay, semLatStay, 'Color', currColor, 'LineWidth', 2);
title('stays')

figure(switchFig); hold on;
errorbar(-trialB:1:trialF, meanLatSw, semLatSw, 'Color', currColor, 'LineWidth', 2);
title('switches')
%%


%% prepare figures
rawRateFig = figure2;
stayRateFig = figure2;
switchRateFig = figure2;


%% all stay lats/all switch lats
rwdInd = 0;
for i = 1:(trialF+trialB+1)
    meanRate(i) = mean(periStiRate(periStiRwds(:,i)==rwdInd,i), 'omitnan');
    meanRateSw(i) = mean(periStiRate(periStiSwitchesNrwd(:,i)==1&periStiRwds(:,i)==rwdInd,i), 'omitnan');
    meanRateStay(i) = mean(periStiRate(periStiSwitchesNrwd(:,i)==0&periStiRwds(:,i)==rwdInd,i), 'omitnan');
    semRate(i) = sem(periStiRate(periStiRwds(:,i)==rwdInd,i));
    semRateSw(i) = sem(periStiRate(periStiSwitchesNrwd(:,i)==1&periStiRwds(:,i)==rwdInd,i));
    semRateStay(i) = sem(periStiRate(periStiSwitchesNrwd(:,i)==0&periStiRwds(:,i)==rwdInd,i));
end

%% all lick rates
figure(rawRateFig);hold on; 
errorbar(-trialB:1:trialF, meanRate, semRate, 'Color', currColor, 'LineWidth', 2);
title('allTrials rate')

figure(stayRateFig); hold on;
errorbar(-trialB:1:trialF, meanRateStay, semRateStay, 'Color', currColor, 'LineWidth', 2);
title('stays rate')

figure(switchRateFig); hold on;
errorbar(-trialB:1:trialF, meanRateSw, semRateSw, 'Color', currColor, 'LineWidth', 2);
title('switches rate')
%%
