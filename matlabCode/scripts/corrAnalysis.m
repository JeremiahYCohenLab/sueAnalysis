dayList = getDayList(xlFile, sheet, category);
for i = 1:length(dayList)
    session = dayList{i};
    % paths
        [animalName, date] = strtok(session, 'd'); 
        animalName = animalName(2:end);
        pd = parseSessionString_df(session, root, sep);
        sortedFolderLocation = [pd.sortedFolder 'session' sep];
        pupilFile = [sortedFolderLocation session '_pupil.mat'];
    % load pupil
    if exist(pupilFile,'file')
        load(pupilFile)
    else
        fprintf(['no pupi in ' session '\n'])
        continue
    end
    if ~exist('sessionPupilChoice', 'var')
        fprintf(['Re-align in ' session '\n']);
        timeAlign(dayList{i}, 1, 1);
    end
    if errorProp>0.5
        fprintf(['No good align in ' session '\n']);
    end
    clear 'sessionPupilChoice'
    clear 'errorProp'
end
%% find pairs
% load F:\tmpData\catWithQchosen.mat
load F:\tmpData\catWithWFFeatures.mat
% get days with multiple units
[uniqueSessions, ia, ic] = unique(allSessions');
sessionPair = {};
unit1Pair = {};
unit2Pair = {};
unit1group = [];
unit2group = [];
for i = 1:length(uniqueSessions)
    if length(find(ic==i)) >= 2
        unitsTmp = allUnits(ic==i);  
        combTmp = nchoosek(1:length(unitsTmp),2);
        for j = 1:size(combTmp,1)
            sessionPair = [sessionPair; uniqueSessions{i}];
            unit1Pair = [unit1Pair; unitsTmp(combTmp(j,1))];
            unit2Pair = [unit2Pair; unitsTmp(combTmp(j,2))];
            unit1group = [unit1group ind(find(ic==i, 1, 'first')+combTmp(j,1)-1)];
            unit2group = [unit2group ind(find(ic==i, 1, 'first')+combTmp(j,2)-1)];
        end
    end
end
%%
clear allCorr
for i = 1:length(sessionPair)
    allCorr(i) = unitCorr(sessionPair{i}, unit1Pair{i}, unit2Pair{i}, 'binSize', 200, 'binSizePost', 50, 'tf', 0.05, 'tb', 4);
end
%%
color1 = [0 0.8 0.8];
color2 = [1 0 1];
edges = -0.2:0.1:1;
%% compare coeff change in different groups
figure;
subplot(2,4,1); hold on;
histogram([allCorr.hIn], edges, 'FaceColor', color1);
histogram([allCorr.hPre], edges, 'FaceColor', color2);
legend({'in Trial', 'pre Trial'})
title('all')

subplot(2,4,2); hold on;
histogram([allCorr.hAll], edges, 'FaceColor', color1);
histogram([allCorr.hPreSess], edges, 'FaceColor', color2);
legend({'in session', 'pre Session'})
title('all')

subplot(2,4,3); hold on;
histogram([allCorr(unit1group == 1 & unit2group == 1).hIn], edges, 'FaceColor', color1);
histogram([allCorr(unit1group == 1 & unit2group == 1).hPre], edges, 'FaceColor', color2);
legend({'in trial', 'pre trial'})
title('na')

subplot(2,4,4); hold on;
histogram([allCorr(unit1group == 1 & unit2group == 1).hAll], edges, 'FaceColor', color1);
histogram([allCorr(unit1group == 1 & unit2group == 1).hPreSess], edges, 'FaceColor', color2);
legend({'in session', 'pre Session'})
title('na')

subplot(2,4,5); hold on;
histogram([allCorr(unit1group == 2 & unit2group == 2).hIn], edges, 'FaceColor', color1);
histogram([allCorr(unit1group == 2 & unit2group == 2).hPre], edges, 'FaceColor', color2);
legend({'in trial', 'pre trial'})
title('neg')

subplot(2,4,6); hold on;
histogram([allCorr(unit1group == 2 & unit2group == 2).hAll], edges, 'FaceColor', color1);
histogram([allCorr(unit1group == 2 & unit2group == 2).hPreSess], edges, 'FaceColor', color2);
legend({'in session', 'pre Session'})
title('neg')

subplot(2,4,7); hold on;
histogram([allCorr(unit1group ~= unit2group).hIn], edges, 'FaceColor', color1);
histogram([allCorr(unit1group ~= unit2group).hPre], edges, 'FaceColor', color2);
legend({'in trial', 'pre trial'})
title('neg-na')

subplot(2,4,8); hold on;
histogram([allCorr(unit1group ~= unit2group).hAll], edges, 'FaceColor', color1);
histogram([allCorr(unit1group ~= unit2group).hPreSess], edges, 'FaceColor', color2);
legend({'in session', 'pre Session'})
title('neg-na')

sgtitle('compare coeff change')
%% compare coeff diff between different groups
colors = cool(2);    
figure;
subplot(1,4,1); hold on;
histogram([allCorr(unit1group == 1 & unit2group == 1).hIn], edges, 'FaceColor', colors(1,:), 'Normalization', 'probability');
xlabel('corr coef')
title('in trial')
ylabel('na')

subplot(1,4,2); hold on;
histogram([allCorr(unit1group == 1 & unit2group == 1).hPre], edges, 'FaceColor', colors(1,:), 'Normalization', 'probability');
xlabel('corr coef')
title('pre trial')


subplot(1,4,3); hold on;
histogram([allCorr(unit1group == 1 & unit2group == 1).hAll], edges, 'FaceColor', colors(1,:), 'Normalization', 'probability');
xlabel('corr coef')
title('in session')

subplot(1,4,4); hold on;
histogram([allCorr(unit1group == 1 & unit2group == 1).hPreSess], edges, 'FaceColor', colors(1,:), 'Normalization', 'probability');
xlabel('corr coef')
title('pre session')

subplot(1,4,1); hold on;
histogram([allCorr(unit1group == 2 & unit2group == 2).hIn], edges, 'FaceColor', colors(2,:), 'Normalization', 'probability');
title('in trial')
legend({'narrow', 'wide'})
ylabel('')
set(gca, 'TickDir', 'out')

subplot(1,4,2); hold on;
histogram([allCorr(unit1group == 2 & unit2group == 2).hPre], edges, 'FaceColor', colors(2,:), 'Normalization', 'probability');
title('pre trial')
set(gca, 'TickDir', 'out')

subplot(1,4,3); hold on;
histogram([allCorr(unit1group == 2 & unit2group == 2).hAll], edges, 'FaceColor', colors(2,:), 'Normalization', 'probability');
title('in session')
set(gca, 'TickDir', 'out')

subplot(1,4,4); hold on;
histogram([allCorr(unit1group == 2 & unit2group == 2).hPreSess], edges, 'FaceColor', colors(2,:), 'Normalization', 'probability');
title('pre session')
set(gca, 'TickDir', 'out')
sgtitle('compare coeff diff')
%%
figure2;


subplot(1,4,1); hold on;
histogram([allCorr(unit1group~=unit2group & (unit1group+unit2group)==3).hIn], edges, 'FaceColor', colors(1,:), 'Normalization', 'probability');
histogram([allCorr(unit1group==unit2group).hIn], edges, 'FaceColor', colors(2,:), 'Normalization', 'probability');
legend('diff', 'same')
set(gca, 'TickDir', 'out')
title('in trial')

subplot(1,4,2); hold on;
histogram([allCorr(unit1group~=unit2group & (unit1group+unit2group)==3).hPre], edges, 'FaceColor', colors(1,:), 'Normalization', 'probability');
histogram([allCorr(unit1group==unit2group).hPre], edges,  'FaceColor', colors(2,:), 'Normalization', 'probability');
set(gca, 'TickDir', 'out')
title('pre trial')


subplot(1,4,3); hold on;
histogram([allCorr(unit1group~=unit2group & (unit1group+unit2group)==3).hAll], edges, 'FaceColor', colors(1,:), 'Normalization', 'probability');
histogram([allCorr(unit1group==unit2group).hAll], edges, 'FaceColor', colors(2,:), 'Normalization', 'probability');
set(gca, 'TickDir', 'out')
title('in session')

subplot(1,4,4); hold on;
histogram([allCorr(unit1group~=unit2group & (unit1group+unit2group)==3).hPreSess], edges, 'FaceColor', colors(1,:), 'Normalization', 'probability');
histogram([allCorr(unit1group==unit2group).hPreSess], edges, 'FaceColor', colors(2,:), 'Normalization', 'probability');
set(gca, 'TickDir', 'out')
title('pre session')

sgtitle('compare coeff diff')
%% crossCorr
clear allCrossCorr;
[root, sep] = currComputer(); 
savePath = ['F:\allUnits\crossCorr\'];
for i = 1:length(sessionPair)
    tmp = unitAutoCorr(sessionPair{i}, unit1Pair{i}, unit2Pair{i}, 'lag', 10, 'binSize', 1, 'jitter', 20, 'numShuffle', 2000);
    allCrossCorr(i) = tmp; 
    if unit1group(i) == 1 && unit2group(i) == 1
        subFolder = 'narrow';
    else
        if unit1group(i) == 2 && unit2group(i) == 2
            subFolder = 'wide';
        else
            subFolder = 'cross';
        end
    end
    
    saveFigurePDF(gcf, [savePath sep subFolder sep sessionPair{i} '_' unit1Pair{i} '_' unit2Pair{i} '.pdf']);
end
%% pupil correlation
savePath = 'F:\allUnits\pupilCorr\';
clear pupilCorr;
parfor i = 1:length(allSessions)
    pupilCorr(i) = unitCorrPupil(allSessions{i}, allUnits{i}, 'lag', 3000, 'binSize', 200, 'binSizePost', 800, 'tf', 0.8, 'binSizePre', 200, 'plotFlag', 0, 'saveFlag', 0);
%     saveFigurePDF(gcf, [savePath sessionListGood{i} unitListGood{i} '_pupilCorr.pdf']);
end
%%
% delay in allTime
delayAllSig = NaN(1,length(allSessions));
delayAll = NaN(1,length(allSessions));
allMaxSig = NaN(1,length(allSessions));
allMax = NaN(1,length(allSessions));
for i = 1:length(allSessions)
    if ~isempty(pupilCorr(i).allCoeffs)
        coeffTmp = pupilCorr(i).allCoeffs; 
        pTmp = pupilCorr(i).allPs;
        coeffTmp = coeffTmp(3:end);
        pTmp = pTmp(3:end);
        coeffTmpSig = coeffTmp;
        coeffTmpSig(pTmp>=0.05) = NaN;
        [~, delayAllSig(i)] = max(abs(coeffTmpSig));  
        allMaxSig(i) = coeffTmp(delayAllSig(i));
        [~, delayAll(i)] = max(abs(coeffTmp));  
        allMax(i) = coeffTmp(delayAll(i));
    end
end

delayAll = delayAll + 2;
delayAllSig = delayAllSig + 2;
% delay in baseline

delayBlSig = NaN(1,length(allSessions));
delayBl = NaN(1,length(allSessions));
blMaxSig = NaN(1,length(allSessions));
blMax = NaN(1,length(allSessions));
for i = 1:length(allSessions)
    if ~isempty(pupilCorr(i).allCoeffs)
        coeffTmp = pupilCorr(i).baselineCoeffs; 
        pTmp = pupilCorr(i).baselinePs;
        coeffTmp = coeffTmp(4:end);
        pTmp = pTmp(4:end);
        coeffTmpSig = coeffTmp;
        coeffTmpSig(pTmp>=0.05) = NaN;
        [~, delayBlSig(i)] = max(abs(coeffTmpSig));  
        blMaxSig(i) = coeffTmp(delayBlSig(i));
        [~, delayBl(i)] = max(abs(coeffTmp));  
        blMax(i) = coeffTmp(delayBl(i));
    end
end

delayBl = delayBl + 3;
delayBlSig = delayBlSig + 3;
% delay in evoke

delayEvokeSig = NaN(1,length(allSessions));
delayEvoke = NaN(1,length(allSessions));
evokeMaxSig = NaN(1,length(allSessions));
evokeMax = NaN(1,length(allSessions));
for i = 1:length(allSessions)
    if ~isempty(pupilCorr(i).allCoeffs)
        coeffTmp = pupilCorr(i).evokeCoeffs; 
        pTmp = pupilCorr(i).evokePs;
        coeffTmp = coeffTmp(2:end);
        pTmp = pTmp(2:end);
        coeffTmpSig = coeffTmp;
        coeffTmpSig(pTmp>=0.05) = NaN;
        [~, delayEvokeSig(i)] = max(abs(coeffTmpSig));  
        evokeMaxSig(i) = coeffTmp(delayEvokeSig(i));
        [~, delayEvoke(i)] = max(abs(coeffTmp));  
        evokeMax(i) = coeffTmp(delayEvoke(i));
    end
end

delayEvoke = delayEvoke + 1;
delayEvokeSig = delayEvokeSig + 1;
%%
figure2;
subplot(3,1,1); hold on;
histogram(200*delayAll, 'FaceColor', [0.7 0.7 0.7], 'Normalization', 'probability');
histogram(200*delayAllSig, 'FaceColor', 'c', 'Normalization', 'probability');
title('across session')
subplot(3,1,2); hold on;
histogram(200*delayBl, 'FaceColor', [0.7 0.7 0.7], 'Normalization', 'probability');
histogram(200*delayBlSig, 'FaceColor', 'c', 'Normalization', 'probability');
title('baseline')
subplot(3,1,3); hold on;
histogram(200*delayEvoke, 'FaceColor', [0.7 0.7 0.7], 'Normalization', 'probability');
histogram(200*delayEvokeSig, 'FaceColor', 'c', 'Normalization', 'probability');
title('evoke')
sgtitle('delay')
%%
figure2;
subplot(3,1,1); hold on;
histogram(allMax, 'FaceColor', [0.7 0.7 0.7], 'Normalization', 'probability');
histogram(allMaxSig, 'FaceColor', 'c', 'Normalization', 'probability');
title('across session')
subplot(3,1,2); hold on;
histogram(blMax, 'FaceColor', [0.7 0.7 0.7], 'Normalization', 'probability');
histogram(blMaxSig, 'FaceColor', 'c', 'Normalization', 'probability');
title('baseline')
subplot(3,1,3); hold on;
histogram(evokeMax, 'FaceColor', [0.7 0.7 0.7]);
histogram(evokeMaxSig, 'FaceColor', 'c');
title('evoke')
sgtitle('max Coeff')
%%
figure2;
subplot(2,3,1); hold on;
edges = linspace(min(delayAll), max(delayAll)+0.001, 15)*200/1000;
histogram(200*delayAll(ind==1)/1000, edges, 'FaceColor', color1, 'Normalization', 'probability');
histogram(200*delayAll(ind==2)/1000, edges, 'FaceColor', color2, 'Normalization', 'probability');
legend({'narrow','wide'})
title('delay max corr allTime')

subplot(2,3,4); hold on;
edges = linspace(min(allMax), max(allMax)+0.001, 15);
histogram(allMax(ind==1), edges, 'FaceColor', color1, 'Normalization', 'probability');
histogram(allMax(ind==2), edges, 'FaceColor', color2, 'Normalization', 'probability');
legend({'narrow','wide'})
title('coeff max corr allTime')

subplot(2,3,2); hold on;
edges = linspace(min(delayBl), max(delayBl)+0.001, 15)*200/1000;
histogram(200*delayBl(ind==1)/1000, edges, 'FaceColor', color1, 'Normalization', 'probability');
histogram(200*delayBl(ind==2)/1000, edges, 'FaceColor', color2, 'Normalization', 'probability');
legend({'narrow','wide'})
title('delay max corr baseline')

subplot(2,3,5); hold on;
edges = linspace(min(blMax), max(blMax)+0.001, 15);
histogram(blMax(ind==1), edges, 'FaceColor', color1, 'Normalization', 'probability');
histogram(blMax(ind==2), edges, 'FaceColor', color2, 'Normalization', 'probability');
legend({'narrow','wide'})
title('coeff max corr baseline')


subplot(2,3,3); hold on;
edges = linspace(min(delayEvoke), max(delayEvoke)+0.001, 15)*200/1000;
histogram(200*delayEvoke(ind==1)/1000, edges, 'FaceColor', color1, 'Normalization', 'probability');
histogram(200*delayEvoke(ind==2)/1000, edges, 'FaceColor', color2, 'Normalization', 'probability');
legend({'narrow','wide'})
title('delay max corr evoke')

subplot(2,3,6); hold on;
edges = linspace(min(evokeMax), max(evokeMax)+0.001, 15);
histogram(evokeMax(ind==1), edges, 'FaceColor', color1, 'Normalization', 'probability');
histogram(evokeMax(ind==2), edges, 'FaceColor', color2, 'Normalization', 'probability');
legend({'narrow','wide'})
title('coeff max corr evoked')
%%