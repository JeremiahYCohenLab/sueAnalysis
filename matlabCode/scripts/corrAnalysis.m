% find pairs
% load F:\tmpData\catWithQchosen.mat
load F:\tmpData\allUnitAUC.mat
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
%     allCorr(i) = unitCorr(sessionPair{i}, unit1Pair{i}, unit2Pair{i}, 'binSize', 200, 'binSizePost', 50, 'tf', 0.05, 'tb', 4);
    allCorr(i) = unitCorr(sessionPair{i}, unit1Pair{i}, unit2Pair{i}, 'binSize', 200, 'binSizePost', 50, 'tf', 0.05, 'tb', 2);
end
%%
color1 = [0 0.8 0.8]; 
color2 = [1 0 1]; 
edges = -0.2:0.05:1; 
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
title('II')

subplot(2,4,4); hold on;
histogram([allCorr(unit1group == 1 & unit2group == 1).hAll], edges, 'FaceColor', color1);
histogram([allCorr(unit1group == 1 & unit2group == 1).hPreSess], edges, 'FaceColor', color2);
legend({'in session', 'pre Session'})
title('II')

subplot(2,4,5); hold on;
histogram([allCorr(unit1group == 2 & unit2group == 2).hIn], edges, 'FaceColor', color1);
histogram([allCorr(unit1group == 2 & unit2group == 2).hPre], edges, 'FaceColor', color2);
legend({'in trial', 'pre trial'})
title('I')

subplot(2,4,6); hold on;
histogram([allCorr(unit1group == 2 & unit2group == 2).hAll], edges, 'FaceColor', color1);
histogram([allCorr(unit1group == 2 & unit2group == 2).hPreSess], edges, 'FaceColor', color2);
legend({'in session', 'pre Session'})
title('I')

subplot(2,4,7); hold on;
histogram([allCorr(unit1group ~= unit2group).hIn], edges, 'FaceColor', color1);
histogram([allCorr(unit1group ~= unit2group).hPre], edges, 'FaceColor', color2);
legend({'in trial', 'pre trial'})
title('I-II')

subplot(2,4,8); hold on;
histogram([allCorr(unit1group ~= unit2group).hAll], edges, 'FaceColor', color1);
histogram([allCorr(unit1group ~= unit2group).hPreSess], edges, 'FaceColor', color2);
legend({'in session', 'pre Session'})
title('I-II')

sgtitle('compare coeff change')
%% compare coeff diff between different groups
colors = cool(2);    
figure;
subplot(1,4,1); hold on;
histogram([allCorr(unit1group == 1 & unit2group == 1).hIn], edges, 'FaceColor', colors(1,:), 'Normalization', 'probability');
histogram([allCorr(unit1group == 2 & unit2group == 2).hIn], edges, 'FaceColor', colors(2,:), 'Normalization', 'probability');
legend({'narrow', 'wide'})
ylabel('')
set(gca, 'TickDir', 'out')
xlabel('corr coef')
title('in trial')
ylabel('na')

subplot(1,4,2); hold on;
histogram([allCorr(unit1group == 1 & unit2group == 1).hPre], edges, 'FaceColor', colors(1,:), 'Normalization', 'probability');
histogram([allCorr(unit1group == 2 & unit2group == 2).hPre], edges, 'FaceColor', colors(2,:), 'Normalization', 'probability');
title('pre trial')
set(gca, 'TickDir', 'out')
xlabel('corr coef')
title('pre trial')


subplot(1,4,3); hold on;
histogram([allCorr(unit1group == 1 & unit2group == 1).hAll], edges, 'FaceColor', colors(1,:), 'Normalization', 'probability');
histogram([allCorr(unit1group == 2 & unit2group == 2).hAll], edges, 'FaceColor', colors(2,:), 'Normalization', 'probability');
title('in session')
set(gca, 'TickDir', 'out')
xlabel('corr coef')
title('in session')

subplot(1,4,4); hold on;
histogram([allCorr(unit1group == 1 & unit2group == 1).hPreSess], edges, 'FaceColor', colors(1,:), 'Normalization', 'probability');
histogram([allCorr(unit1group == 2 & unit2group == 2).hPreSess], edges, 'FaceColor', colors(2,:), 'Normalization', 'probability');
title('pre session')
set(gca, 'TickDir', 'out')
sgtitle('compare coeff diff')
xlabel('corr coef')
title('pre session')



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
%%
figure2;
subplot(1,3,1);
currCoeffs = [allCorr.hPreCue];
edges = linspace(min(currCoeffs)-0.001, max(currCoeffs)+0.001, 10);
histogram(currCoeffs, edges, 'FaceColor', [0.7 0.7 0.7], 'EdgeColor', 'none', 'Normalization', 'probability');

subplot(1,3,2); hold on;
currCoeffs = [allCorr.hPreCue];
edges = linspace(min(currCoeffs)-0.001, max(currCoeffs)+0.001, 10);
histogram(currCoeffs(unit1group~=unit2group), edges, 'FaceColor', color1, 'EdgeColor', 'none', 'Normalization', 'probability');
histogram(currCoeffs(unit1group==unit2group), edges, 'FaceColor', color2, 'EdgeColor', 'none', 'Normalization', 'probability');
legend({'diff',  'same'});

subplot(1,3,3); hold on;
currCoeffs = [allCorr.hPreCue];
edges = linspace(min(currCoeffs)-0.001, max(currCoeffs)+0.001, 10);
histogram(currCoeffs(unit1group==1 & unit2group==1), edges, 'FaceColor', color1, 'EdgeColor', 'none', 'Normalization', 'probability');
histogram(currCoeffs(unit1group==2 & unit2group==2), edges, 'FaceColor', color2, 'EdgeColor', 'none', 'Normalization', 'probability');
legend({'II',  'I'});
%%
% plot traces to show correlation
tb = 2;
tf = 0;
stepSize = 100;
binSize = 100;
myKernel = ones(1,10);
myKernel = myKernel./sum(myKernel);
myKernelShort = ones(1,5);
myKernelShort = myKernelShort./sum(myKernelShort);
numBins = 6;
numBinsRaw = 4;
for i = 1:length(sessionPair)
    [~, matCue] = getUnitMatCue(sessionPair{i}, unit1Pair{i}, tb, tf, stepSize, binSize);
    unit1PreCue = sum(matCue, 2)/tb;
    [~, matCue] = getUnitMatCue(sessionPair{i}, unit2Pair{i}, tb, tf, stepSize, binSize);
    unit2PreCue = sum(matCue, 2)/tb;   
    unit1Smoothed = conv(unit1PreCue, myKernel);
    unit1Smoothed = unit1Smoothed(round(0.5*length(myKernel)):end-round(0.5*length(myKernel)));
    unit2Smoothed = conv(unit2PreCue, myKernel);
    unit2Smoothed = unit2Smoothed(round(0.5*length(myKernel)):end-round(0.5*length(myKernel)));
    unit1SmoothedShort = conv(unit1PreCue, myKernelShort);
    unit1SmoothedShort = unit1SmoothedShort(round(0.5*length(myKernelShort)):end-round(0.5*length(myKernelShort)));
    unit2SmoothedShort = conv(unit2PreCue, myKernelShort);
    unit2SmoothedShort = unit2SmoothedShort(round(0.5*length(myKernelShort)):end-round(0.5*length(myKernelShort)));
    figure; 
    subplot(3,1,1);hold on;
    plot(unit1PreCue, 'LineWidth', 1, 'Color', color1);
    plot(unit2PreCue, 'LineWidth', 1, 'Color', color2);
    legend({num2str(3-unit1group(i)), num2str(3-unit2group(i))})
    
    subplot(3,1,2);hold on;
    plot(unit1Smoothed, 'LineWidth', 1, 'Color', color1);
    plot(unit2Smoothed, 'LineWidth', 1, 'Color', color2);
    title('long kernel')
    
    subplot(3,1,3);hold on;
    plot(unit1SmoothedShort, 'LineWidth', 1, 'Color', color1);
    plot(unit2SmoothedShort, 'LineWidth', 1, 'Color', color2);
    title('short kernel')
    
    sgtitle([num2str(i) '' sessionPair{i} unit1Pair{i} unit2Pair{i} ' ' num2str(allCorr(i).hPreCue)], 'Interpreter', 'none');
    screen = get(0,'Screensize');
    screen(4) = screen(4) - 100;
    set(gcf, 'Position', screen)
    
    
    edges = linspace(min(unit1Smoothed)-0.01, max(unit1Smoothed)+0.01, numBins+1);
    edges = quantile(unit1Smoothed, linspace(0, 1, numBins+1));
    edges(1) = edges(1) - 0.001;
    edges(end) = edges(end) + 0.001;
    mean1 = zeros(1, numBins);
    sem1 = zeros(1, numBins);
    mean2 = zeros(1, numBins);
    sem2 = zeros(1, numBins);
    for j = 1:numBins
        mean1(j) = mean(unit1Smoothed(unit1Smoothed>=edges(j) & unit1Smoothed<edges(j+1)));
        sem1(j) = sem(unit1Smoothed(unit1Smoothed>=edges(j) & unit1Smoothed<edges(j+1)));        
        mean2(j) = mean(unit2Smoothed(unit1Smoothed>=edges(j) & unit1Smoothed<edges(j+1)));
        sem2(j) = sem(unit2Smoothed(unit1Smoothed>=edges(j) & unit1Smoothed<edges(j+1)));
    end
    
    figure2Wide;
    subplot(1,2,1);
    hold on;
    scatter(unit1Smoothed, unit2Smoothed, 10, [0.4 0.4 0.4], 'filled','MarkerEdgeColor', 'none', 'MarkerFaceAlpha', 0.25);
    errorbar(mean1, mean2, sem2, sem2, sem1, sem1, 'LineWidth', 2, 'Color', [0 0.8 0.8]); 
    xlabel(num2str(3-unit1group(i)));
    ylabel(num2str(3-unit2group(i)));
    title('long kernel');
    
    edges = quantile(unit1PreCue, linspace(0, 1, numBinsRaw+1));
    edges(1) = edges(1) - 0.001;
    edges(end) = edges(end) + 0.001;
    mean1 = zeros(1, numBinsRaw);
    sem1 = zeros(1, numBinsRaw);
    mean2 = zeros(1, numBinsRaw);
    sem2 = zeros(1, numBinsRaw);
    for j = 1:numBinsRaw
        mean1(j) = mean(unit1PreCue(unit1PreCue>=edges(j) & unit1PreCue<edges(j+1)));
        sem1(j) = sem(unit1PreCue(unit1PreCue>=edges(j) & unit1PreCue<edges(j+1)));        
        mean2(j) = mean(unit2PreCue(unit1PreCue>=edges(j) & unit1PreCue<edges(j+1)));
        sem2(j) = sem(unit2PreCue(unit1PreCue>=edges(j) & unit1PreCue<edges(j+1)));
    end
    subplot(1,2,2);
    hold on;
    edgeX = linspace(min(unit1PreCue)-0.01, max(unit1PreCue)+0.01, length(unique(unit1PreCue))+1);
    edgeY = linspace(min(unit2PreCue)-0.01, max(unit2PreCue)+0.01, length(unique(unit2PreCue))+1);
   
%     N = histcounts2(unit1PreCue,unit2PreCue, edgeX, edgeY);
%     imagesc(x, y);
    scatter(unit1PreCue, unit2PreCue, 15, [0.4 0.4 0.4], 'filled','MarkerEdgeColor', 'none', 'MarkerFaceAlpha', 0.10);
    errorbar(mean1, mean2, sem2, sem2, sem1, sem1, 'LineWidth', 2, 'Color', [0 0.8 0.8]); 
    xlabel(num2str(3-unit1group(i)));
    ylabel(num2str(3-unit2group(i)));
    title('raw')
    sgtitle([num2str(i) '' sessionPair{i} unit1Pair{i} unit2Pair{i} ' ' num2str(allCorr(i).hPreCue)], 'Interpreter', 'none');
end

%%

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
% pupil align
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
%%
load F:\tmpData\allUnitAUC.mat
savePath = 'F:\allUnits\pupilCorr\';
clear pupilCorr;
for i = 1:length(allSessions)
    pupilCorr(i) = unitCorrPupil(allSessions{i}, allUnits{i}, 'lag', 3000, 'binSize', 200, 'binSizePost', 500, 'tf', 0.5, 'binSizePre', 200, 'plotFlag', 1, 'saveFlag', 0);
%     sgtitle([allSessions{i} ' ' allUnits{i} ' ' num2str(ind(i))], 'Interpreter', 'none');
%     saveFigurePDF(gcf, [savePath allSessions{i} allUnits{i} '_pupilCorr.pdf']);
%     close gcf
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
color1 = [0.2 0.2 0.2];
color2 = [0.7 0.7 0.7];
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

%% baseline
figure2; hold on;
ecdf([pupilCorr(ind==1).baselineCorr], 'Bounds','on');
ecdf([pupilCorr(ind==2).baselineCorr], 'Bounds','on');
legend({'II', '', '', 'I'})
title('baselineCorrTrialbyTrial')
xlabel('corrCoeff')
%%
figure2; hold on;
edges = linspace(min([pupilCorr.baselineCorr]), max([pupilCorr.baselineCorr]), 20);
histogram([pupilCorr(ind==1).baselineCorr], edges, 'FaceColor', color1, 'EdgeColor', 'none');
histogram([pupilCorr(ind==2).baselineCorr], edges, 'FaceColor', color2, 'EdgeColor', 'none');
legend({'II', 'I'})
title('baselineCorrTrialbyTrial')
xlabel('corrCoeff')
[h, p, ~, stats] = ttest2([pupilCorr(ind==1).baselineCorr], [pupilCorr(ind==2).baselineCorr]);
title(['p=' num2str(p) 'tstat=' num2str(stats.tstat), 'df=' num2str(stats.df)])
[h, p, ~, stats] = ttest([pupilCorr.baselineCorr]);
ylabel(['p=' num2str(p) 'tstat=' num2str(stats.tstat), 'df=' num2str(stats.df)])
%%
maxCorr = NaN(size(allSessions));
maxLag = NaN(size(allSessions));
minCorr = NaN(size(allSessions));
for i = 1:length(allSessions)
    if ~isempty(pupilCorr(i).allCoeffs)
        [maxCorr(i), maxLag(i)] = max([pupilCorr(i).allCoeffs]);
        minCorr(i) = pupilCorr(i).allLm(1);
    end
end
maxLag = 0.2*maxLag;
figure2Wide; 
subplot(1,3,1); hold on;
edges = linspace(min(maxCorr), max(maxCorr), 15);
histogram(maxCorr(ind==1), edges, 'FaceColor', color1, 'EdgeColor', 'none');
histogram(maxCorr(ind==2), edges, 'FaceColor', color2, 'EdgeColor', 'none');
legend({'II', 'I'})
title('allCorrWholeSession')
xlabel('maxCorrCoeff')
[h, p, ~, stats] = ttest(maxCorr);
ylabel(['p=' num2str(p) 'tstat=' num2str(stats.tstat), 'df=' num2str(stats.df)])
[h, p, ~, stats] = ttest2(maxCorr(ind==1),maxCorr(ind==2));
title(['p=' num2str(p) 'tstat=' num2str(stats.tstat), 'df=' num2str(stats.df)])

subplot(1,3,2); hold on;
edges = linspace(min(maxLag), max(maxLag), 10);
histogram(maxLag(ind==1), edges, 'FaceColor', color1, 'EdgeColor', 'none');
histogram(maxLag(ind==2), edges, 'FaceColor', color2, 'EdgeColor', 'none');
legend({'II', 'I'})
title('allCorrWholeSession')
xlabel('maxLag (s)')
[h, p, ~, stats] = ttest2(maxLag(ind==1),maxLag(ind==2));
title(['type I vs type II' 'p=' num2str(p) 'tstat=' num2str(stats.tstat), 'df=' num2str(stats.df)])


subplot(1,3,3); hold on;
edges = linspace(min(minCorr), max(minCorr), 10);
histogram(minCorr(ind==1), edges, 'FaceColor', color1, 'EdgeColor', 'none');
histogram(minCorr(ind==2), edges, 'FaceColor', color2, 'EdgeColor', 'none');
legend({'II', 'I'})
title('allCorrWholeSession')
xlabel('minCorrCoeff')
[h, p, ~, stats] = ttest2(minCorr(ind==1),minCorr(ind==2));
title(['type I vs type II' 'p=' num2str(p) 'tstat=' num2str(stats.tstat), 'df=' num2str(stats.df)])

%%

maxCorr = NaN(size(allSessions));
maxLag = NaN(size(allSessions));

for i = 1:length(allSessions)
    if ~isempty(pupilCorr(i).allCoeffs)
        [maxCorr(i), maxLag(i)] = max([pupilCorr(i).evokeCoeffs]);
%         maxCorr(i) = pupilCorr(i).evokeCoeffs(4);
%         minCorr(i) = pupilCorr(i).allLm(1);
    end
end
maxLag = 0.5*maxLag;
figure2Wide; 
subplot(1,2,1); hold on;
edges = linspace(min(maxCorr), max(maxCorr), 13);
histogram(maxCorr(ind==1), edges, 'FaceColor', color1, 'EdgeColor', 'none');
histogram(maxCorr(ind==2), edges, 'FaceColor', color2, 'EdgeColor', 'none');
legend({'II', 'I'})
title('evoked')

xlabel('maxCorrCoeff')
[h, p, ~, stats] = ttest2(maxCorr(ind==1),maxCorr(ind==2));
title(['p=' num2str(p) 'tstat=' num2str(stats.tstat), 'df=' num2str(stats.df)])
[h, p, ~, stats] = ttest(maxCorr);
ylabel(['p=' num2str(p) 'tstat=' num2str(stats.tstat), 'df=' num2str(stats.df)])



subplot(1,2,2); hold on;
edges = linspace(min(maxLag), max(maxLag), 9);
histogram(maxLag(ind==1), edges, 'FaceColor', color1, 'EdgeColor', 'none');
histogram(maxLag(ind==2), edges, 'FaceColor', color2, 'EdgeColor', 'none');
legend({'II', 'I'})
title('evoked')
xlabel('maxLag (s)')
[p, h, stats] = ranksum(maxLag(ind==1),maxLag(ind==2));
title(['p=' num2str(p)])
%%