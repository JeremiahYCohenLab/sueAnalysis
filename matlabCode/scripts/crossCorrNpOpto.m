dayList = getDayList('npOpto', 'optoTagging', 'withUnits');
waveformRange = {[1:7], [1:60]}; % 2:6, 5:15
%% autocorreltion combine session by session data
root = 'F:\npOptoRecordings\withOpto\';
allWFs = [];
loc = [];
allAutoCorr = [];
allMetrics = [];
allSessions = [];
allUnits = [];
frSpont = [];
% lag = 30000; % in ms
% binSize = 200; % in ms
lag = 30; % in ms
binSize = 2; % in ms
for i = 1:length(dayList)
    % load data
    session = dayList{i};
    wfPath = [root session '\sorted\wfLoc.mat'];
    load(wfPath, 'maxChannelOpto','waveforms2DAligned', 'waveformsAligned', 'waveformsOpto', 'allUnitsOpto', 'optoMetrics');
    % load waveforms
    currWFLin = waveforms2DAligned(:, waveformRange{1}, waveformRange{2});
    currWFLin = permute(currWFLin, [3,2,1]);
    baseline = mean(currWFLin(1:3,:,:), 1);
    currWFLin = currWFLin - baseline;
    currWFLin = reshape(currWFLin, [], size(waveforms2DAligned, 1))'; 
    currWFLin = currWFLin./max(squeeze(currWFLin(:, 181:240)), [], 2);
    allWFs = [allWFs; currWFLin];
    loc = [loc; [maxChannelOpto - mean(maxChannelOpto)]/(max(maxChannelOpto) - min(maxChannelOpto))];
    % auto correlation 
    for j = 1:length(allUnitsOpto)
        s = unitAutoCorrNp(session, allUnitsOpto{j}, allUnitsOpto{j}, 'lag', lag, 'binSize', binSize, 'plot', 0);
        allAutoCorr = [allAutoCorr; s.corrCoAll];
        frSpont = [frSpont, s.fr];
    end

    allMetrics = [allMetrics; optoMetrics];
    
    currSession = cell(size(allUnitsOpto, 1), 1);
    currSession(:) = {session};
    allSessions = [allSessions; currSession];

    allUnits = [allUnits; allUnitsOpto];
   
end
%%
load('F:\tmpData\npOptoWFC4.mat', 'ind', 'score')
load('F:\tmpData\npOptoCrossCorrLong.mat', 'noDrift')
%% threshold
allAutoCorr = allAutoCorr(:,1:(lag/binSize));
% longCorrThresh = 0.1;
% shortCorrThresh = 0.4;
% noDrift = mean(allAutoCorr(:, 1:2), 2) < longCorrThresh & mean(allAutoCorr(:, end-1:end), 2) < shortCorrThresh;
% %% normalize


%%
allAutoCorrFocus = allAutoCorr(:, end-10:end) - mean(allAutoCorr(:,1:100), 2);
noDrift = noDrift & allAutoCorrFocus(:,1)<0.1;
allAutoCorrFocusZ = zscore(allAutoCorrFocus(noDrift,:), 0 , 1);
allAutoCorrFocusZAll = NaN(size(allAutoCorrFocus));
allAutoCorrFocusZAll(noDrift, :) = allAutoCorrFocusZ;
[coeff,score,latent, ~, explained, mu] = pca(allAutoCorrFocus);
%%
figure2;
gscatter(score(:,1), score(:,2), ind)
%%
figure2;
for i = 1:4
    tempScore = allAutoCorrFocusZAll(:,end+1-i);
    subplot(4, 1, i);
    edges = linspace(min(tempScore), max(tempScore), 25);
    hold on;
    for j = 1:4
        histogram(tempScore(ind==j & noDrift), edges)
    end
    legend

end
%%
figure2;
for i = 1:8
    tempScore = score(:, i);
    subplot(8, 1, i);
    edges = linspace(min(tempScore(noDrift)), max(tempScore(noDrift)), 25);
    hold on;
        histogram(tempScore((ind==1) & noDrift), edges, 'Normalization', 'probability')
        histogram(tempScore((ind~=1) & noDrift), edges, 'Normalization', 'probability')

    legend

end
%%
[minCorr, minInd] = min(allAutoCorrFocus, [], 2);
figure2;
gscatter(minCorr(noDrift), minInd(noDrift), ind(noDrift));



%% crosscorreltion combine session by session data
root = 'F:\npOptoRecordings\withOpto\';
allCrossCorr = [];
% lag = 30000; % in ms
% binSize = 200; % in ms
lag = 20000; % in ms
binSize = 200; % in ms
allCmbs = [];
allCI = [];
allSessionsCmb = [];
allIndsCmb = [];
allSurrMean = [];
for i = 1:length(dayList)
    % load data
    session = dayList{i};
    wfPath = [root session '\sorted\wfLoc.mat'];
    load(wfPath, 'maxChannelOpto','waveforms2DAligned', 'waveformsAligned', 'waveformsOpto', 'allUnitsOpto', 'optoMetrics');
    
    C = nchoosek(1:length(allUnitsOpto),2); 
    unitC = cell(size(C));
    unitC(:,1) = allUnitsOpto(C(:,1));
    unitC(:,2) = allUnitsOpto(C(:,2));
    currInds = ind(cellfun(@(x) strcmp(x, session), allSessions));
    currIndCmbs = zeros(size(C,1), 2);
    currIndCmbs(:,1) = currInds(C(:,1));
    currIndCmbs(:,2) = currInds(C(:,2));
    % cross correlation 
    for j = 1:size(C,1)
        s = unitAutoCorrNp(session, unitC{j,1}, unitC{j,2}, 'lag', lag, 'binSize', binSize, 'plot', 0);
        allCrossCorr = [allCrossCorr; s.corrCoAll];
        allCI = cat(3, allCI, s.CI);
        allSurrMean = [allSurrMean; s.mean];
    end
    allCmbs = [allCmbs; unitC];
    currSession = cell(size(C, 1), 1);
    currSession(:) = {session};
    allSessionsCmb = [allSessionsCmb; currSession];
    allIndsCmb = [allIndsCmb; currIndCmbs];
   
end
%%
unitName1 = strcat(allSessionsCmb, allCmbs(:,1));
unitName2 = strcat(allSessionsCmb, allCmbs(:,2));
driftNames = strcat(allSessions(~noDrift), allUnits(~noDrift));

driftCmbs = zeros(length(allSessionsCmb),1);
for i = 1:length(driftNames)
    tempInd = cellfun(@(x) strcmp(driftNames{i}, x), unitName1) | cellfun(@(x) strcmp(driftNames{i}, x), unitName2);
    driftCmbs(tempInd) = 1;
end
%% detect significance. location and sign
deviation = zeros(size(allCI));
deviation(1, :, :) = - squeeze(allCI(1,:,:)) + allCrossCorr';
deviation(2, :, :) = squeeze(allCI(2,:,:)) - allCrossCorr';
% deviation(:,151,:) = 0;

figure2;
imagescHeat(squeeze(deviation(1,:,~driftCmbs))', 0);
title('lower')

figure2
imagescHeat(squeeze(deviation(2,:,~driftCmbs))', 0);
title('upper')


%% detect level and location of break points
breakPoints = squeeze(deviation(1,:,:))'<0 | squeeze(deviation(2,:,:))'<0;
allCrossCorrSig = allCrossCorr;
allCrossCorrSig(~breakPoints) = 0;
figure2;
imagescHeat(allCrossCorrSig(~driftCmbs,:), 0)
%% 
figure2;

for i = 1:4
    for j = 1:4
       currInd = allIndsCmb(:,1) == i & allIndsCmb(:,2) == j;
       if sum(currInd & ~driftCmbs) > 0
           subplot(4,4,(i-1)*4 + j)
           imagescHeat(allCrossCorrSig(currInd & ~driftCmbs, 151-20:151+20), 0);
           title([num2str(i) '---' num2str(j)])
       end
    end
end
%% flip by pair
figure2;

for i = 1:4
    for j = i:4
       currInd = allIndsCmb(:,1) == i & allIndsCmb(:,2) == j;
       currCorr = allCrossCorrSig(currInd & ~driftCmbs, :);
       if i ~= j
            currIndFlip = allIndsCmb(:,1) == j & allIndsCmb(:,2) == i;
            currCorrFlip = allCrossCorrSig(currIndFlip & ~driftCmbs, :);
            currCorr = [currCorr; flip(currCorrFlip, 2)];
       end


           subplot(4,4,(i-1)*4 + j)
           imagescHeat(currCorr(:, 151-20:151+20), 0);
           title([num2str(i) '---' num2str(j)])

    end
end
%% find breaking points existence for short
breakInd = squeeze(deviation(2,11-2:11+2,:)<0);
breakIndSum = sum(breakInd,1)>0;

breakRatio = zeros(4, 4);
for i = 1:4
    for j = i:4
       currInd = (allIndsCmb(:,1) == i & allIndsCmb(:,2) == j)|(allIndsCmb(:,1) == j & allIndsCmb(:,2) == i);
       currBreaks = breakIndSum(currInd & ~driftCmbs);
       breakRatio(i, j) = sum(currBreaks)/length(currBreaks);

    end
end

figure2;
imagescHeat(breakRatio, 0.5);

%% find breaking points existence for short
breakInd = squeeze(deviation(2,11-2:11+2,:)<0);
breakIndSum = sum(breakInd,1)>0;

breakRatio = zeros(4, 4);
for i = 1:4
    for j = i:4
       currInd = (allIndsCmb(:,1) == i & allIndsCmb(:,2) == j)|(allIndsCmb(:,1) == j & allIndsCmb(:,2) == i);
       currBreaks = breakIndSum(currInd & ~driftCmbs);
       breakRatio(i, j) = sum(currBreaks)/length(currBreaks);

    end
end

figure2;
imagescHeat(breakRatio, 0.001);

% type I and type II
breakInd = squeeze(deviation(2,11-2:11+2,:)<0);
breakIndSum = sum(breakInd,1)>0;

breakRatio = zeros(2, 2);

indType = allIndsCmb;
indType(allIndsCmb==2|allIndsCmb==3) = 1;
indType(allIndsCmb==1|allIndsCmb==4) = 2;

for i = 1:2
    for j = i:2
       currInd = (indType(:,1) == i & indType(:,2) == j)|(indType(:,1) == j & indType(:,2) == i);
       currBreaks = breakIndSum(currInd & ~driftCmbs);
       breakRatio(i, j) = sum(currBreaks)/length(currBreaks);

    end
end

figure2;
imagescHeat(breakRatio, 0.001);
%% find breaking points existence for long
breakInd = squeeze(deviation(2,151-2:151+2,:)<0);
breakIndSum = sum(breakInd,1)>0;

breakRatio = zeros(4, 4);
for i = 1:4
    for j = i:4
       currInd = (allIndsCmb(:,1) == i & allIndsCmb(:,2) == j)|(allIndsCmb(:,1) == j & allIndsCmb(:,2) == i);
       currBreaks = breakIndSum(currInd & ~driftCmbs);
       breakRatio(i, j) = sum(currBreaks)/length(currBreaks);

    end
end

figure2;
imagescHeat(breakRatio, 0.5);
%% correlation between groups, long
figure2;

corrGroup = cell(4,4);
for i = 1:4
    for j = i:4
       currInd = (allIndsCmb(:,1) == i & allIndsCmb(:,2) == j)|(allIndsCmb(:,1) == j & allIndsCmb(:,2) == i);
       currCorr = allCrossCorr(currInd & ~driftCmbs & ~lowCmbs, 101);
       corrGroup{i, j} = currCorr;
       subplot(4, 4, 4*(i-1)+j)
       histogram(corrGroup{i, j});
    end
end

figure2;
imagescHeat(cellfun(@mean, corrGroup), 0)
%% by specific group
figure2
hold on
temp = [corrGroup{2, 2}; corrGroup{3, 3}];
edges = linspace(min(temp)-0.001, max(temp)+0.001, 20);
histogram(corrGroup{2, 2}, edges, 'Normalization', 'probability', 'FaceColor','b', 'EdgeColor','none');
histogram(corrGroup{3, 3}, edges, 'Normalization', 'probability', 'FaceColor','r', 'EdgeColor','none');
legend({'Type I-a', 'Type I-b'})


%% correlation between types
figure2;
indType = allIndsCmb;
indType(allIndsCmb==2|allIndsCmb==3) = 1;
indType(allIndsCmb==1|allIndsCmb==4) = 2;
corrAll = allCrossCorr(:, lag/binSize + 1);
edges = linspace(min(corrAll)-0.01, max(corrAll)-0.01, 30);
hold on;
temp1 = corrAll(indType(:,1) == 1 & indType(:,2) == 1 & ~driftCmbs & ~lowCmbs);
histogram(temp1, edges, 'FaceColor', [0.3 0.3 0.3], 'EdgeColor', 'none', 'Normalization', 'probability');
temp2 = corrAll(indType(:,1) == 2 & indType(:,2) == 2 & ~driftCmbs & ~lowCmbs);
histogram(temp2, edges, 'FaceColor', [0.7 0.7 0.7], 'EdgeColor', 'none', 'Normalization', 'probability');
legend({'typeI', 'typeII'})
%% low fr units
unitName1 = strcat(allSessionsCmb, allCmbs(:,1));
unitName2 = strcat(allSessionsCmb, allCmbs(:,2));
lowNames = strcat(allSessions(frSpont(1, :)<=1.5), allUnits(frSpont(1, :)<=1.5));

lowCmbs = zeros(length(allSessionsCmb),1);
for i = 1:length(lowNames)
    tempInd = cellfun(@(x) strcmp(lowNames{i}, x), unitName1) | cellfun(@(x) strcmp(lowNames{i}, x), unitName2);
    lowCmbs(tempInd) = 1;
end



%% short
figure2;

% type I and type II
breakInd = squeeze(deviation(2,11-2:11+2,:)<0);
breakIndSum = sum(breakInd,1)>0;

breakRatio = zeros(2, 2);

indType = allIndsCmb;
indType(allIndsCmb==2|allIndsCmb==3) = 1;
indType(allIndsCmb==1|allIndsCmb==4) = 2;

for i = 1:2
    for j = i:2
       currInd = (indType(:,1) == i & indType(:,2) == j)|(indType(:,1) == j & indType(:,2) == i);
       currBreaks = breakIndSum(currInd & ~driftCmbs);
       breakRatio(i, j) = sum(currBreaks)/length(currBreaks);

    end
end


figure2;
imagescHeat(breakRatio, 0.001);
%% long area under curve by type and group

% type I and type II
breakSum = sum(allCrossCorrCorrected(:, [101-4:99, 103:101+4]), 2);

breakSumG = cell(4, 4);

indType = allIndsCmb;
indType(allIndsCmb==2|allIndsCmb==3) = 1;
indType(allIndsCmb==1|allIndsCmb==4) = 2;
figure2;
for i = 1:2
    for j = i:2
        subplot(2, 2, 2*(i-1) + j)
       currInd = (indType(:,1) == i & indType(:,2) == j)|(indType(:,1) == j & indType(:,2) == i);
       currBreaks = breakSum(currInd & ~driftCmbs & ~lowCmbs);
       histogram(currBreaks);
    end
end

% groups
figure2;
for i = 1:4
    for j = i:4
       subplot(4, 4, 4*(i-1) + j)
       currInd = (allIndsCmb(:,1) == i & allIndsCmb(:,2) == j)|(allIndsCmb(:,1) == j & allIndsCmb(:,2) == i);
       currBreaks = breakSum(currInd & ~driftCmbs & ~lowCmbs);
       histogram(currBreaks, 10)
       breakSumG{i, j} = currBreaks;
    end
end

figure2;
imagescHeat(cellfun(@mean, breakSumG), 0);
%%
figure2
hold on
temp = [breakSumG{2, 2}; breakSumG{3, 3}];
edges = linspace(min(temp)-0.001, max(temp)+0.001, 20);
histogram(breakSumG{2, 2}, edges, 'Normalization', 'probability', 'FaceColor','b', 'EdgeColor','none');
histogram(breakSumG{3, 3}, edges, 'Normalization', 'probability', 'FaceColor','r', 'EdgeColor','none');
legend({'Type I-a', 'Type I-b'})
%% plot waveforms
figure2
hold on
time = [(1:60) - 21]/30000 * 1000;
plot(time, mean(allWFs(ind==1 | ind==4, 181:240)), 'Color', [0.6 0.6 0.6], 'LineWidth', 2)
plot(time, mean(allWFs(ind==2 | ind==3, 181:240)), 'Color', [0.2 0.2 0.2], 'LineWidth', 2)
xlabel('time from peak (ms)')
%% relative sig
allCrossCorrCorrected = allCrossCorr - allSurrMean;
breakSum = sum(allCrossCorrCorrected(:, [101-5:99, 103:101+5]), 2);
figure;
indType = allIndsCmb;
indType(allIndsCmb==2|allIndsCmb==3) = 1;
indType(allIndsCmb==1|allIndsCmb==4) = 2;
edges = linspace(min(breakSum)-0.01, max(breakSum)-0.01, 30);
hold on;
temp1 = breakSum(indType(:,1)==1 & indType(:,2) == 1 & ~driftCmbs & ~lowCmbs);
histogram(temp1, edges, 'FaceColor', [0.3 0.3 0.3], 'EdgeColor', 'none', 'Normalization', 'probability');
temp2 = breakSum(indType(:,1) == 2 & indType(:,2) == 2 & ~driftCmbs & ~lowCmbs);
histogram(temp2, edges, 'FaceColor', [0.7 0.7 0.7], 'EdgeColor', 'none', 'Normalization', 'probability');
legend({'typeI', 'typeII'})
%% barplot for example neurons (long)
for i = 1340
    focusInd = i;
    time = (-lag:binSize:lag)/1000;   
    figure2;
    subplot(1,2,1)
    hold on;
    bar(time, allCrossCorr(focusInd, :), 'k')
    plot(time, allCI(1,:,focusInd), 'LineStyle', '--', 'Color', [0.4 0.4 0.4])
    plot(time, allCI(2,:,focusInd), 'LineStyle', '--', 'Color', [0.4 0.4 0.4])
    plot(time, allSurrMean(focusInd,:), 'Color', 'k')
    xlim([-2, 2])
    title(num2str(i))
    subplot(1,2,2)
    hold on;
    bar(time, allCrossCorrCorrected(focusInd, :), 'k')
    % plot(time, allSurrMean(focusInd,:), 'Color', 'k')
    xlim([-2, 2])
    title(num2str(i))
end
figure2;
hold on;
bar(time, allCrossCorrCorrected(focusInd, :), 'FaceColor', 'none', 'EdgeColor', [0.4 0.4 0.4])
bar(time([101-5:99, 103:101+5]), allCrossCorrCorrected(focusInd, [101-5:99, 103:101+5]), 'FaceColor', 'r')
xlim([-2, 2])
title(num2str(i))
%% barplot for example neurons (short)
for i = 1342
    focusInd = i;
    time = (-lag:binSize:lag);   
    figure2;
    subplot(1,2,1)
    hold on;
    bar(time, allCrossCorr(focusInd, :), 'k')
    plot(time, allCI(1,:,focusInd), 'LineStyle', '--', 'Color', [0.4 0.4 0.4])
    plot(time, allCI(2,:,focusInd), 'LineStyle', '--', 'Color', [0.4 0.4 0.4])
    plot(time, allSurrMean(focusInd,:), 'Color', 'k')
    xlim([-10, 10])
    title(num2str(i))
    subplot(1,2,2)
    hold on;
    bar(time, allCrossCorrCorrected(focusInd, :), 'k')
    % plot(time, allSurrMean(focusInd,:), 'Color', 'k')
    xlim([-10, 10])
    title(num2str(i))
end
% figure2;
% hold on;
% bar(time, allCrossCorrCorrected(focusInd, :), 'FaceColor', 'none', 'EdgeColor', [0.4 0.4 0.4])
% bar(time([101-5:99, 103:101+5]), allCrossCorrCorrected(focusInd, [101-5:99, 103:101+5]), 'FaceColor', 'r')
% xlim([-2, 2])
% title(num2str(i))
% 

%%



