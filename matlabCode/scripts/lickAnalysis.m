clear all
session = 'mZS082d20220607';
[root, sep] = currComputer();
pd = parseSessionString_df(session, root, sep);
load([pd.sortedFolder 'lickSession.mat'], 'lickSession');
s = behAnalysisNoPlot_opMD(session, 'simpleFlag', 1);

% lick detection
tresh = 0.7;
upTresh = 0.7;
downTresh = 0.2; % default 0.2
oddTresh = 0.1;
focusTime = 1.5;
minILI = 15; % in ms
erosionTresh = 0.6;
erosionDisTresh = 50;
% tounge denoising
maxDeltaX = 10;
maxDeltaY = 10;


% get frameRate
frameRate = max([lickSession.videoFRs]);
% low pass
fc = 100;
sampFreq = frameRate;
[b, a] = butter(2,fc/(sampFreq/2),'low');
% get time
focusFrames = round(focusTime * frameRate);
timeCue = NaN(length(lickSession), 2);
% figure2; hold on;
for j = 1:length(lickSession)
    if s.CSplus(j)
        currTime = lickSession(j).time;
        currTime = currTime(1:focusFrames,3);
        timeCue(j,:) = findchangepts(currTime, 'MaxNumChanges', 2);
%         plot(currTime);
    end
end
timeCue = timeCue(:,1);
timeCue(isnan(timeCue)) = round(mean(timeCue, 'omitnan'));
% plot all first licks
% fc = 50;
% [b, a] = butter(2,fc/(frameRate/2),'low');
clear allWindows;
clear lickLick
clear allLicks
allWindows = cell(length(s.responseInds), 1);
lickLick = cell(length(s.responseInds), 1);
colNum = ceil(length(s.responseInds)/15);
% timeFig = figure;
% screen = get(0,'Screensize');
% screen(4) = screen(4) - 100;
% set(gcf, 'Position', screen);
% trajFig = figure;
% screen = get(0,'Screensize');
% screen(4) = screen(4) - 100;
% set(gcf, 'Position', screen);
%%
for j = 1:length(s.responseInds)
%%
    currTongue = lickSession(s.responseInds(j)).tongue;
    time = 1:length(currTongue);
    time = 1000/frameRate* (time - timeCue(s.responseInds(j)));
    
    x = currTongue(:,1);
    y = currTongue(:,2);
    conf = currTongue(:,3);
    
%     figure(timeFig);
%     subplot(15,colNum,j); hold on; 
    figure2;
    subplot(2,1,1); hold on;
    plot(time, conf);
    confDiff = diff(conf);
    oddFrames = find((abs(confDiff(1:end-1))>oddTresh | abs(confDiff(2:end))>oddTresh) & confDiff(1:end-1).*confDiff(2:end)<0) + 1;
    while ~isempty(oddFrames)
        conf(oddFrames) = 0.5*(conf(oddFrames-1) + conf(oddFrames+1));
        confDiff = diff(conf);
        oddFrames = find((abs(confDiff(1:end-1))>oddTresh | abs(confDiff(2:end))>oddTresh) & confDiff(1:end-1).*confDiff(2:end)<0) + 1;               
    end
    % denoise conf
    plot(time, conf);
    xlim([0 s.lickLat(j)+200])
    plot([s.lickLat(j) s.lickLat(j)], [0 1], 'Color', 'r', 'LineStyle', '--', 'LineWidth', 2);
    plot([s.lickLat(j)+ s.rwdDelay + 50 s.lickLat(j) + s.rwdDelay + 50], [0 1], 'Color', 'r', 'LineStyle', '--', 'LineWidth', 2);
%     scatter(tStarts, ones(size(tStarts)))
%     scatter(tEnds, ones(size(tEnds)))
    fStarts = conf - upTresh>0;
    fEnds = conf - downTresh>0;
    
    fStarts = find(fStarts(2:end-1) > fStarts(1:end-2) & fStarts(3:end) > fStarts(1:end-2))+1;
    tStarts = time(fStarts);
    fEnds = find(fEnds(2:end-1) < fEnds(1:end-2) & fEnds(3:end) < fEnds(1:end-2));
    fEnds = fEnds(fEnds>0);
    tEnds = time(fEnds);    
    
    
%     shortInd = find(diff(tStarts)< minILI) + 1;
%     tStarts(shortInd) = NaN;
%     tStarts = tStarts(~isnan(tStarts));
%     
%     shortInd = find(diff(tEnds)< minILI) + 1;
%     tEnds(shortInd) = NaN;
%     tEnds = tEnds(~isnan(tEnds));

    % merge close and same edges
    [sorted, sortedInds] = sort([fStarts; fEnds]);
%     [sorted, sortedInds] = sort([tStarts, tEnds]);
    IDsorted = [ones(size(fStarts)); zeros(size(fEnds))];
    IDsorted = IDsorted(sortedInds);
    shortID = diff(sorted)<frameRate*minILI/1000 & diff(IDsorted)==0; % same type and short latency
    shortSID = find(shortID & IDsorted(1:end-1) == 1)+1; % both start and short latency
    shortEID = find(shortID & IDsorted(1:end-1) == 0); % both end and short latency
    
    sorted(shortSID) = NaN;
    sorted(shortEID) = NaN;
    
    fStarts = sorted(IDsorted == 1);
    fEnds = sorted(IDsorted == 0);
    
    fStarts = fStarts(~isnan(fStarts));
    fEnds = fEnds(~isnan(fEnds));
    % only take those before outcome delivery
    scatter(time(fStarts), ones(size(fStarts)), 12, 'b')
    scatter(time(fEnds), ones(size(fEnds)), 12, 'r')
    title([num2str(j) '  choice' num2str(s.allChoices(j))]);
    xlim([0 s.lickLat(j)+s.rwdDelay+100])
    
    % plot port retraction
    if s.allChoices(j) > 0
        port = lickSession(s.responseInds(j)).porR;
    else
        port = lickSession(s.responseInds(j)).portL;
    end
    
    yChange = port(:,2);
    yChange = (yChange-min(yChange))/(max(yChange) - min(yChange));
    plot(time, yChange, 'Color', [0.7 0.7 0.7], 'LineWidth', 2)
    
    % find out all start-end intervals to see if mean higher than tresh,
    % decide all licks here
    [sorted, sortedInds] = sort([fStarts; fEnds]);
    IDsorted = [ones(size(fStarts)); zeros(size(fEnds))];
    IDsorted = IDsorted(sortedInds);
    
    windows = find(IDsorted(1:end-1)==1 & IDsorted(2:end)==0);
    windows = [sorted(windows), sorted(windows+1)];
    windows = windows((windows(:,2)-windows(:,1))>=frameRate*minILI/1000, :);
    
    windowsConf = NaN(size(windows,1), 1);
    allCurrLicks = cell(size(windows, 1), 1);
    allCurrLicksFilt = cell(size(windows, 1), 1);
    for w = 1:size(windows,1)
        currConf = conf(windows(w,1):windows(w,2));
        % window erosion
        % confidence 
        lastInd = find(currConf >= erosionTresh, 1, 'last');
        windows(w,2) = windows(w,1) + lastInd - 1;
        if isempty(lastInd)
            windows(w,2) = windows(w,1) + 1;
        else
            windows(w,2) = windows(w,1) + lastInd - 1;
        end
        % distance erosion
        tempLick = currTongue(windows(w,1):windows(w,2),:);
        diffTongueX = diff(tempLick(:,1));
        diffTongueY = diff(tempLick(:,2));
        disTemp = sqrt(diffTongueX.^2 + diffTongueY.^2);
        lastInd = find(disTemp > erosionDisTresh);
        if ~isempty(lastInd)
            windows(w,2) = windows(w,1) + lastInd(1) - 1;
        end
        %new lick
        tempLick = currTongue(windows(w,1):windows(w,2),:);
        diffTongueX = diff(tempLick(:,1));
        diffTongueY = diff(tempLick(:,2));
        % calcualte window conf
        windowsConf(w) = mean(conf(windows(w,1):windows(w,2)));    
%         subplot(2,1,2); hold on;
%         plot(tempLick(:,1), tempLick(:,2), 'Color', [0.6 0.6 0.6], 'LineWidth', 2, 'LineStyle', '--');
        % calculate consecutive vector angles
        if (windows(w,2)-windows(w,1))<frameRate*minILI/1000
            continue
        end
        tongueVec = diff(tempLick(:,1:2));
        cosSign = diag(tongueVec(1:end-1,:)*tongueVec(2:end,:)');
        oddFrames = find((abs(diffTongueX(1:end-1,1))>maxDeltaX | abs(diffTongueX(2:end,1))>maxDeltaX ... 
                        | abs(diffTongueY(1:end-1,1))>maxDeltaY | abs(diffTongueY(2:end,1))>maxDeltaY) ...
                        & cosSign<0) + 1;
           
        
        while ~isempty(oddFrames)
            tempLick(oddFrames,1:2) = 0.5*(tempLick(oddFrames-1, 1:2) + tempLick(oddFrames+1, 1:2));
            % calculate consecutive vector angles
            diffTongueX = diff(tempLick(:,1));
            diffTongueY = diff(tempLick(:,2));
            tongueVec = diff(tempLick(:,1:2));
            cosSign = diag(tongueVec(1:end-1,:)*tongueVec(2:end,:)');
            oddFrames = find((abs(diffTongueX(1:end-1,1))>maxDeltaX | abs(diffTongueX(2:end,1))>maxDeltaX ... 
                            | abs(diffTongueY(1:end-1,1))>maxDeltaY | abs(diffTongueY(2:end,1))>maxDeltaY) ...
                            & cosSign<0) + 1;
        end
        
        text(mean(time(windows(w,:))), 0.5, num2str(windowsConf(w),2), 'FontSize', 12, 'HorizontalAlignment', 'center');
        set(gca, 'YColor', 'none')
%         set(gca, 'XColor', 'none')
        allCurrLicks{w} = tempLick;       
        tempLick(:,1) = filtfilt(b, a, tempLick(:,1));
        tempLick(:,2) = filtfilt(b, a, tempLick(:,2));
        allCurrLicksFilt{w} = tempLick;
    end
   
    allCurrLicks = allCurrLicks(windowsConf>=tresh & (windows(:,2)-windows(:,1))>=frameRate*minILI/1000 & (windows(:,2)-windows(:,1))<frameRate*200/1000);
    allCurrLicksFilt = allCurrLicksFilt(windowsConf>=tresh & (windows(:,2)-windows(:,1))>=frameRate*minILI/1000 & (windows(:,2)-windows(:,1))<frameRate*200/1000);
    windows = windows(windowsConf>=tresh & (windows(:,2)-windows(:,1))>=frameRate*minILI/1000 & (windows(:,2)-windows(:,1))<frameRate*200/1000, :);
    % focus on first lick
    % find out which is the trigger lick 
    decisionID = find(time(windows(:,1))<=s.lickLat(j)+50 & time(windows(:,2))>=s.lickLat(j));
    allLicks(j).decisionID = decisionID;
    allLicks(j).windows = windows;
    allLicks(j).allLicks = allCurrLicks;
    allLicks(j).allLicksFilt = allCurrLicksFilt;
    allLicks(j).time = time;
    allLicks(j).conf = conf;
    
    % misc but useful positions
    portR = lickSession(s.responseInds(j)).porR;
    portR = mean(portR(portR(:,3)>0.9 & time'<0,1:2));
    portL = lickSession(s.responseInds(j)).portL;
    portL = mean(portL(portL(:,3)>0.9 & time'<0,1:2));
    
    midPoint = 0.5*(portL + portR);
    
    allLicks(j).portL = portL;
    allLicks(j).portR = portR;
    allLicks(j).midPoint = midPoint;
    
    if ~isempty(decisionID)
%         figure(trajFig);
%         subplot(15,colNum,j); hold on; 
        
        subplot(2,1,2); hold on;
        tempLick = currTongue(windows(decisionID,1):windows(decisionID,2),1:2)- midPoint;
        scatter(tempLick(1,2), tempLick(1,1), 25, 'k', 'filled')
        plot(tempLick(:,2), tempLick(:,1), 'Color', [0.6 0.6 0.6], 'LineWidth', 2, 'LineStyle', '--');
        tempLick = allCurrLicks{decisionID};
        tempLick = tempLick(:,1:2) - midPoint;
        plot(tempLick(:,2), tempLick(:,1), 'Color', 'r', 'LineWidth', 1);
        tempLick = allCurrLicksFilt{decisionID};
        tempLick = tempLick(:,1:2) - midPoint;
        plot(tempLick(:,2), tempLick(:,1), 'Color', 'B', 'LineWidth', 1, 'LineStyle', '--');   
        set(gca, 'YColor', 'none')
        set(gca, 'XColor', 'none')
        
%         xlabel('front       vertical       back')
%         ylabel('right      horizontal      left')
    end
end
%%
savepath = [pd.sortedFolder session '_tongue.mat'];
save(savepath, 'allLicks', 'allLicks');
%%
savepath = [pd.sortedFolder session '_tongue.mat'];
load(savepath);
%% focus on firstLick
% plot first trials on all licks
firstLicksX = [];
firstLicksY = [];
firstLicks = [];
firstTime = [];
tb = 5; % in ms
figure2;
for j = 1:length(allLicks)
    if allLicks(j).decisionID
        decisionID = allLicks(j).decisionID;
        disRL = norm(allLicks(j).portL - allLicks(j).portR);
        temp = allLicks(j).allLicksFilt{decisionID}*4/disRL;
        tempX = temp(:,1);
        tempY = temp(:,2);
        firstLicks = [firstLicks; temp];
        time = allLicks(j).time(allLicks(j).windows(decisionID,1): allLicks(j).windows(decisionID,2));
        horiLocation = 2*length(firstLicks);
        if s.allChoices(j) == 1
            subplot(1,4,1); hold on;
            plot(time-s.lickLat(j), tempX + horiLocation, 'Color', 'k');
            set(gca, 'YTick', [])
            subplot(1,4,2); hold on;
            plot(time-s.lickLat(j), tempY + horiLocation, 'Color', 'k');
            set(gca, 'YTick', [])
        else
            subplot(1,4,3); hold on;
            plot(time-s.lickLat(j), tempX + horiLocation, 'Color', 'k');
            set(gca, 'YTick', [])
            subplot(1,4,4); hold on;
            plot(time-s.lickLat(j), tempY + horiLocation, 'Color', 'k');
            set(gca, 'YTick', [])
        end
    end
end
%%
figure2;
for j = 1:length(allLicks)
    if allLicks(j).decisionID
        decisionID = allLicks(j).decisionID;
        temp = allLicks(j).allLicksFilt{decisionID};
        tempX = temp(:,1);
        tempY = temp(:,2);
        firstLicks = [firstLicks; temp];
        time = allLicks(j).time(allLicks(j).windows(decisionID,1): allLicks(j).windows(decisionID,2));
        horiLocation = 30*length(firstLicks);
        hold on; 
        plot(tempX-allLicks(j).midPoint(1), tempY-allLicks(j).midPoint(2), 'Color', [0.6 0.6 0.6], 'LineStyle', ':');
        
%             subplot(1,2,1); hold on;
%             plot(tempX-allLicks(j).midPoint(1), tempY-allLicks(j).midPoint(2), 'Color', [0.6 0.6 0.6], 'LineStyle', ':');
%             set(gca, 'YTick', [])
%         else
%             subplot(1,2,2); hold on;
%             plot(tempX-allLicks(j).midPoint(1), tempY-allLicks(j).midPoint(2), 'Color', [0.6 0.6 0.6], 'LineStyle', ':');
%             set(gca, 'YTick', [])
%         end
    end
end
%% load model variables;
modelName = '5params'; 
col = 'cueOnGood';
sampNum = 2000;
[root,sep] = currComputer();
pd = parseSessionString_df(session, root, sep);
params = getStanModelParams_sampsOnly(pd.aniName, col, modelName, sampNum, 'sessionName', session, 'biasFlag',1);
t = inferModelVar(session, params, modelName, 'perturb', 0);
s = behAnalysisNoPlot_opMD(session, 'simpleFlag', 1);
%% value extraction

Qdiff = t.Q(:,2) - t.Q(:,1);
QdiffChosen = Qdiff;
QdiffChosen(s.allChoices==-1) = - QdiffChosen(s.allChoices==-1);
pRight = t.probChoice;
pRight(s.allChoices==1) = 1 - t.probChoice(s.allChoices==1);

ind  = 1:length(allLicks);
meanX = NaN(1, length(ind));
meanY = NaN(1, length(ind));
meanXPre = NaN(1, length(ind));
meanXMax = NaN(1, length(ind));
meanYMax = NaN(1, length(ind));
speedAtPort = NaN(1, length(ind));
vXAtPort = NaN(1, length(ind));
vYAtPort = NaN(1, length(ind));
vXPrePeak = NaN(1, length(ind));
vYPrePeak = NaN(1, length(ind));
vXFixInter = NaN(1, length(ind));
vYFixInter = NaN(1, length(ind));
speedPrePeak = NaN(1, length(ind));
crossInd = NaN(1, length(ind));
maxIndAll = NaN(1, length(ind));
crossX = NaN(1, length(ind));
crossY = NaN(1, length(ind));
thresh = 0;
preBins = 6;

for j = 1:length(ind)
    if allLicks(ind(j)).decisionID
        decisionID = allLicks(ind(j)).decisionID;
        interval = allLicks(ind(j)).time(2) - allLicks(ind(j)).time(1);
        disRL = norm(allLicks(ind(j)).portL - allLicks(ind(j)).portR);
        temp = allLicks(ind(j)).allLicksFilt{decisionID};
        tempX = (temp(:,1) - allLicks(ind(j)).midPoint(1)) * 4/disRL;
        tempY = -(temp(:,2) - allLicks(ind(j)).midPoint(2)) * 4/disRL;
        [meanYMax(ind(j)), maxInd] = max(tempY);
        maxIndAll(ind(j)) = maxInd;
        meanXMax(ind(j)) = tempX(maxInd);
        meanXPre(ind(j)) = mean(tempX(1:maxInd));
        
        meanX(ind(j)) = mean(tempX);
        meanY(ind(j)) = mean(tempY);
        time = allLicks(ind(j)).time(allLicks(ind(j)).windows(decisionID,1): allLicks(ind(j)).windows(decisionID,2));

        lickDisTmpX = diff(tempX);
        lickDisTmpY = diff(tempY);
        lickSpeed = sqrt(lickDisTmpX.^2 + lickDisTmpY.^2)/(interval);
        lickVX = lickDisTmpX/(allLicks(ind(j)).time(2) - allLicks(ind(j)).time(1));
        lickVY = lickDisTmpY/(allLicks(ind(j)).time(2) - allLicks(ind(j)).time(1));
       
        if maxInd >= preBins+1 && maxInd < length(tempX)-2
            vXPrePeak(ind(j)) = (tempX(maxInd) - tempX(maxInd-preBins))/(preBins*interval);
            vYPrePeak(ind(j)) = (tempY(maxInd) - tempY(maxInd-preBins))/(preBins*interval);
            speedPrePeak(ind(j)) = mean(lickSpeed(maxInd-1:maxInd));
            hold on;
        end
        % find threshold crossing
        cross = tempY >thresh;
        cross = find(cross(1:end-1)==0 & cross(2:end)==1, 1);
        if ~isempty(cross)
            crossInd(ind(j)) = cross;
            speedAtPort(ind(j)) = lickSpeed(cross);
            vXAtPort(ind(j)) = lickVX(cross);
            vYAtPort(ind(j)) = lickVY(cross);
            crossX(ind(j)) = mean(tempX(cross));
            crossY(ind(j)) = mean(tempY(cross));
        end
        
        % use fixed inter
        if ~isempty(cross) && maxInd >= preBins+1 && maxInd < length(tempX)-2
            vXFixInter(ind(j)) = (tempX(maxInd) - tempX(cross))/((maxInd-cross)*interval);
            vYFixInter(ind(j)) = (tempY(maxInd) - tempY(cross))/((maxInd-cross)*interval);
        end

    end
end
%%
% speedAtPort(s.allChoices==1 & ~isnan(speedAtPort)) = zscore(speedAtPort(s.allChoices==1 & ~isnan(speedAtPort)));
% speedAtPort(s.allChoices==-1 & ~isnan(speedAtPort)) = zscore(speedAtPort(s.allChoices==-1 & ~isnan(speedAtPort)));

vXFixInter(s.allChoices==1 & ~isnan(vXFixInter)) = zscore(vXFixInter(s.allChoices==1 & ~isnan(vXFixInter)));
vXFixInter(s.allChoices==-1 & ~isnan(vXFixInter)) = zscore(vXFixInter(s.allChoices==-1 & ~isnan(vXFixInter)));
vXFixInterAbs = vXFixInter;
vXFixInterAbs(s.allChoices==1) = -vXFixInterAbs(s.allChoices==1);
%% extract value and plot
Qdiff = t.Q(:,2) - t.Q(:,1);
QdiffChosen = Qdiff;
QdiffChosen(s.allChoices==-1) = - QdiffChosen(s.allChoices==-1);
pRight = t.probChoice;
pRight(s.allChoices==1) = 1 - t.probChoice(s.allChoices==1);

[~, ind] = sort(Qdiff);
colors = cool(length(s.allChoices));

figure2; hold on;
meanX = NaN(1, length(ind));
meanY = NaN(1, length(ind));
meanXPre = NaN(1, length(ind));
meanXMax = NaN(1, length(ind));
meanYMax = NaN(1, length(ind));
speedAtPort = NaN(1, length(ind));
vXAtPort = NaN(1, length(ind));
vYAtPort = NaN(1, length(ind));
vXPrePeak = NaN(1, length(ind));
vYPrePeak = NaN(1, length(ind));
vXFixInter = NaN(1, length(ind));
vYFixInter = NaN(1, length(ind));
speedPrePeak = NaN(1, length(ind));
crossInd = NaN(1, length(ind));
maxIndAll = NaN(1, length(ind));
crossX = NaN(1, length(ind));
crossY = NaN(1, length(ind));
thresh = 0;
preBins = 6;
figure2;
for j = 1:length(ind)
    if allLicks(ind(j)).decisionID
        decisionID = allLicks(ind(j)).decisionID;
        interval = allLicks(ind(j)).time(2) - allLicks(ind(j)).time(1);
        disRL = norm(allLicks(ind(j)).portL - allLicks(ind(j)).portR);
        temp = allLicks(ind(j)).allLicksFilt{decisionID};
        tempX = (temp(:,1) - allLicks(ind(j)).midPoint(1)) * 4/disRL;
        tempY = -(temp(:,2) - allLicks(ind(j)).midPoint(2)) * 4/disRL;
        [meanYMax(ind(j)), maxInd] = max(tempY);
        maxIndAll(ind(j)) = maxInd;
        meanXMax(ind(j)) = tempX(maxInd);
        meanXPre(ind(j)) = mean(tempX(1:maxInd));
        
        s.meanX(ind(j)) = mean(tempX);
        meanY(ind(j)) = mean(tempY);
%         firstLicks = [firstLicks; temp];
        time = allLicks(ind(j)).time(allLicks(ind(j)).windows(decisionID,1): allLicks(ind(j)).windows(decisionID,2));
%         horiLocation = 30*length(firstLicks);
        hold on; 
        plot(tempX, tempY, 'LineWidth', 0.5, 'LineStyle', ':', 'Color', colors(j,:));
%         patchline(tempX, tempY, 'linestyle', ':', 'edgecolor', colors(j,:), 'linewidth', 0.5, 'edgealpha', 0.4);
        % scatter(meanXMax(ind(j)), meanYMax(ind(j)), 15, colors(j,:), 'filled', 'MarkerFaceAlpha', 0.75);

        lickDisTmpX = diff(tempX);
        lickDisTmpY = diff(tempY);
        lickSpeed = sqrt(lickDisTmpX.^2 + lickDisTmpY.^2)/(interval);
        lickVX = lickDisTmpX/(allLicks(ind(j)).time(2) - allLicks(ind(j)).time(1));
        lickVY = lickDisTmpY/(allLicks(ind(j)).time(2) - allLicks(ind(j)).time(1));
       
        % find peak       
        if maxInd >= preBins+1 && maxInd < length(tempX)-2
            vXPrePeak(ind(j)) = (tempX(maxInd) - tempX(maxInd-preBins))/(preBins*interval);
            vYPrePeak(ind(j)) = (tempY(maxInd) - tempY(maxInd-preBins))/(preBins*interval);
            speedPrePeak(ind(j)) = mean(lickSpeed(maxInd-1:maxInd));
            hold on;
%             if s.allChoices(ind(j)) == 1
%                 scatter(tempX(maxInd) - tempX(maxInd-2), rand(1), 10, 'r');
%             else
%                 scatter(tempX(maxInd) - tempX(maxInd-2), rand(1), 10, 'b');
%             end
%             figure2; hold on;
%             plot(tempX, tempY, 'Color', 'k');
%             scatter(tempX(1), tempY(1), 15, 'r');
%             scatter(tempX(maxInd), tempY(maxInd), 20, 'k');
%             scatter(tempX(maxInd-preBins), tempY(maxInd-preBins), 20, 'g');
        end
        % find threshold crossing
        cross = tempY >thresh;
        cross = find(cross(1:end-1)==0 & cross(2:end)==1, 1);
        if ~isempty(cross)
            crossInd(ind(j)) = cross;
            speedAtPort(ind(j)) = lickSpeed(cross);
            vXAtPort(ind(j)) = lickVX(cross);
            vYAtPort(ind(j)) = lickVY(cross);
            crossX(ind(j)) = mean(tempX(cross));
            crossY(ind(j)) = mean(tempY(cross));
        end
        
        % use fixed inter
        if ~isempty(cross) && maxInd >= preBins+1 && maxInd < length(tempX)-2
            vXFixInter(ind(j)) = (tempX(maxInd) - tempX(cross))/((maxInd-cross)*interval);
            vYFixInter(ind(j)) = (tempY(maxInd) - tempY(cross))/((maxInd-cross)*interval);
        end

    end
end

speedAtPort(s.allChoices==1 & ~isnan(speedAtPort)) = zscore(speedAtPort(s.allChoices==1 & ~isnan(speedAtPort)));
speedAtPort(s.allChoices==-1 & ~isnan(speedAtPort)) = zscore(speedAtPort(s.allChoices==-1 & ~isnan(speedAtPort)));

% vYAtPort(s.allChoices==1 & ~isnan(vYAtPort)) = zscore(speedAtPort(s.allChoices==1 & ~isnan(vYAtPort)));
% vYAtPort(s.allChoices==-1 & ~isnan(vYAtPort)) = zscore(speedAtPort(s.allChoices==-1 & ~isnan(vYAtPort)));

vXFixInter(s.allChoices==1 & ~isnan(vXFixInter)) = zscore(vXFixInter(s.allChoices==1 & ~isnan(vXFixInter)));
vXFixInter(s.allChoices==-1 & ~isnan(vXFixInter)) = zscore(vXFixInter(s.allChoices==-1 & ~isnan(vXFixInter)));
%% colorCodeQdiffC on maxPosition
figure2;
[~, ind] = sort(abs(QdiffChosen));
colors = cool(length(s.allChoices));
scatter(meanXMax(ind), meanYMax(ind), 15, colors,  'filled', 'MarkerFaceAlpha', 0.75);
%% tuning with location
numBins = 3;
target = QdiffChosen;
m = meanXMax;
% edges = quantile(target, linspace(0, 1, numBins+1));
edges = linspace(min(target)-0.001, max(target)+0.001, numBins+1);
edges(1) = edges(1) - 0.0001;   
edges(end) = edges(end) + 0.0001;
meanTarget = NaN(1, numBins);
meanXLoc = NaN(1, numBins);
meanYLoc = NaN(1, numBins);
semXLoc = NaN(1, numBins);
semYLoc = NaN(1, numBins);

for j = 1:numBins
    meanTarget(j) = mean(target(target>=edges(j) & target<edges(j+1)), 'omitnan');
    meanXLoc(j) = mean(m(target>=edges(j) & target<edges(j+1)), 'omitnan');
    semXLoc(j) = sem(m(target>=edges(j) & target<edges(j+1)));
    meanYLoc(j) = mean(meanYMax(target>=edges(j) & target<edges(j+1)), 'omitnan');
    semYLoc(j) = sem(meanYMax(target>=edges(j) & target<edges(j+1)));
end

figure2;
% subplot(1,2,1);
hold on;
plot(meanTarget, meanXLoc, 'Color', 'k', 'LineWidth', 2);
patch([meanTarget, flip(meanTarget)], [meanXLoc-semXLoc, flip(meanXLoc+semXLoc)], [0.5 0.5 0.5], 'FaceAlpha', 0.25,  'EdgeColor','none')
% x = target;
% y = m;
% x = x(~isnan(y));
% y = y(~isnan(y));
% [R, P] = corrcoef(x, y);
lmTmp = fitlm(target, m);
xlabel('Qdiff')
ylabel('maxX')
title(['p=' num2str(lmTmp.Coefficients.pValue(end)) 'Coeff=' num2str(num2str(lmTmp.Coefficients.Estimate(end))), 'tState=' num2str(num2str(lmTmp.Coefficients.tStat(end)))]);
%%
figure;
% subplot(1,2,1);
hold on;
plot(meanTarget, meanYLoc, 'Color', 'k', 'LineWidth', 2);
patch([meanTarget, flip(meanTarget)], [meanYLoc-semYLoc, flip(meanYLoc+semYLoc)], [0.5 0.5 0.5], 'FaceAlpha', 0.25,  'EdgeColor','none')
x = target;
y = meanYMax;
x = x(~isnan(y));
y = y(~isnan(y));
[R, P] = corrcoef(x, y);
xlabel('QdiffChosen')
ylabel('Ymax')
title(['p=' num2str(P(2,1)) 'Coeff=' num2str(R(2,1)), 'n=' num2str(length(x))]);
% lmTmp = fitlm(target, meanYMax);
% xlabel('QdiffChosen')
% ylabel('maxY')
% title(['p=' num2str(lmTmp.Coefficients.pValue(end)) 'Coeff=' num2str(num2str(lmTmp.Coefficients.Estimate(end))), 'tState=' num2str(num2str(lmTmp.Coefficients.tStat(end)))]);

%% scatters
% colors = [linspace(0, 1, length(ind))', zeros(size(ind)), zeros(size(ind))];
% vXFixInter needs to be centered and flipped before being used
% vXFixInter(s.allChoices == -1) = -vXFixInter(s.allChoices == -1);
[sortedTarget, indSpeeds] = sort(vXFixInterAbs);
colors = cool(sum(~isnan(vYFixInter)));
figure2; hold on;
scatter(meanXMax(indSpeeds(1:sum(~isnan(vXFixInter)))), meanYMax(indSpeeds(1:sum(~isnan(vXFixInter)))), 12, colors, 'filled', 'MarkerFaceAlpha', 0.7);
%% tuning with speed
numBins = 3;
target = QdiffChosen;
m = vXFixInterAbs;
% m = m + 0.01;
% m(s.allChoices==-1) = -m(s.allChoices==-1);
edges = quantile(target, linspace(0, 1, numBins+1));
% edges = linspace(min(target)-0.001, max(target)+0.001, numBins+1);
edges(1) = edges(1) - 0.0001;
edges(end) = edges(end) + 0.0001;
meanTarget = NaN(1, numBins);
meanXLoc = NaN(1, numBins);
meanYLoc = NaN(1, numBins);
semXLoc = NaN(1, numBins);
semYLoc = NaN(1, numBins);

for j = 1:numBins
    meanTarget(j) = mean(target(target>=edges(j) & target<edges(j+1)), 'omitnan');
    meanXLoc(j) = mean(m(target>=edges(j) & target<edges(j+1)), 'omitnan');
    semXLoc(j) = sem(m(target>=edges(j) & target<edges(j+1)));
    meanYLoc(j) = mean(meanYMax(target>=edges(j) & target<edges(j+1)), 'omitnan');
    semYLoc(j) = sem(meanYMax(target>=edges(j) & target<edges(j+1)));
end

figure;
% subplot(1,2,1);
hold on;
plot(meanTarget, meanXLoc, 'Color', 'k', 'LineWidth', 2);
patch([meanTarget, flip(meanTarget)], [meanXLoc-semXLoc, flip(meanXLoc+semXLoc)], [0.5 0.5 0.5], 'FaceAlpha', 0.25,  'EdgeColor','none')
% x = target;
% y = m;
% x = x(~isnan(y));
% y = y(~isnan(y));
% [R, P] = corrcoef(x, y);
lmTmp = fitlm(target, m);
xlabel('Qdiff')
ylabel('VX(mm/ms)')
title(['p=' num2str(lmTmp.Coefficients.pValue(end)) 'Coeff=' num2str(num2str(lmTmp.Coefficients.Estimate(end))), 'tState=' num2str(num2str(lmTmp.Coefficients.tStat(end)))]);
%% compare stay and switch
figure2; hold on;
scatter(meanXMax(s.stayChoice_Inds), meanYMax(s.stayChoice_Inds));
scatter(meanXMax(s.changeChoice_Inds), meanYMax(s.changeChoice_Inds));
legend({'stay', 'switch'})
figure2; hold on;
edges = linspace(min(vXFixInterAbs)-0.0001, max(vXFixInterAbs)+0.0001, 10);
histogram(vXFixInterAbs(s.stayChoice_Inds), edges, 'Normalization', 'probability')
histogram(vXFixInterAbs(s.changeChoice_Inds), edges, 'Normalization', 'probability')
legend('stay', 'switch')
[h, p, ~, stats] = ttest2(vXFixInterAbs(s.stayChoice_Inds), vXFixInterAbs(s.changeChoice_Inds));
title(['svsVX' ' p=' num2str(p)]);
figure2; hold on;
edges = linspace(min(abs(meanXMax)-0.0001), max(abs(meanXMax))+0.0001, 10);
histogram(abs(meanXMax(s.stayChoice_Inds)), edges, 'Normalization', 'probability')
histogram(abs(meanXMax(s.changeChoice_Inds)), edges, 'Normalization', 'probability')
legend('stay', 'switch')
[h, p, ~, stats] = ttest2(abs(meanXMax(s.stayChoice_Inds)), abs(meanXMax(s.changeChoice_Inds)));
title(['maxX' ' p=' num2str(p)]);
%% compare lasered or not, current time
figure2; hold on;
scatter(meanXMax(s.laser==1), meanYMax(s.laser==1));
scatter(meanXMax(s.laser==0), meanYMax(s.laser==0));
legend({'laser', 'no laser'})
figure2; hold on;
edges = linspace(min(vXFixInterAbs)-0.0001, max(vXFixInterAbs)+0.0001, 10);
histogram(vXFixInterAbs(s.laser==1), edges, 'Normalization', 'probability')
histogram(vXFixInterAbs(s.laser==0), edges, 'Normalization', 'probability')
legend('laser', 'no laser')
[h, p, ~, stats] = ttest2(vXFixInterAbs(s.laser==1), vXFixInterAbs(s.laser==0));
title(['laserVX' ' p=' num2str(p)]);
figure2; hold on;
edges = linspace(min(abs(meanXMax))-0.0001, max(abs(meanXMax))+0.0001, 10);
histogram(abs(meanXMax(s.laser==1)), edges, 'Normalization', 'probability')
histogram(abs(meanXMax(s.laser==0)), edges, 'Normalization', 'probability')
legend('laser', 'no laser')
[h, p, ~, stats] = ttest2(abs(meanXMax(s.laser==1)), abs(meanXMax(s.laser==0)));
title(['laserxMax' ' p=' num2str(p)]);
%% compare lasered or not, next trial stay
figure2; hold on;
scatter(meanXMax(intersect(find(s.laser==1), s.stayChoice_Inds)), meanYMax(intersect(find(s.laser==1), s.stayChoice_Inds)));
scatter(meanXMax(intersect(find(s.laser==0), s.stayChoice_Inds)), meanYMax(intersect(find(s.laser==0), s.stayChoice_Inds)));
legend({'laser', 'no laser'})
figure2; hold on;
edges = linspace(min(vXFixInterAbs)-0.0001, max(vXFixInterAbs)+0.0001, 10);
histogram(vXFixInterAbs(intersect(find(s.laser==1), s.stayChoice_Inds)), edges, 'Normalization', 'probability')
histogram(vXFixInterAbs(intersect(find(s.laser==0), s.stayChoice_Inds)), edges, 'Normalization', 'probability')
legend('laser', 'no laser')
[h, p, ~, stats] = ttest2(vXFixInterAbs(intersect(find(s.laser==1), s.stayChoice_Inds)), vXFixInterAbs(intersect(find(s.laser==0), s.stayChoice_Inds)));
title(['laserVX' ' p=' num2str(p)]);
figure2; hold on;
edges = linspace(min(abs(meanXMax)), max(abs(meanXMax)), 10);
histogram(abs(meanXMax(intersect(find(s.laser==1), s.stayChoice_Inds))), edges, 'Normalization', 'probability')
histogram(abs(meanXMax(intersect(find(s.laser==0), s.stayChoice_Inds))), edges, 'Normalization', 'probability')
legend('laser', 'no laser')
[h, p, ~, stats] = ttest2(abs(meanXMax(intersect(find(s.laser==1), s.stayChoice_Inds))), abs(meanXMax(intersect(find(s.laser==0), s.stayChoice_Inds))));
title(['laserxMax' ' p=' num2str(p)]);
%%
