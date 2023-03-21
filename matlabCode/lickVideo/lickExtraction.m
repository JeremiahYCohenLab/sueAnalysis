function lickExtraction(session, plotFlag)
[root, sep] = currComputer();
pd = parseSessionString_df(session, root, sep);
load([pd.sortedFolder 'lickSession.mat']);
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
    if plotFlag
        figure2;
        subplot(2,1,1); hold on;
        plot(time, conf);
    end
    confDiff = diff(conf);
    oddFrames = find((abs(confDiff(1:end-1))>oddTresh | abs(confDiff(2:end))>oddTresh) & confDiff(1:end-1).*confDiff(2:end)<0) + 1;
    while ~isempty(oddFrames)
        conf(oddFrames) = 0.5*(conf(oddFrames-1) + conf(oddFrames+1));
        confDiff = diff(conf);
        oddFrames = find((abs(confDiff(1:end-1))>oddTresh | abs(confDiff(2:end))>oddTresh) & confDiff(1:end-1).*confDiff(2:end)<0) + 1;               
    end
    % denoise conf
    if plotFlag
        plot(time, conf);
        xlim([0 s.lickLat(j)+200])
        plot([s.lickLat(j) s.lickLat(j)], [0 1], 'Color', 'r', 'LineStyle', '--', 'LineWidth', 2);
        plot([s.lickLat(j)+ s.rwdDelay + 50 s.lickLat(j) + s.rwdDelay + 50], [0 1], 'Color', 'r', 'LineStyle', '--', 'LineWidth', 2);
    %     scatter(tStarts, ones(size(tStarts)))
    %     scatter(tEnds, ones(size(tEnds)))
    end
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
    if plotFlag
        scatter(time(fStarts), ones(size(fStarts)), 12, 'b')
        scatter(time(fEnds), ones(size(fEnds)), 12, 'r')
        title([num2str(j) '  choice' num2str(s.allChoices(j))]);
        xlim([0 s.lickLat(j)+s.rwdDelay+100])
    end
    % plot port retraction
    if s.allChoices(j) > 0
        port = lickSession(s.responseInds(j)).porR;
    else
        port = lickSession(s.responseInds(j)).portL;
    end
    
    yChange = port(:,2);
    yChange = (yChange-min(yChange))/(max(yChange) - min(yChange));
    if plotFlag
        plot(time, yChange, 'Color', [0.7 0.7 0.7], 'LineWidth', 2)
    end
    
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
        if plotFlag
            text(mean(time(windows(w,:))), 0.5, num2str(windowsConf(w),2), 'FontSize', 12, 'HorizontalAlignment', 'center');
            set(gca, 'YColor', 'none')
        end
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
        if plotFlag
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
        end
%         xlabel('front       vertical       back')
%         ylabel('right      horizontal      left')
    end
end

savepath = [pd.sortedFolder session '_tongue.mat'];
save(savepath, 'allLicks', 'allLicks');
end