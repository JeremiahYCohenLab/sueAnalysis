function lickPreprocessing(session)
%% get model settings
cdnn = 'tongueTrackingStraightDec1';
shuffle = 1;
iter = 200000;
%% get all video names
[root, sep] = currComputerVideo();
pd = parseSessionString_df(session, root, sep);
allFiles = dir(pd.lickPath);
allFiles = {allFiles.name};
% aviExpression = ['^' session 'DLC' '\w*' '200000.csv' '$'];
aviExpression = ['.avi' '$'];

videoInd = cellfun(@(x) ~isempty(regexp(x, aviExpression, 'once')), allFiles);
videoNames = allFiles(videoInd)';
%% get information of all videos create tabel
lickS = struct;
for i = 1:length(videoNames)
    v = VideoReader([pd.lickPath videoNames{i}]);
    lickS(i).name = v.Name;
    lickS(i).videoFRs = v.FrameRate;
    lickS(i).videoLength = v.Duration;
    time = split(videoNames{i}, {'-', '.'});
    time = cellfun(@str2num, time(1:4));
    if time(1)<=8
        time(1) = time(1)+24;
    end
    lickS(i).videoTime = time(1)*60*60*1000 + time(2)*60*1000 + time(3)*1000 + time(4);
end
%% clean up short videos in ITI
lickS = lickS([lickS.videoLength]>3);
[~, sortID] = sort([lickS.videoTime]);
lickS = lickS([sortID]);
%% time align
s = behAnalysisNoPlot_opMD(session, 'simpleFlag', 1);
startTimes = [s.behSessionData.CSon] - s.behSessionData(1).CSon;
videoTime = [lickS.videoTime]-lickS(1).videoTime;
% linear regression
wellAligned = false;
humanCheck = true;
if length(startTimes)== length(videoTime)
    
    lm = fitlm(startTimes,videoTime);
    transformedTime = lm.Fitted;
    if lm.MSE < 200^2
        wellAligned = true;
        humanCheck = false;
        lickSession = lickS;
    end
else
    if length(startTimes)<length(videoTime)
        b = abs(startTimes' - videoTime); % align by first 100 trials
        [minDis, bestVideo] = min(b, [], 2);
         if max(minDis) < 1000*length(lickS)/200
            lickSession = lickS(bestVideo);
            wellAligned = true;
            humanCheck = false;   
         end
    else
        b = abs(startTimes' - videoTime);
        [minDis, bestTrial] = min(b, [], 1);
        if max(minDis) < 1000
            lickSession(bestTrial) = lickS;
            wellAligned = true;
            humanCheck = false;   
        end
    end
end


%  fprintf([session ' potential mismatch double check \n'])
% alignment check
 xBeh = zeros(1,max(startTimes)+10);
 xBeh(startTimes + 1) = 1;
 figure2;hold on;
 plot(xBeh)
 scatter(videoTime, ones(size(videoTime)));
 ylim([0.95 1.05])
 xlim([0 200000])
 labels = num2cell(1:length(startTimes));
 labels = cellfun(@num2str, labels, 'UniformOutput', false);
 text(startTimes, 1.005*ones(size(startTimes)), labels, 'HorizontalAlignment', 'center');
 sgtitle([session ' ' num2str(double(wellAligned))]);
 
 
%% load all traced body parts
suffix = ['DLC_resnet50_' cdnn 'shuffle' num2str(shuffle) '_' num2str(iter) '.csv'];
bodyParts = {'midUp', 'midBottom', 'portL', 'porR', 'lip', 'tongue', 'time'};
for i = 1:length(lickSession)
%     [~, dayList, ~] = xlsread([root xlFile '.xlsx'], sheet);
    nameCurr = split(lickSession(i).name, '.avi');
    nameCurr = nameCurr{1};
    M = readmatrix([pd.lickPath nameCurr suffix]);
    for j = 1:length(bodyParts)
        currInds = (3*(j-1)+2):(3*j+1);
        currLocs = M(:,currInds);
        lickSession(i).(bodyParts{j}) = currLocs;
    end
end
[root, sep] = currComputer();
pd = parseSessionString_df(session, root, sep);
save([pd.sortedFolder 'lickSession.mat'], 'lickSession')


