% load from session folder: 
videoName = '19-11-11.958';
rootVideo = 'D:\video\';
[root, sep] = currComputer(); 
dlcName = 'tongueTrackingStraightDec1shuffle1_250000';
% play raw video
pdVideo = parseSessionString_df(session, rootVideo, sep);
lickVideo = [pdVideo.lickPath videoName '.avi'];
% implay(lickVideo);
% load trajectories
pd = parseSessionString_df(session, root, sep);
s = behAnalysisNoPlot_opMD(session, 'simpleFlag', 1);
trajPath = [pd.sortedFolder 'lickSession.mat'];
load(trajPath);
tonguePath = [pd.sortedFolder session '_tongue.mat'];
load(tonguePath);
%% find target lick
ind = contains({lickSession.name}, videoName);
indResp = ind(s.responseInds);

tongueWindows = allLicks(indResp).windows;
tongueTrajs = allLicks(indResp).allLicksFilt;
time = allLicks(indResp).time;
rawTrace = lickSession(ind).tongue;
LLoc = allLicks(indResp).portL;
RLoc = allLicks(indResp).portR;

ratio = norm(LLoc - RLoc)/4;
tongueStates = zeros(size(time));
tongueTrajsWhole = NaN(size(time,2), 2);
for w = 1:length(tongueWindows)
    tongueStates(tongueWindows(w,1):tongueWindows(w,2)) = 1;
    tongueTrajsWhole(tongueWindows(w,1):tongueWindows(w,2),:) = tongueTrajs{w}(:,1:2);
end
%% trial structure
cueStates = zeros(size(time));
cueStates(time>0 & time<500) = 1;

responseStates  = zeros(size(time));
responseStates(time>s.lickLat(indResp) & time<s.lickLat(indResp)+20) = 1;

rewardStates = zeros(size(time));
rewardStates(time>s.lickLat(indResp)+s.rwdDelay & time<s.lickLat(indResp) + 20 + s.rwdDelay) = 1;

%% make tongue on video whole video
Vobj = VideoReader(lickVideo);
n = Vobj.NumFrames;


v = VideoWriter([lickVideo(1:end-4) '_test.avi']);
v.FrameRate = 90;

open(v);

startF = find(time>-200, 1);
endF = find(time>1100, 1);
figure('Position', [0 0 1500 750]);

for i = startF:endF
    matrix_image = read(Vobj, i);
    subplot(1,2,1)
    imshow(matrix_image);
    hold on;

    if tongueStates(i)>0
        scatter(tongueTrajsWhole(i, 1), tongueTrajsWhole(i, 2), 70, 'red', 'filled');
    end
    plot([630 630-2*ratio], [460 460], 'LineWidth', 2, 'Color', 'w')
    text(630, 430, sprintf('%d mm', 2), 'FontSize', 20, 'Color', [1, 1, 1], 'HorizontalAlignment','right')    
    text(315, 510, sprintf('%d ms', round(time(i))), 'FontSize', 20, 'Color', [0 0 0], 'HorizontalAlignment','center')
    hold off
    subplot(10,2,6)
    plot(time(startF:i), cueStates(startF:i), 'LineWidth', 3, 'Color', 'k');
    xlim([time(startF) time(endF)])
    ylim([0 1])
    set(gca, 'XColor', 'none')
    set(gca, 'YColor', 'none')
    text(-100+time(startF), 0.5, 'Cue', 'FontSize', 20, 'HorizontalAlignment', 'right')


    subplot(10,2,8)
    plot(time(startF:i), responseStates(startF:i), 'LineWidth', 3, 'Color', 'k');
    xlim([time(startF) time(endF)])
    ylim([0 1])
    set(gca, 'XColor', 'none')
    set(gca, 'YColor', 'none')
    text(-50+time(startF), 0.5, 'Choice', 'FontSize', 20, 'HorizontalAlignment', 'right')

    subplot(10,2,10)
    plot(time(startF:i), rewardStates(startF:i), 'LineWidth', 3, 'Color', 'k');
    xlim([time(startF) time(endF)])
    ylim([0 1])
    set(gca, 'XColor', 'none')
    set(gca, 'YColor', 'none')
    text(-50+time(startF), 0.5, 'Reward', 'FontSize', 20, 'HorizontalAlignment', 'right')

    subplot(10,2,[12, 14, 16])
    plot(time(startF:i), -tongueTrajsWhole(startF:i, 2), 'LineWidth', 3, 'Color', 'k');
    xlim([time(startF) time(endF)])
    ylim([-350 -200])
    set(gca, 'XColor', 'none')
    set(gca, 'YColor', 'none')
    text(-50+time(startF), -275, 'Tongue', 'FontSize', 20, 'HorizontalAlignment', 'right')

    frame = getframe(gcf);
    
 
  writeVideo(v,frame);  
 
end
close(v)
%% single lick
Vobj = VideoReader(lickVideo);
n = Vobj.NumFrames;


v = VideoWriter([lickVideo(1:end-4) '_testFirst.avi']);
v.FrameRate = 22;

open(v);

startF = tongueWindows(2, 1)-3;
endF = tongueWindows(2, 2)+3;
figure('Position', [0 0 1000 1000]);
crop = [100, 480; 1, 400];
[~, peakFrame] = min(tongueTrajsWhole(tongueWindows(2, 1):tongueWindows(2, 2), 2)); 
peakFrame = peakFrame + tongueWindows(2, 1) - 1;
for i = startF:endF
    matrix_image = read(Vobj, i);
    % subplot(4,4,[5:7, 9:11, 13:15])
    imshow(matrix_image(crop(1,1):crop(1,2), crop(2,1):crop(2,2), :));
    hold on;

    if tongueStates(i)>0
        plot(tongueTrajsWhole(startF:i, 1)+crop(2,1), tongueTrajsWhole(startF:i, 2)-crop(1,1), 'red', 'LineWidth', 3);
        if i >= peakFrame
           scatter(tongueTrajsWhole(peakFrame, 1)+crop(2,1), tongueTrajsWhole(peakFrame, 2)-crop(1,1), 60, 'k', 'LineWidth', 3);
        end
    end
    plot([crop(1,2)-crop(1,1) crop(1,2)-crop(1,1)-1*ratio], [360 360], 'LineWidth', 2, 'Color', 'w')
    text(crop(1,2)-crop(1,1)-5, crop(1,2)-crop(1,1)-50, sprintf('%d mm', 1), 'FontSize', 20, 'Color', [1, 1, 1], 'HorizontalAlignment','right') 
    text(0.5*(crop(1,2)-crop(1,1)), crop(2,2)-crop(2,1)+20, sprintf('%d ms', round(time(i)-time(startF))), 'FontSize', 20, 'Color', [0 0 0], 'HorizontalAlignment','center')
    hold off
    % subplot(10,2,6)
    % plot(time(startF:i), cueStates(startF:i), 'LineWidth', 3, 'Color', 'k');
    % xlim([time(startF) time(endF)])
    % ylim([0 1])
    % set(gca, 'XColor', 'none')
    % set(gca, 'YColor', 'none')
    % text(-100+time(startF), 0.5, 'Cue', 'FontSize', 20, 'HorizontalAlignment', 'right')
    % 
    % 
    % subplot(10,2,8)
    % plot(time(startF:i), responseStates(startF:i), 'LineWidth', 3, 'Color', 'k');
    % xlim([time(startF) time(endF)])
    % ylim([0 1])
    % set(gca, 'XColor', 'none')
    % set(gca, 'YColor', 'none')
    % text(-50+time(startF), 0.5, 'Choice', 'FontSize', 20, 'HorizontalAlignment', 'right')
    % 
    % subplot(10,2,10)
    % plot(time(startF:i), rewardStates(startF:i), 'LineWidth', 3, 'Color', 'k');
    % xlim([time(startF) time(endF)])
    % ylim([0 1])
    % set(gca, 'XColor', 'none')
    % set(gca, 'YColor', 'none')
    % text(-50+time(startF), 0.5, 'Reward', 'FontSize', 20, 'HorizontalAlignment', 'right')
    % 
    % subplot(10,2,[12, 14, 16])
    % plot(time(startF:i), -tongueTrajsWhole(startF:i, 2), 'LineWidth', 3, 'Color', 'k');
    % xlim([time(startF) time(endF)])
    % ylim([-350 -200])
    % set(gca, 'XColor', 'none')
    % set(gca, 'YColor', 'none')
    % text(-50+time(startF), -275, 'Tongue', 'FontSize', 20, 'HorizontalAlignment', 'right')

    frame = getframe(gcf);
    
 
  writeVideo(v,frame);  
 
end
close(v)
%%









%%





