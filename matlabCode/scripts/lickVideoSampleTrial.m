videoName = 'F:\lickSampleVidoes\lick1\21-20-46.328';
dlcName = 'DLC_resnet50_lick1Apr28shuffle1_300000';
rawPosition = csvread([videoName dlcName '.csv'],3,0);
FR = 400; 
thresh = 0.8;
time = -500 + 1000/FR*[1:size(rawPosition,1)];
timeLED = rawPosition(:, 34);
confTongue = rawPosition(:, 31);
xTongue = rawPosition(:,29);
xTongue(confTongue<thresh) = NaN; 
yTongue = rawPosition(:,30);
yTongue(confTongue<thresh) = NaN; 
xPortL = rawPosition(:, 17);
yPortL = rawPosition(:, 18);
xPortR = rawPosition(:, 20);
yPortR = rawPosition(:, 21);
%% time LED ll
figure2;
plot(time, timeLED);
ylim([0 1.2])
title('time')
xlabel('time (ms)')
%% correct for time
ledThresh = 0.95;
ledOn = find(timeLED>=ledThresh, 1);
timeCorrected = 1000/FR*([1:size(rawPosition,1)] - ledOn);
figure2;
plot(timeCorrected, timeLED);
ylim([0 1.2])
title('corrected time')
xlabel('time (ms)')
%% calculate scale bar
% LR port distance
L = [mean(xPortL(timeCorrected<0)), mean(yPortL(timeCorrected<0))];
R = [mean(xPortR(timeCorrected<0)), mean(yPortR(timeCorrected<0))];
disRL = norm(L-R);
%% plot events
figure2;
subplot(4,1,1);
plot(timeCorrected, timeLED>=ledThresh, 'Color', 'k', 'LineWidth', 2);
set(gca, 'TickDir', 'out')
set(gca, 'box', 'off')

ylabel('cue')
subplot(4,1,2);
plot(timeCorrected, confTongue>=thresh, 'Color', 'k', 'LineWidth', 2);
ylabel('tongue')
set(gca, 'TickDir', 'out')
set(gca, 'box', 'off')

Lconf = rawPosition(:,19);
yPortL = interp1(timeCorrected(Lconf>ledThresh), yPortL(Lconf>ledThresh), timeCorrected);
xPortL = interp1(timeCorrected(Lconf>ledThresh), xPortL(Lconf>ledThresh), timeCorrected);
Rconf = rawPosition(:,22);
yPortR = interp1(timeCorrected(Rconf>ledThresh), yPortR(Rconf>ledThresh), timeCorrected);
xPortR = interp1(timeCorrected(Rconf>ledThresh), yPortR(Rconf>ledThresh), timeCorrected);

subplot(4,1,3)
plot(timeCorrected, -4/disRL*(yPortL-L(2)), 'Color', 'k', 'LineWidth', 2);
ylabel('Lport (mm)')
ylim([-0.2 5])
set(gca, 'TickDir', 'out')
set(gca, 'box', 'off')

subplot(4,1,4);
plot(timeCorrected, -4/disRL*(yPortR-R(2)), 'Color', 'k', 'LineWidth', 2);
ylabel('Rport (mm)')
ylim([-0.2 5])
set(gca, 'TickDir', 'out')
set(gca, 'box', 'off')

%% conf licks x axis
Lconf = rawPosition(:,19);
L = [mean(rawPosition(timeCorrected'<=0 & Lconf>=thresh, 17)), mean(rawPosition(timeCorrected'<=0 & Lconf>=thresh, 18))];
Rconf = rawPosition(:,22);
R = [mean(rawPosition(timeCorrected'<=0 & Rconf>=thresh, 20)), mean(rawPosition(timeCorrected'<=0 & Rconf>=thresh, 21))];
midPoint = 0.5*(L + R);
figure2;
subplot(1,3,1); hold on;
plot(timeCorrected, xTongue);
plot(timeCorrected, mean([L(1), R(1)])* ones(size(timeCorrected)), 'LineStyle', '--', 'Color','r' )
title('x position')
xlabel('time (ms)')
subplot(1,3,2); hold on;
plot(timeCorrected, yTongue);
plot(timeCorrected, mean([L(2), R(2)])* ones(size(timeCorrected)), 'LineStyle', '--', 'Color', 'r' )
title('y position')
xlabel('time (ms)')
subplot(1,3,3);
plot(timeCorrected, confTongue);
title('lick confidence')
xlabel('time (ms)')
%% plot all licks
% split tongue movement into licks 
nanInd = ~isnan(xTongue);
startPoints = find(nanInd(1:end-1)== 0 & nanInd(2:end)==1) + 1;
endPoints = find(nanInd(1:end-1)== 1 & nanInd(2:end)==0);
lickX = {};
lickY = {};
lickTime = {};
for i = 1:length(startPoints)
    if (endPoints(i)-startPoints(i))*1000/FR >= 20
        lickX = [lickX, xTongue(startPoints(i):endPoints(i))];
        lickY = [lickY, yTongue(startPoints(i):endPoints(i))];
        lickTime = [lickTime, timeCorrected(startPoints(i):endPoints(i))];
    end
end
figure2; hold on;
colors = cool(length(lickX));
for i = 1:length(lickX)
    plot(lickTime{i}, lickY{i}, 'color', colors(i,:), 'LineWidth', 2);
    
end
plot([min(lickTime{1}) max(lickTime{end})+50*(length(lickX)-1)], [mean([L(2), R(2)]) mean([L(2), R(2)])], 'LineStyle', '--', 'Color', 'r' )
set(gca, 'YDir','reverse')
%%
figure2; hold on;
colors = cool(length(lickX));
for i = 1:length(lickX)
    plot(lickX{i} + 50*(i -1), lickY{i}, 'color', colors(i,:), 'LineWidth', 2);
    
end
plot([min(lickX{1}) max(lickX{1})], [mean([L(2), R(2)]) mean([L(2), R(2)])], 'LineStyle', '--', 'Color', 'r' )
plot([mean([L(1), R(1)]) mean([L(1), R(1)])], [min(lickY{1}) max(lickY{1})], 'LineStyle', '--', 'Color', 'r' )

%%
v = VideoReader([videoName '.avi']);
f = read(v, mean([startPoints(find(startPoints>ledOn, 1)), endPoints(find(startPoints>ledOn, 1))]-3));
figure2; 
imshow(f);hold on;
plot([100 100+2*disRL/4], [50 50], 'Color', [1 1 1], 'LineWidth', 2)
text(100, 25, '2mm', 'Color', [1 1 1], 'FontSize', 12)
%%
hold on;
focusInd = 8;
for i = 1:length(lickX)
    if i ~= focusInd
        plot(lickX{i}, lickY{i}, 'color', [0 0 0], 'LineWidth', 1);
    end
end
    plot(lickX{focusInd}, lickY{focusInd}, 'color', colors(focusInd,:), 'LineWidth', 2);

%% filter location data
fc = 50; %cutoff frequency in Hz (anything lower than this passes)
fs = FR; %sampling rate in Hz
[b,a] = butter(2,fc/(fs/2),'low');
lickXFilt = [];
lickYFilt = [];

for i = 1:length(lickX)
    lickXFilt{i} = 4/disRL * (filtfilt(b,a, lickX{i}) - midPoint(1));
    lickYFilt{i} = -4/disRL * (filtfilt(b,a, lickY{i}) - midPoint(2));
end
%% calculate speed data
lickSpeed = [];
lickTimeSpeed = [];

for i = 1:length(lickX)
    % time
    timeTemp = lickTime{i};
    lickTimeSpeed{i} = 0.5 * (timeTemp(1:end-1) + timeTemp(2:end));
    
    lickDisTmpX = diff(lickXFilt{i});
    lickDisTmpY = diff(lickYFilt{i});
    lickSpeed{i} = sqrt(lickDisTmpX.^2 + lickDisTmpY.^2)/(1000/FR);
end
%% calculate acc
lickAccX = [];
lickAccY = [];
lickVY = {};
lickTimeAcc = [];

for i = 1:length(lickX)
    % time
    timeTemp = lickTimeSpeed{i};
    lickTimeAcc{i} = 0.5 * (timeTemp(1:end-1) + timeTemp(2:end));
    
    lickDisTmpX = diff(lickXFilt{i});
    lickVX = lickDisTmpX/(1000/FR); % mm/ms
    lickAccX{i} = diff(lickVX)/(1000/FR); % mm/ms^2
    
    
    lickDisTmpY = diff(lickYFilt{i});
    lickVY{i} = lickDisTmpY/(1000/FR); % mm/ms
    lickAccY{i} = diff(lickVY{i})/(1000/FR); % mm/ms^2
end
%% plot y location
figure; hold on;
colors = cool(length(lickXFilt));
for i = 1:length(lickXFilt)
    plot(lickTime{i}, lickYFilt{i}, 'color', colors(i,:), 'LineWidth', 2);
end
plot([min(lickTime{1}) max(lickTime{end})], [0 0], 'LineStyle', '--', 'Color', 'r' )
% set(gca, 'YDir','reverse')
title('Filter y location') 
ylabel('distance')
%% plot speed
figure; hold on;
colors = cool(length(lickSpeed));
for i = 1:length(lickSpeed)
    plot(lickTimeSpeed{i}, lickSpeed{i}, 'color', colors(i,:), 'LineWidth', 2);
end
plot([min(lickTimeSpeed{1}) max(lickTimeSpeed{end})], [0 0], 'LineStyle', '--', 'Color', 'r' )
% set(gca, 'YDir','reverse')
title('Speed mm/ms') 
ylabel('Speed mm/ms')
xlabel('Time - go cue (ms)')
%% plot velocity on Y
figure; hold on;
colors = cool(length(lickVY));
for i = 1:length(lickVY)
    plot(lickTimeSpeed{i}, lickVY{i}, 'color', colors(i,:), 'LineWidth', 2);
end
plot([min(lickTimeSpeed{1}) max(lickTimeSpeed{end})], [0 0], 'LineStyle', '--', 'Color', 'r' )
% set(gca, 'YDir','reverse')
title('velocity mm/ms') 
ylabel('velocity mm/ms')
xlabel('Time - go cue (ms)')
%% plot acc Y
figure; hold on;
colors = cool(length(lickAccY));
for i = 1:length(lickAccY)
    plot(lickTimeAcc{i}, lickAccY{i}, 'color', colors(i,:), 'LineWidth', 2);
end
plot([min(lickTimeAcc{1}) max(lickTimeAcc{end})], [0 0], 'LineStyle', '--', 'Color', 'r' )
% set(gca, 'YDir','reverse')
title('AccY mm/ms^2') 
ylabel('AccY mm/ms^2')
xlabel('Time - go cue (ms)')
 %%
figure2; hold on;
colors = cool(length(lickX));
for i = 1:length(lickX)
    plot(lickXFilt{i} + 0.7*(i-1), lickYFilt{i}, 'color', colors(i,:), 'LineWidth', 2);
    scatter(lickXFilt{i}(1) + 0.7*(i-1), lickYFilt{i}(1), 12, 'r', 'filled')
end
plot([min(lickXFilt{1}) max(lickXFilt{end})+0.5*(length(lickX)-1)], [0 0], 'LineStyle', '--', 'Color', 'r' )
% set(gca, 'YDir','reverse')
title('Filter xy location') 
% plot([0 1], [-1,5 -1.5], 'k')
% plot([250 250], [300 300-1*disRL/4], 'k')
% % text(250, 270, '1 mm', 'FontSize', 12)
xlabel('xlocation')
ylabel('ylocation')
%%
figure2; hold on;
colors = cool(length(lickX));
for i = length(lickX):-1:1
    plot(lickXFilt{i}, lickYFilt{i}, 'color', colors(i,:), 'LineWidth', 2);
   scatter(lickXFilt{i}(1), lickYFilt{i}(1), 20, 'r', 'filled')
end
plot([-2 2], [0 0], 'LineStyle', '--', 'Color', 'r' )
plot([0 0], [-2 1], 'LineStyle', '--', 'Color', 'r' )
%%
figure2; hold on;
xLoc = 4/disRL * ((1:size(f, 2)) - midPoint(1));
yLoc = -4/disRL * ((1:size(f, 1)) - midPoint(2));
image(xLoc, yLoc, f);
colormap gray
hold on;
focusInd = 1;
for i = 1:length(lickX)
    if i ~= focusInd
        plot(lickXFilt{i}, lickYFilt{i}, 'color',  colors(i,:), 'LineWidth', 1.5, 'LineStyle', ':');
    end
end
plot(lickXFilt{focusInd}, lickYFilt{focusInd}, 'color', colors(focusInd,:), 'LineWidth', 2);
scatter(lickXFilt{focusInd}(1), lickYFilt{focusInd}(1), 20, 'r', 'filled')
plot([100 100+2*disRL/4], [50 50], 'Color', [1 1 1], 'LineWidth', 2)
text(100, 25, '2mm', 'Color', [1 1 1], 'FontSize', 12)
%%