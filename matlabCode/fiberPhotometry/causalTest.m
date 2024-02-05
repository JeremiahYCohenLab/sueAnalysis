%% pupil dynamics
clear all
animal = '684890';
test = 'test2';
[root, sep] = currComputer();
%%
fileDir = [root animal '\' animal 'optoTest\20231222\' test];
allFiles = dir(fileDir);
ind = contains({allFiles.name}, '000.csv');

fileDir = [root animal '\' animal 'optoTest\20231222\' test sep];
allFiles = dir(fileDir);
optoFolder = {allFiles([allFiles.isdir]==1).name};
optoFolder = optoFolder{cellfun(@length, optoFolder)>5 & ~contains(optoFolder, 'spon')};
fpDir = [fileDir optoFolder sep];
%%
a = readmatrix([fileDir sep allFiles(ind).name]);
% convolution
LL = a(:,16);
figure2; hold on;
plot(LL)
kernel = ones(10,1);
accuSig = conv(LL, kernel);
accuSig = accuSig(length(kernel):end);
% denoise (erosion)
LL(accuSig<=8) = 0;
plot(LL)
%%
% pupil timing
upThresh = 0.95;
downThresh = 0.3;
upFrames = find(LL(2:end)>=upThresh & LL(1:end-1)<upThresh)+1; 
%%
downFrames = find(LL(2:end)<=downThresh & LL(1:end-1)>downThresh)+1; 
frameRate = mean(diff(upFrames))/20;   % frames/s
%% signal
dia = sqrt((a(:,2) - a(:,5)).^2 + (a(:,3) - a(:,6)).^2);
% filtering low pass
fcut = 2;
[b, a] = butter(2, fcut/(frameRate/2), 'low');
dia = filtfilt(b, a, dia);
diaZ = zscore(dia);
% cut into trials
time = [-2 10];
frameNum = round(frameRate*(time(2) - time(1)));
frameNumPre = round(frameRate*(0 - time(1)));
timeStamps = linspace(time(1), time(2), frameNum);

diaMat = zeros(length(upFrames), frameNum);
diaMatZ = zeros(length(upFrames), frameNum);
for i = 1:length(upFrames)
    diaMat(i,:) = dia(upFrames(i)-frameNumPre: upFrames(i)-frameNumPre+frameNum-1)';
    diaMatZ(i,:) = diaZ(upFrames(i)-frameNumPre: upFrames(i)-frameNumPre+frameNum-1)';
end

diaBl = mean(diaMat(:, 1:frameNumPre), 2);
dilation = diaMat./diaBl - 1;

save([fileDir sep animal 'optoTestPupil.mat'], 'dia', 'LL', 'diaMat', 'diaMatZ', 'dilation', 'diaBl', 'timeStamps');
%% plot pupil response
figPupil = figure2;
sgtitle([animal ' ' test ' pupil'])
subplot(1,3,1)
plotFilled(timeStamps, diaMat, 'b');
title([animal ' ' test ' raw'])
  
subplot(1,3,2)
plotFilled(timeStamps, diaMatZ, 'b');
title([animal ' ' test ' z'])

subplot(1,3,3)
plotFilled(timeStamps, dilation, 'b');   
title([animal ' ' test ' dilation'])

screen = get(0,'Screensize');
screen(4) = screen(4) - 100;
set(figPupil, 'Position', screen);
saveFigurePDF(figPupil, [fpDir animal test ' pupil.pdf']);

%% fiberphotometry
tb = 2;
tf = 5;
samplingFreq = 20;


allFiles = dir(fpDir);
allFiles = {allFiles([allFiles.isdir]==0).name}';

% load signal
GInd = contains(allFiles, 'FIP_DataG');
IsoInd = contains(allFiles, 'FIP_DataIso');
RInd =  contains(allFiles, 'FIP_DataR');
numChannels = 3;

GSig = readmatrix([fpDir allFiles{GInd}]);
timeStampsG = GSig(:,1);
GSig = GSig(:, 2:2+numChannels-1);

IsoSig = readmatrix([fpDir allFiles{IsoInd}]);
timeStampsIso = IsoSig(:,1);
IsoSig = IsoSig(:, 2:2+numChannels-1);

RSig = readmatrix([fpDir allFiles{RInd}]);
timeStampsR = RSig(:,1);
RSig = RSig(:, 3);
%% detect trial starts
timeFIP = timeStampsG;

thresh = 200;
upInds = find(RSig(1:end-1)<thresh & RSig(2:end)>=thresh)+1;
upInds = upInds([1; find(diff(upInds) > 100)+1]);
upTimes = timeFIP(upInds);

GSig = GSig(timeFIP >= (timeFIP(upInds(1)) - (tb*1000+2000)) & timeFIP < (timeFIP(upInds(end)) + tf*1000+2000),:);
IsoSig = IsoSig(timeFIP >= (timeFIP(upInds(1)) - (tb*1000+2000)) & timeFIP < (timeFIP(upInds(end)) + tf*1000+2000),:);
RSig = RSig(timeFIP >= (timeFIP(upInds(1)) - (tb*1000+2000)) & timeFIP < (timeFIP(upInds(end)) + tf*1000+2000),:);
timeFIP = timeFIP(timeFIP >= (timeFIP(upInds(1)) - (tb*1000+2000)) & timeFIP < (timeFIP(upInds(end)) + tf*1000+2000),:);

timeFIP = timeFIP - upTimes(1);
upTimes = upTimes - upTimes(1);

upInds = find(RSig(1:end-1)<thresh & RSig(2:end)>=thresh)+1;
upInds = upInds([1; find(diff(upInds) > 100)+1]);

L = min(RSig);
H = max(RSig);

figure2; hold on;
plot(timeFIP, RSig);
plot([upTimes upTimes]', [L*ones(length(upTimes),1) H*ones(length(upTimes),1)]', 'Color', 'r', 'LineWidth', 2, 'LineStyle', '-', 'Marker', 'none')
%% preprocessing
% truncate
% use timeStampsG as default time

dFF = cell(3,1);
winLeft = samplingFreq * 20; % 20s before sample
p = 10; % percentage
winRight = 0;
for i = 1:numChannels
    GSig(:,i) = denoising(GSig(:,i), samplingFreq, 9);
    IsoSig(:,i) = denoising(IsoSig(:,i), samplingFreq, 9);
    fit = fitlm(zscore(IsoSig(:,i)), zscore(GSig(:,i)));
    dFF{i} = zscore(GSig(:,i)) - (zscore(IsoSig(:,i))*fit.Coefficients.Estimate(2) +  fit.Coefficients.Estimate(1));
%     baseline = running_percentile_filter(dFF{i}, winLeft, winRight, p); 
%     dFF{i} = dFF{i} - baseline;
end
%% align signal to triggers
dFFTrigger = cell(numChannels, 1);

for i = 1:numChannels
    currSig = dFF{i};
    sigMatTmp = [];
    for j = 1:length(upTimes)
        tmp = currSig((upInds(j) - tb*samplingFreq) : (upInds(j)+tf*samplingFreq))';
        sigMatTmp = [sigMatTmp; tmp];
    end
    dFFTrigger{i} = sigMatTmp;
end
timeStamps = linspace(-tb*1000, tf*1000, (tb+tf)*samplingFreq+1);
%%
figure2;
sgtitle([animal ' ' test])
for i = 1:2
    subplot(4, 1, i);
    hold on;
    yyaxis left;
    set(gca, 'YColor', 'r');
    plot(timeFIP, RSig, 'LineWidth', 0.5, 'Color', [1 0.5 0.5]);
    L = min(RSig);
    H = max(RSig);
    plot([upTimes upTimes]', [L*ones(length(upTimes),1) H*ones(length(upTimes),1)]', 'Color', 'r', 'LineWidth', 2, 'LineStyle', '-', 'Marker', 'none')
    yyaxis right;
    plot(timeFIP, dFF{i}, 'LineWidth', 1.5, 'Color', 'b');
    title(['Fiber' num2str(i)])
    set(gca, 'YColor', 'b');
    ylabel('dF/F')
    xlabel('time (ms)')
    xlim([-tb*1000 upTimes(end) + 15000]);
end

for i = 1:2
    subplot(2, 4, 2*(i-1)+1+4);
    hold on;
    % plot(timeStamps, dFFTrigger{i}, 'Color', [0.5 0.5 1]);
    plotFilled(timeStamps, dFFTrigger{i}, 'b');
    patch([0 2000 2000 0], [min(dFF{i}) min(dFF{i}) max(dFF{i}) max(dFF{i})], 'r', 'edgeColor', 'none', 'FaceAlpha', 0.5);
    title(['Fiber' num2str(i)])
    ylabel('dF/F')
    xlabel('time (ms)')
    
    subplot(2, 4, 2*(i-1)+2+4);
    hold on;
    bl = dFFTrigger{i};
    bl = bl(:,timeStamps<0);
    bl = mean(bl, 2);
    sig = dFFTrigger{i} - bl;
    plot(timeStamps, sig, 'Color', [0.5 0.5 1]);
    plotFilled(timeStamps, sig, 'b');
    patch([0 2000 2000 0], [min(sig, [], 'all') min(sig, [], 'all') max(sig, [], 'all') max(sig, [], 'all')], 'r', 'edgeColor', 'none', 'FaceAlpha', 0.5);
    title(['Fiber' num2str(i)])
    ylabel('dF/F - bl')
    xlabel('time (ms)')
end
screen = get(0,'Screensize');
screen(4) = screen(4) - 100;
set(gcf, 'Position', screen);
saveFigurePDF(gcf, [fpDir animal test 'photometry.pdf']);
%%





