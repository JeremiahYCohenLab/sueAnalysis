function photometryPreprocessing(session, regions, frameRate, varargin)
%% load data fromm file
% all time in s. 
p = inputParser;
p.addParameter('binSize', 100/1000);
p.addParameter('stepSize', 100/1000);
p.addParameter('tb', 2);
p.addParameter('tf', 5);
p.parse(varargin{:});

[root, sep] = currComputer();
pd = parseSessionString_df(session, root, sep);
dataPath = [pd.baseFolder 'photometry'  sep session 'neuron.csv'];
trialStartPath = [pd.baseFolder 'photometry'  sep session 'trialStart.csv'];
photometry_data = readmatrix(dataPath);
trialStart = readtable(trialStartPath);
trialStart = trialStart.Seconds;
% Hopkins frameRate = 40; Allen frameRate = 20;
% channel names
binSize = p.Results.binSize; % stepSize in ms
stepSize = p.Results.stepSize; % stepSize in ms
tb = p.Results.tb; % time before cue in second
tf = p.Results.tf; % time following cue in second 
midPoints = [0.5*binSize:stepSize:((tb+tf)-0.5*binSize)] - tb;
%% truncate after plugging in
truncateTimeStart = trialStart(1) - 2; % 2 seconds before the first trial 
truncateTimeEnd = trialStart(end) + 10;
photometry_data = photometry_data(photometry_data(:,2)>= truncateTimeStart & photometry_data(:,2) <= truncateTimeEnd,:);
photometry_data = photometry_data(1:2*floor(size(photometry_data,1)/2),:); % make sure frame number is even 
%% Photometry Data Extraction and Normalization
sigCode = [18, 274, 530, 786];
isoCode = [17, 273, 529, 785];
frameType = photometry_data(:,3);
if rem(length(frameType), 2)==1
    frameType(end)=0;
end
raw_gcamptimestamps = photometry_data(:,2);
raw_gcampdataR1 = photometry_data(:,4);
raw_gcampdataR2 = photometry_data(:,5);
raw_R1_sig = raw_gcampdataR1(ismember(frameType, sigCode));
raw_R1_ref = raw_gcampdataR1(ismember(frameType, isoCode));
raw_R2_sig = raw_gcampdataR2(ismember(frameType, sigCode));
raw_R2_ref = raw_gcampdataR2(ismember(frameType, isoCode));
if (length(raw_R1_sig) + length(raw_R2_sig)) ~= length(raw_gcampdataR1)
    fprintf('mistake in frame alignment');
    return
end
% Timestamps for gcampdata
timestampsSig = raw_gcamptimestamps(ismember(frameType, sigCode));
timestampsRef = raw_gcamptimestamps(ismember(frameType, isoCode));
% plot raw data
figure;
subplot(4,4,1:3); hold on;
plot(timestampsSig, raw_R1_sig);
plot(timestampsRef, raw_R1_ref);
legend({'sig', 'ref'})
title([regions{1} 'raw'])
subplot(4,4,5:7); hold on;
plot(timestampsSig, raw_R2_sig);
plot(timestampsRef, raw_R2_ref);
title([regions{2}, 'raw'])
xlabel('time in s')

%% Denoising the signal 
dn2_R1_sig = denoising(raw_R1_sig, frameRate);
dn2_R1_ref = denoising(raw_R1_ref, frameRate);
dn2_R2_sig = denoising(raw_R2_sig, frameRate);
dn2_R2_ref = denoising(raw_R2_ref, frameRate);
%% plotting denoised
% plot raw data
% figure;
% subplot(2,1,1); hold on;
% plot(timestampsSig, dn2_R1_sig);
% plot(timestampsRef, dn2_R1_ref);
% title(regions{1})
% subplot(2,1,2); hold on;
% plot(timestampsSig, dn2_R2_sig);
% plot(timestampsRef, dn2_R2_ref);
% title(regions{2})
% xlabel('time in s')
% sgtitle('denoised data')
%% montion correction
fit = fitlm(zscore(dn2_R1_ref), zscore(dn2_R1_sig));
mc_R1_sig = zscore(dn2_R1_sig) - (zscore(dn2_R1_ref)*fit.Coefficients.Estimate(2) +  fit.Coefficients.Estimate(1));
fit = fitlm(zscore(dn2_R2_ref), zscore(dn2_R2_sig));
mc_R2_sig = zscore(dn2_R2_sig) - (zscore(dn2_R2_ref)*fit.Coefficients.Estimate(2) +  fit.Coefficients.Estimate(1));
% %% no correction
% mc_R1_sig = dn2_R1_sig;
% mc_R2_sig = dn2_R2_sig;
%% plot montion corrected data
% figure;
% subplot(2,1,1); hold on;
% plot(timestampsSig, mc_R1_sig);
% title(regions{1})
% subplot(2,1,2); hold on;
% plot(timestampsSig, mc_R2_sig);
% title(regions{2})
% xlabel('time in s')
% sgtitle('montion corrected data')


%% Obtain Baseline Trace
%estimate 10th percentile of F using 15 sec. sliding window 
%this is to generate baseline fluorescence estimate
SR = 30;
win = SR * 20;
p = 10;
winright = 0;
[ ys ] = running_percentile_filter(mc_R1_sig, win, winright, p); 
delta_F_R1 = mc_R1_sig - ys;
[ ys ] = running_percentile_filter(mc_R2_sig, win, winright, p); 
delta_F_R2 = mc_R2_sig - ys;
%%
subplot(4,4,9:11); hold on;
plot(timestampsSig, delta_F_R1, 'k', 'LineWidth',1.5);
title([regions{1} 'deltaF'])
subplot(4,4,13:15); hold on;
plot(timestampsSig, delta_F_R2, 'LineWidth',1.5);
title([regions{2}  'deltaF'])
xlabel('time in s')
% set(gca, 'XColor', 'none');
% set(gca, 'YColor', 'none');
%% find prominent LC peaks
% fc = 0.2; %cutoff frequency in Hz (anything lower than this passes)
% fs = 30; %sampling rate in Hz
% [b,a] = butter(2,fc/(fs/2),'low');
% sm_R1_sig = filtfilt(b,a,delta_F_R1);
%% detect Peaks in LC
if strcmp(regions{1}, 'LC')
    peakInds = peakseek(delta_F_R1,200,2.5);
else
    peakInds = peakseek(delta_F_R2,200,2.5);
end

timeB = 4;  % in s
timeF = 5;  % in s
G0max = NaN(length(peakInds), (timeB + timeF)*SR+1);
G1max = NaN(length(peakInds), (timeB + timeF)*SR+1);
for i = 1:length(peakInds)
    startTS = max([1, peakInds(i) - timeB * SR]);
    endTS = min([peakInds(i) + timeF * SR, length(delta_F_R1)]);
    G0max(i,(timeB*SR-(peakInds(i)-startTS)+1):end-(timeF*SR-(endTS-peakInds(i)))) = delta_F_R1(startTS:endTS);
    G1max(i,(timeB*SR-(peakInds(i)-startTS)+1):end-(timeF*SR-(endTS-peakInds(i)))) = delta_F_R2(startTS:endTS);
end
%%
subplot(4,4,[4, 8]);hold on; 
time = linspace(-timeB, timeF, (timeB + timeF)*SR+1);
plotFilled(time, G0max, 'k');
plot([0 0], [0 4*10^-4], 'LineStyle', '--', 'Color', [0.7 0.7 0.7]);
% ylim([0 4*10^-4]);
% plot([-4 -2]', [0.3*10^-4 0.3*10^-4]', 'k', 'LineWidth', 2)
% text(-3.5, 0.2*10^-4, '2s', 'FontSize', 12)
% set(gca, 'XColor', 'none');
% set(gca, 'YColor', 'none');
title(regions{1}, 'FontSize', 15)
subplot(4,4,[12 16]);hold on
time = linspace(-timeB, timeF, (timeB + timeF)*SR+1);
plotFilled(time, G1max, [0 0.2 0.9]);
plot([0 0], [1*10^-5 6*10^-5], 'LineStyle', '--', 'Color', [0.7 0.7 0.7]);
% set(gca, 'XColor', 'none');
% set(gca, 'YColor', 'none');
title(regions{2}, 'FontSize', 15)
sgtitle([session 'preprocessing'])

screen = get(0,'Screensize');
screen(4) = screen(4) - 100;
set(gcf, 'Position', screen)

savePath = pd.saveFigFolder;
if isempty(dir(savePath))
    mkdir(savePath)
end
saveFigurePDF(gcf,[savePath session '_photometryPre.pdf'])
%% Extract behavior timestamps from analog_inputs
mPFCmat = zeros(length(trialStart), length(midPoints));
LCmat = zeros(length(trialStart), length(midPoints));

for i = 1:length(trialStart)
    tmpTime = timestampsSig(timestampsSig > (trialStart(i)-tb) & timestampsSig <= (trialStart(i)+tf));
    tmp1 = delta_F_R1(timestampsSig > (trialStart(i)-tb) & timestampsSig <= (trialStart(i)+tf));
    tmp2 = delta_F_R2(timestampsSig > (trialStart(i)-tb) & timestampsSig <= (trialStart(i)+tf));
    for j = 1:length(midPoints)
        if strcmp(regions{1}, 'mPFC')
            mPFCmat(i,j) = mean(tmp1(tmpTime>trialStart(i)+midPoints(j)-0.5*binSize & tmpTime<=trialStart(i)+midPoints(j)+0.5*binSize));
            LCmat(i,j) = mean(tmp2(tmpTime>trialStart(i)+midPoints(j)-0.5*binSize & tmpTime<=trialStart(i)+midPoints(j)+0.5*binSize));
        else
            LCmat(i,j) = mean(tmp1(tmpTime>trialStart(i)+midPoints(j)-0.5*binSize & tmpTime<=trialStart(i)+midPoints(j)+0.5*binSize));
            mPFCmat(i,j) = mean(tmp2(tmpTime>trialStart(i)+midPoints(j)-0.5*binSize & tmpTime<=trialStart(i)+midPoints(j)+0.5*binSize));
        end
    end
end
%% 
os = behAnalysisNoPlot_opMD(session);
mPFCmatChoice = zeros(length(os.responseInds), length(midPoints));
LCmatChoice = zeros(length(os.responseInds), length(midPoints));

for i = 1:length(os.responseInds)
    currTime = trialStart(os.responseInds(i)) + 1/1000 * os.lickLat(i);
    tmpTime = timestampsSig(timestampsSig > (currTime-tb) & timestampsSig <= (currTime+tf));
    tmp1 = delta_F_R1(timestampsSig > (currTime-tb) & timestampsSig <= (currTime+tf));
    tmp2 = delta_F_R2(timestampsSig > (currTime-tb) & timestampsSig <= (currTime+tf));
    for j = 1:length(midPoints)
        if strcmp(regions{1}, 'mPFC')
            mPFCmatChoice(i,j) = mean(tmp1(tmpTime>currTime+midPoints(j)-0.5*binSize & tmpTime<=currTime+midPoints(j)+0.5*binSize));
            LCmatChoice(i,j) = mean(tmp2(tmpTime>currTime+midPoints(j)-0.5*binSize & tmpTime<=currTime+midPoints(j)+0.5*binSize));
        else
            LCmatChoice(i,j) = mean(tmp1(tmpTime>currTime+midPoints(j)-0.5*binSize & tmpTime<=currTime+midPoints(j)+0.5*binSize));
            mPFCmatChoice(i,j) = mean(tmp2(tmpTime>currTime+midPoints(j)-0.5*binSize & tmpTime<=currTime+midPoints(j)+0.5*binSize));
        end
    end
end

%% save
if ~exist(pd.sortedFolder, 'dir')
    mkdir(pd.sortedFolder)
end
save([pd.sortedFolder session '_photometry.mat'], 'mPFCmat', 'LCmat', 'mPFCmatChoice', 'LCmatChoice', 'tb', 'tf', 'frameRate', 'binSize', 'stepSize', 'midPoints');




