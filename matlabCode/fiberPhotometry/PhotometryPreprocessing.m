%% load data fromm file
session = 'mZS079d20220307';
[root, sep] = currComputer();
pd = parseSessionString_df(session, root, sep);
dataPath = [pd.baseFolder 'photometry'  sep session];
photometry_data = readmatrix(dataPath);
%% truncate after plugging in
truncateFrame = 60 * (3169- photometry_data(1,2))+ 1; % need to be odd number to keep all frames intact
photometry_data = photometry_data(photometry_data(:,1)>= truncateFrame,:);
photometry_data = photometry_data(1:2*floor(size(photometry_data,1)/2),:); % make sure everything is even 
%% Photometry Data Extraction and Normalization
sigCode = 18;
isoCode = 17;
frameType = photometry_data(:,3);
if rem(length(frameType), 2)==1
    frameType(end)=0;
end
raw_gcamptimestamps = photometry_data(:,2);
raw_gcampdataR1 = photometry_data(:,4);
raw_gcampdataR2 = photometry_data(:,5);
raw_R1_sig = raw_gcampdataR1(frameType == sigCode);
raw_R1_ref = raw_gcampdataR1(frameType == isoCode);
raw_R2_sig = raw_gcampdataR2(frameType == sigCode);
raw_R2_ref = raw_gcampdataR2(frameType == isoCode);


% Timestamps for gcampdata
timestampsSig = raw_gcamptimestamps(frameType == sigCode);
timestampsRef = raw_gcamptimestamps(frameType == isoCode);
%% plot raw data
figure;
subplot(2,1,1); hold on;
plot(timestampsSig, raw_R1_sig);
plot(timestampsRef, raw_R1_ref);
title('LC')
subplot(2,1,2); hold on;
plot(timestampsSig, raw_R2_sig);
plot(timestampsRef, raw_R2_ref);
title('mPFC')
xlabel('time in s')
sgtitle('raw data')

%% Denoising the signal 
dn2_R1_sig = denoising(raw_R1_sig);
dn2_R1_ref = denoising(raw_R1_ref);
dn2_R2_sig = denoising(raw_R2_sig);
dn2_R2_ref = denoising(raw_R2_ref);
%% plotting denoised
%% plot raw data
figure;
subplot(2,1,1); hold on;
plot(timestampsSig, dn2_R1_sig);
plot(timestampsRef, dn2_R1_ref);
title('LC')
subplot(2,1,2); hold on;
plot(timestampsSig, dn2_R2_sig);
plot(timestampsRef, dn2_R2_ref);
title('mPFC')
xlabel('time in s')
sgtitle('denoised data')
%% montion correction
fit = fitlm(dn2_R1_ref, dn2_R1_sig);
mc_R1_sig = dn2_R1_sig - (dn2_R1_ref*fit.Coefficients.Estimate(2) +  fit.Coefficients.Estimate(1));
fit = fitlm(dn2_R2_ref, dn2_R2_sig);
mc_R2_sig = dn2_R2_sig - (dn2_R2_ref*fit.Coefficients.Estimate(2) +  fit.Coefficients.Estimate(1));
%% no correction
mc_R1_sig = dn2_R1_sig;
mc_R2_sig = dn2_R2_sig;
%% plot montion corrected data
figure;
subplot(2,1,1); hold on;
plot(timestampsSig, mc_R1_sig);
title('LC')
subplot(2,1,2); hold on;
plot(timestampsSig, mc_R2_sig);
title('mPFC')
xlabel('time in s')
sgtitle('montion corrected data')


%% Obtain Baseline Trace
%estimate 10th percentile of F using 15 sec. sliding window 
%this is to generate baseline fluorescence estimate
SR = 30;
win = SR * 20;
p = 10;
% [ ys ] = running_percentile(gcampdata, win, p);
%if you want to run a left-sided sliding window 
winright = 0;
[ ys ] = running_percentile_filter(mc_R1_sig, win, winright, p); 
delta_F_R1 = mc_R1_sig - ys;
delta_FoverF_R1 = delta_F_R1 ./ ys ;
[ ys ] = running_percentile_filter(mc_R2_sig, win, winright, p); 
delta_F_R2 = mc_R2_sig - ys;
delta_FoverF_R2 = delta_F_R2 ./ ys ;
%exclusion criterion: only sessions where 97.5% of DF/F0 across the entire session
%exceeded 1% are included
delta_F_p = delta_FoverF .* 100; %percent change from baseline fluorescence
exc = prctile(delta_F_p, 97.5);
pass = (exc > 1);
%delta_F_z = zscore(delta_F,0,'all'); %%WHY NEED THIS??
%%
figure;
subplot(2,1,1); hold on;
plot(timestampsSig, delta_F_R1, 'k', 'LineWidth',1.5);
% plot([2300 2350], [-0.00015 -0.00015], 'k', 'LineWidth',1.5)
% text(2320, -0.0002, '50 s')
% set(gca, 'XColor', 'none');
% set(gca, 'YColor', 'none');
title('LC')
subplot(2,1,2); hold on;
plot(timestampsSig, delta_F_R2, 'LineWidth',1.5);
title('mPFC')
xlabel('time in s')
% set(gca, 'XColor', 'none');
% set(gca, 'YColor', 'none');
sgtitle('deltaF')
%%
figure;
subplot(2,1,1); hold on;
plot(timestampsSig, delta_FoverF_R1);
title('LC')
subplot(2,1,2); hold on;
plot(timestampsSig, delta_FoverF_R2);
title('mPFC')
xlabel('time in s')
sgtitle('deltaFoverF')
%% find prominent LC peaks
fc = 0.2; %cutoff frequency in Hz (anything lower than this passes)
fs = 30; %sampling rate in Hz
[b,a] = butter(2,fc/(fs/2),'low');
sm_R1_sig = filtfilt(b,a,delta_F_R1);
%% detect Peaks in LC
peakInds = peakseek(delta_F_R1,200,0.5*10^-4);

timeB = 4;  % in s
timeF = 5;  % in s
PFCmax = NaN(length(peakInds), (timeB + timeF)*SR+1);
LCmax = NaN(length(peakInds), (timeB + timeF)*SR+1);
for i = 1:length(peakInds)
    startTS = max([1, peakInds(i) - timeB * SR]);
    endTS = min([peakInds(i) + timeF * SR, length(delta_F_R1)]);
    PFCmax(i,(timeB*SR-(peakInds(i)-startTS)+1):end-(timeF*SR-(endTS-peakInds(i)))) = delta_F_R2(startTS:endTS);
    LCmax(i,(timeB*SR-(peakInds(i)-startTS)+1):end-(timeF*SR-(endTS-peakInds(i)))) = delta_F_R1(startTS:endTS);
end
%%
figure2;
subplot(1,2,1);hold on
time = linspace(-timeB, timeF, (timeB + timeF)*SR+1);
plotFilled(time, LCmax, 'k');
plot([0 0], [0 4*10^-4], 'LineStyle', '--', 'Color', [0.7 0.7 0.7]);
ylim([0 4*10^-4]);
% plot([-4 -2]', [0.3*10^-4 0.3*10^-4]', 'k', 'LineWidth', 2)
% text(-3.5, 0.2*10^-4, '2s', 'FontSize', 12)
set(gca, 'XColor', 'none');
set(gca, 'YColor', 'none');
title('LC', 'FontSize', 15)
subplot(1,2,2);hold on
time = linspace(-timeB, timeF, (timeB + timeF)*SR+1);
plotFilled(time, PFCmax, [0 0.2 0.9]);
plot([0 0], [1*10^-5 6*10^-5], 'LineStyle', '--', 'Color', [0.7 0.7 0.7]);
set(gca, 'XColor', 'none');
set(gca, 'YColor', 'none');
title('PFC', 'FontSize', 15)

%% Extract behavior timestamps from analog_inputs
% %%
% % Subtract logicals from each other in LLP column
% switches_LLP = diff(beh_data(:,2)); %2
% % Subtract logicals from each other in RLP column
% switches_RLP = diff(beh_data(:,4)); %4
% % Subtract logicals from each other in HE column
% switches_HE = diff(beh_data(:,3)); %3
% % Subtract logicals from each other in lick column
% switches_lick = diff(beh_data(:,5)); %5
% 
% % When the difference is 1, that's event ON, find the indices for ON
% LLP_ON_index = find(switches_LLP == 1)+1;
% RLP_ON_index = find(switches_RLP == 1)+1;
% HE_ON_index = find(switches_HE == 1)+1;
% lick_ON_index = find(switches_lick == 1)+1;
% 
% % When the difference is -1, that's event OFF, find indices for OFF
% LLP_OFF_index = find(switches_LLP == -1)+1;
% RLP_OFF_index = find(switches_RLP == -1)+1;
% HE_OFF_index = find(switches_HE == -1)+1;
% lick_OFF_index = find(switches_lick == -1)+1;
% 
% % Get the time stamps from those event indices
% LLP_ON_timestamps = beh_data(LLP_ON_index, 1);
% LLP_OFF_timestamps = beh_data(LLP_OFF_index, 1);
% RLP_ON_timestamps = beh_data(RLP_ON_index, 1);
% RLP_OFF_timestamps = beh_data(RLP_OFF_index, 1);
% HE_ON_timestamps = beh_data(HE_ON_index, 1);
% HE_OFF_timestamps = beh_data(HE_OFF_index,1);
% lick_ON_timestamps = beh_data(lick_ON_index, 1);
% lick_OFF_timestamps = beh_data(lick_OFF_index,1);
% 
% %% Plot raw GCAMP data
% % Make some pretty colors for later plotting
% % http://math.loyola.edu/~loberbro/matlab/html/colorsInMatlab.html
% red = [0.8500, 0.3250, 0.0980];
% green = [0.4660, 0.6740, 0.1880];
% cyan = [0.3010, 0.7450, 0.9330];
% yellow = [0.9290, 0.6940, 0.1250];
% purple = [0.4940, 0.1840, 0.5560];
% gray1 = [.7 .7 .7];
% gray2 = [.8 .8 .8];
% 
% figure('Position',[100, 100, 800, 400])
% hold on;
% p1 = plot(time, LED1,'color',red,'LineWidth',1);
% p2 = plot(time, gcampdata,'color',cyan,'LineWidth',1);
% p3 = plot(time, ys,'color',green,'LineWidth',1);
% %p3 = plot(time, Y_exp_fit_all(time),'color',green,'LineWidth',10);
% p4 = plot(time,delta_F_p,'color',gray1,'LineWidth',1);
% ylabel('F','fontsize',16);
% axis tight;
% legend([p1 p2 p3 p4], {'Raw Signal', 'Denoised and Decay Corrected','Baseline','df/F'});
