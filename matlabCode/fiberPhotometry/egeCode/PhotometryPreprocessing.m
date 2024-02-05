%% Photometry Data Extraction and Normalization

raw_gcamptimestamps = photometry_data(:,1);
raw_gcampdata = photometry_data(:,2);
% if there are an odd # of rows after removing headers, delete the last row
% in both wave and wavetime
A = raw_gcampdata(1:2:end,:);
B = raw_gcampdata(2:2:end,:);

%truncate the data from the LED that has more datapoints to match data
%dimension
if length(A) > length(B)
  A = A(1:end-1,:);
elseif length(A) < length(B)
  B = B(1:end-1,:);
end

%the raw fluorescence values from LED that are consistently greater belong
%to 470 nm (i.e. LED1)
if sum(A>B) > (length(A)*.99) %% 0.99 seems arbitrary here?
  LED1 = A;
  LED2 = B;
else
  LED1 = B;
  LED2 = A;
end

% Timestamps for gcampdata
if LED1 == A
    k = 1;
    j = 2;
elseif LED1 == B
    k = 2;
    j = 1;
end
gcampdata_timestamps = raw_gcamptimestamps(k:2:end,:);
gcampdata_timestamps2 = raw_gcamptimestamps(j:2:end,:);

%% Denoising the signal 
LED1dn = medfilt1(LED1,5); %median filter gets rid of random/ salt and pepper noise 
%If you have rare outliers an order of 5 should be good. 
%Otherwise you can increase to 7, 9 or even 11 but not more because you`ll
%distort.
fc = 10; %cutoff frequency in Hz (anything lower than this passes)
fs = 20; %sampling rate in Hz
%below returns transfer function coefficients of an nth-order lowpass digital filter 
%Butterworth filter with normalized cutoff frequency Wn
%For digital filters, the cutoff frequencies must lie between 0 and 1
%where 1 corresponds to the Nyquist frequency
%rad/sample)
%[b,a] = butter(2,fc/(fs/2),'low');
[b,a] = butter(2,0.99999999,'low');
LED1dn2 = filtfilt(b,a,LED1dn); %filtfilt is zero-phase so no temporal distortion

%% High-pass filter method of accounting for slow decay 
fc2 = 0.001; %cutoff frequency in Hz (anything higher than this f passes)
[b2,a2] = butter(2,fc2/(fs/2),'high'); %this accounts for slow time course changes on the order of ~16 minutes
LED1c = filtfilt(b2,a2,LED1dn2); %filtfilt is zero-phase
gcampdata = LED1c+50000; %add a constant so the trace doesn't go down to zero on y-axis
SR = 20;
time = (1:length(gcampdata))/SR; % in seconds

%% Double Exponential Method for accounting for decay across session (3.25.19)
%% PICK ONE OR THE OTHER
% SR = 20;
%time = (1:length(LED1))/SR; % in seconds
% time = (1:length(LED1dn2))/SR; % in seconds
% %Y_exp_fit_all = fit(time',LED1,'exp2');
% %Y_exp_fit_all = fit(time',LED1dn2,'exp2');
% f_0 = Y_exp_fit_all(0);
% normDat = (LED1dn2 .* f_0) ./ Y_exp_fit_all(time);
% gcampdata = normDat;

%% Obtain Baseline Trace
%estimate 10th percentile of F using 15 sec. sliding window 
%this is to generate baseline fluorescence estimate
win = SR * 15;
p = 10;
% [ ys ] = running_percentile(gcampdata, win, p);
%if you want to run a left-sided sliding window 
winright = 0;
[ ys ] = running_percentile_filter(gcampdata, win, winright, p); 
delta_F = gcampdata - ys;
delta_FoverF = delta_F ./ ys ;
%exclusion criterion: only sessions where 97.5% of DF/F0 across the entire session
%exceeded 1% are included
delta_F_p = delta_FoverF .* 100; %percent change from baseline fluorescence
exc = prctile(delta_F_p, 97.5);
pass = (exc > 1) 
%delta_F_z = zscore(delta_F,0,'all'); %%WHY NEED THIS??

%% Extract behavior timestamps from analog_inputs
%%
% Subtract logicals from each other in LLP column
switches_LLP = diff(beh_data(:,2)); %2
% Subtract logicals from each other in RLP column
switches_RLP = diff(beh_data(:,4)); %4
% Subtract logicals from each other in HE column
switches_HE = diff(beh_data(:,3)); %3
% Subtract logicals from each other in lick column
switches_lick = diff(beh_data(:,5)); %5

% When the difference is 1, that's event ON, find the indices for ON
LLP_ON_index = find(switches_LLP == 1)+1;
RLP_ON_index = find(switches_RLP == 1)+1;
HE_ON_index = find(switches_HE == 1)+1;
lick_ON_index = find(switches_lick == 1)+1;

% When the difference is -1, that's event OFF, find indices for OFF
LLP_OFF_index = find(switches_LLP == -1)+1;
RLP_OFF_index = find(switches_RLP == -1)+1;
HE_OFF_index = find(switches_HE == -1)+1;
lick_OFF_index = find(switches_lick == -1)+1;

% Get the time stamps from those event indices
LLP_ON_timestamps = beh_data(LLP_ON_index, 1);
LLP_OFF_timestamps = beh_data(LLP_OFF_index, 1);
RLP_ON_timestamps = beh_data(RLP_ON_index, 1);
RLP_OFF_timestamps = beh_data(RLP_OFF_index, 1);
HE_ON_timestamps = beh_data(HE_ON_index, 1);
HE_OFF_timestamps = beh_data(HE_OFF_index,1);
lick_ON_timestamps = beh_data(lick_ON_index, 1);
lick_OFF_timestamps = beh_data(lick_OFF_index,1);

%% Plot raw GCAMP data
% Make some pretty colors for later plotting
% http://math.loyola.edu/~loberbro/matlab/html/colorsInMatlab.html
red = [0.8500, 0.3250, 0.0980];
green = [0.4660, 0.6740, 0.1880];
cyan = [0.3010, 0.7450, 0.9330];
yellow = [0.9290, 0.6940, 0.1250];
purple = [0.4940, 0.1840, 0.5560];
gray1 = [.7 .7 .7];
gray2 = [.8 .8 .8];

figure('Position',[100, 100, 800, 400])
hold on;
p1 = plot(time, LED1,'color',red,'LineWidth',1);
p2 = plot(time, gcampdata,'color',cyan,'LineWidth',1);
p3 = plot(time, ys,'color',green,'LineWidth',1);
%p3 = plot(time, Y_exp_fit_all(time),'color',green,'LineWidth',10);
p4 = plot(time,delta_F_p,'color',gray1,'LineWidth',1);
ylabel('F','fontsize',16);
axis tight;
legend([p1 p2 p3 p4], {'Raw Signal', 'Denoised and Decay Corrected','Baseline','df/F'});
