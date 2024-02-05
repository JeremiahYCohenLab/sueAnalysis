folder = 'F:\michiganProbe\679458_070723_leftLC\679458_2023-07-07_13-51-24_d3500';
unit = 'TT4_SS_01';
%% set parameters
Wn = [600, 6000];

%%
session = Session(folder);
streamKey = session.recordNodes{1,1}.recordings{1,1}.continuous.keys();
streamKey = streamKey{1};
samples = session.recordNodes{1,1}.recordings{1,1}.continuous(streamKey).samples;
timeStamps =  session.recordNodes{1,1}.recordings{1,1}.continuous(streamKey).timestamps;
samplingFreq = session.recordNodes{1,1}.recordings{1,1}.continuous(streamKey).metadata.sampleRate;
numChannels = session.recordNodes{1,1}.recordings{1,1}.continuous(streamKey).metadata.numChannels;
%%
Wn = Wn /(samplingFreq/2);
[b, a] = butter(2, Wn, 'bandpass');
referenceC = mean(samples(1:13,:), 1);

% curr channel
chan = split(unit, '_');
chan = chan{1};
currC = str2double(chan(3:end));
currSamp = double(samples(currC, :)) - referenceC;
bit2volts = session.recordNodes{1,1}.recordings{1,1}.info.continuous.channels(currC).bit_volts;

% filtering
currSamp = filtfilt(b, a, currSamp)*bit2volts;
%%
% plot raster
figure2Wide;
hold on;
plot(timeStamps - timeStamps(1), currSamp);
xlabel('time (s)')
ylabel('uV')
title('LFP')
%% load single units
sortingFolder = [folder '\spikes\' ];
spikeTimes = readmatrix([sortingFolder unit '.txt'])/1000;
plot([spikeTimes - timeStamps(1), spikeTimes - timeStamps(1)]', [-250*ones(size(spikeTimes)), -270*ones(size(spikeTimes))]', 'k', 'LineWidth', 1)
set(gca, 'TickDir', 'out')
%%
% plot waveform
tmp_TTname = ['TT' num2str(currC) '.ntt'];
TTdir = fullfile(sortingFolder, tmp_TTname);
[tt_ts, tt_sig] = Nlx2MatSpike(TTdir, [1 0 0 0 1], 0, 1, 1);
waveform = squeeze(tt_sig(:, 1, ismember(tt_ts, spikeTimes*1000)))' * bit2volts/10;

figure2;
hold on;
sampleTime = 1:32;
sampleTime = sampleTime/samplingFreq * 1000;
plotFilledStd(sampleTime, waveform, 'k')

xlabel('time in ms')
ylabel('uV')
%%