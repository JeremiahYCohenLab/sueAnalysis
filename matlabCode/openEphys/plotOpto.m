% plot neurons against opto
fileName = '2023-06-07_12-28-27';
folder = ['F:\allenData\683497_leftLC_060723\' fileName];
%% prepare neuron signal
% signal unit by themselves
unitTbl = readtable([folder '\spikes\' fileName '.xlsx']);
unitName = unitTbl.cellID;
unitName = unitName(unitTbl.Lratio<0.1);
allUnits = cell(1, length(unitName));
for i = 1:length(unitName)
    allUnits{i} = csvread([folder '\spikes\' unitName{i} '.txt'], 0, 0);
end
% multi-unit by each channel
allMultiUnits = cell(16, 1);
for i = 1:16
    ttName = ['TT' num2str(i) '.ntt'];
    [Timestamps,Samples] = Nlx2MatSpike([folder '\spikes\' ttName],[1 0 0 0 1],0,1,[]);
    allMultiUnits{i} = Timestamps;
end
% LFP by channel
session = Session(folder);
streamKey = session.recordNodes{1,1}.recordings{1,1}.continuous.keys();
streamKey = streamKey{1};
samples = session.recordNodes{1,1}.recordings{1,1}.continuous(streamKey).samples;
timeStamps =  session.recordNodes{1,1}.recordings{1,1}.continuous(streamKey).timestamps;
samplingFreq = session.recordNodes{1,1}.recordings{1,1}.continuous(streamKey).metadata.sampleRate;
numChannels = session.recordNodes{1,1}.recordings{1,1}.continuous(streamKey).metadata.numChannels;
lowpass = 300;
Wn = lowpass/(samplingFreq/2);
[b, a] = butter(2, Wn, 'low');  
% median subtranction
referenceC = mean(samples, 1);
samples = double(samples) - referenceC;
lfps = -filtfilt(b, a, samples')';


% plot by opto file
preLen = 500; % in ms
postLen = 1000; % in ms
preSamp = round(500/1000 * samplingFreq);
postSamp = round(samplingFreq);

load([folder '\opto\opto.mat']);

% specify led file
for i = 1:length(led)
    currLed = unique(led{i});
    for j = 1:length(currLed)
        currLedTime = ttlTime{i}(led{i} == currLed(j));
        % create cell for raster
        rasterCell = cell(length(currLedTime), length(allUnits));
        for s = 1:length(currLedTime)
            for u = 1:length(allUnits)
                temp = allUnits{u};
                rasterCell{s, u} = temp(temp > (1000*currLedTime(s)-preLen) & temp <= (1000*currLedTime(s)+postLen)) - 1000*currLedTime(s);
            end
        end
        % create cell for lfp
        lfpCell = cell(1, numChannels);
        for c = 1:numChannels
            lfptemp = [];
            for s = 1:length(currLedTime)
                [~, indS] = min(abs(timeStamps - currLedTime(s)));
                lfptemp = [lfptemp; lfps(c, indS-preSamp:indS+postSamp)];
            end
            lfpCell{c} = lfptemp;
        end
        % create cell for multi-unit
        multiCell = cell(length(currLedTime), numChannels);
        for c = 1:numChannels
            temp = allMultiUnits{c};
            for s = 1:length(currLedTime)               
                multiCell{s, c} = temp(temp > (1000*currLedTime(s)-preLen) & temp <= (1000*currLedTime(s)+postLen)) - 1000*currLedTime(s);
            end
        end        
        
        
        ledResp(j).raster = rasterCell;
        ledResp(j).lfp = lfpCell;
        ledResp(j).multi = multiCell;
    end 
    optoResp(i).resp = ledResp;
    optoResp(i).unit = unitName;
    optoResp(i).preLen = preLen;
    optoResp(i).postLen = postLen;
    optoResp(i).preSamp = preSamp;
    optoResp(i).postSamp = postSamp;
end

save([folder '\opto\optoResp.mat'], 'optoResp');
%% plot aligned files and save
fileName = '2023-06-07_12-28-27';
folder = ['F:\allenData\683497_leftLC_060723\' fileName '\opto\'];
fileName = [folder 'optoResp.mat'];
load(fileName);
% plot single units 

