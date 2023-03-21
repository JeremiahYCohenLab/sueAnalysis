function optoIDTT(session, varargin)
shutterOffset = 0.8;
p = inputParser;
% default parameters if none given
p.addParameter('Session', 'opto')
p.addParameter('subFolder', '')
p.addParameter('Pulses', 10)
p.addParameter('Trains', 10)
p.addParameter('PulseWidth', 10) %ms
p.addParameter('ResponseWindow', 20000) %us
p.addParameter('MedianRemoval', true)  
p.addParameter('HighPassCutoffInHz', 300);
p.addParameter('SamplingFreq', 32000);
p.parse(varargin{:});


%get session info
[root, sep] = currComputer();
pd = parseSessionString_df(session, root, sep);
sortedPath = [pd.nLynxFolder p.Results.Session sep p.Results.subFolder sep];
dataPath = [pd.nLynxFolder p.Results.Session sep p.Results.subFolder sep];
savePath = [pd.saveFigFolder p.Results.Session sep p.Results.subFolder];
if ~exist(savePath, 'dir')
    mkdir(savePath)
end

%get sorted and raw ephys data
% optoFiles = dir(fullfile(dataPath,'*.rhd'));
sortedFiles = dir(fullfile(sortedPath,'*TT*.txt'));

% load(fullfile(sortedPath, 'settings.mat'), 'scaling_for_int16', 'frequency_parameters');
sampFreq = p.Results.SamplingFreq;

%events data
[tsEv, Ev] = Nlx2MatEV([sortedPath 'Events.nev'], [1 0 1 0 0], 0, 1);

% TTL events
laser = 4;
laserRaw = 1024; %16384; % raw data file uses this as laser TTL 

% isolate just laser on times
Ev(Ev == laserRaw) = laser;

biEv = de2bi(Ev);
biEv = biEv(:,3);
laserInd = (biEv(2:end)==1 & biEv(1:end-1)==0); %for some reason, _0_ is the laser on in this rig

laserOnTimes = tsEv(laserInd)/1000;
laserOnTimes = laserOnTimes + shutterOffset;
laser = laserOnTimes*1000;
PulseFreq = round(1000000/min(diff(laser)));   

% %set params for butterworth filter
% Wn = p.Results.HighPassCutoffInHz / (sampFreq/2);
% [b, a] = butter(2, Wn, 'high');


%get list of TTs with sorted units
for i = 1:length(sortedFiles)
    tmpInds = strfind(sortedFiles(i).name, '_');
    ttList(i) = str2double(sortedFiles(i).name(tmpInds(1)-1));
end
TTprev = '';
%set window and stim paramaters
tB = 500000;% in us
tA = 500000;
pulseInds = (1:p.Results.Pulses:p.Results.Pulses*p.Results.Trains);
respWin = p.Results.ResponseWindow;
rasterLength = length(-1*tB:(p.Results.Pulses*(1000000/PulseFreq)+tA));

header = Nlx2MatCSC([sortedPath 'CSC1.ncs'], [0 0 0 0 0], 1, 1, []);

AD2uV = split(header{contains(header, '-ADBitVolts')}, 'Volts');
AD2uV = str2double(AD2uV{2})*10^6;

%%

for i = 1:length(sortedFiles)
    [cellName, ~] = strtok(sortedFiles(i).name, '.');
    spikeTimes = [load(strcat(sortedPath, sortedFiles(i).name))]';

    % baseline freq
    sponsFreq = 1000000*sum(spikeTimes > laser(1) - 5000000 & spikeTimes < laser(1))/...
    (min(5000000, laser(1)-spikeTimes(1)));

    spikeRast = [];
    for j = 1:p.Results.Trains
        spikeRast{j} = spikeTimes(spikeTimes > (laser(pulseInds(j)) - tB) &...
            spikeTimes < (laser(pulseInds(j)+ p.Results.Pulses - 1) + p.Results.PulseWidth + tA));
        spontSpikeRast{j} = spikeTimes((spikeTimes > (laser(pulseInds(j)) - tB) & spikeTimes < laser(pulseInds(j))) |...
            (spikeTimes > (laser(pulseInds(j)+ p.Results.Pulses - 1) + p.Results.PulseWidth + respWin) & spikeTimes < (laser(pulseInds(j)+9) + tA)));
        if ~isempty(spikeRast{j})
            spikeRast{j} = spikeRast{j} - laser(pulseInds(j)); %puts in time relative to first light pulse    
        end
    end
    
    %find times when there is no light for control comparison
    laserSham = [linspace(-3000000, -tB, p.Results.Pulses/2) ...
        linspace(p.Results.Pulses*1000000/PulseFreq, rasterLength+3000000, p.Results.Pulses/2)];
    
    spikeLat = nan(p.Results.Trains,p.Results.Pulses);
    spikeLatSham = nan(p.Results.Trains,p.Results.Pulses);
    spikeNum = zeros(p.Results.Trains,p.Results.Pulses);
    spikeNumSham = zeros(p.Results.Trains,p.Results.Pulses);
    lightSpikeTimes = [];
    for j = 1:p.Results.Trains              %for all pulses in all trains, find spikes within the response window
        for k = 1:p.Results.Pulses
            spikeRespTmp = spikeTimes(spikeTimes > laser(pulseInds(j)+k-1) & ...
                spikeTimes < laser(pulseInds(j)+k-1) + respWin);
            spikeRespTmpSham = spikeTimes(spikeTimes > laserSham(k) + laser(pulseInds(j)) & ...
                spikeTimes < laserSham(k) + laser(pulseInds(j)) + respWin);
            if ~isempty(spikeRespTmp)
                spikeLat(j,k) = spikeRespTmp(1) - laser(pulseInds(j)+k-1);
                spikeNum(j,k) = length(spikeRespTmp);
                lightSpikeTimes = [lightSpikeTimes spikeRespTmp];
            end
            if ~isempty(spikeRespTmpSham)
                spikeLatSham(j,k) = spikeRespTmpSham(1) - (laserSham(k) + laser(pulseInds(j)));
                spikeNumSham(j,k) = length(spikeRespTmpSham);
            end           
        end
    end
    spontSpikeTimes = spikeTimes;
    spontSpikeTimes(ismember(spontSpikeTimes, lightSpikeTimes)) = [];
    avgSpikeLat = nanmean(spikeLat);    avgSpikeLatSham = nanmean(spikeLatSham);        %find average spikeLat and P(spike)
    semSpikeLat = nanstd(spikeLat)/sqrt(p.Results.Trains);      semSpikeLatSham = nanstd(spikeLatSham)/sqrt(p.Results.Trains); 
    spikeProb = mean(~isnan(spikeLat)); spikeProbSham = mean(~isnan(spikeLatSham));
    %get waveforms for light evoked and spontaneous spikes
    [TTname, unitNum] = strtok(cellName, 'SS');
    TTname = TTname(1:end-1);
    unitNum = unitNum(end);
    if strcmp(TTname, TTprev) == false % if the tetrode has changed, load a new one
        tmp_TTname = [TTname '.ntt'];
        TTdir = fullfile(sortedPath, tmp_TTname);
        [tt_ts, tt_sig] = Nlx2MatSpike(TTdir, [1 0 0 0 1], 0, 1, 1);
        TTprev = TTname;
    end
    
    % low pass waveform shape
    fc = 6000;
    [b, a] = butter(2,fc/(sampFreq/2),'low');
    tt_sig = filtfilt(b, a, tt_sig);
    
    for j = 1:4
        lightWaveForm{j} = AD2uV*squeeze(tt_sig(:, j, ismember(tt_ts, lightSpikeTimes)))';
        spontWaveForm{j} = AD2uV*squeeze(tt_sig(:, j, ismember(tt_ts, spontSpikeTimes)))';
    end
    % find channel with max peak
    % no.TT
    CSCnum = 4*str2double(TTname(end));
    CSCnum = CSCnum-3 : CSCnum;
    % no.CSC
    meanWaveform = cellfun(@(x) mean(x,1), spontWaveForm, 'UniformOutput', false);
    for k = 1:4
        maxChannel(k) = max(meanWaveform{k});
    end
    [~,maxChannel] = max(maxChannel);
    maxChannel = CSCnum(maxChannel);
    % get raw trace for trials
    [ts, samp] = Nlx2MatCSC([sortedPath 'CSC' num2str(maxChannel) '.ncs'], [1 0 0 0 1], 0, 1, []);
    samp = reshape(samp,[],1);
    tSamp = 1/p.Results.SamplingFreq * 1e6; % time for each sample
    if length(unique(diff(ts))) == 1 % no pausing
        ts_interp = ts(1):tSamp:ts(1) + tSamp*(size(samp,1) - 1);
    else % if pausing/skip due to data loss, use the proper for loop
        ts_interp = NaN(1, size(samp,1)); 
        for k = 1:length(ts)
            ts_interp(512*(k - 1) + 1:512*(k)) = ts(k):tSamp:ts(k) + tSamp*511;
        end
    end
    
    % find traces for trials
    rawTraces = cell(p.Results.Trains,1);
    time = cell(p.Results.Trains,1);
    for j = 1:p.Results.Trains
        rawTraces{j} = AD2uV*samp(ts_interp>laser(pulseInds(j))- tB & ts_interp<laser(pulseInds(j)+p.Results.Pulses-1) + tA); 
        time{j} = (ts_interp(ts_interp>laser(pulseInds(j))- tB & ts_interp<laser(pulseInds(j)+p.Results.Pulses-1) + tA)-laser(pulseInds(j)))/1000000;
    end
    
    
    %% plot everything
    
    rasters = figure; subplot(4,3,[1:3]); hold on; title(strcat(session, '_', cellName),'Interpreter','none')
    xlabel('Time (us)'); ylabel('Trials')
    LineFormat.Color = 'k'; LineFormat.LineWidth = 1;
    plotSpikeRaster(spikeRast,'PlotType','vertline','XLimForCell',[-1*tB rasterLength-tB],'LineFormat',LineFormat);
    hold on;
    x = linspace(0, ((p.Results.Pulses-1)*1000/PulseFreq), p.Results.Pulses);
    xx = x + p.Results.PulseWidth;
    for j = 1:length(x)
        plotShaded(1000*[x(j) xx(j)],[0 0; 1+p.Results.Trains 1+p.Results.Trains],'b');
    end
    
    text(0, 0.25, ['baseline activity ', num2str(sponsFreq) 'Hz']);
    
    subplot(4,3,[4:6]); hold on; 
    rowH = 1.2*range(rawTraces{1});
    x = x/1000;
    xx = xx/1000;
    for j = 1:p.Results.Trains
        plot(time{j}, rawTraces{j} + rowH*(j-1), 'color', 'k');
        for k = 1:p.Results.Pulses
            line([x(k) xx(k)], [max(rawTraces{j})+rowH*(j-1) max(rawTraces{j})+rowH*(j-1)], 'color', [0 0 1], 'LineWidth', 2);
        end
    end
    
    plot([-0.25 -0.25], [-700 -200],'color', 'k', 'lineWidth',2);
    plot([-0.25 0], [-700 -700],'color', 'k', 'lineWidth',2);
    text(-0.5,-350, '500 uV', 'HorizontalAlignment','left')
    text(-0.1,-900, '0.25 s', 'HorizontalAlignment','center')
%     for j = 1:length(x)
%         plotShaded([x(j) xx(j)]/1000,[min(rawTraces{1}) min(rawTraces{1}); max(rawTraces{1})+rowH*(p.Results.Trains-1) max(rawTraces{1})+rowH*(p.Results.Trains-1)],'b');
%     end
    xlim(minmax(time{1}));
    
    subplot(4,3,7); hold on;
    xlabel('Pulse'); ylabel('Latency (ms)'); ylim([0 respWin/1000]); xlim([0 p.Results.Pulses+1])
    errorbar(avgSpikeLat/1000, semSpikeLat/1000, 'b', 'LineWidth', 2);
%     errorbar(avgSpikeLatSham/1000, semSpikeLatSham/1000, 'k', 'LineWidth', 2);
%     legend('laser','control');
    ylim([0 90])
    
    subplot(4,3,8); hold on;
    xlabel('Pulse'); ylabel('P(spike)'); ylim([-0.1 1.1]); xlim([0 p.Results.Pulses+1])
    plot(spikeProb, 'b', 'LineWidth', 2);
    plot(spikeProbSham, 'k', 'LineWidth', 2);

    subplot(4,3,9); hold on;
    xlabel('Pulse'); ylabel('spikeNum'); ylim([-0.1 max([spikeNum, spikeNumSham],[],'all')]); xlim([0 p.Results.Pulses+1])
    plot(mean(spikeNum), 'b', 'LineWidth', 2);
    plot(mean(spikeNumSham), 'k', 'LineWidth', 2);
    

    subplot(4,3,10); hold on;
    ylabel('Amplitude (\muV)');
    if ~isempty(spontWaveForm{1})
        if size(spontWaveForm{1}, 1) > 1
            for j = 1:4
                plotFilled([1:32]+32*(j-1), spontWaveForm{j}, 'k');
            end
        else
            for j = 1:4
                plot([1:32]+32*(j-1), spontWaveForm{j}, 'k');
            end
        end
    end
    
    subplot(4,3,11); hold on;
    ylabel('Amplitude (\muV)');
    if ~isempty(lightWaveForm{1})
        if size(lightWaveForm{1}, 1) > 1
            for j = 1:4
                plotFilled([1:32]+32*(j-1), lightWaveForm{j}, 'b');
            end
        else
            for j = 1:4
                plot([1:32]+32*(j-1), lightWaveForm{j}, 'b');
            end
        end
    end
    
    subplot(4,3,12); hold on
    ylabel('Amplitude (\muV)')
    if ~isempty(lightWaveForm{1})
        if size(lightWaveForm{1}, 1) > 1
            for j = 1:4
                plotFilled([1:32]+32*(j-1), lightWaveForm{j}, 'b');
            end
        else
            for j = 1:4
                plot([1:32]+32*(j-1), lightWaveForm{j}, 'b');
            end
        end
    end
    if ~isempty(spontWaveForm{1})
        if size(spontWaveForm{1}, 1) > 1
            for j = 1:4
                plotFilled([1:32]+32*(j-1), spontWaveForm{j}, 'k');
            end
        else
            for j = 1:4
                plot([1:32]+32*(j-1), spontWaveForm{j}, 'k');
            end
        end
    end
    
    plot([5 5], [70 120],'color', 'k', 'lineWidth',2);
    plot([5 21], [70 70],'color', 'k', 'lineWidth',2);
    text(10,100, '50 uV', 'HorizontalAlignment','left')
    text(15,60, '0.5 ms', 'HorizontalAlignment','center')
    screen = get(0,'Screensize');
    screen(4) = screen(4) - 100;
    set(rasters, 'Position', screen)
    saveFigurePDF(rasters,[savePath sep session '_' cellName '_' p.Results.Session 'ID'])
end
  