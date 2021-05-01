function rasters = qualityCheck(session, saveFlag, varargin)
shutterOffset = 0.8;
p = inputParser;
% default parameters if none given
p.addParameter('Session', 'opto')
p.addParameter('Pulses', 10)
p.addParameter('ResponseWindow', 20000) %us
p.addParameter('MedianRemoval', true)  
p.addParameter('HighPassCutoffInHz', 300);
p.addParameter('SamplingFreq', 32000);
p.addParameter('unit', '');
p.parse(varargin{:});

%get session info
[root, sep] = currComputer();
pd = parseSessionString_df(session, root, sep);
savePath = [pd.saveFigFolder sep 'qualityCheck'];
xlFile = [pd.animalName '.xlsx'];
unit = p.Results.unit;
% get corresponding unit files
[nums, unitsInfo,~] = xlsread([root xlFile], 'neurons');
[row,~] = find(contains(unitsInfo, session));
if unit
   [rowU,~] = find(contains(unitsInfo(:,2), unit));
   row = intersect(row,rowU);
end

if ~exist(savePath, 'dir')
    mkdir(savePath)
end

% load(fullfile(sortedPath, 'settings.mat'), 'scaling_for_int16', 'frequency_parameters');
sampFreq = p.Results.SamplingFreq;

% TTL events
laser = 4;
laserRaw = 1024; %16384; % raw data file uses this as laser TTL 


TTprev = '';

%%

for i = 1:length(row)
    Trains = 10;
    subFolder = unitsInfo{row(i), 9};
    pulseWidth = str2double(strtok(subFolder,'ms'));
    sortedPath = [pd.nLynxFolder p.Results.Session sep subFolder sep];
    cellName = [unitsInfo{row(i),8} '.txt'];
    spikeTimes = (load(strcat(sortedPath, cellName)))';
    %Event data
    %events data
    [tsEv, Ev] = Nlx2MatEV([sortedPath 'Events.nev'], [1 0 1 0 0], 0, 1);
    laser = 4;
    % isolate just laser on times
    Ev(Ev == laserRaw) = laser;
    
    biEv = de2bi(Ev);
    biEv = biEv(:,3);
    laserInd = (biEv(2:end)==1 & biEv(1:end-1)==0); %for some reason, _0_ is the laser on in this rig
   % laserInd = [false laserInd'];
    
    laserOnTimes = tsEv(laserInd)/1000;
    laserOnTimes = laserOnTimes + shutterOffset; 
    
    if length(laserOnTimes) ~= 100
        figure2;
        plot(laserOnTimes);
        adjustTime = true;
    end
    laser = laserOnTimes*1000;
    PulseFreq = round(1000000/min(diff(laser)));  

    %set window and stim paramaters
    tB = 500000;% in us
    tA = 500000;
    pulseInds = (1:p.Results.Pulses:p.Results.Pulses*Trains);
    respWin = p.Results.ResponseWindow;
    rasterLength = length(-1*tB:(p.Results.Pulses*(1000000/PulseFreq)+tA));

    header = Nlx2MatCSC([sortedPath 'CSC1.ncs'], [0 0 0 0 0], 1, 1, []);

    AD2uV = split(header{contains(header, '-ADBitVolts')}, 'Volts');
    AD2uV = str2double(AD2uV{2})*10^6;

    % baseline freq
    sponsFreq = 1000000*sum(spikeTimes > laser(1) - 5000000 & spikeTimes < laser(1))/...
    (min(5000000, laser(1)-spikeTimes(1)));

    spikeRast = [];
    for j = 1:Trains
        spikeRast{j} = spikeTimes(spikeTimes > (laser(pulseInds(j)) - tB) &...
            spikeTimes < (laser(pulseInds(j)+ p.Results.Pulses - 1) + pulseWidth + tA));
        spontSpikeRast{j} = spikeTimes((spikeTimes > (laser(pulseInds(j)) - tB) & spikeTimes < laser(pulseInds(j))) |...
            (spikeTimes > (laser(pulseInds(j)+ p.Results.Pulses - 1) + pulseWidth + respWin) & spikeTimes < (laser(pulseInds(j)+9) + tA)));
        if ~isempty(spikeRast{j})
            spikeRast{j} = spikeRast{j} - laser(pulseInds(j)); %puts in time relative to first light pulse    
        end
    end
    
    %find times when there is no light for control comparison
    laserSham = [linspace(-3000000, -tB, p.Results.Pulses/2) ...
        linspace(p.Results.Pulses*1000000/PulseFreq, rasterLength+3000000, p.Results.Pulses/2)];
    
    spikeLat = nan(Trains,p.Results.Pulses);
    spikeLatSham = nan(Trains,p.Results.Pulses);
    spikeNum = zeros(Trains,p.Results.Pulses);
    spikeNumSham = zeros(Trains,p.Results.Pulses);
    lightSpikeTimes = [];
    for j = 1:Trains              %for all pulses in all trains, find spikes within the response window
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
    semSpikeLat = nanstd(spikeLat)/sqrt(Trains);      semSpikeLatSham = nanstd(spikeLatSham)/sqrt(Trains); 
    spikeProb = mean(~isnan(spikeLat)); spikeProbSham = mean(~isnan(spikeLatSham));
    %get waveforms for light evoked and spontaneous spikes
    [TTname, unitNum] = strtok(cellName, 'SS');
    TTname = TTname(1:end-1);
    unitNum = unitNum(end);
    tmp_TTname = [TTname '.ntt'];
    TTdir = fullfile(sortedPath, tmp_TTname);
    [tt_ts, tt_sig] = Nlx2MatSpike(TTdir, [1 0 0 0 1], 0, 1, 1);

    for j = 1:4
        lightWaveForm{j} = AD2uV*squeeze(tt_sig(:, j, ismember(tt_ts, lightSpikeTimes)))';
        spontWaveForm{j} = AD2uV*squeeze(tt_sig(:, j, ismember(tt_ts, spontSpikeTimes)))';
    end
    
    %session data
    cellNameSession = [unitsInfo{row(i),2} '.txt'];
    sortedPathSession = [pd.nLynxFolder 'session' sep];
    spikeTimesSession = (load(strcat(sortedPathSession, cellNameSession)))';
    [~, unitNum] = strtok(cellNameSession, 'SS');
    unitNum = unitNum(end);
    TTdir = fullfile(sortedPathSession, tmp_TTname);
    [tt_ts, tt_sig] = Nlx2MatSpike(TTdir, [1 0 0 0 1], 0, 1, 1);

    for j = 1:4
        WaveForm{j} = AD2uV*squeeze(tt_sig(:, j, ismember(tt_ts, spikeTimesSession)))';
    end
    
    %% plot everything
    
    rasters = figure; subplot(4,3,[1:6]); hold on; title(strcat(session, '_', cellNameSession),'Interpreter','none')
    xlabel('Time (us)'); ylabel('Trials')
    LineFormat.Color = 'k'; LineFormat.LineWidth = 1;
    plotSpikeRaster(spikeRast,'PlotType','vertline','XLimForCell',[-1*tB rasterLength-tB],'LineFormat',LineFormat);
    hold on;
    x = linspace(0, ((p.Results.Pulses-1)*1000/PulseFreq), p.Results.Pulses);
    xx = x + pulseWidth;
    for j = 1:length(x)
        plotShaded(1000*[x(j) xx(j)],[0 0; 1+Trains 1+Trains],'b');
    end
    
    text(0, 0.25, ['baseline activity ', num2str(sponsFreq) 'Hz']);
    
    subplot(4,3,7); hold on;
    xlabel('Pulse'); ylabel('Latency (ms)'); ylim([0 respWin/1000]); xlim([0 p.Results.Pulses+1])
    errorbar(avgSpikeLat/1000, semSpikeLat/1000, 'b', 'LineWidth', 2);
%     errorbar(avgSpikeLatSham/1000, semSpikeLatSham/1000, 'k', 'LineWidth', 2);
%     legend('laser','control');
    ylim([0 30])
    
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
    for j = 1:4
        plotFilled([1:32]+32*(j-1), WaveForm{j}, 'k');
    end
    line([0 128], [0 0], 'color', [0.7, 0.7, 0.7]);
    text(10, 150, sprintf('spikeNumber %d' , length(spikeTimesSession)));
    text(10, -100, sprintf('Lratio %d' , nums(row(i)-1,1)));
    ylim([-200 200]);
    
    
    subplot(4,3,11); hold on;
    ylabel('density');
    h = histogram(diff(spikeTimesSession),'BinWidth',1000,'Normalization','probability');
    set(gca, 'XScale', 'log')
    line([2000 2000],[0 0.002],'color',[1 0 0]);
    formatSpec = ['RP violation = %.2f ' '%'];
    text(1000, 0.002, [sprintf(formatSpec, 100*sum(diff(spikeTimesSession)<2000)/(length(spikeTimesSession)-1)) '%']);
    title('ISI')
    
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
    if saveFlag
        saveFigurePDF(rasters,[savePath sep session '_' cellNameSession '_opto_' cellName '.pdf'])
    end
end
  