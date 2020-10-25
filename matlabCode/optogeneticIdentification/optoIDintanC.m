function optoIDintanC(session, varargin)

p = inputParser;
% default parameters if none given
p.addParameter('Session', 'opto')
p.addParameter('Pulses', 10)
p.addParameter('Trains', 10)
p.addParameter('PulseWidth', 10)
p.addParameter('PulseFreq', 10)
p.addParameter('ResponseWindow', 30)
p.addParameter('MedianRemoval', true)
p.addParameter('HighPassCutoffInHz', 300);
p.parse(varargin{:});


%get session info
[root, sep] = currComputer();
pd = parseSessionString_oM(session, root, sep);
sortedPath = [pd.ephysPath 'spikeSort_' p.Results.Session sep];
dataPath = [pd.ephysPath p.Results.Session sep];
savePath = [pd.saveFigFolder 'optoID'];
if ~exist(savePath, 'dir')
    mkdir(savePath)
end

%get sorted and raw ephys data
optoFiles = dir(fullfile(dataPath,'*.rhd'));
sortedFiles = dir(fullfile(sortedPath,'*.txt'));
load(fullfile(sortedPath, 'settings.mat'), 'scaling_for_int16', 'frequency_parameters');
sampFreq = frequency_parameters.amplifier_sample_rate;

load(fullfile(sortedPath, 'events.mat'))
laser = events.DIN07;
laser = laser + 0.8; % 800us offset

%set params for butterworth filter
Wn = p.Results.HighPassCutoffInHz / (sampFreq/2);
[b, a] = butter(2, Wn, 'high');


%get list of channels with sorted units
for i = 1:length(sortedFiles)
    tmpInds = strfind(sortedFiles(i).name, '_');
    ttList(i) = str2double(sortedFiles(i).name(tmpInds(1)+1:tmpInds(2)-1));
end
Cprev = '';

%set window and stim paramaters
tB = 500;
tA = 500;
pulseInds = (1:p.Results.Pulses:p.Results.Pulses*p.Results.Trains);
respWin = p.Results.ResponseWindow;
rasterLength = length(-1*tB:(p.Results.Pulses*(1000/p.Results.PulseFreq)+tA));
    

%% 

for i = 1:length(sortedFiles)
    [cellName, ~] = strtok(sortedFiles(i).name, '.');
    spikeTimes = [load(strcat(sortedPath, sortedFiles(i).name))]';
    
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
    laserSham = [linspace(-tB, 0-1000/p.Results.PulseFreq, p.Results.Pulses/2) ...
        linspace(p.Results.Pulses*1000/p.Results.PulseFreq, rasterLength-1000/p.Results.PulseFreq, p.Results.Pulses/2)];
    
    spikeLat = nan(p.Results.Trains,p.Results.Pulses);
    spikeLatSham = nan(p.Results.Trains,p.Results.Pulses);
    lightSpikeTimes = [];
    for j = 1:p.Results.Trains              %for all pulses in all trains, find spikes within the response window
        for k = 1:p.Results.Pulses
            spikeRespTmp = spikeTimes(spikeTimes > laser(pulseInds(j)+k-1) & ...
                spikeTimes < laser(pulseInds(j)+k-1) + respWin);
            spikeRespTmpSham = spikeTimes(spikeTimes > laserSham(k) + laser(pulseInds(j)) & ...
                spikeTimes < laserSham(k) + laser(pulseInds(j)) + respWin);
            if ~isempty(spikeRespTmp)
                spikeLat(j,k) = spikeRespTmp(1) - laser(pulseInds(j)+k-1);
                lightSpikeTimes = [lightSpikeTimes spikeRespTmp(1)];
            end
             if ~isempty(spikeRespTmpSham)
                spikeLatSham(j,k) = spikeRespTmpSham(1) - (laserSham(k) + laser(pulseInds(j)));
            end           
        end
    end
    spontSpikeTimes = spikeTimes;
    spontSpikeTimes(ismember(spontSpikeTimes, lightSpikeTimes)) = [];
    avgSpikeLat = nanmean(spikeLat);    avgSpikeLatSham = nanmean(spikeLatSham);        %find average spikeLat and P(spike)
    semSpikeLat = nanstd(spikeLat)/sqrt(p.Results.Trains);      semSpikeLatSham = nanstd(spikeLatSham)/sqrt(p.Results.Trains); 
    spikeProb = mean(~isnan(spikeLat)); spikeProbSham = mean(~isnan(spikeLatSham));
    
    %get waveforms for light evoked and spontaneous spikes
    [Cname, unitNum] = strtok(cellName, 'SS');
    Cname = Cname(1:end-1);
    unitNum = unitNum(end);
    if strcmp(Cname, Cprev) == false % if the tetrode has changed, load a new one
        tmp_Cname = [Cname '.ntt'];
        Cdir = fullfile(sortedPath, tmp_Cname);
        [c_ts, c_cn, c_sig] = Nlx2MatSpike(Cdir, [1 0 1 0 1], 0, 1, 1);
        Cprev = Cname;
    end
    lightWaveForm = squeeze(c_sig(:, 1, ismember(c_ts, lightSpikeTimes)))' / scaling_for_int16;
    spontWaveForm = squeeze(c_sig(:, 1, ismember(c_ts, spontSpikeTimes)))' / scaling_for_int16;
 
  
    
    %% plot everything
    
    rasters = figure; subplot(4,3,[1:6]); hold on; title(strcat(session, '_', cellName),'Interpreter','none')
    xlabel('Time (ms)'); ylabel('Trials')
    LineFormat.Color = 'k'; LineFormat.LineWidth = 1;
    plotSpikeRaster(spikeRast,'PlotType','vertline','XLimForCell',[-1*tB rasterLength-tB],'LineFormat',LineFormat);
    hold on;
    x = linspace(0, ((p.Results.Pulses-1)*1000/p.Results.PulseFreq), p.Results.Pulses);
    xx = x + p.Results.PulseWidth;
    for j = 1:length(x)
        plotShaded([x(j) xx(j)],[0 0; 1+p.Results.Trains 1+p.Results.Trains],'b');
    end
    
    subplot(4,3,[7]); hold on;
    xlabel('Pulse'); ylabel('Latency (ms)'); ylim([0 30]); xlim([0 p.Results.Pulses+1])
    errorbar(avgSpikeLat, semSpikeLat, 'b', 'LineWidth', 2);
    errorbar(avgSpikeLatSham, semSpikeLatSham, 'k', 'LineWidth', 2);
    legend('laser','control');
    
    subplot(4,3,10); hold on;
    xlabel('Pulse'); ylabel('P(spike)'); ylim([-0.1 1.1]); xlim([0 p.Results.Pulses+1])
    plot(spikeProb, 'b', 'LineWidth', 2);
    plot(spikeProbSham, 'k', 'LineWidth', 2);

    subplot(4,3,8); hold on;
    ylabel('Amplitude (\muV)');
    plotFilled([1:32], spontWaveForm, 'k');
    
    subplot(4,3,11); hold on;
    ylabel('Amplitude (\muV)');
    if ~isempty(lightWaveForm)
        if size(lightWaveForm, 1) > 1
            plotFilled([1:32], lightWaveForm, 'b');
        else
            plot([1:32], lightWaveForm, 'b');
        end
    end
    
    subplot(4,3,[9 12]); hold on
    xlabel('Time (ms)'); ylabel('Amplitude (\muV)')
    if ~isempty(lightWaveForm)
        if size(lightWaveForm, 1) > 1
            plotFilled([1:32], lightWaveForm, 'b');
        else
            plot([1:32], lightWaveForm, 'b');
        end
    end
    plotFilled([1:32], spontWaveForm, 'k');
 
    set(rasters, 'Position', get(0,'Screensize'))
    saveFigurePDF(rasters,[savePath sep session '_' cellName '_' p.Results.Session 'ID'])
end
  