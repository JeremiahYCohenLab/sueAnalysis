function met = getClusterMetric(session, unit, plotFlag, saveFlag, varargin)
shutterOffset = 0.8;
p = inputParser;
% default parameters if none given
p.addParameter('Session', 'opto')
p.addParameter('Pulses', 10)
p.addParameter('ResponseWindow', 20000) %us
p.addParameter('MedianRemoval', true)  
p.addParameter('HighPassCutoffInHz', 300);
p.addParameter('SamplingFreq', 32000);
p.parse(varargin{:});

%get session info
[root, sep] = currComputer();
pd = parseSessionString_df(session, root, sep);

sampFreq = p.Results.SamplingFreq;

% load broken channel info
[animalName, date] = strtok(session, 'd'); 
animalName = animalName(2:end);
xlFile = [animalName '.xlsx'];
brokenChannels = xlsread([root xlFile], 'brokenChannels');
% get corresponding unit files
[nums, unitsInfo,~] = xlsread([root xlFile], 'neurons');
[row,~] = find(contains(unitsInfo, session));
[rowU,~] = find(contains(unitsInfo(:,2), unit));
row = intersect(row,rowU);
subFolder = unitsInfo{row, 9};
optoUnit = unitsInfo{row,8};
% pathInfo
    Trains = 10;
    pulseWidth = str2double(strtok(subFolder,'ms'));
    sortedPath = [pd.nLynxFolder 'opto' sep subFolder sep];
    savePath = [pd.saveFigFolder 'optoMet' sep];
    cellName = [optoUnit '.txt'];
    %% nev and ntt file
    spikeTimes = (load(strcat(sortedPath, cellName)))';
    %Event data
    %events data
    [tsEv, Ev] = Nlx2MatEV([sortedPath 'Events.nev'], [1 0 1 0 0], 0, 1);
    laser = 4;
    laserRaw = 1024;
    % isolate just laser on times
    Ev(Ev == laserRaw) = laser;
    
    biEv = de2bi(Ev);
    biEv = biEv(:,3);
    laserInd = (biEv(2:end)==1 & biEv(1:end-1)==0); %for some reason, _0_ is the laser on in this rig
   % laserInd = [false laserInd'];
    
    laserOnTimes = tsEv(laserInd)/1000;
    laserOnTimes = laserOnTimes + shutterOffset; 
    
    if length(laserOnTimes) ~= 100
        if rem(length(laserOnTimes),10)==0
            Trains = length(laserOnTimes)/10;
        else
            figure2;
            plot(laserOnTimes);
            adjustTime = true;
        end
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
    
    % correct spiketime to first point different from baseline
    % get TT file waveforms, find peak times on tallest channel
    [TTname, optoUnitNum] = strtok(cellName, 'SS');
    TTname = TTname(1:end-1);
    optoUnitNum = optoUnitNum(end);
    tmp_TTname = [TTname '.ntt'];
    TTdir = fullfile(sortedPath, tmp_TTname);
    [tt_ts, tt_sig] = Nlx2MatSpike(TTdir, [1 0 0 0 1], 0, 1, 1);
    allWaveForm = squeeze(tt_sig(:, :, ismember(tt_ts, spikeTimes)));
    meanWaveForm = mean(allWaveForm, 3);
    [~, peakChannel] = max(max(meanWaveForm,[], 1));
    lowestChannel = find(max(meanWaveForm, [], 1)== 0);
    [~, peakTimes] = max(allWaveForm(:,peakChannel,:), [], 1);
    % realign spikeTime
    tSamp = 1/p.Results.SamplingFreq * 1e6; % time per sample in microseconds
    if length(spikeTimes) == length(squeeze(peakTimes))
     spikeTimes = spikeTimes + tSamp*(squeeze(peakTimes)' - 11); % change depend on how different were the peak from 11th sample
    end
    %% get waveform from CSC and find the earliest different time point
    ttNum = str2double(strtok(TTname,'TT'));
    chan = (ttNum*4 -3):ttNum*4;
    [ts, samp0] = Nlx2MatCSC([sortedPath 'CSC' num2str(chan(1)) '.ncs'], [1 0 0 0 1], 0, 1, []);
    [samp1] = Nlx2MatCSC([sortedPath 'CSC' num2str(chan(2)) '.ncs'], [0 0 0 0 1], 0, 1, []);
    [samp2] = Nlx2MatCSC([sortedPath 'CSC' num2str(chan(3)) '.ncs'], [0 0 0 0 1], 0, 1, []);
    [samp3] = Nlx2MatCSC([sortedPath 'CSC' num2str(chan(4)) '.ncs'], [0 0 0 0 1], 0, 1, []);
    samp = cat(3, samp0, samp1, samp2, samp3);
    samp = reshape(samp, [] ,4); 
    % find reference
    if ~isempty(lowestChannel)
        if length(lowestChannel)==1 % with one broken or one ref
            ref = samp(:,lowestChannel);  
        else
            if length(lowestChannel)>1 % with both broken and ref
                ref = setdiff(chan(lowestChannel),brokenChannels);  % find the unbroken ref
                ref = samp(:,chan==ref); 
            end

        end
        % change reference to the smallest channel if ntt file did that.
        for j = 1:4
            samp(:,j) = samp(:,j) - ref;
        end
    end
        % remove broken channel
    if ~isempty(ismember(chan,brokenChannels))
        samp(:,ismember(chan,brokenChannels)) = 0;
    end
    % calculate time in CSC
    if length(unique(diff(ts))) == 1 % no pausing
        ts_interp = ts(1):tSamp:ts(1) + tSamp*(size(samp,1) - 1);
    else % if pausing/skip due to data loss, use the proper for loop
        ts_interp = NaN(1, size(samp,1)); 
        for i = 1:length(ts)
            ts_interp(512*(i - 1) + 1:512*(i)) = ts(i):tSamp:ts(i) + tSamp*511;
        end
    end

    leg = -50:50;
    peakInds = zeros(length(spikeTimes),1);
    for k = 1:length(spikeTimes)
        [~, tmpInd] = min(abs(ts_interp - spikeTimes(k)));
        peakInds(k) = tmpInd;
    end
    intactInds = (peakInds+leg(1))>0 & (peakInds+leg(end))<size(samp,1);
    peakInds = peakInds(intactInds);
    spikeTimes = spikeTimes(intactInds);
    waveforms = zeros(length(leg),4,length(peakInds));
    % waveform from CSC
    for w = 1:length(leg)
        waveforms(w,:,:) = samp(peakInds + leg(w), :)';
    end
    % find samples not affected by spikes
    nospikeSamp = [];
    for w = 1:20
        tmp = samp(peakInds(2:end)- 40 - w, :)';
        nospikeSamp = [nospikeSamp, tmp];
    end
    peakLagAll = zeros(4,1);
    peakEndAll = zeros(4,1);
    tmpH = zeros(length(leg),4);

    baseline = mean(nospikeSamp,2);
    for w = 1:4
        for k = 1:length(leg)
            h = ttest(waveforms(k,w,:), baseline(w),'Alpha',0.001);
            tmpH(k,w) = h;
        end
        tmpHCov = conv(tmpH(:,w), ones(1,8));
        tmpHCov = tmpHCov(8:end);
        if ~isempty(find(tmpHCov(1:end-2) == 7 & tmpHCov(2:end-1)>5 & tmpHCov(3:end)>5))
          peakLagAll(w) = min(find(tmpHCov(1:end-2) == 7 & tmpHCov(2:end-1)>5 & tmpHCov(3:end)>5)); % find first continued 5 sig points
          peakEndAll(w) = max(find(tmpHCov == 5))+4; % find the end of the spike
        else
            peakLagAll(w) = 51;
            peakEndAll(w) = 51;
        end
        
    end
    
    peakLag = min(peakLagAll) - 51; 
    peakEnd = max(peakEndAll) - 51; % sample number relative to peak time. 
    spikeTimes = spikeTimes + peakLag*tSamp; % correct spike time to deviation time
%% baseline freq
    baselineSpikes = 0;
    for j = 1:Trains
        if j == 1
            baselineSpikes = sum(spikeTimes < laser(1));
        else
            tmp = sum(spikeTimes > (laser(1 + p.Results.Pulses*(j-1)) - 1000000) & spikeTimes < laser(1 + p.Results.Pulses*(j-1)));
            baselineSpikes = baselineSpikes + tmp;
        end
    end
    baselineTime = 1000000*(Trains - 1) + laser(1) - spikeTimes(1);
    sponsFreq = 1000000*baselineSpikes/baselineTime;
    % get rasters
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
    % baseline respond num
    spikeNumSham = sponsFreq*pulseWidth/1000;
    spikePropSham = sponsFreq*pulseWidth/1000;
    
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
    
    if ~isempty(lightSpikeTimes)
        lightWaveForm = waveforms(:, :, ismember(spikeTimes, lightSpikeTimes));
    else
        lightWaveForm = [];
    end
    if ~isempty(lightSpikeTimes)
        spontWaveForm = waveforms(:, :, ismember(spikeTimes, spontSpikeTimes));
    else
        spontWaveForm = [];
    end
%% summarize metrices
    % waveform differece
    mainInds = 40:70;
    mainWaveform = mean(squeeze(waveforms(mainInds,peakChannel,:)),2);
    energy = norm(mainWaveform)/sqrt(length(mainWaveform));
    distance = mean(squeeze(lightWaveForm(mainInds,peakChannel,:)),2)-mean(squeeze(spontWaveForm(mainInds,peakChannel,:)),2);
    distance = norm(distance)/sqrt(length(distance));
    distance = distance/energy;
    met.distance = distance;
    % spikeLat
    met.spikeLat = nanmean(spikeLat);
    % spikeProp
    met.spikeProp = spikeProb;
    % spikePropSham (also spikeNumSham)
    met.spikePropSham = spikePropSham;
    % baselineFreq
    met.baseline = sponsFreq;
    % width (trough to peak)
    meanWaveform = mean(squeeze(waveforms(peakLag+51:peakEnd+51,peakChannel,:)),2);
    peakSign = sign(meanWaveform(-peakLag+1));
    if ~isempty(lightSpikeTimes)
        if peakSign > 0
            valleys = find(diff(meanWaveform(1:end-1))<0 & diff(meanWaveform(2:end))>0)+1;
        else
            valleys = find(diff(meanWaveform(1:end-1))>0 & diff(meanWaveform(2:end))<0)+1;
        end
        valleyInd = min(find(valleys>-peakLag));
        valleyInd = valleys(valleyInd);
        width = valleyInd - ( -peakLag + 1);
        met.width = width;
    else
        met.width = NaN;
    end
    % Lratio
    met.Lratio = nums(row-1,1);
    % pulseWidth
    met.pulseWidth = pulseWidth;
    % spikeNum
    met.spikeNum = spikeNum;
    %% waveforms
    met.optoWaveform = AD2uV.*mean(waveforms(peakLag+51:peakEnd+51,:,:),3);
    %% session data
    cellNameSession = [unit '.txt'];
    sortedPathSession = [pd.nLynxFolder 'session' sep];
    spikeTimesSession = (load(strcat(sortedPathSession, cellNameSession)))';
    [~, unitNum] = strtok(cellNameSession, 'SS');
    unitNum = unitNum(end);
    TTdir = fullfile(sortedPathSession, tmp_TTname);
    [tt_ts, tt_sig] = Nlx2MatSpike(TTdir, [1 0 0 0 1], 0, 1, 1);

    for j = 1:4
        WaveForm{j} = AD2uV.*squeeze(tt_sig(:, j, ismember(tt_ts, spikeTimesSession)))';
    end
    
    met.isiV = sum(diff(spikeTimesSession)<2000)/(length(spikeTimesSession)-1);
    met.waveform = AD2uV.*mean(tt_sig(:,:,ismember(tt_ts, spikeTimesSession)),3);
    %% plots
if plotFlag
    rasters = figure;
    subplot(4,3,1:3); hold on; title(strcat(session, '_', unit),'Interpreter','none')
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
    
    subplot(4,3,4); hold on;
    xlabel('Pulse'); ylabel('Latency (ms)'); ylim([0 respWin/1000]); xlim([0 p.Results.Pulses+1])
    errorbar(avgSpikeLat/1000, semSpikeLat/1000, 'b', 'LineWidth', 2);
%     errorbar(avgSpikeLatSham/1000, semSpikeLatSham/1000, 'k', 'LineWidth', 2);
%     legend('laser','control');
    ylim([0 90])
    
    subplot(4,3,5); hold on;
    xlabel('Pulse'); ylabel('P(spike)'); ylim([-0.1 1.1]); xlim([0 p.Results.Pulses+1])
    plot(spikeProb, 'b', 'LineWidth', 2);
    plot(spikeProbSham, 'k', 'LineWidth', 2);

    subplot(4,3,6); hold on;
    xlabel('Pulse'); ylabel('spikeNum'); ylim([-0.1 max([spikeNum, spikeNumSham],[],'all')]); xlim([0 p.Results.Pulses+1])
    plot(mean(spikeNum), 'b', 'LineWidth', 2);
    plot(mean(spikeNumSham), 'k', 'LineWidth', 2); 
    
    subplot(4,3, 7:8); hold on
    ylabel('Amplitude (\muV)');
    for j = 1:4
        plotFilled([1:32]+32*(j-1), WaveForm{j}, 'k');
    end
    line([0 128], [0 0], 'color', [0.7, 0.7, 0.7]);
    text(10, 150, sprintf('spikeNumber %d' , length(spikeTimesSession)));
    text(10, -100, sprintf('Lratio %d' , nums(row-1,1)));
    ylim([min(-150, AD2uV*min(mainWaveform)) max(180, AD2uV*max(mainWaveform))])
    
    subplot(4,3,10:11); hold on
    ylabel('Amplitude (\muV)')
    if ~isempty(lightWaveForm)
        if size(lightWaveForm, 3) > 1
            for j = 1:4
                plotFilled([1:length(leg)]+length(leg)*(j-1), AD2uV.*squeeze(lightWaveForm(:,j,:))', 'b');
            end
        else
            for j = 1:4
                plot([1:length(leg)]+length(leg)*(j-1), AD2uV.*squeeze(lightWaveForm(:,j,:))', 'b');
            end
        end
    end
    if ~isempty(spontWaveForm)
        if size(spontWaveForm, 3) > 1
            for j = 1:4
                plotFilled([1:length(leg)]+length(leg)*(j-1), AD2uV.*squeeze(spontWaveForm(:,j,:))', 'k');
            end
        else
            for j = 1:4
                plot([1:length(leg)]+length(leg)*(j-1), AD2uV.*squeeze(spontWaveForm(:,j,:))', 'k');
            end
        end
    end
    for j = 1:4
        plot([1:length(leg)]+length(leg)*(j-1), 50*tmpH(:,j), 'm');
        line([peakEndAll(j)+length(leg)*(j-1) peakEndAll(j)+length(leg)*(j-1)], [-100 100], 'color', 'c')
        line([peakLagAll(j)+length(leg)*(j-1) peakLagAll(j)+length(leg)*(j-1)], [-100 100], 'color', 'm')
        if j == peakChannel
           fill([peakLag+51+length(leg)*(j-1) peakLag+51+length(leg)*(j-1) peakEnd+51+length(leg)*(j-1) peakEnd+51+length(leg)*(j-1)],...
               [50 -50 -50 50], [0.5 1 1], 'LineStyle','none', 'FaceAlpha', 0.5)
            line([valleyInd+(51+peakLag-1)+length(leg)*(j-1) valleyInd+51+peakLag+length(leg)*(j-1)], [-100 100], 'color', 'b', 'LineWidth', 2) 
        end
        line([1+length(leg)*(j-1) length(leg)+length(leg)*(j-1)], [AD2uV*baseline(j) AD2uV*baseline(j)], 'color', 'r', 'LineStyle','--')
    end

    plot([5 5], [70 120],'color', 'k', 'lineWidth',2);
    plot([5 21], [70 70],'color', 'k', 'lineWidth',2);
    text(10,100, '50 uV', 'HorizontalAlignment','left')
    text(15,60, '0.5 ms', 'HorizontalAlignment','center')
    text(10,-120, ['Energy:' num2str(AD2uV*energy) ' Distance:' num2str(distance) ' Width:' num2str(width)], 'HorizontalAlignment','left')
    ylim([min(-150, AD2uV*min(mainWaveform)) max(180, AD2uV*max(mainWaveform))])
    screen = get(0,'Screensize');
    screen(4) = screen(4) - 100;
    set(rasters, 'Position', screen)
    
    subplot(4,3,[9,12]); hold on
    ylabel('density');
    h = histogram(diff(spikeTimesSession),'BinWidth',1000,'Normalization','probability');
    set(gca, 'XScale', 'log')
    line([2000 2000],[0 0.002],'color',[1 0 0]);
    formatSpec = ['RP violation = %.2f ' '%'];
    text(1000, 0.002, [sprintf(formatSpec, 100*sum(diff(spikeTimesSession)<2000)/(length(spikeTimesSession)-1)) '%']);
    title('ISI')
end
if saveFlag
        save([sortedPathSession sep session '_' unit '_met.mat'], 'met')
        if plotFlag
            if exist(savePath)
                saveFigurePDF(rasters,[savePath session '_' unit '.pdf'])
            else
                mkdir(savePath)
                saveFigurePDF(rasters,[savePath session '_' unit '.pdf'])
            end
        end
end


  