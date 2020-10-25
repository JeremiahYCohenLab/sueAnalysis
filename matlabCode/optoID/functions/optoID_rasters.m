function optoStruct = optoID_rasters(folder, varargin)

% set appropriate settings
[root, sep] = currComputer();
p = inputParser;
% default parameters if none given
p.addParameter('SaveFigFlag', true);
p.addParameter('SaveStructFlag', true);
p.addParameter('Root', root)
p.addParameter('Separator', sep)
p.addParameter('Sequences', [10 25 50])
p.addParameter('TrainsPerSequence', 10)
p.addParameter('PulsesPerTrain', 10)
p.addParameter('HighPassCutoffInHz', 300)
p.addParameter('LowPassCutoffInHz', []);
p.addParameter('SamplingFreq', 32000)
p.addParameter('LightEvokedMaxLat', 20);
p.addParameter('SpontSpike_LightRefracPeriod', 100);
p.addParameter('OptoFolder','')
p.addParameter('WaveformSamples', 192);
p.addParameter('Raster_TBack', -500);
p.addParameter('Raster_TFor', 1500);
p.addParameter('FilterTrace_Flag', true);

pathData = parseSessionString_df(folder, root, sep);
if str2double(pathData.date) >= 20180322 % 5ms pulse width before this, 3ms pulse width on and after
    p.addParameter('PulseDuration', 10);
else
    p.addParameter('PulseDuration', 5);
end
p.parse(varargin{:});

shutterOffset = 0.8; % 295H neuralynx rig off by 0.8ms how to measure this?

wfSamp_left = floor(1/2*p.Results.WaveformSamples);
wfSamp_right = ceil(1/2*p.Results.WaveformSamples) - 1;
% wfMax_expected = wfSamp_left - 2;
wfMax_expected = wfSamp_left + 1;



nLynxFolder = [pathData.baseFolder p.Results.OptoFolder 'neuralynx\'];
sortedFolder = [pathData.baseFolder p.Results.OptoFolder 'sorted\opto\'];
saveFolder = [pathData.baseFolder p.Results.OptoFolder 'figures\'];
if ~isdir(saveFolder)
    mkdir(saveFolder);
end
sortedFiles = dir(sortedFolder);
% load tetrode data
for i = find(contains({sortedFiles.name}, 'TT') & contains({sortedFiles.name}, 'txt'))
    load([sortedFolder sortedFiles(i).name]);
end
allVar = who;
allClust = find(contains(allVar,'SS') & contains(allVar,'TT'));
TTstruct = [];
for m = allClust'
    eval(sprintf('TTstruct.%s = %s/1000;', allVar{m}, allVar{m}))
end
    
[tsEv, Ev] = Nlx2MatEV([nLynxFolder 'Events.nev'], [1 0 1 0 0], 0, 1);
% TTL events
waterR = 1; % pin 10
waterL = 2; % pin 11
lickR = 32; % pin 15
lickL = 16; % pin 14
laser = 4;
laserRaw = 1024; %16384; % raw data file uses this as laser TTL 

% isolate just laser on times
Ev(Ev == laser+waterR | Ev == laser+waterL | Ev == laser+lickL | Ev == laser+lickR) = laser;
Ev(Ev == waterR | Ev == waterL | Ev == lickR | Ev == lickL) = 0;
laserEventInds = find(Ev == laser | Ev == laserRaw);
consecLaserInds = diff(laserEventInds) == 1; % laserOn TTLs without laserOff TTLs
Ev(laserEventInds(consecLaserInds) + 1) = 0; % set all repeats to 0
laserOn_mask = ismember(Ev, laser) | ismember(Ev, laserRaw);
% only for rig 295F
laserOn_mask = [laserOn_mask(2:end) false]; % offset by an index because, for some reason, _0_ is the laser on in this rig
laserOnTimes = tsEv(laserOn_mask);
laserOnTimes = laserOnTimes / 1e3; % to ms
laserOnTimes = laserOnTimes + shutterOffset;

if isempty(p.Results.LowPassCutoffInHz)
    Wn = p.Results.HighPassCutoffInHz / (p.Results.SamplingFreq/2);
    [b, a] = butter(2, Wn, 'high');
else
    Wn = [p.Results.HighPassCutoffInHz p.Results.LowPassCutoffInHz] / (p.Results.SamplingFreq/2);
    [b, a] = butter(2, Wn, 'bandpass');
end

TT_last = 0;
for currN = fields(TTstruct)'
    currN = currN{:};
    TTnum = currN(3:regexp(currN, '_SS_') - 1);
    TTnum = str2double(TTnum);
    
    chan0 = TTnum + (3*(TTnum - 1));
    chan1 = chan0 + 1;
    chan2 = chan0 + 2;
    chan3 = chan0 + 3;
    
    fprintf('Currently on %s.\n', currN)
    
    if TT_last ~= TTnum % only redo all of this if on a different TT
    
        [ts, samp0, head0] = Nlx2MatCSC([nLynxFolder 'CSC' num2str(chan0) '.ncs'], [1 0 0 0 1], 1, 1, []);
        [samp1, head1] = Nlx2MatCSC([nLynxFolder 'CSC' num2str(chan1) '.ncs'], [0 0 0 0 1], 1, 1, []);
        [samp2, head2] = Nlx2MatCSC([nLynxFolder 'CSC' num2str(chan2) '.ncs'], [0 0 0 0 1], 1, 1, []);
        [samp3, head3] = Nlx2MatCSC([nLynxFolder 'CSC' num2str(chan3) '.ncs'], [0 0 0 0 1], 1, 1, []);
    
        if p.Results.FilterTrace_Flag == true
            samp0 = filtfilt(b, a, samp0(:));
            samp1 = filtfilt(b, a, samp1(:));
            samp2 = filtfilt(b, a, samp2(:));
            samp3 = filtfilt(b, a, samp3(:));
        else
            samp0 = samp0(:);
            samp1 = samp1(:);
            samp2 = samp2(:);
            samp3 = samp3(:);
        end

        tSamp = 1/p.Results.SamplingFreq * 1e6; % time per sample in microseconds

        if length(unique(diff(ts))) == 1
            ts_interp = ts(1):tSamp:ts(1) + tSamp*(length(samp0) - 1);
        else
            ts_interp = NaN(1, length(samp0));
            for i = 1:length(ts)
                ts_interp(512*(i - 1) + 1:512*(i)) = ts(i):tSamp:ts(i) + tSamp*511;
            end
        end
        ts_interp = ts_interp / 1e3; % time in milliseconds
        t_indivWF = (-wfSamp_left:wfSamp_right)*tSamp / 1e3; % time for each waveform; 0 is the spike time
    end
    
    % spontaneous waveform - exclude from laserOn to laserOn + 100ms
    spont_excludeInds = [];
    for lightOn_time = laserOnTimes
        spont_excludeInds = [spont_excludeInds find(ts_interp - lightOn_time >= 0, 1):find(ts_interp - (lightOn_time + p.Results.SpontSpike_LightRefracPeriod) >= 0, 1)];
    end
    spont_includeInds = [];    
    
    wf0 = NaN(length(TTstruct.(currN)), p.Results.WaveformSamples);
    wf1 = NaN(length(TTstruct.(currN)), p.Results.WaveformSamples);
    wf2 = NaN(length(TTstruct.(currN)), p.Results.WaveformSamples);
    wf3 = NaN(length(TTstruct.(currN)), p.Results.WaveformSamples);
    t_spikeAll = NaN(length(TTstruct.(currN)), 1);
    for currSpike_ind = 1:length(TTstruct.(currN))
        t_spike = TTstruct.(currN)(currSpike_ind);
        ind_closest = find(ts_interp - t_spike >= 0, 1); % t_spike and ts_interp are off slightly; this finds the nearest index
        
        if ind_closest + wfSamp_right <= length(ts_interp) && ind_closest - wfSamp_left >= wfSamp_left % spike is within range of data
            if ~any(ind_closest == spont_excludeInds) % if this spike was probably not light-evoked, include it as a spontaneous spike
                spont_includeInds = [spont_includeInds currSpike_ind];
            end
            
            t_wf0 = samp0(ind_closest - wfSamp_left:ind_closest + wfSamp_right);
            t_wf1 = samp1(ind_closest - wfSamp_left:ind_closest + wfSamp_right);
            t_wf2 = samp2(ind_closest - wfSamp_left:ind_closest + wfSamp_right);
            t_wf3 = samp3(ind_closest - wfSamp_left:ind_closest + wfSamp_right);

            [~, maxInd] = max(max([t_wf0'; t_wf1'; t_wf2'; t_wf3'])); % find index of peak; in case it's jittered slightly
            if maxInd ~= wfMax_expected
                ind_closest = ind_closest - (wfMax_expected - maxInd);
                
                t_wf0 = samp0(ind_closest - wfSamp_left:ind_closest + wfSamp_right);
                t_wf1 = samp1(ind_closest - wfSamp_left:ind_closest + wfSamp_right);
                t_wf2 = samp2(ind_closest - wfSamp_left:ind_closest + wfSamp_right);
                t_wf3 = samp3(ind_closest - wfSamp_left:ind_closest + wfSamp_right);
            end
            t_spikeAll(currSpike_ind) = ts_interp(ind_closest);
            wf0(currSpike_ind, :) = t_wf0;
            wf1(currSpike_ind, :) = t_wf1;
            wf2(currSpike_ind, :) = t_wf2;
            wf3(currSpike_ind, :) = t_wf3;
        end
    end
    
    wf0_bkg = wf0(spont_includeInds, :);
    wf1_bkg = wf1(spont_includeInds, :);
    wf2_bkg = wf2(spont_includeInds, :);
    wf3_bkg = wf3(spont_includeInds, :);
    
    % find where waveform significantly deviates from 0
    wfAll_bkg = [mean(wf0_bkg)' mean(wf1_bkg)' mean(wf2_bkg)' mean(wf3_bkg)'];
    [~, maxWfTet_ind] = max(max(wfAll_bkg));
    maxWfTet_ind = maxWfTet_ind - 1; % same indices as tetrodes
    switch maxWfTet_ind
        case 0
            wf_relev = wf0_bkg;
        case 1
            wf_relev = wf1_bkg;
        case 2
            wf_relev = wf2_bkg;
        case 3
            wf_relev = wf3_bkg;
    end
%     sigThresh = 5;
%     sigFound = false;
%     sigCount = 0;
    [~, sigInd] = max(mean(wf_relev));
    
    zeroCrossings = find(diff(mean(wf_relev) > 0) == 1); % all crossings from negative to positive
    sigInd = zeroCrossings(find(zeroCrossings < sigInd, 1, 'last'));
%     while sigFound == false && sigInd > 1
% %         for currSamp = 1:size(wf_relev, 2)
% %             if ttest(wf_relev(:, currSamp), 0, 'Alpha', 1e-4', 'Tail', 'right') == 1 % increment by 1 if significant
% %                 sigCount = sigCount + 1;
% %             else
% %                 sigCount = 0;
% %             end
% %             if sigCount == sigThresh
% %                 sigFound = true;
% %                 break
% %             end
% %         end
%         if ttest(wf_relev(:, sigInd), 0, 'Alpha', 1e-4, 'Tail', 'right') == 1
%             sigInd = sigInd - 1;
%         else
%             sigFound = true;
%         end
% %         if sigFound == false
% %             sigFound = true;
% %         end
%     end
%     tOffset_sig = (find(t_indivWF == 0) - (currSamp - sigThresh + 1)) * (tSamp / 1e3);
    tOffset_sig = (find(t_indivWF == 0) - (sigInd)) * (tSamp / 1e3);
%     tOffset_sig = 0;
    t_spikeAll = t_spikeAll - tOffset_sig;
   
    for currSeq = 1:length(p.Results.Sequences)
        latencyToSpike = NaN(p.Results.TrainsPerSequence, p.Results.PulsesPerTrain);
        responseFraction = NaN(p.Results.TrainsPerSequence, p.Results.PulsesPerTrain);
        wf0_lht = [];
        wf1_lht = [];
        wf2_lht = [];
        wf3_lht = [];
        rasterInfo = cell(1, p.Results.TrainsPerSequence);
        
        for currTrain = 1:p.Results.TrainsPerSequence
            relevLaserInds = (1:p.Results.PulsesPerTrain) + (currTrain - 1)*p.Results.PulsesPerTrain + (currSeq - 1)*p.Results.TrainsPerSequence*p.Results.PulsesPerTrain;
            ref = laserOnTimes(relevLaserInds(1));
            rasterInfo(currTrain) = {(t_spikeAll(t_spikeAll - ref > p.Results.Raster_TBack & t_spikeAll - ref < p.Results.Raster_TFor) - ref)'};
            
            pulseOnTimes = laserOnTimes(relevLaserInds);
            if currTrain == 1
                firstPulseTime = pulseOnTimes(1);
            end
            for currPulse = 1:length(pulseOnTimes)
                nearestSpikeInd = find(t_spikeAll > pulseOnTimes(currPulse), 1);
                if currTrain == 1 && currPulse == 1
                    firstSpikeInd = nearestSpikeInd;
                elseif currTrain == p.Results.TrainsPerSequence && currPulse == p.Results.PulsesPerTrain
                    lastSpikeInd = nearestSpikeInd;
                end
                nearestSpike = t_spikeAll(nearestSpikeInd);
                if nearestSpike - pulseOnTimes(currPulse) <= p.Results.LightEvokedMaxLat
                    latencyToSpike(currTrain, currPulse) = nearestSpike - pulseOnTimes(currPulse);
                    responseFraction(currTrain, currPulse) = true;
                    wf0_lht = [wf0_lht; wf0(nearestSpikeInd, :)];
                    wf1_lht = [wf1_lht; wf1(nearestSpikeInd, :)];
                    wf2_lht = [wf2_lht; wf2(nearestSpikeInd, :)];
                    wf3_lht = [wf3_lht; wf3(nearestSpikeInd, :)];
                else
                    responseFraction(currTrain, currPulse) = false;
                end
            end
        end
        
        % plot figure        
        optoID_fig = figure;
        rasterPlot = subplot(2,3,[1 2 3]); hold on
        latencyPlot = subplot(2,3,4); hold on
        responseFracPlot = subplot(2,3,5); hold on
        waveformPlot = subplot(2,3,6); hold on
        
        subplot(rasterPlot)
        title(sprintf('Session: %s. Neuron: %s. Freq: %i Hz', pathData.sessionFolder, currN, p.Results.Sequences(currSeq)), ...
            'interpreter', 'none')
        rasterInfo(cellfun(@isempty,rasterInfo)) = {zeros(1,0)};
        if ~all(cellfun(@isempty, rasterInfo))
            plotSpikeRaster(rasterInfo, 'PlotType','vertline'); hold on
        end
        xlim([p.Results.Raster_TBack p.Results.Raster_TFor])
        xlabel('Time - Train Onset (ms)')
        ylabel('Train Number')
        for laserLine = 1:p.Results.PulsesPerTrain
            x = (laserLine - 1)*(1/p.Results.Sequences(currSeq))*1e3;
            x = [x x+p.Results.PulseDuration];
            fill([x fliplr(x)], [0 0 p.Results.TrainsPerSequence p.Results.TrainsPerSequence] + 0.5, 'b', 'edgecolor', 'none', 'facealpha', 0.1);
        end
        ylim([0 p.Results.TrainsPerSequence] + 0.5)
        
        subplot(latencyPlot)
        plot(nanmean(latencyToSpike)); ylim([0 p.Results.LightEvokedMaxLat]);
        xlabel('Pulse Number')
        ylabel('Latency to Light-Evoked Spike (ms)')
        xlim([1 p.Results.PulsesPerTrain])
        
        subplot(responseFracPlot)
        plot(mean(responseFraction)); ylim([0 1])
        xlabel('Pulse Number')
        ylabel('P(Spike)')
        
        scaleInd = contains(head0,'-ADBitVolts');
        wf0_uV = str2double(head0{scaleInd}(13:end))*1e6;
        wf1_uV = str2double(head1{scaleInd}(13:end))*1e6;
        wf2_uV = str2double(head2{scaleInd}(13:end))*1e6;
        wf3_uV = str2double(head3{scaleInd}(13:end))*1e6;
        
        if size(wf0_lht, 1) == 1 % only 1 spike
            wf0_lht_mean = wf0_lht*wf0_uV;
            wf1_lht_mean = wf1_lht*wf1_uV;
            wf2_lht_mean = wf2_lht*wf2_uV;
            wf3_lht_mean = wf3_lht*wf3_uV;
            wf0_lht_sem = zeros(1, size(wf0_lht, 2)); 
            wf1_lht_sem = zeros(1, size(wf1_lht, 2)); 
            wf2_lht_sem = zeros(1, size(wf2_lht, 2)); 
            wf3_lht_sem = zeros(1, size(wf3_lht, 2)); 
        else
            wf0_lht_mean = nanmean(wf0_lht*wf0_uV);
            wf1_lht_mean = nanmean(wf1_lht*wf1_uV);
            wf2_lht_mean = nanmean(wf2_lht*wf2_uV);
            wf3_lht_mean = nanmean(wf3_lht*wf3_uV);
            wf0_lht_sem = sem(wf0_lht*wf0_uV);
            wf1_lht_sem = sem(wf1_lht*wf1_uV);
            wf2_lht_sem = sem(wf2_lht*wf2_uV);
            wf3_lht_sem = sem(wf3_lht*wf3_uV);
        end
        if size(wf0_bkg, 1) == 1 % only 1 spike
            wf0_bkg_mean = wf0_bkg*wf0_uV;
            wf1_bkg_mean = wf1_bkg*wf1_uV;
            wf2_bkg_mean = wf2_bkg*wf2_uV;
            wf3_bkg_mean = wf3_bkg*wf3_uV;
            wf0_bkg_sem = zeros(1, size(wf0_bkg, 2));
            wf1_bkg_sem = zeros(1, size(wf1_bkg, 2));
            wf2_bkg_sem = zeros(1, size(wf2_bkg, 2));
            wf3_bkg_sem = zeros(1, size(wf3_bkg, 2));
        else
            wf0_bkg_mean = nanmean(wf0_bkg*wf0_uV);
            wf1_bkg_mean = nanmean(wf1_bkg*wf1_uV);
            wf2_bkg_mean = nanmean(wf2_bkg*wf2_uV);
            wf3_bkg_mean = nanmean(wf3_bkg*wf3_uV);
            wf0_bkg_sem = sem(wf0_bkg*wf0_uV);
            wf1_bkg_sem = sem(wf1_bkg*wf1_uV);
            wf2_bkg_sem = sem(wf2_bkg*wf2_uV);
            wf3_bkg_sem = sem(wf3_bkg*wf3_uV);
        end
        
        subplot(waveformPlot)
        ylabel('Voltage (uV)')
        spacing = p.Results.WaveformSamples + 8;
        t_wf0 = (1:p.Results.WaveformSamples) + spacing*0;
        t_wf1 = (1:p.Results.WaveformSamples) + spacing*1;
        t_wf2 = (1:p.Results.WaveformSamples) + spacing*2;
        t_wf3 = (1:p.Results.WaveformSamples) + spacing*3;
        if all(~isnan(wf0_lht_mean)) % any spikes
            plot(t_wf0, wf0_lht_mean, 'b')
            fill([t_wf0 fliplr(t_wf0)], [wf0_lht_mean + wf0_lht_sem fliplr([wf0_lht_mean - wf0_lht_sem])], 'b', 'edgecolor', 'none', 'facealpha', 0.35);
            
            plot(t_wf1, wf1_lht_mean, 'b')
            fill([t_wf1 fliplr(t_wf1)], [wf1_lht_mean + wf1_lht_sem fliplr([wf1_lht_mean - wf1_lht_sem])], 'b', 'edgecolor', 'none', 'facealpha', 0.35);
            
            plot(t_wf2, wf2_lht_mean, 'b')
            fill([t_wf2 fliplr(t_wf2)], [wf2_lht_mean + wf2_lht_sem fliplr([wf2_lht_mean - wf2_lht_sem])], 'b', 'edgecolor', 'none', 'facealpha', 0.35);
            
            plot(t_wf3, wf3_lht_mean, 'b')
            fill([t_wf3 fliplr(t_wf3)], [wf3_lht_mean + wf3_lht_sem fliplr([wf3_lht_mean - wf3_lht_sem])], 'b', 'edgecolor', 'none', 'facealpha', 0.35);
        end
        if all(~isnan(wf0_bkg_mean)) % any spikes
            plot(t_wf0, wf0_bkg_mean, 'k')
            fill([t_wf0 fliplr(t_wf0)], [wf0_bkg_mean + wf0_bkg_sem fliplr([wf0_bkg_mean - wf0_bkg_sem])], 'k', 'edgecolor', 'none', 'facealpha', 0.35);
            
            plot(t_wf1, wf1_bkg_mean, 'k')
            fill([t_wf1 fliplr(t_wf1)], [wf1_bkg_mean + wf1_bkg_sem fliplr([wf1_bkg_mean - wf1_bkg_sem])], 'k', 'edgecolor', 'none', 'facealpha', 0.35);
            
            plot(t_wf2, wf2_bkg_mean, 'k')
            fill([t_wf2 fliplr(t_wf2)], [wf2_bkg_mean + wf2_bkg_sem fliplr([wf2_bkg_mean - wf2_bkg_sem])], 'k', 'edgecolor', 'none', 'facealpha', 0.35);
            
            plot(t_wf3, wf3_bkg_mean, 'k')
            fill([t_wf3 fliplr(t_wf3)], [wf3_bkg_mean + wf3_bkg_sem fliplr([wf3_bkg_mean - wf3_bkg_sem])], 'k', 'edgecolor', 'none', 'facealpha', 0.35);
            switch maxWfTet_ind
                case 0
%                     plot(t_wf0(currSamp - sigThresh + 1), wf0_bkg_mean(currSamp - sigThresh + 1), 'ro');
                    plot(t_wf0(sigInd), wf0_bkg_mean(sigInd), 'ro');
                case 1
                    plot(t_wf1(sigInd), wf1_bkg_mean(sigInd), 'ro');
                case 2
                    plot(t_wf2(sigInd), wf2_bkg_mean(sigInd), 'ro');
                case 3
                    plot(t_wf3(sigInd), wf3_bkg_mean(sigInd), 'ro');
            end
        end

        xlim([t_wf0(1) t_wf3(end)])
        
        set(gcf, 'Position', get(0,'Screensize'))
        if p.Results.SaveFigFlag == true
            saveFigurePDF(optoID_fig,[saveFolder '\' strrep(pathData.sessionFolder, '\', '') '_' currN '_' int2str(p.Results.Sequences(currSeq)) 'Hz'])
        end
        pause(0.1)
        
        optoStruct.(currN).bkg.wf0 = wf0_bkg;
        optoStruct.(currN).bkg.wf1 = wf1_bkg;
        optoStruct.(currN).bkg.wf2 = wf2_bkg;
        optoStruct.(currN).bkg.wf3 = wf3_bkg;
        optoStruct.(currN).(['freq' num2str(p.Results.Sequences(currSeq)) 'Hz']).lat = latencyToSpike;
        optoStruct.(currN).(['freq' num2str(p.Results.Sequences(currSeq)) 'Hz']).rf = responseFraction;
        optoStruct.(currN).(['freq' num2str(p.Results.Sequences(currSeq)) 'Hz']).lht.wf0 = wf0_lht;
        optoStruct.(currN).(['freq' num2str(p.Results.Sequences(currSeq)) 'Hz']).lht.wf1 = wf1_lht;
        optoStruct.(currN).(['freq' num2str(p.Results.Sequences(currSeq)) 'Hz']).lht.wf2 = wf2_lht;
        optoStruct.(currN).(['freq' num2str(p.Results.Sequences(currSeq)) 'Hz']).lht.wf3 = wf3_lht;
    end
    optoStruct.(currN).bkg.wf0 = wf0_bkg;
    optoStruct.(currN).bkg.wf1 = wf1_bkg;
    optoStruct.(currN).bkg.wf2 = wf2_bkg;
    optoStruct.(currN).bkg.wf3 = wf3_bkg;
    optoStruct.(currN).tOffset_sig = tOffset_sig;
    
    TT_last = TTnum;
    
    close all
end

optoStruct.metaData.shutterOffset = shutterOffset;
optoStruct.metaData.parameterSet = p.Results;
optoStruct.metaData.wfSamp_left = wfSamp_left;
optoStruct.metaData.wfSamp_right = wfSamp_right;
optoStruct.metaData.wfMax_expected = wfMax_expected;

if p.Results.SaveStructFlag == true
    save([sortedFolder strrep(pathData.sessionFolder, '\', '') '_optoStruct.mat'], 'optoStruct');
end

fprintf('Finished generating optoID rasters\n');


