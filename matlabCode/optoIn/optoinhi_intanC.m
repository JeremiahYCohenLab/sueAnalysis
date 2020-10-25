function optoinhi_intanC(filename, varargin)

p = inputParser;
% default parameters if none given
p.addParameter('Session', 'opto')
p.addParameter('SamplingFreq', 30000)
p.addParameter('Trains', 10)
p.addParameter('PulseWidth', 500)
p.addParameter('MedianRemoval', true)
p.addParameter('scaling_for_int16', 100)
p.addParameter('windowSize', 50)

p.parse(varargin{:});

%get session info
[root, sep] = currComputer();
[animalName, date] = strtok(filename, 'd'); 
animalName = animalName(2:end);
date = date(1:9);
pen = date(6:7);
power = str2double(date(8:9))/10;
sessionFolder = ['m' animalName date];

%specify and make directories
sortedPath = [root animalName sep sessionFolder sep 'sorted' sep p.Results.Session sep];
unsortedPath = [root animalName sep sessionFolder sep 'ephys' sep p.Results.Session sep];
cPath = [root animalName sep sessionFolder sep 'ephys' sep 'spikeSort_' p.Results.Session sep];
saveDir = [root animalName sep sessionFolder sep 'figures'];
if ~exist(saveDir, 'dir')
    mkdir(saveDir)
end

%get sorted and raw ephys data
optoFiles = dir(fullfile(unsortedPath,'*.rhd'));
sortedFiles = dir(fullfile(sortedPath,'*.txt'));

if exist([cPath 'settings.mat'])
    load([cPath 'settings.mat']);
    p.parse('SamplingFreq', frequency_parameters.board_adc_sample_rate);
    p.parse('scaling_for_int16', scaling_for_int16); 
end

%get list of channels with sorted units
for i = 1:length(sortedFiles)
    chanList(i) = sscanf(sortedFiles(i).name,'C%d');
end

%combine traces and DI data from raw files, only from relvant channels
laser = [];
traces = [];
for i = 1:length(optoFiles)
    [digInTmp, tracesTmp, ~] = readIntan(unsortedPath, optoFiles(i).name);
    laser = [laser digInTmp(8,:)];
    traces = [traces tracesTmp];
end
if p.Results.MedianRemoval
    traces = traces - median(traces);
end
traces = traces(chanList,:);
laserIn = find(laser(1:end-1) == 0 &  laser(2:end) > 0) + 1 + 24; %800us offseet,
laser = laserIn/30; %convert to ms. 



%set window and stim paramaters
tB = 500;
tA = 1500;

rasterLength = length(-1*tB:p.Results.PulseWidth+tA);
    

%% 

for i = 1:length(sortedFiles)
    % get spikesIns for each pulse
    [cellName, ~] = strtok(sortedFiles(i).name, '.');
    cell = cellName(end);
    spikeTimes = load(strcat(sortedPath, sortedFiles(i).name))';
    spikeRast = [];
    spikeSlide = zeros(p.Results.Trains,rasterLength);
    window = ones(1,p.Results.windowSize);
    binnedFiring = [];
    for j = 1:p.Results.Trains
        spikeRast{j} = spikeTimes(spikeTimes > (laser(j) - tB) &...
            spikeTimes < (laser(j) + p.Results.PulseWidth + tA));

         if ~isempty(spikeRast{j})
            spikeRast{j} = spikeRast{j} - laser(j); %puts in index relative to laser start 
            spikeSlide(j,round(spikeRast{j} + tB)) = 1;  
            binnedF = conv(spikeSlide(j,:),window);
            binnedFiring(j,:) = binnedF(p.Results.windowSize:end-p.Results.windowSize+1);
        end
    end
    
    
    %get raw traces for inhibition
    rawTrace = [];
    traceTmp = traces(i,:);
    for j = 1:p.Results.Trains
        rawTrace{j} = traceTmp(laserIn(j)-30*tB:laserIn(j) + 30*p.Results.PulseWidth + 30*tA);
    end
    
    %get spike waveforms
    [Cname, unitNum] = strtok(cellName, 'SS');
    Cname = Cname(1:end-1);
    unitNum = unitNum(end);
    Cprev = '';
    if strcmp(Cname, Cprev) == false % if the tetrode has changed, load a new one
        tmp_Cname = [Cname '.ntt'];
        Cdir = fullfile(cPath, tmp_Cname);
        [c_ts, c_cn, c_sig] = Nlx2MatSpike(Cdir, [1 0 1 0 1], 0, 1, 1);
        Cprev = Cname;
    end
    for j = 1:4
        waveform{j} = squeeze(c_sig(:, j, ismember(c_ts, spikeTimes)))' / p.Results.scaling_for_int16;
    end
  
    %% plot everything
    
    rasters = figure; 
    screen = get(0,'Screensize');
    screen(3) = screen(3)/2;
    screen(4) = screen(4)-100;
    set(rasters, 'Position', screen);
    suptitle(strcat(animalName, 'pen', pen, Cname, 'cell', cell, 'power', num2str(power)));
    subplot(4,4,1:4); hold on; 
    xlabel('Time (ms)'); ylabel('Trials')
    LineFormat.Color = 'k'; LineFormat.LineWidth = 1;
    plotSpikeRaster(spikeRast,'PlotType','vertline','XLimForCell',[-1*tB rasterLength-tB],'LineFormat',LineFormat);
    hold on;
    plotShaded([0 p.Results.PulseWidth],[0 0; 1+p.Results.Trains 1+p.Results.Trains],'b');
    
    subplot(4,4,5:8); hold on;
    title('Raw Trace')
     xlabel('Time (ms)'); ylabel('Amplitude (\muV)'); 
     xlim([-tB rasterLength-tB]); ylim([-200 250]); 
     for j = 1:p.Results.Trains
         traceTmp = rawTrace{j};
         plot(-tB:rasterLength-tB-1,-traceTmp(1:30:end));
         hold on;
     end
     
    subplot(4,4,9:12); hold on;
    title(['Binned Firings ' num2str(p.Results.windowSize) 'ms']);
    xlabel('Time (ms)'); ylabel('Firing No.');
    plot(-tB+p.Results.windowSize/2:rasterLength-tB-p.Results.windowSize/2,sum(binnedFiring));
    plotShaded([0 p.Results.PulseWidth],[0 0; max(sum(binnedFiring)) max(sum(binnedFiring))],'b');
 
     
    for j = 1:4
        subplot(4,4,12+j); hold on;
        title('Waveform C1')
        ylabel('Amplitude (\muV)');
        if size(waveform{j}, 1) > 1
            plotFilledStd(1:32, waveform{j}, 'b');
        else
            plot(1:32, waveform{j}, 'b');
        end
    end


     saveFigurePDF(rasters,[saveDir sep animalName '_' cellName '_' p.Results.Session 'ID'])
 end
    
    
    
    
    
    
    