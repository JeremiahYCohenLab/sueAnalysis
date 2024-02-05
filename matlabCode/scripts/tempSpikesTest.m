channels = 1:8;
opto = readmatrix('F:\npOptoRecordings\687697_2023-09-13_12-04-06\230913121414_687697.opto.csv');
freq = sum(opto(2:end, [6, 9]), 2);
freq = 1000./freq;
allFreq = unique(freq);
sites = opto(2:end, 2);
allSites = unique(sites);
preTime = 500;
postTime = 1000;
%%
timeBins = round(-preTime/1000 * samplingFreq : (postTime)/1000 * samplingFreq);
timeStamps = timeBins/samplingFreq * 1000;
samples = samples - mean(samples, 1);
samples = double(samples) - mean(samples, 1);
channels = 1:8;
preTime = 500;
postTime = 1000;
timeBins = round(-preTime/1000 * samplingFreq) : round((postTime)/1000 * samplingFreq));
time = timeBins/samplingFreq * 1000;

%%
sortedDir = 'F:\npOptoRecordings\687697_2023-09-13_12-04-06\spikes';
sortedFiles = dir(sortedDir);
sortedFiles = {sortedFiles.name}';
unitInd = cellfun(@(x) contains(x, '.txt'), sortedFiles) & cellfun(@(x) contains(x, 'SS'), sortedFiles);
units = sortedFiles(unitInd);
focusSite = 2;
focusFreq = 5;
currTime = ttl2onTime(freq == focusFreq & sites == focusSite);
for u = 1:length(units)
    currUnit = load([sortedDir '\' units{u}]);
    currUnit = currUnit/1000;
    rasterCell = cell(length(currTime), 1);
    for t = 1:length(currTime)
        rasterCell{t} = currUnit(currUnit > currTime(t) - preTime/1000 & currUnit < currTime(t) + postTime/1000) - currTime(t);
    end
    figure2;
    plotSpikeRaster(rasterCell, 'PlotType','vertline');
    sgtitle(['unit' num2str(units{u}) ' site' num2str(focusSite)], 'Interpreter', 'none')
end
%%



