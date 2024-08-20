function PhotometryPreprocessingAllen_Gi_Th(session)
samplingFreq = 20;
fc = 9;
tb = 5000;
tf = 10000;
binSize = 100;
stepSize = 100;
[root, sep] = currComputer();
pd = parseSessionString_df(session, root, sep);
% % prevent it from over-writing previous data
% if exist([pd.sortedFolder session '_photometry.mat'], "file")
%     fprintf([session ' photometry data preprocessed before. \n']);
%     return
% end
% %
s = behAnalysisNoPlot_opMD(session, 'simpleFlag', 1);  
p = parseSessionString_df(session, root, sep);
fpDir = [p.baseFolder sep 'photometry' sep];
if ~exist(fpDir, 'dir')
    fprintf([session ' no photometry' '\n']);
    return
end
allFiles = dir(fpDir);
allFiles = {allFiles([allFiles.isdir]==0).name}';
% load ttl file
ttlInd = contains(allFiles, 'TTL') & contains(allFiles, 'csv');
ttlFile = allFiles{ttlInd};
ttlSigInd = contains(allFiles, 'TTL') & ~contains(allFiles, 'csv');
ttlSignalFile = allFiles{ttlSigInd};
ttlTime = readmatrix([fpDir ttlFile]); % time in ms, 1 Hz, updates each second
fileID = fopen([fpDir ttlSignalFile]);
ttlSignal = fread(fileID, 'float64'); % sampling rate 1kHz
fclose(fileID);
% load signal
GInd = contains(allFiles, 'FIP_DataG');
IsoInd = contains(allFiles, 'FIP_DataIso');
numChannels = 3;

GSig = readmatrix([fpDir allFiles{GInd}]);
timeStampsG = GSig(:,1);
GSig = GSig(:, 2:2+numChannels-1);

IsoSig = readmatrix([fpDir allFiles{IsoInd}]);
timeStampsIso = IsoSig(:,1);
IsoSig = IsoSig(:, 2:2+numChannels-1);
%% detect trial starts
upThresh = 4;
downThresh = 0.5;
upInds = find(ttlSignal(1:end-1)<upThresh & ttlSignal(2:end)>=upThresh)+1;
downInds = find(ttlSignal(1:end-1)>downThresh & ttlSignal(2:end)<=downThresh)+1;
time = linspace(ttlTime(1), ttlTime(end)+1000, length(ttlSignal)+1) - 1000; % correction for signal lag
time = time(1:length(ttlSignal));
trialStarts = time(upInds);
trialStartsBeh = [s.behSessionData.CSon];


if length(trialStarts) == length(trialStartsBeh)
    fprintf([session ' perfect match with behavior data \n']);
else
    if length(trialStarts) > length(trialStartsBeh)
        testInds = 1:((length(trialStarts) - length(trialStartsBeh)) + 1);
        testCs = zeros(size(testInds));
        for j = 1:length(testInds)
            currStarts = trialStarts(testInds(j): testInds(j)+length(trialStartsBeh)-1);
            [c, ~] = corrcoef(diff(currStarts), diff(trialStartsBeh));
            testCs(j) = c(1, 2);
        end
        [maxC, bestInd] = max(testCs);
        trialStarts = trialStarts(testInds(bestInd): testInds(bestInd)+length(trialStartsBeh)-1);
        fprintf([session ' re-aligned with behavior data with ' num2str(maxC) '\n']);
    end
end
%% preprocessing
% truncate
% use timeStampsG as default time
timeFIP = timeStampsG; % in ms 
% detect when recording started and ended
diffSig = max(GSig, [], 2) - min(GSig, [], 2);
startInd = find(diffSig>30,1)+1;

GSig = GSig(startInd:end,:);
IsoSig = IsoSig(startInd:end,:);
timeFIP = timeFIP(startInd:end,:);
dFF = cell(3,1);
winLeft = samplingFreq * 60; % used to be 20s before sample
p = 10; % percentage
winRight = samplingFreq * 0;
for i = 1:numChannels
    GSig(:,i) = denoising(GSig(:,i), samplingFreq, fc);
    IsoSig(:,i) = denoising(IsoSig(:,i), samplingFreq, fc);
    [b, a] = butter(2, 0.01/(fc/2), 'low');
    isoLow = filtfilt(b, a, IsoSig(:,i));
    GLow = filtfilt(b, a, zscore(GSig(:,i)));
    fit = fitlm(IsoSig(:,i), zscore(GSig(:,i))); 
    dFF{i} = zscore(GSig(:,i)) - fit.predict;
    bl{i} = running_percentile_filter(dFF{i}, winLeft, winRight, p); 
    dFF{i} = dFF{i} - bl{i};
    isoBl{i} = fit.predict; 
end

%% truncate into trials
% aligned to go cue
midPoints = (-tb + 0.5*binSize):stepSize:(tf - 0.5*binSize);
mPFCmatG = zeros(length(trialStarts), length(midPoints));
THmatG = zeros(length(trialStarts), length(midPoints));
GimatG = zeros(length(trialStarts), length(midPoints));

for i = 1:length(trialStarts)
    tmpTime = timeFIP(timeFIP > (trialStarts(i)-tb) & timeFIP <= (trialStarts(i)+tf));
    tmp1 = dFF{1}(timeFIP > (trialStarts(i)-tb) & timeFIP <= (trialStarts(i)+tf));
    tmp2 = dFF{2}(timeFIP > (trialStarts(i)-tb) & timeFIP <= (trialStarts(i)+tf));
    tmp3 = dFF{3}(timeFIP > (trialStarts(i)-tb) & timeFIP <= (trialStarts(i)+tf));
    for j = 1:length(midPoints)
        mPFCmatG(i,j) = mean(tmp1(tmpTime>trialStarts(i)+midPoints(j)-0.5*binSize & tmpTime<=trialStarts(i)+midPoints(j)+0.5*binSize));
        THmatG(i,j) = mean(tmp2(tmpTime>trialStarts(i)+midPoints(j)-0.5*binSize & tmpTime<=trialStarts(i)+midPoints(j)+0.5*binSize));
        GimatG(i,j) = mean(tmp3(tmpTime>trialStarts(i)+midPoints(j)-0.5*binSize & tmpTime<=trialStarts(i)+midPoints(j)+0.5*binSize));
    end
end
% aligned to choice
mPFCmatChoiceG = zeros(length(s.responseInds), length(midPoints));
THmatChoiceG = zeros(length(s.responseInds), length(midPoints));
GimatChoiceG = zeros(length(s.responseInds), length(midPoints));
for i = 1:length(s.responseInds)
    currTime = trialStarts(s.responseInds(i)) + s.lickLat(i);
    tmpTime = timeFIP(timeFIP > (currTime-tb) & timeFIP <= (currTime+tf));
    tmp1 = dFF{1}(timeFIP > (currTime-tb) & timeFIP <= (currTime+tf));
    tmp2 = dFF{2}(timeFIP > (currTime-tb) & timeFIP <= (currTime+tf));
    tmp3 = dFF{3}(timeFIP > (currTime-tb) & timeFIP <= (currTime+tf));
    for j = 1:length(midPoints)     
        mPFCmatChoiceG(i,j) = mean(tmp1(tmpTime>currTime+midPoints(j)-0.5*binSize & tmpTime<=currTime+midPoints(j)+0.5*binSize));
        THmatChoiceG(i,j) = mean(tmp2(tmpTime>currTime+midPoints(j)-0.5*binSize & tmpTime<=currTime+midPoints(j)+0.5*binSize));
        GimatChoiceG(i,j) = mean(tmp3(tmpTime>currTime+midPoints(j)-0.5*binSize & tmpTime<=currTime+midPoints(j)+0.5*binSize));
    end
end
frameRate = samplingFreq;
%% save
if ~exist(pd.sortedFolder, 'dir')
    mkdir(pd.sortedFolder)
end

midPoints = midPoints/1000;

figure;
num_channels_plot = 3;
for i = 1:num_channels_plot
    subplot(4,num_channels_plot,i);
    hold on;
    plot(GSig(:,i), 'Color','b')
    plot(bl{i})
    plot(isoBl{i})
    legend({'signal', 'bl', 'isobl'})
    title(['G-channel ' num2str(i)])

    subplot(4,num_channels_plot,num_channels_plot+i)
    plot(IsoSig(:,i), 'Color','k')
    title('Iso')

    subplot(4,num_channels_plot,2*num_channels_plot+i)
    plot(GSig(:,i), 'Color', 'b')
    title('iso removal')

    subplot(4,num_channels_plot,3*num_channels_plot+i)
    plot(dFF{i}, 'Color', 'k')
    title('baseline removal')    


end
sgtitle(session)
screen = get(0,'Screensize');
screen(4) = screen(4) - 100;
set(gcf, 'Position', screen);

savePath = pd.saveFigFolder;
saveFigurePDF(gcf,[savePath session 'preCompare.pdf']);

save([pd.sortedFolder session '_photometryGi.mat'], 'mPFCmatG', 'THmatG', 'GimatG', 'mPFCmatChoiceG', 'THmatChoiceG', 'GimatChoiceG','tb', 'tf', 'frameRate', 'binSize', 'stepSize', 'midPoints', 'GSig', 'dFF', 'timeFIP', 'trialStarts');

