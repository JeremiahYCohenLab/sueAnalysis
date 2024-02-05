function PhotometryPreprocessingAllen(session)
tb = 2000;
tf = 5000;
binSize = 100;
stepSize = 100;
[root, sep] = currComputer();
samplingFreq = 20;
s = behAnalysisNoPlot_opMD(session, 'simpleFlag', 1);  
p = parseSessionString_df(session, root, sep);
fpDir = [p.baseFolder sep 'photometry' sep];
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

if length(trialStarts) ~= length(s.behSessionData)
    fprintf([session ' mis-match with behavior data \n']);
end
%% preprocessing
% truncate
% use timeStampsG as default time
timeFIP = timeStampsG; % in ms 
% detect when recording started and ended
diffSig = max(GSig, [], 2) - min(GSig, [], 2);
startInd = find(diffSig>50, 1);

GSig = GSig(startInd:end,:);
IsoSig = IsoSig(startInd:end,:);
timeFIP = timeFIP(startInd:end,:);
dFF = cell(3,1);
winLeft = samplingFreq * 20; % 20s before sample
p = 10; % percentage
winRight = 0;
for i = 1:numChannels
    GSig(:,i) = denoising(GSig(:,i), samplingFreq);
    IsoSig(:,i) = denoising(IsoSig(:,i), samplingFreq);
    fit = fitlm(zscore(IsoSig(:,i)), zscore(GSig(:,i))); 
    dFF{i} = zscore(GSig(:,i)) - (zscore(IsoSig(:,i))*fit.Coefficients.Estimate(2) +  fit.Coefficients.Estimate(1));
    % baseline = running_percentile_filter(dFF{i}, winLeft, winRight, p); 
    % dFF{i} = dFF{i} - baseline;
end

%% truncate into trials
% aligned to go cue
midPoints = (-tb + 0.5*binSize):stepSize:(tf - 0.5*binSize);
mPFCmat = zeros(length(trialStarts), length(midPoints));
LCNmat = zeros(length(trialStarts), length(midPoints));
LCmat = zeros(length(trialStarts), length(midPoints));

for i = 1:length(trialStarts)
    tmpTime = timeFIP(timeFIP > (trialStarts(i)-tb) & timeFIP <= (trialStarts(i)+tf));
    tmp1 = dFF{1}(timeFIP > (trialStarts(i)-tb) & timeFIP <= (trialStarts(i)+tf));
    tmp2 = dFF{2}(timeFIP > (trialStarts(i)-tb) & timeFIP <= (trialStarts(i)+tf));
    tmp3 = dFF{3}(timeFIP > (trialStarts(i)-tb) & timeFIP <= (trialStarts(i)+tf));
    for j = 1:length(midPoints)
        mPFCmat(i,j) = mean(tmp1(tmpTime>trialStarts(i)+midPoints(j)-0.5*binSize & tmpTime<=trialStarts(i)+midPoints(j)+0.5*binSize));
        LCNmat(i,j) = mean(tmp2(tmpTime>trialStarts(i)+midPoints(j)-0.5*binSize & tmpTime<=trialStarts(i)+midPoints(j)+0.5*binSize));
        LCmat(i,j) = mean(tmp3(tmpTime>trialStarts(i)+midPoints(j)-0.5*binSize & tmpTime<=trialStarts(i)+midPoints(j)+0.5*binSize));
    end
end
% aligned to choice
mPFCmatChoice = zeros(length(s.responseInds), length(midPoints));
LCNmatChoice = zeros(length(s.responseInds), length(midPoints));
LCmatChoice = zeros(length(s.responseInds), length(midPoints));
for i = 1:length(s.responseInds)
    currTime = trialStarts(s.responseInds(i)) + s.lickLat(i);
    tmpTime = timeFIP(timeFIP > (currTime-tb) & timeFIP <= (currTime+tf));
    tmp1 = dFF{1}(timeFIP > (currTime-tb) & timeFIP <= (currTime+tf));
    tmp2 = dFF{2}(timeFIP > (currTime-tb) & timeFIP <= (currTime+tf));
    tmp3 = dFF{3}(timeFIP > (currTime-tb) & timeFIP <= (currTime+tf));
    for j = 1:length(midPoints)     
        mPFCmatChoice(i,j) = mean(tmp1(tmpTime>currTime+midPoints(j)-0.5*binSize & tmpTime<=currTime+midPoints(j)+0.5*binSize));
        LCNmatChoice(i,j) = mean(tmp2(tmpTime>currTime+midPoints(j)-0.5*binSize & tmpTime<=currTime+midPoints(j)+0.5*binSize));
        LCmatChoice(i,j) = mean(tmp3(tmpTime>currTime+midPoints(j)-0.5*binSize & tmpTime<=currTime+midPoints(j)+0.5*binSize));
    end
end

%% save
if ~exist(pd.sortedFolder, 'dir')
    mkdir(pd.sortedFolder)
end
save([pd.sortedFolder session '_photometry.mat'], 'mPFCmat', 'LCNmat', 'LCmat', 'mPFCmatChoice', 'LCNmatChoice', 'LCmatChoice','tb', 'tf', 'frameRate', 'binSize', 'stepSize', 'midPoints');

