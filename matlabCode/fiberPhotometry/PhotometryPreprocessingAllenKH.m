function PhotometryPreprocessingAllenKH(session)
samplingFreq = 20;
tb = 5000;
tf = 10000;
binSize = 100;
stepSize = 100; 
fc = 5;
[root, sep] = currComputer();
pd = parseSessionString_df(session, root, sep);
%% prevent it from over-writing previous data
% if exist([pd.sortedFolder session '_photometry.mat'], "file")
%     fprintf([session ' photometry data preprocessed before. \n']);
%     return
% end
%%
 
p = parseSessionString_df(session, root, sep);
fpDir = [p.baseFolder sep 'photometry' sep];
s = behAnalysisNoPlot_opMD(session, 'simpleFlag', 1);
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

ttlSig1 = ttlSignal(1:3:end);
ttlSig2 = ttlSignal(2:3:end);
ttlSig3 = ttlSignal(3:3:end);



% load signal
GInd = contains(allFiles, 'FIP_DataG');
IsoInd = contains(allFiles, 'FIP_DataIso');
RInd = contains(allFiles, 'FIP_DataR');
numChannels = 3;

GSig = readmatrix([fpDir allFiles{GInd}]);
timeStampsG = GSig(:,1);
GSig = GSig(:, 2:2+numChannels-1);

IsoSig = readmatrix([fpDir allFiles{IsoInd}]);
timeStampsIso = IsoSig(:,1);
IsoSig = IsoSig(:, 2:2+numChannels-1);

RSig = readmatrix([fpDir allFiles{RInd}]);
timeStampsR = RSig(:,1);
RSig = RSig(:, 2:2+numChannels-1);

%% detect trial starts
load([pd.sortedFolder session '_sessionData_behav.mat']);
trialStarts = [behSessionData.CSon];
%% preprocessing
% truncate
% use timeStampsG as default time
timeFIP = timeStampsG; % in ms 
% detect when recording started and ended

startInd = 50;

GSig = GSig(startInd:end,:);
IsoSig = IsoSig(startInd:end,:);
RSig = RSig(startInd:end,:);
timeFIP = timeFIP(startInd:end,:);
dFF = cell(3,1);
winLeft = samplingFreq * 20; % 20s before sample
p = 10; % percentage
winRight = samplingFreq * 20;
allfitsBeta = zeros(3,2);
for i = 1:numChannels
    GSig(:,i) = denoisingKH(GSig(:,i), samplingFreq, fc);
    RSig(:,i) = denoisingKH(RSig(:,i), samplingFreq, fc);
    IsoSig(:,i) = denoisingKH(IsoSig(:,i), samplingFreq, fc);
    fit = fitlm(IsoSig(:,i), zscore(GSig(:,i)));
    dFF{i} = zscore(GSig(:,i)) - fit.predict;
    baseline{i} = running_percentile_filter(dFF{i}, winLeft, winRight, p); 
    dFFbl{i} = dFF{i} - baseline{i};

    [b, a] = butter(2, 1/(samplingFreq/2), "low");
    isoBl = filtfilt(b,a, IsoSig(:,i));   
    isoFit{i} = fit.predict;  
    fitR = fitlm(IsoSig(:,i), zscore(RSig(:,i)));
    dFFR{i} = zscore(RSig(:,i)) - fitR.predict;

    % fitFast = fitlm(IsoSig(:,i)-isoBl, dFFIso{i});
    % dFFIso{i} = dFFIso{i} - fitFast.predict;
    % isoFit{i} = isoFit{i} + fitFast.predict;
    % dFFori{i} = dFF{i};
    if i == 1
        [b, a] = butter(2, 0.01/(samplingFreq/2), "low");
        bl{1} = filtfilt(b,a, IsoSig(:,i));
        sig{1} = filtfilt(b,a, GSig(:,i));
        [b, a] = butter(2, [0.01 0.1]/(samplingFreq/2), "bandpass");
        bl{2} = filtfilt(b,a, IsoSig(:,i));
        sig{2} = filtfilt(b,a, GSig(:,i));
        [b, a] = butter(2, [0.1 1]/(samplingFreq/2), "bandpass");
        bl{3} = filtfilt(b,a, IsoSig(:,i));
        sig{3} = filtfilt(b,a, GSig(:,i));
        [b, a] = butter(2, 1/(samplingFreq/2), "high");
        bl{4} = filtfilt(b,a, IsoSig(:,i));
        sig{4} = filtfilt(b,a, GSig(:,i));
        for j = 1:4
            fit = fitlm(bl{j}, sig{j});
            allfitsBeta(j,:) = fit.Coefficients.Estimate;
            allfitsT(j,:) = fit.Coefficients.tStat;
        end
        figure2;
        subplot(2,2,1);
        plot(allfitsBeta(:,1));
        title('beta0')
        subplot(2,2,2)
        plot(allfitsBeta(:,2));
        title('beta1')
        subplot(2,2,3);
        plot(allfitsT(:,1));
        title('Tstats0')
        subplot(2,2,4)
        plot(allfitsT(:,2));
        title('Tstats1')
        
        sgtitle(session)
    end


end

%% truncate into trials

% aligned to go cue
midPoints = (-tb + 0.5*binSize):stepSize:(tf - 0.5*binSize);
mPFCmatG = zeros(length(trialStarts), length(midPoints));
NACmatG = zeros(length(trialStarts), length(midPoints));
LCmatG = zeros(length(trialStarts), length(midPoints));
mPFCmatR = zeros(length(trialStarts), length(midPoints));
NACmatR = zeros(length(trialStarts), length(midPoints));
LCmatR = zeros(length(trialStarts), length(midPoints));

for i = 1:length(trialStarts)
    tmpTime = timeFIP(timeFIP > (trialStarts(i)-tb) & timeFIP <= (trialStarts(i)+tf));
    tmp1 = dFF{1}(timeFIP > (trialStarts(i)-tb) & timeFIP <= (trialStarts(i)+tf));
    tmp2 = dFF{2}(timeFIP > (trialStarts(i)-tb) & timeFIP <= (trialStarts(i)+tf));
    tmp3 = dFF{3}(timeFIP > (trialStarts(i)-tb) & timeFIP <= (trialStarts(i)+tf));
    for j = 1:length(midPoints)
        mPFCmatG(i,j) = mean(tmp1(tmpTime>trialStarts(i)+midPoints(j)-0.5*binSize & tmpTime<=trialStarts(i)+midPoints(j)+0.5*binSize));
        NACmatG(i,j) = mean(tmp2(tmpTime>trialStarts(i)+midPoints(j)-0.5*binSize & tmpTime<=trialStarts(i)+midPoints(j)+0.5*binSize));
        LCmatG(i,j) = mean(tmp3(tmpTime>trialStarts(i)+midPoints(j)-0.5*binSize & tmpTime<=trialStarts(i)+midPoints(j)+0.5*binSize));
    end
    tmp1 = dFFR{1}(timeFIP > (trialStarts(i)-tb) & timeFIP <= (trialStarts(i)+tf));
    tmp2 = dFFR{2}(timeFIP > (trialStarts(i)-tb) & timeFIP <= (trialStarts(i)+tf));
    tmp3 = dFFR{3}(timeFIP > (trialStarts(i)-tb) & timeFIP <= (trialStarts(i)+tf));
    for j = 1:length(midPoints)
        mPFCmatR(i,j) = mean(tmp1(tmpTime>trialStarts(i)+midPoints(j)-0.5*binSize & tmpTime<=trialStarts(i)+midPoints(j)+0.5*binSize));
        NACmatR(i,j) = mean(tmp2(tmpTime>trialStarts(i)+midPoints(j)-0.5*binSize & tmpTime<=trialStarts(i)+midPoints(j)+0.5*binSize));
        LCmatR(i,j) = mean(tmp3(tmpTime>trialStarts(i)+midPoints(j)-0.5*binSize & tmpTime<=trialStarts(i)+midPoints(j)+0.5*binSize));
    end
end
% aligned to choice
mPFCmatChoiceG = zeros(length(s.responseInds), length(midPoints));
NACmatChoiceG = zeros(length(s.responseInds), length(midPoints));
LCmatChoiceG = zeros(length(s.responseInds), length(midPoints));
mPFCmatChoiceR = zeros(length(s.responseInds), length(midPoints));
NACmatChoiceR = zeros(length(s.responseInds), length(midPoints));
LCmatChoiceR = zeros(length(s.responseInds), length(midPoints));
for i = 1:length(s.responseInds)
    currTime = trialStarts(s.responseInds(i)) + s.lickLat(i);
    tmpTime = timeFIP(timeFIP > (currTime-tb) & timeFIP <= (currTime+tf));
    tmp1 = dFF{1}(timeFIP > (currTime-tb) & timeFIP <= (currTime+tf));
    tmp2 = dFF{2}(timeFIP > (currTime-tb) & timeFIP <= (currTime+tf));
    tmp3 = dFF{3}(timeFIP > (currTime-tb) & timeFIP <= (currTime+tf));
    for j = 1:length(midPoints)     
        mPFCmatChoiceG(i,j) = mean(tmp1(tmpTime>currTime+midPoints(j)-0.5*binSize & tmpTime<=currTime+midPoints(j)+0.5*binSize));
        NACmatChoiceG(i,j) = mean(tmp2(tmpTime>currTime+midPoints(j)-0.5*binSize & tmpTime<=currTime+midPoints(j)+0.5*binSize));
        LCmatChoiceG(i,j) = mean(tmp3(tmpTime>currTime+midPoints(j)-0.5*binSize & tmpTime<=currTime+midPoints(j)+0.5*binSize));
    end
    tmp1 = dFFR{1}(timeFIP > (currTime-tb) & timeFIP <= (currTime+tf));
    tmp2 = dFFR{2}(timeFIP > (currTime-tb) & timeFIP <= (currTime+tf));
    tmp3 = dFFR{3}(timeFIP > (currTime-tb) & timeFIP <= (currTime+tf));
    for j = 1:length(midPoints)     
        mPFCmatChoiceR(i,j) = mean(tmp1(tmpTime>currTime+midPoints(j)-0.5*binSize & tmpTime<=currTime+midPoints(j)+0.5*binSize));
        NACmatChoiceR(i,j) = mean(tmp2(tmpTime>currTime+midPoints(j)-0.5*binSize & tmpTime<=currTime+midPoints(j)+0.5*binSize));
        LCmatChoiceR(i,j) = mean(tmp3(tmpTime>currTime+midPoints(j)-0.5*binSize & tmpTime<=currTime+midPoints(j)+0.5*binSize));
    end
end
frameRate = samplingFreq;
%% save
if ~exist(pd.sortedFolder, 'dir')
    mkdir(pd.sortedFolder)
end

midPoints = midPoints/1000;

figure;
for i = 1:2
    subplot(4,2,i); hold on;
    plot(zscore(GSig(:,i)), 'Color','b')
    plot(isoFit{i})
    plot(zscore(baseline{i}))
    legend({'signal', 'isoBl', 'baseline'})
    title('G')
    subplot(4,2,2+i)
    plot(IsoSig(:,i), 'Color','k')
    title('Iso')
    subplot(4,2,4+i)
    plot(dFFbl{i}, 'Color', 'b')
    title('baseline removal')

    subplot(4,2,6+i)
    plot(dFF{i}, 'Color', 'k')
    title('iso removal')    


end
sgtitle(session)
screen = get(0,'Screensize');
screen(4) = screen(4) - 100;
set(gcf, 'Position', screen);

savePath = pd.saveFigFolder;
saveFigurePDF(gcf,[savePath session 'preCompare.pdf']);
save([pd.sortedFolder session '_photometryCombinewithKH.mat'], 'mPFCmatG', 'NACmatG', 'LCmatG', 'mPFCmatChoiceG', 'NACmatChoiceG', 'LCmatChoiceG', 'mPFCmatR', 'NACmatR', 'LCmatR', 'mPFCmatChoiceR', 'NACmatChoiceR', 'LCmatChoiceR', 'tb', 'tf', 'frameRate', 'binSize', 'stepSize', 'midPoints', 'GSig', 'dFF', 'timeFIP', 'trialStarts');

