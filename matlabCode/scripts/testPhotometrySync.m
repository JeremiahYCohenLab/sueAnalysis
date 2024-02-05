folder = 'F:\photometrySyncTest\test2\';
load([folder '\beh.mat']);

fpDir = folder;
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

GSig = readmatrix([fpDir allFiles{GInd}]);
timeStampsG = GSig(:,1);
GSig = GSig(:, 2);
%% detect trial starts
upThresh = 4;
downThresh = 0.5;
upInds = find(ttlSignal(1:end-1)<upThresh & ttlSignal(2:end)>=upThresh)+1;
downInds = find(ttlSignal(1:end-1)>downThresh & ttlSignal(2:end)<=downThresh)+1;
time = linspace(ttlTime(1), ttlTime(end)+1000, length(ttlSignal)+1) - 1000;
time = time(1:length(ttlSignal));
trialStarts = time(upInds);

if length(trialStarts) ~= length(behSessionData)
    fprintf([session ' mis-match with behavior data \n']);
end
%%