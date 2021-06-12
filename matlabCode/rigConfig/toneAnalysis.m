filepath = 'F:\toneTest\20210518\2021-05-18_20-46-27\';
[timestamps, eventID, TTL, Evstring] = Nlx2MatEV([filepath 'Events.nev'],[1 1 1 0 1], 0, 1);
biTTL = de2bi(TTL);
CSplusCol = [4, 8];
CSminusCol = [6, 8];
toneCol = 5;

CSplusInds = biTTL(:,CSplusCol(1))==1 & biTTL(:,CSplusCol(2))==1;
CSplusonInds = [false; ~CSplusInds(1:end-1) & CSplusInds(2:end)];
CSplusoffInds = [false; CSplusInds(1:end-1) & ~CSplusInds(2:end)];
CSplusonTimes = timestamps(CSplusonInds);
CSplusoffTimes = timestamps(CSplusoffInds);
CSduration = CSplusoffTimes-CSplusonTimes;

CSminusInds = biTTL(:,CSminusCol(1))==1 & biTTL(:,CSminusCol(2))==1;
CSminusonInds = [false; ~CSminusInds(1:end-1) & CSminusInds(2:end)];
CSminusonTimes = timestamps(CSminusonInds);

toneOnTime = timestamps(biTTL(:,toneCol)==1);
toneSwitchOnTime = timestamps([false; biTTL(1:end-1,toneCol)~=1 & biTTL(2:end,toneCol)==1]);
%%
toneLatCSplus = zeros(size(CSplusonTimes));
toneLengthCSplus = zeros(size(CSplusonTimes));
toneFreqCSplus = zeros(size(CSplusonTimes));
toneLatCSminus = zeros(size(CSminusonTimes));
toneLengthCSminus = zeros(size(CSminusonTimes));
toneFreqCSminus = zeros(size(CSminusonTimes));

for i = 1:length(CSplusonTimes)
    toneTemp = toneOnTime(toneOnTime>=CSplusonTimes(i)-500 & toneOnTime < CSplusonTimes(i) + 1000000);
    toneSwitchTemp = toneSwitchOnTime(toneSwitchOnTime>=CSplusonTimes(i)-500 & toneSwitchOnTime < CSplusonTimes(i) + 1000000);
    toneLatCSplus(i) = min(toneTemp - CSplusonTimes(i)); % in us
    toneLengthCSplus(i) = max(toneTemp)- min(toneTemp); % in us
    toneFreqCSplus(i) = 1000*length(toneSwitchTemp)/toneLengthCSplus(i); % in kHz
end

for i = 1:length(CSminusonTimes)
    toneTemp = toneOnTime(toneOnTime>=CSminusonTimes(i)-500 & toneOnTime < CSminusonTimes(i) + 1000000);
    toneLatCSminus(i) = min(toneTemp - CSminusonTimes(i)); % in us
    toneLengthCSminus(i) = max(toneTemp)-min(toneTemp); % in us
    toneFreqCSminus(i) = 1000*length(toneTemp)/toneLengthCSminus(i); % in kHz
end
%%
figure2('Position',[0 0 300 900]);
subplot(3,1,1); hold on;
histogram(toneLengthCSplus/1000, 0:25:550, 'Facecolor', 'c', 'Normalization', 'Probability');
histogram(toneLengthCSminus/1000, 0:25:550, 'Facecolor', 'm', 'Normalization', 'Probability');
xlabel('ms')
subplot(3,1,2); hold on;
histogram(toneLatCSplus/1000, 0:0.01:0.2, 'Facecolor', 'c', 'Normalization', 'Probability');
histogram(toneLatCSminus/1000, 0:0.01:0.2, 'Facecolor', 'm', 'Normalization', 'Probability');
xlabel('ms')
subplot(3,1,3); hold on;
histogram(toneFreqCSplus, 7:0.5:16, 'Facecolor', 'c', 'Normalization', 'Probability');
histogram(toneFreqCSminus, 7:0.5:16, 'Facecolor', 'm', 'Normalization', 'Probability');
xlabel('kHz')
%%
figure2;
scatter3(toneLengthCSplus/1000, toneLatCSplus/1000, toneFreqCSplus);
xlabel('len')
ylabel('lat')
zlabel('freq')
%% rasters
allEventsCSplus = cell(length(CSplusonTimes),1);
for i = 1:length(CSplusonTimes)
    timeTmp = timestamps(timestamps>=(CSplusonTimes(i)-500)&timestamps<=(CSplusonTimes(i)+1000000))- CSplusonTimes(i);
    allEventsCSplus{i}=timeTmp;
end
figure;
plotSpikeRaster(allEventsCSplus(toneLengthCSplus<400000,1),'PlotType','vertline');
figure;
plotSpikeRaster(allEventsCSplus(toneLengthCSplus>400000,1),'PlotType','vertline');

%%
allEventsCSminus = cell(length(CSminusonTimes),1);
for i = 1:length(CSminusonTimes)
    timeTmp = timestamps(timestamps>=(CSminusonTimes(i)-500)&timestamps<=(CSminusonTimes(i)+1000000))- CSminusonTimes(i);
    allEventsCSminus{i}=timeTmp;
end
plotSpikeRaster(allEventsCSminus(toneLengthCSminus<400000,1),'PlotType','vertline');
figure;
plotSpikeRaster(allEventsCSminus(toneLengthCSminus>400000,1),'PlotType','vertline');
%%
filePath = 'F:\toneTest\20210521\';
fileName = 'ToneTest.m4a';
[y,Fs] = audioread([filePath fileName]);
sampleTime = 1000000/Fs; % in us 
time = 0:sampleTime:sampleTime*(length(y)-1);

CSon = [behSessionData.CSon]*1000;
CSInd = find(y(6:end)>0.3 & abs(y(3:end-3))<0.2  & abs(y(2:end-4))<0.2  & abs(y(1:end-5))<0.2)+5;
CSon = CSon-CSon(1)+time(CSInd(1));

tb = 1000000;
tf = 1000000;
tempFrame = -tb:sampleTime:tf-sampleTime;
toneCell = cell(length(CSon),1);
for i = 1:length(CSon)
    temp = y(time>=CSon(i)-tb & time <= CSon(i)+tf);
    toneCell{i} = temp(1:length(tempFrame));
end
%%
figure2;hold on;
CSplusInds = find(contains({behSessionData.trialType}, 'CSplus'));
for i = 1:length(CSplusInds)
    plot(tempFrame, i+2.5*toneCell{CSplusInds(i)}, 'k','LineWidth',0.1);
end
ylim([0 length(CSplusInds)+1])

figure2; hold on;
CSminusInds = find(contains({behSessionData.trialType}, 'CSminus'));
for i = 1:length(CSminusInds)
    plot(tempFrame, i+2.5*toneCell{CSminusInds(i)}, 'k','LineWidth',0.1);
end
ylim([0 length(CSminusInds)+1])
%%
    

