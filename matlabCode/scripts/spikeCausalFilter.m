load('F:\tmpData\allUnitAUC.mat');

%% fit all rpes
[root, sep] = currComputer();
col = 'good';
modelName = '5params';

allPe = cell(length(allSessions),1);
for i = 1:length(allSessions)
    session = allSessions{i};
    unit = allUnits{i};
    s = behAnalysisNoPlot_opMD(session, 'simpleFlag', 1);
    pd = parseSessionString_df(session, root, sep);
    sampFile = [pd.animalName col '_' modelName];
    path = [root pd.animalName sep pd.animalName 'sorted' sep 'stan' sep 'bernoulli' sep modelName sep col sep];
    [t,~,noSession] = getStanModelParams_samps(modelName, [path sampFile '.mat'], 4000, 'sessionName', session);
    target = t.pe;
    allPe{i} = t.pe;
    
end

%% use the new window to do new LM fitting
[root, sep] = currComputer();

numBinsPSTH = 6;
tbPSTH = 4;
tfPSTH = 4;
binSizePSTH = 200;
stepSizePSTH = 100;
allPSTH = [];
% filter
tauD = 1000; % ms
tauR = 50; % ms
len = 2500;
delay = 0; % ms

myFilter = causalFilter(len, tauR, tauD, delay);

% activation
cHalf = 2;
k = 2;
peak = 3;
for i = 1:length(allSessions)
    session = allSessions{i};
    unit = allUnits{i};
    s = behAnalysisNoPlot_opMD(session, 'simpleFlag', 1);
    pd = parseSessionString_df(session, root, sep);
    sampFile = [pd.animalName col '_' modelName];
    path = [root pd.animalName sep pd.animalName 'sorted' sep 'stan' sep 'bernoulli' sep modelName sep col sep];
    
    [~, matChoice, ~, slideTime] = getUnitMatChoice(session, unit, tbPSTH, tfPSTH, stepSizePSTH, binSizePSTH);
    % filt and activation
    matChoiceFilt = [];
    for j = 1:length(s.allChoices)
        matChoiceFilt(j,:) = conv(matChoice(j,:), myFilter);
    end
    matChoiceFilt = caActivation(matChoiceFilt, cHalf, k, peak);
    matChoiceFiltDS = zeros(size(matChoiceFilt,1), size(matChoiceFilt,2)/10);
    for j = 1:size(matChoiceFilt,1)
        temp = matChoiceFilt(j,:);
        temp = reshape(temp, 10, []);
        matChoiceFiltDS(j,:) = mean(temp,1);
    end
    target = allPe{i};
    tempPSTH = NaN(numBinsPSTH, size(matChoiceFiltDS, 2));
   
    % PSTH
    edges = [linspace(min(target)- 0.01, 0, 0.5*numBinsPSTH+1) linspace(0, max(target)+ 0.01,0.5*numBinsPSTH+1)];
    edges = [edges(1:0.5*numBinsPSTH+1), edges(0.5*numBinsPSTH+3:end)];
    for j = 1:numBinsPSTH
        if ~isempty(find(target >= edges(j) & target < edges(j+1), 1))
            tempPSTH(j,:) = mean(matChoiceFiltDS(target >= edges(j) & target < edges(j+1),:));
        end
    end
    zsInd = sum(isnan(tempPSTH),2)==0;
    zsTemp = zscore(tempPSTH(zsInd,:), [], 'all');
    tempPSTH(zsInd,:) = zsTemp;
    allPSTH = cat(3, allPSTH, tempPSTH);
    
end
%% causalFilter
% axonDelay = 20; % in ms
% halfLife = 500; % in ms
% len = 2000; % in ms
% a = -log(2)/halfLife;
% decayTime = 0:(len-axonDelay-1);
% decayCurve = exp(a*decayTime);
% myFilter = zeros(1, len);
% myFilter(1:axonDelay) = 0;
% myFilter(axonDelay+1:end) = decayCurve;
% myFilter = myFilter/sum(myFilter);
% myFilter = flip(myFilter)/sum(myFilter);
% new filter 


% filter all traces
% dim = size(allPSTH);
% allPSTHfilt = [];
% binSize = 100;
% time = (1000*-tbPSTH + 0.5*binSize):binSize:1000*tfPSTH;
% for i = 1:dim(1)
%     for j = 1:dim(3)
%         temp = conv(allPSTH(i,:,j), myFilter);
% %         smooth and simplify
%         temp = temp(1:1000*(tbPSTH+tfPSTH));
%         tempSmooth = reshape(temp, [binSize, length(temp)/binSize]);
%         tempSmooth = mean(tempSmooth);
%         allPSTHfilt(i,:,j) = tempSmooth;
%     end
% end
%%
coeffsMax = [lmMaxAll.coeffs]';
tStatsMax = [lmMaxAll.tStats]';
sigMax = [lmMaxAll.ps]'<0.05;
colors = [1 0 0;
          1 0.3 0.3;
          1 0.6 0.6;
%           0.6 0.6 1;
%           1 0.6 0.6;
          0.6 0.6 1;
          0.3 0.3 1;
          0 0 1];
time = linspace(0, size(matChoiceFilt, 2), size(allPSTH,2)) - 1000*tbPSTH;
% plot without activation function
figure2;
subplot(1,2,1);

for i = 1:numBinsPSTH
    plotFilled(time, squeeze(allPSTH(i,:,sigMax(:,1) & tStatsMax(:,1)<0 & ind==1))', colors(i,:));
end
% plot([s.rwdDelay s.rwdDelay], [-0.2 0.4], 'Color', [0.6 0.6 0.6], 'LineWidth', 2, 'LineStyle', '--');
set(gca, 'Box', 'off')
set(gca,'tickdir', 'out')
% set(gca, 'XTick', [-1000:1000:3000])
% set(gca, 'YTick', [-0.2:0.1:0.4])
xlabel('time from choice (ms)', 'FontSize', 18)
ylabel('simulated signal (zscored)', 'FontSize', 18)
xlim([-1000 2500])
% ylim([-0.15 0.25])

title('Type II', 'FontSize', 18)

subplot(1,2,2);
for i = 1:numBinsPSTH
    plotFilled(time, squeeze(allPSTH(i,:,sigMax(:,1) & tStatsMax(:,1)>0 & ind==2))', colors(i,:));
end
% plot([s.rwdDelay s.rwdDelay], [-0.2 0.4], 'Color', [0.6 0.6 0.6], 'LineWidth', 2, 'LineStyle', '--');
set(gca, 'Box', 'off')
set(gca,'tickdir', 'out')
% set(gca, 'XTick', [-1000:1000:3000])
% set(gca, 'YTick', [-0.2:0.1:0.4])
xlabel('time from choice (ms)', 'FontSize', 18)
ylabel('simulated signal (zscored)', 'FontSize', 18)
xlim([-1000 2500])
% ylim([-0.15 0.25])

title('Type I', 'FontSize', 18)
%% plot with activation function
cHalf = 5;
k = 3;
peak = 3;

% filter all traces
dim = size(allPSTH);
allPSTHfilt = [];
binSize = 100;
time = (1000*-tbPSTH):1000*tfPSTH;
for i = 1:dim(1)
    for j = 1:dim(3)
        temp = conv(allPSTH(i,:,j), myFilter);
%         smooth and simplify
        temp = temp(1:1000*(tbPSTH+tfPSTH)+1);
        allPSTHfilt(i,:,j) = temp;
    end
end

% allPSTHfilt = caActivation(allPSTHfilt, cHalf, k, peak);


figure2;
subplot(1,2,1);

for i = 1:numBinsPSTH
    plotFilled(time, squeeze(allPSTHfilt(i,:,sigMax(:,1) & tStatsMax(:,1)<0 & ind==1))', colors(i,:));
end
% plot([s.rwdDelay s.rwdDelay], [-0.2 0.4], 'Color', [0.6 0.6 0.6], 'LineWidth', 2, 'LineStyle', '--');
set(gca, 'Box', 'off')
set(gca,'tickdir', 'out')
set(gca, 'XTick', [-1000:1000:3000])
% set(gca, 'YTick', [-0.2:0.1:0.4])
xlabel('time from choice (ms)', 'FontSize', 18)
ylabel('spikes/s (zscored)', 'FontSize', 18)
xlim([-1000 2500])
% ylim([-0.15 0.25])

title('Type II', 'FontSize', 18)

subplot(1,2,2);
for i = 1:numBinsPSTH
    plotFilled(time, squeeze(allPSTHfilt(i,:,sigMax(:,1) & tStatsMax(:,1)>0 & ind==2))', colors(i,:));
end
% plot([s.rwdDelay s.rwdDelay], [-0.2 0.4], 'Color', [0.6 0.6 0.6], 'LineWidth', 2, 'LineStyle', '--');
set(gca, 'Box', 'off')
set(gca,'tickdir', 'out')
set(gca, 'XTick', [-1000:1000:3000])
% set(gca, 'YTick', [-0.2:0.1:0.4])
xlabel('time from choice (ms)', 'FontSize', 18)
ylabel('spikes/s (zscored)', 'FontSize', 18)
xlim([-1000 2500])
% ylim([-0.15 0.25])

title('Type I', 'FontSize', 18)
%%



%% new filter
tauD = 1500; % ms
tauR = 50; % ms
len = 2000;
delay = 20; % ms

decay = exp(-(1:len)/tauD);
rise = 1 - exp(-(1:len)/tauR);
kernel = [zeros(1, delay) decay.*rise];
kernel = kernel/sum(kernel);










