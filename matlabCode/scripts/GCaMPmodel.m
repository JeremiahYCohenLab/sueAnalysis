clear all;
load 'F:\tmpData\GCaMP.mat';
%% prepare data for fitting
N = length(allOutcome);
Tsesh = cellfun(@length, allOutcome);
T = max(Tsesh);

rpeMat = zeros(N, T);
nRwdMat = zeros(N, T);
signalMat = zeros(N ,T);

for i = 1:N
    rpeMat(i, 1:Tsesh(i)) = allPe{i}';
    nRwd = -allOutcome{i}';
    nRwd(nRwd==0) = 1;
    nRwdMat(i, 1:Tsesh(i)) = nRwd;
    signalMat(i, 1:Tsesh(i)) = allSignalFocus{i}';
end
%%
x = -1:0.02:1;
rpe = x;
nrwd = -sign(x);
figure2; hold on;
plot(x, rpe);
plot(x, nrwd);
%%
w = 0:0.05:1;
colors = cool(length(w));
figure2; hold on;
for i = 1:length(w)
    currY = w(i)*rpe + (1-w(i))*nrwd;
    plot(x, currY, 'Color', colors(i,:));
end
yMid = 0.5*rpe + 0.5*nrwd;
plot(x, yMid, 'Color', 'r');
%%
numBins = 4;
for i = 1:length(allOutcome)
    figure2;
    meanPe = NaN(1, numBins);
    meanSig = NaN(1, numBins);
    edges = linspace(min(allPe{i})-0.01, max(allPe{i})+0.01, numBins+1);
    for j = 1:numBins
        meanPe(j) = mean(allPe{i}(allPe{i}>=edges(j) & allPe{i}<edges(j+1)), 'omitnan');
        semPe(j) = sem(allPe{i}(allPe{i}>=edges(j) & allPe{i}<edges(j+1)));
        meanSig(j) = mean(allSignalFocus{i}(allPe{i}>=edges(j) & allPe{i}<edges(j+1)), 'omitnan');
    end
    errorbar(meanPe, meanSig, semPe)
    title(dayList{i});
end
%%
[root, sep] = currComputer();
xlFile = 'photometry'; 
sheet = 'all'; 
category = 'goodFP';
region = 'mPFC';
modelName = 'gcampNoPrior';
session_dat = struct('N',N,'T',T, 'Tsesh', Tsesh, 'rpeMat', rpeMat, 'nRwdMat', nRwdMat, 'signalMat', signalMat);
savePath = [root xlFile sep sheet 'sorted' sep 'stan' sep modelName sep category sep];
iter = 2000;
warmup = 0.5*iter;
numChains = 8;

if ~exist(savePath)
    mkdir(savePath);      
end

filePath = 'C:\Users\zhixi\Documents\gitRepositories\sueAnalysis\matlabCode\operantMatching\learningModels\stan\GCaMP\';
fullName = [modelName '.stan'];
fit = stan('file',[filePath fullName],'data',session_dat,'verbose',true,...
            'iter', iter, 'warmup', warmup, 'working_dir', savePath, 'chains', numChains, 'refresh', 200);
%%
[~, summary] = fit.print();
samples = fit.extract('permuted',true);
%% calculate all estimated params
paramNames = {'half', 'slope', 'maxF', 'intercept', 'w'};
allParams = zeros(length(paramNames), length(dayList));

for i = 1:length(dayList)
    for j = 1:length(paramNames)
        currParams = samples.(paramNames{j})(:,i);
        edges = linspace(min(currParams)-0.0001, max(currParams)+0.0001, 50);
        counts = histcounts(currParams, edges);
        [~, ind] = max(counts);
        allParams(j, i) = mean(currParams(currParams>=edges(ind) & currParams<edges(ind+1)));
    end
end
%% plot all the params
figure2;
colors = cool(length(paramNames));
for i = 1:length(paramNames)
    subplot(1,length(paramNames),i);
    histogram(allParams(i,:), 5, 'FaceColor', colors(i,:), 'EdgeColor','none');
    title(paramNames{i});
end
%%
figure; hold on;
edgesAll = linspace(0, 1, 51);
edgesAll(1) = edgesAll(1)-0.001;
edgesAll(end) = edgesAll(end)+0.001;
for i = 1:length(dayList)
    [counts, edges] = histcounts(samples.w(:,i), edgesAll);
    counts = counts/(sum(counts));
    mid = 0.5*(edges(1:end-1) + edges(2:end));
    patch([mid flip(mid)], [counts, zeros(size(mid))]+i*0.1, 'k', 'FaceAlpha', 0.5, 'EdgeColor', 'none');
    p5 = summary{{[paramNames{5} '[' num2str(i) ']']}, {'p5_'}};
    p95 = summary{{[paramNames{5} '[' num2str(i) ']']}, {'p95_'}};
    scatter([p95], [i*0.1], 15, 'r', 'filled');
    scatter([p5], [i*0.1], 15, 'b', 'filled');
end
plot([0.5 0.5], [0 0.1*(length(dayList)+1)], 'Color', [0.6 0.6 0.6], 'LineStyle', '--')
xlabel('w')
ylabel('posterior')
set(gca, 'TickDir', 'out')
%% plot actual data with model fitting no combine
numBins = 4;
for i = 1:length(allOutcome)
    figure; hold on;
    meanPe = NaN(1, numBins);
    meanSig = NaN(1, numBins);
    edges = linspace(min(allPe{i})-0.01, max(allPe{i})+0.01, numBins+1);
    for j = 1:numBins
        meanPe(j) = mean(allPe{i}(allPe{i}>=edges(j) & allPe{i}<edges(j+1)), 'omitnan');
        meanSig(j) = mean(allSignalFocus{i}(allPe{i}>=edges(j) & allPe{i}<edges(j+1)), 'omitnan');
    end
    plot(meanPe, meanSig)
    title(dayList{i});
    prediction = caActivation(allPe{i}, allParams(1,i), allParams(2,i), allParams(3,i)) + allParams(4,i);
    scatter(allPe{i}, prediction);
end
%% plot actual data with model fitting combine
numBins = 6;
for i = 1:length(allOutcome)
    figure; hold on;
    meanPe = NaN(1, numBins);
    meanSig = NaN(1, numBins);
    edges = quantile(allPe{i}, numBins);
    edges(1) = edges(1)-0.001;
    edges(end) = edges(end)+0.001;
    edges = linspace(min(allPe{i})-0.001, max(allPe{i})+0.001, numBins+1);
    p5 = summary{{[paramNames{5} '[' num2str(i) ']']}, {'p5_'}};
    p95 = summary{{[paramNames{5} '[' num2str(i) ']']}, {'p95_'}};
    for j = 1:numBins
        meanPe(j) = mean(allPe{i}(allPe{i}>=edges(j) & allPe{i}<edges(j+1)), 'omitnan');
        semPe(j) = sem(allSignalFocus{i}(allPe{i}>=edges(j) & allPe{i}<edges(j+1)));
        meanSig(j) = mean(allSignalFocus{i}(allPe{i}>=edges(j) & allPe{i}<edges(j+1)), 'omitnan');
    end
    nRwd = -allOutcome{i};
    nRwd(nRwd==0) = 1;
    errorbar(meanPe, meanSig, semPe,'LineWidth', 2, 'Color', [0.4 0.4 0.4]);
%     title([num2str(allParams(1, i)) ' ' num2str(allParams(2, i)) ' ' num2str(allParams(2, i)) ' ' num2str(allParams(1, i))]);
    title(['w = ' num2str(allParams(5,i)) ' ' num2str(p5) ' ' num2str(p95)])
    scatter([allParams(1, i)], [allParams(4, i)], 10, 'r', 'filled');
    plot([-1 1], [allParams(3, i)+allParams(4, i)]*[1 1], 'LineStyle', '--', 'Color', [0.7 0.7 0.7]);
    prediction = caActivation(allParams(5,i)*allPe{i}+(1-allParams(5,i))*nRwd, allParams(1,i), allParams(2,i), allParams(3,i)) + allParams(4,i);
%     scatter(allPe{i}, allSignalFocus{i}, 10, [0.7 0.7 0.7]);
    scatter(allPe{i}, prediction, 10, 'k');
end
%%
