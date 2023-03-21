%% settings
xlFile = 'inhibitionAll.xlsx';
sheet = 'ZS071';
% col1 = 'control2PairCtrl';
% col2 = 'control1PairCtrl';
col1 = 'inhibitionNrwd';
col2 = 'controlNrwd';
modelName = '5params';
iter = 10000;
warmup = iter * 0.5;
numChains = 6;
control = struct('delta', 0.8);

paramNames = getParamNames_dF(modelName,0);
[root, sep] = currComputer();
dayList1 = getDayList(xlFile, sheet, col1);
dayList2 = getDayList(xlFile, sheet, col2);
pathCmp = [root sheet sep sheet 'sorted' sep 'stan' sep 'bernoulli' sep 'comparePaired' sep modelName sep];
sampFileCmp = [sheet col1 'vs' col2 '_' modelName];
%% load beh data
if length(dayList1)~=length(dayList2)
    print('sessions not paired \n')
    return
end
dayList = [dayList1; dayList2];
Tsesh = [];
    for i = 1:length(dayList)
        sessionName = dayList{i};
        filename = [sessionName '.asc'];
        %fprintf([sessionName '\n']);
        behSessionData = loadBehavioralData(filename, 0);
        behavStruct = parseBehavioralData(behSessionData, 1000);

        choiceTmp{i} = behavStruct.allChoices;
        choiceTmp{i}(choiceTmp{i} == -1) = 0;

        outcomeTmp{i} = abs(behavStruct.allRewards); 
        ITItemp{i} = behavStruct.timeBtwn;
        Tsesh(i,1) = length(outcomeTmp{i});
    end

    T = max(Tsesh);
    N = length(dayList1);
    choice = zeros(N, T);
    outcome = zeros(N, T);
    ITI = zeros(N,T);

    for i = 1:2*N
        choice(i, 1:Tsesh(i)) = choiceTmp{i};
        outcome(i, 1:Tsesh(i)) = outcomeTmp{i};
        ITI(i, 1:Tsesh(i)) = ITItemp{i};
    end
%%
allFits = struct;
summaries = struct;
%% fitting
for i = 1:length(paramNames)
    cmpParam = i;
    savePath = [root sheet sep sheet 'sorted' sep 'stan' sep 'bernoulli' sep 'comparePaired' sep modelName sep col1 'vs' col2 sep paramNames{cmpParam} sep];
    if ~exist(savePath)
        mkdir(savePath);  
    end
    %create data structure to feed into stan model
    session_dat = struct('N',N,'T',T, 'Tsesh', Tsesh, 'choice', choice, 'outcome', outcome, 'ITI', ITI, 'cmpParam', cmpParam);
    filePath = 'C:\Users\zhixi\Documents\gitRepositories\sueAnalysis\matlabCode\operantMatching\learningModels\stan\bernoulli\';
    fullName = ['stan_qLearning_' modelName 'CmpPair.stan'];
    % run the stan model

    allFits.(paramNames{i}) = stan('file',[filePath fullName],'data',session_dat,'verbose',true,...
                'iter', iter, 'warmup', warmup, 'working_dir', savePath, 'chains', numChains, 'refresh', 200, 'control', control);
    % [~, summaries.(paramNames{i})] = allFits.(paramNames{i}).print();
end
   %% summary
for i = 1:length(paramNames)
[~, summaries.(paramNames{i})] = allFits.(paramNames{i}).print();
end 
save([pathCmp col1 'vs' col2 sep 'summaries.mat'], 'summaries');
%% read samples from file and save
saveFlag = 1;
biasFlag = 1;
for currParam = 1:length(paramNames)
      samples = struct;
   fieldNames = [paramNames strcat('mu_', paramNames) strcat('mu_', paramNames, '_diff') ];
   fieldNames = [fieldNames 'log_lik'];
   fieldNames = [fieldNames 'log_likMean'];
   fieldNames = [fieldNames 'sigma'];
   fieldNames = [fieldNames 'mu_diff'];
   fieldNames = [fieldNames 'sigmaDiff'];
   fieldNames = [fieldNames 'divergent__'];
   if biasFlag
       fieldNames = [fieldNames 'bias'];   
   end
   for f = 1:length(fieldNames)
       samples.(fieldNames{f}) = [];
   end
    % read cvs files
   modelFiles = dir([pathCmp col1 'vs' col2 sep paramNames{currParam} sep ]);
   outputExpression = ['^' 'output'];
   outputFiles = {modelFiles(~cellfun(@isempty, cellfun(@(x) regexp(x, outputExpression), {modelFiles.name}, 'UniformOutput', false))).name};
   for i = 1:length(outputFiles)
       currOutput = outputFiles{i};
       [nums, texts, all] = xlsread([pathCmp col1 'vs' col2 sep paramNames{currParam} sep currOutput], currOutput(1:end-4));
       titles = texts(39,:);
       titles = cellfun(@(x) strtok(x, '.'), titles, 'UniformOutput', false);
        for f = 1:length(fieldNames)
            tmpInd = cellfun(@(x) strcmp(x, fieldNames(f)), titles, 'UniformOutput', false);
            tmpInd = cell2mat(tmpInd);
            tmp = nums(2:end, tmpInd);
            samples.(fieldNames{f}) = [[samples.(fieldNames{f})]; tmp];
        end
   end
   
   if saveFlag
       save([pathCmp sep col1 'vs' col2 sep paramNames{currParam}  sep 'samples.mat'], 'samples', 'dayList1', 'dayList2')
   end
end
%% load all samples
allSamps = struct;
load([pathCmp col1 'vs' col2 sep 'summaries.mat']);
for currCmp = 1:length(paramNames)
    load([pathCmp sep col1 'vs' col2 sep paramNames{currCmp} sep 'samples.mat'])
    allSamps.(paramNames{currCmp}) = samples;
end
%% plot all posteriors
colors = cool(length(paramNames));
allP = figure;
numParams = length(paramNames);
for currCmp = 1:numParams
    tmpSamp = allSamps.(paramNames{currCmp});
    for currP = 1:numParams
        subplot(numParams,numParams, sub2ind([numParams, numParams], currP, currCmp)); hold on;
        mu = tmpSamp.(['mu_' paramNames{currP}]);
        mu_diff = tmpSamp.(['mu_' paramNames{currP} '_diff']);
        if currP == currCmp
            edges = linspace(min([mu; mu_diff]), max([mu; mu_diff]), 50);
            histogram(mu, edges, 'FaceColor', colors(currP,:), 'EdgeColor', 'none', 'Normalization', 'probability');
            histogram(mu_diff, edges, 'FaceColor', [0 1 0], 'FaceAlpha', 0.4, 'EdgeColor', 'none', 'Normalization', 'probability');
        else
            edges = linspace(min(mu), max(mu), 50);
            histogram(mu, edges, 'FaceColor', colors(currP,:), 'EdgeColor', 'none', 'Normalization', 'probability')
        end
        if currCmp == 1
            title(paramNames{currP});
        end
    end
end
suptitle([sheet ' ' col1 'vs' col2 ' ' modelName]);
screenSize = get(0,'Screensize');
screenSize(4) = screenSize(4) - 100;
set(allP, 'renderer', 'painters', 'position', screenSize)
saveFigurePDF(allP,[pathCmp col1 'vs' col2 sep modelName '_' col1 'vs' col2 'allParamPosterior.pdf'])

%% plot compared posteriors
cmpP = figure;
for currCmp = 1:numParams
    tmpSamp = allSamps.(paramNames{currCmp});
    muDiff = [tmpSamp.(['mu_' paramNames{currCmp} '_diff'])] - [tmpSamp.(['mu_' paramNames{currCmp}])];
    diff = tmpSamp.mu_diff;
    subplot(3, numParams, sub2ind([numParams,3], currCmp, 1)); hold on;
    edges = linspace(min(muDiff), max(muDiff), 50);
    histogram(muDiff, edges, 'FaceColor', colors(currCmp,:), 'EdgeColor', 'none', 'Normalization', 'probability');
    line([0 0], [0 0.1], 'LineStyle', '--', 'Color', [0.5 0.5 0.5],'LineWidth',2)
    title([paramNames{currCmp} ' mu2-mu1'])
    subplot(3, numParams, sub2ind([numParams,3], currCmp, 2)); hold on;
    edges = linspace(min(diff), max(diff), 50);
    histogram(diff, edges, 'FaceColor', colors(currCmp,:), 'EdgeColor', 'none', 'Normalization', 'probability');
    p5__ = summaries.(paramNames{currCmp}){{'mu_diff'}, {'p5_'}};
    p95__ = summaries.(paramNames{currCmp}){{'mu_diff'}, {'p95_'}};
    line([p5__ p5__], [0 0.1], 'LineStyle', '--', 'Color', [1 0.5 0.5],'LineWidth',2)
    line([p95__ p95__], [0 0.1], 'LineStyle', '--', 'Color', [0.5 0.5 1],'LineWidth',2)
    title([paramNames{currCmp} ' muDiff'])
    subplot(3,numParams,sub2ind([numParams,3], currCmp, 3)); hold on;
    loglikMean = mean(tmpSamp.log_likMean,2);
    edges = linspace(min(loglikMean), max(loglikMean), 50);
    histogram(loglikMean, edges, 'FaceColor', colors(currCmp,:), 'EdgeColor', 'none', 'Normalization', 'probability');
    title('logLL/trial')
end
suptitle([sheet ' ' col1 'vs' col2 ' ' modelName]);
set(cmpP, 'renderer', 'painters', 'position', screenSize)
saveFigurePDF(cmpP,[pathCmp col1 'vs' col2 sep modelName '_' col1 'vs' col2 'cmpParamPosterior.pdf'])


%% compute priors
samps = 1000000;
mu_common_pr = normrnd(0,sqrt(0.5),samps,1);
diff_pr = normrnd(0,sqrt(0.5), samps,1);

mu_1 = normcdf(mu_common_pr - diff_pr);
mu_2 = normcdf(mu_common_pr + diff_pr);
diff = mu_2 - mu_1;

figure; 
subplot(1,2,1);
histogram(mu_common_pr, 'FaceColor', [0.6 0.6 0.6], 'EdgeColor', 'none', 'Normalization', 'probability');
title('mu-common');
subplot(1,2,2);
histogram(diff_pr, 'FaceColor', [1 0 0], 'EdgeColor', 'none', 'Normalization', 'probability');
title('mu-diff');

figure;
subplot(1,2,1); hold on;
histogram(mu_1, 0:0.05:1, 'Normalization', 'probability', 'EdgeColor', 'none');
histogram(mu_2, 0:0.05:1, 'Normalization', 'probability', 'EdgeColor', 'none');
legend({'mu-Ctrl','mu-Inhi'})
title('mu');

subplot(1,2,2);
histogram(diff)
title('mu-Inhi-mu-Ctrl')

%% session prior 
var = 1;
sigma = abs(cauchyrnd(0, var, samps, 1));
mu = normrnd(0, 1, samps, 1);
sessPr = normrnd(0, 1, samps, 1);
sessParam = normcdf(mu + sigma.*sessPr);
figure;
histogram(sessParam, 'Normalization', 'probability');
%%
sigma1 = abs(cauchyrnd(0, var, samps, 1));
sigma2 = abs(cauchyrnd(0, 0.1, samps, 1));
sessionPr1 = normrnd(0, 1, samps, 1);
sessionPr2 = normrnd(0, 1, samps, 1);
sessParam_1 = normcdf(mu_common_pr + sigma1.*sessionPr1 - diff_pr - sigma2.*sessionPr2);
sessParam_2 = normcdf(mu_common_pr + sigma1.*sessionPr1 + diff_pr + sigma2.*sessionPr2);

figure;
subplot(1,2,1);
histogram(sessParam_1, 'Normalization', 'probability');
subplot(1,2,2);
histogram(sessParam_2, 'Normalization', 'probability');
%%
pFig = figure2;
numParams = length(paramNames);
blue = [0 1 1];
purp = [0.7 0 1];
colors = cool(numParams);
for i = 1:numParams
    subplot(1,numParams,i); hold on;
    histogram(eval(['samples.mu_' paramNames{i}]) , 100,...
        'Normalization', 'Probability', 'FaceColor', colors(i,:), 'EdgeColor', 'none')
    set(gca,'tickdir', 'out') 
    title(paramNames{i})
end

%% save
saveFile = [sampFile '.mat'];
eval([sampFile,  ' = samples;']);
save([path saveFile], sampFile);
saveFigurePDF(pFig,[path model  '_posteriors.pdf'])
%%


%extract samples from the stan fit object
samples = [];
while isempty(samples)
    samples = fit.extract('permuted',true);
end
if isfield(samples, 'y_pred')
    samples = rmfield(samples, 'y_pred');
end
[~, summary] = fit.print();

%generate best estimates of parameters
 paramEsts = [];
if p.Results.nonfixedParams
    tmp = eval(['samples.mu_' p.Results.nonfixedParams]);
    [n,e] = histcounts(tmp, 50);
    [~, maxInd] = max(n);
    paramEsts = median(tmp(tmp > e(maxInd) & tmp < e(maxInd+1)));
else
    allSamples = [];    
    edges = cell(1,length(paramInds));
    for i = 1:length(paramInds)
        tmp = eval(['samples.mu_' paramNames{i}]);
        allSamples = [allSamples tmp];
        edges{i} = linspace(min(tmp), max(tmp),40);
    end
    n = histcnd(allSamples,edges); %bin samples by multiple dimensions
    [~, inds] = myMaxAll(n); %find the bin with max num in bin
    for i = 1:length(paramInds) %use median in bin as best estimate
        tmp = allSamples(:,i);
        edgeTmp = edges{i};
        if inds(i) < 50
            paramEsts(i) = median(tmp(tmp >= edgeTmp(inds(i)) & tmp < edgeTmp(inds(i)+1)));
        else
            paramEsts(i) = edgeTmp(inds(i));
        end
    end
end


%plot the distributions of the mouse-level parameters
pFig = figure2('position', [0 0 800 400]); 
if p.Results.nonfixedParams
    histogram(eval(['samples.mu_' p.Results.nonfixedParams]), 100,...
            'Normalization', 'Probability', 'FaceColor', 'k')
        set(gca,'tickdir', 'out')
        xlabel(p.Results.nonfixedParams)
else
    numParams = length(paramInds);
    blue = [0 1 1];
    purp = [0.7 0 1];
    colors = [linspace(blue(1),purp(1),numParams)', linspace(blue(2),purp(2),numParams)', linspace(blue(3),purp(3),numParams)'];
    for i = 1:numParams
        subplot(1,numParams,i); hold on;
        histogram(eval(['samples.mu_' paramNames{paramInds(i)}]) , 100,...
            'Normalization', 'Probability', 'FaceColor', colors(i,:), 'EdgeColor', 'none')
%         line([paramEsts(i) paramEsts(i)], [0 0.05], 'color', [0 0 0]);
        set(gca,'tickdir', 'out') 
        title(paramNames{paramInds(i)})
    end
end
titleTxt = strrep([sheet ' - ' p.Results.modelName], '_', ' ');
suptitle(titleTxt);
set(gcf,'Renderer', 'Painters')

dFig = figure;
for currPy = 1:numParams
    tmpY = eval(['samples.mu_' paramNames{currPy}]);
    tmpY_d = tmpY(logical(samples.divergent__));
    tmpY = tmpY(~logical(samples.divergent__));
    for currPx = 1:numParams
        subplot(numParams,numParams,[(currPy-1)*numParams + currPx]); hold on;
        
        if currPy == currPx
            h = histogram(tmpY, 30, 'FaceColor', 'c', 'normalization', 'probability');
            histogram(tmpY_d, h.BinEdges, 'FaceColor', 'm', 'normalization', 'probability')
        else       
            tmpX = eval(['samples.mu_' paramNames{currPx}]);
            tmpX_d = tmpX(logical(samples.divergent__));
            tmpX = tmpX(~logical(samples.divergent__));
            scatter(tmpX, tmpY, [], 'c')
            scatter(tmpX_d, tmpY_d, [], 'm')
        end
        
        if currPx == 1
            ylabel(paramNames{currPy})
        end
        if currPy == numParams
            xlabel(paramNames{currPx})
        end
        
    end
end
titleTxt = [titleTxt ' (divergence rate = ' num2str(sum(samples.divergent__)/length(samples.divergent__)) ')'];
suptitle(titleTxt);
set(gcf, 'renderer', 'painters', 'position', [-1919 41 1920 963])


if p.Results.saveFlag
    %save the full samples
    if p.Results.simFlag
        sampFile = ['sim_', p.Results.modelName];
        saveFile = [sampFile '.mat'];
        eval([sampFile,  ' = samples;']);
        saveFigurePDF(pFig,[savePath p.Results.modelName  '_posteriors'])
        saveFigurePDF(dFig,[savePath p.Results.modelName  '_divergence'])
        save([savePath saveFile], sampFile, 'paramEsts', 'params', 'outcome', 'choice');
    else
        if p.Results.nonfixedParams
            sampFile = [sheet category, '_', p.Results.modelName, '_', p.Results.nonfixedParams];
            saveFile = [sampFile '.mat'];
            eval([sampFile,  ' = samples;']);
        else
            sampFile = [sheet category '_', p.Results.modelName];
            saveFile = [sampFile '.mat'];
            eval([sampFile,  ' = samples;']);
        end
        saveFigurePDF(pFig,[savePath sheet category '_' p.Results.modelName  '_posteriors'])
        saveFigurePDF(dFig,[savePath sheet category '_' p.Results.modelName  '_divergence'])
        save([savePath saveFile], sampFile, 'paramEsts', 'dayList', 'summary');
        
    end
    %save([savePath saveFile], sampFile, 'paramEsts', 'dayList', 'tbl');
    
end
    
