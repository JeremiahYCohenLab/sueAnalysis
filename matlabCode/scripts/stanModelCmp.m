% stan_qLearningFit(xlFile, sheet, category, varargin)
xlFile = 'inhibitionAll.xlsx';
sheet = 'ZS066';
col1 = 'inhibition';
col2 = 'control2';
modelName = '5params';
cmpParam = 4; 
iter = 10000;
warmup = iter * 0.5;
numChains = 6;
control = struct('delta', 0.8);
%%
paramNames = getParamNames_dF(modelName,0);
[root, sep] = currComputer();
savePath = [root sheet sep sheet 'sorted' sep 'stan' sep 'bernoulli' sep 'compare' sep modelName sep paramNames{cmpParam} sep col1 'vs' col2 sep];
if ~exist(savePath)
    mkdir(savePath);
end
dayList1 = getDayList(xlFile, sheet, col1);
dayList2 = getDayList(xlFile, sheet, col2);
dayList = [dayList1; dayList2];
N1 = length(dayList);
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
    N = length(dayList);
    N1 = length(dayList1);
    choice = zeros(N, T);
    outcome = zeros(N, T);
    ITI = zeros(N,T);

    for i = 1:N
        choice(i, 1:Tsesh(i)) = choiceTmp{i};
        outcome(i, 1:Tsesh(i)) = outcomeTmp{i};
        ITI(i, 1:Tsesh(i)) = ITItemp{i};
    end
%%
allFits = struct;
summaries = struct;
%%
for i = 1:4
cmpParam = i;
savePath = [root sheet sep sheet 'sorted' sep 'stan' sep 'bernoulli' sep 'compare' sep modelName sep paramNames{cmpParam} sep col1 'vs' col2 sep];
if ~exist(savePath)
    mkdir(savePath);  
end
%create data structure to feed into stan model
session_dat = struct('N',N,'T',T, 'Tsesh', Tsesh, 'choice', choice, 'outcome', outcome, 'ITI', ITI, 'N1', N1, 'cmpParam', cmpParam);
filePath = 'C:\Users\zhixi\Documents\gitRepositories\sueAnalysis\matlabCode\operantMatching\learningModels\stan\bernoulli\';
fullName = ['stan_qLearning_' modelName 'Cmp.stan'];
% run the stan model

allFits.(paramNames{i}) = stan('file',[filePath fullName],'data',session_dat,'verbose',true,...
            'iter', iter, 'warmup', warmup, 'working_dir', savePath, 'chains', numChains, 'refresh', 200, 'control', control);
% [~, summaries.(paramNames{i})] = allFits.(paramNames{i}).print();
end
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
    
