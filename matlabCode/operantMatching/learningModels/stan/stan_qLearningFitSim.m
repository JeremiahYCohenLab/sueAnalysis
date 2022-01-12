function [paramEsts] = stan_qLearningFitSim(modelName, sessionNum, varargin)

p = inputParser;
% default parameters if none given
p.addParameter('revForFlag', 0);
p.addParameter('bernFlag', 1);
p.addParameter('iter', 10000);
p.addParameter('warmup', []);
p.addParameter('saveFlag', 1);
p.addParameter('maxTrial', 500);
p.addParameter('numChains', 6);
p.addParameter('simFlag', 0);
p.addParameter('control', struct('delta', 0.85))
p.parse(varargin{:});

paramNames = getParamNames_dF(modelName,1);
fullName = ['stan_qLearning_' modelName '.stan'];

if isempty(p.Results.warmup)
    warmup = max(floor(p.Results.iter/2),1);
else
    warmup = p.Results.warmup;
end

[root, sep] = currComputer();

if p.Results.bernFlag
    savePath = [root 'sim' sep 'stan' sep 'bernoulli' sep modelName sep];
else
    savePath = [root 'sim' sep 'stan' sep modelName sep];
end

if ~exist(savePath)
    mkdir(savePath);
end
% generate random paramters
params = zeros(sessionNum, length(paramNames));
allParams = struct;
allParams.a = betarnd(3,4, sessionNum, 1); % a
allParams.aN  = betarnd(3, 5, sessionNum, 1); % aN 
allParams.aP = betarnd(5, 2, sessionNum, 1); % aP
allParams.aF = betarnd(4, 2, sessionNum, 1); % aF
allParams.beta = normrnd(5, 2, sessionNum, 1); % beta
allParams.beta(allParams.beta<1) = 2; % beta
allParams.k = 2*betarnd(1, 3, sessionNum, 1); % kappa
allParams.bias = 2*(rand(sessionNum, 1)-0.5); % bias

% parameters
for i = 1:length(paramNames)
    params(:,i) = allParams.(paramNames{i});
end
% simulation
choice = zeros(sessionNum, p.Results.maxTrial);
outcome = zeros(sessionNum, p.Results.maxTrial);
for sim = 1:sessionNum
    %simulation
%     if contains(p.Results.modelName, 'aF')
%         [~, outcomeSim, choiceSim] = vkfSim_aF('params', params(sim,:),'randomSeed', sim,'maxTrials', p.Results.maxTrial, 'plotFlag', 0);
%     else
%         if contains(p.Results.modelName, 'kappa')
%            [~, outcomeSim, choiceSim] = vkfSim_kappa('params', params(sim,:),'randomSeed', sim,'maxTrials', p.Results.maxTrial, 'plotFlag', 0); 
%         else
%             [~, outcomeSim, choiceSim] = vkfSim('params', params(sim,:),'randomSeed', sim,'maxTrials', p.Results.maxTrial, 'plotFlag', 0);
%         end
%     end
    currParams = params(sim,:);
    expression = ['[~, outcomeSim, choiceSim] = qLearningModel_', modelName, '_simNoPlot(currParams, p.Results.maxTrial, sim);'];
    eval(expression)
    choiceSim(choiceSim<0) = 0;
    outcomeSim = abs(outcomeSim);
    choice(sim,:) = choiceSim;
    outcome(sim,:) = outcomeSim;
end

T = p.Results.maxTrial;
N = sessionNum;
Tsesh = p.Results.maxTrial*ones(sessionNum,1);

session_dat = struct('N',N,'T',T, 'Tsesh', Tsesh, 'choice', choice, 'outcome', outcome);

%run the stan model
if p.Results.bernFlag
    filePath = 'C:\Users\zhixi\Documents\gitRepositories\sueAnalysis\matlabCode\operantMatching\learningModels\stan\bernoulli\';
else
    filePath = 'C:\Users\zhixi\Documents\gitRepositories\sueAnalysis\matlabCode\operantMatching\learningModels\stan\';
end
fit = stan('file',[filePath fullName],'data',session_dat,'verbose',true,...
            'iter', p.Results.iter, 'warmup', warmup, 'working_dir', savePath, 'chains', p.Results.numChains, 'refresh', 200, 'control', p.Results.control);
%read command line output to stall matlab until stan is finished processing
doneFlag = 0;
diary([savePath 'diaryTmp.txt']); diary off;
fid = fopen([savePath 'diaryTmp.txt'],'rt');
tmp = textscan(fid,'%s','Delimiter','\n');
fclose(fid);
baseCount = find(~cellfun(@isempty,strfind(tmp{1}, '[100%]')) == 1);
while doneFlag == 0
    diary([savePath 'diaryTmp.txt']); pause(5); diary off;
    fid = fopen([savePath 'diaryTmp.txt'],'rt');
    tmp = textscan(fid,'%s','Delimiter','\n');
    fclose(fid);
    tmpCount = find(~cellfun(@isempty,strfind(tmp{1}, '[100%]')) == 1);
    if ~isempty(tmpCount)
        if length(tmpCount) == length(baseCount) + p.Results.numChains
            doneFlag = 1;
        end
    end
end
    delete([savePath 'diaryTmp.txt'])
fit.block();
pause(30);
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

allSamples = [];
edges = cell(1,length(paramNames)-1);
for i = 1:length(paramNames)-1
    tmp = eval(['samples.mu_' paramNames{i}]);
    allSamples = [allSamples tmp];
    edges{i} = linspace(min(tmp), max(tmp),40);
end
n = histcnd(allSamples,edges); %bin samples by multiple dimensions
[~, inds] = myMaxAll(n); %find the bin with max num in bin
for i = 1:length(paramNames)-1 %use median in bin as best estimate
    tmp = allSamples(:,i);
    edgeTmp = edges{i};
    if inds(i) < 50
        paramEsts(i) = median(tmp(tmp >= edgeTmp(inds(i)) & tmp < edgeTmp(inds(i)+1)));
    else
        paramEsts(i) = edgeTmp(inds(i));
    end
end


%plot the distributions of the mouse-level parameters
pFig = figure2('position', [0 0 800 400]); 
numParams = length(paramNames)-1;
blue = [0 1 1];
purp = [0.7 0 1];
colors = [linspace(blue(1),purp(1),numParams)', linspace(blue(2),purp(2),numParams)', linspace(blue(3),purp(3),numParams)'];
for i = 1:numParams
    subplot(1,numParams,i); hold on;
    histogram(eval(['samples.mu_' paramNames{i}]) , 100,...
        'Normalization', 'Probability', 'FaceColor', colors(i,:), 'EdgeColor', 'none')
%         line([paramEsts(i) paramEsts(i)], [0 0.05], 'color', [0 0 0]);
    set(gca,'tickdir', 'out') 
    title(paramNames{i})
end
titleTxt = strrep(['sim' ' - ' modelName], '_', ' ');
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
    sampFile = ['sim_', modelName];
    saveFile = [sampFile '.mat'];
    eval([sampFile,  ' = samples;']);
    saveFigurePDF(pFig,[savePath modelName  '_posteriors'])
    saveFigurePDF(dFig,[savePath modelName  '_divergence'])
    save([savePath saveFile], sampFile, 'paramEsts', 'params', 'outcome', 'choice')
    
end
    
