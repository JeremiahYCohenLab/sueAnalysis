function [paramEsts] = stan_qLearningFit_laserRPE(xlFile, sheet, category, varargin)

p = inputParser;
% default parameters if none given
p.addParameter('revForFlag', 0);
p.addParameter('bernFlag', 1);
p.addParameter('nonfixedParams', 0);
p.addParameter('fixedParams', []);
p.addParameter('modelName', '5paramsLaserNegRPE');
% p.addParameter('modelName', 'vkf_fixV_kappa');
p.addParameter('iter', 10000);
p.addParameter('warmup', []);
p.addParameter('saveFlag', 1);
p.addParameter('maxTrial', 1000);
p.addParameter('numChains', 8);
p.addParameter('simFlag', 0);
p.addParameter('shuffle', 0);
p.addParameter('control', struct('delta', 0.85))
p.parse(varargin{:});

paramNames = getParamNames_dF(p.Results.modelName,0);
if ~p.Results.nonfixedParams
    paramInds = 1:length(paramNames);
    fullName = ['stan_qLearning_' p.Results.modelName '.stan'];
else
    paramInds = find(~contains(paramNames, p.Results.nonfixedParams));
    fullName = ['stan_qLearning_' p.Results.modelName '_' p.Results.nonfixedParams '.stan'];
end

if contains(p.Results.modelName, 'vkf')
    fullName = ['stan_' p.Results.modelName '.stan'];
end

if isempty(p.Results.warmup)
    warmup = max(floor(p.Results.iter/2),1);
else
    warmup = p.Results.warmup;
end

[root, sep] = currComputer();



if p.Results.bernFlag
    if p.Results.shuffle
        savePath = [root sheet sep sheet 'sorted' sep 'stan' sep 'bernoulli' sep p.Results.modelName sep category sep 'shuffle' sep];
    else
        savePath = [root sheet sep sheet 'sorted' sep 'stan' sep 'bernoulli' sep p.Results.modelName sep category sep];
    end

else
    savePath = [root sheet sep sheet 'sorted' sep 'stan' sep p.Results.modelName sep category sep];
end
if ~exist(savePath)
    mkdir(savePath);
end
dayList = getDayList(xlFile, sheet, category);

for i = 1:length(dayList)
    sessionName = dayList{i};
    behavStruct = behAnalysisNoPlot_opMD(sessionName, 'simpleFlag', 1);

    choiceTmp{i} = behavStruct.allChoices;
    if p.Results.bernFlag
        choiceTmp{i}(choiceTmp{i} == -1) = 0;
    else
        choiceTmp{i} = choiceTmp{i} + 1;
        choiceTmp{i}(choiceTmp{i} == 0) = 1;
    end
    outcomeTmp{i} = abs(behavStruct.allRewards); 
    ITItemp{i} = behavStruct.timeBtwn;
    Tsesh(i,1) = length(outcomeTmp{i});
    if p.Results.shuffle
        tmp = rand(size(behavStruct.laser));
        tmp(tmp <= 1-mean(behavStruct.laser)) = 0;
        tmp(tmp > 1-mean(behavStruct.laser)) = 1;
        laserTemp{i} = tmp;
    else
        laserTemp{i} = [behavStruct.laser];
    end
end

T = max(Tsesh);
N = length(dayList);
choice = zeros(N, T);
outcome = zeros(N, T);
laser = zeros(N, T);
ITI = zeros(N,T);

for i = 1:N
    choice(i, 1:Tsesh(i)) = choiceTmp{i};
    outcome(i, 1:Tsesh(i)) = outcomeTmp{i};
    laser(i, 1:Tsesh(i)) = laserTemp{i};
    ITI(i, 1:Tsesh(i)) = ITItemp{i}/1000; % convert from ms to s
end
choice = choice(:,1:min([T p.Results.maxTrial]));
outcome = outcome(:,1:min([T p.Results.maxTrial]));
ITI = ITI(:,1:min([T p.Results.maxTrial]));
Tsesh(Tsesh>p.Results.maxTrial) = p.Results.maxTrial;
T = min([T p.Results.maxTrial]);

session_dat = struct('N',N,'T',T, 'Tsesh', Tsesh, 'choice', choice, 'outcome', outcome, 'ITI', ITI, 'laser', laser);

%%
%create data structure to feed into stan model

if p.Results.nonfixedParams
    for j = 1:length(paramInds)
        session_dat.(paramNames{paramInds(j)}) = p.Results.fixedParams(:,paramInds(j));
    end
end

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
%% plot results
posterior = figure;
colors = [0, 0.8, 0.8;
          1, 0.3, 1];
tmpSamp = samples;
numParams = length(paramNames);
for currP = 1:length(paramNames)
    subplot(1,numParams, currP); hold on;
    mu = tmpSamp.(['mu_' paramNames{currP}]);
    edges = linspace(min(mu), max(mu), 50);
    histogram(mu, edges, 'FaceColor', colors(1,:), 'EdgeColor', 'none', 'Normalization', 'probability')
    title(paramNames{currP});
    if strcmp(paramNames{currP}, 'diff')
        p5__ = summary{{'mu_diff'}, {'p5_'}};
        p95__ = summary{{'mu_diff'}, {'p95_'}};
        line([p5__ p5__], [0 0.1], 'LineStyle', '--', 'Color', [1 0.5 0.5],'LineWidth',2)
        line([p95__ p95__], [0 0.1], 'LineStyle', '--', 'Color', [0.5 0.5 1],'LineWidth',2)
    end
end


screen = get(0,'Screensize');
screen(4) = screen(4) - 500;
set(posterior, 'Position', screen)

samps = tmpSamp;
dFig = figure;
for currPy = 1:numParams + size(tmpSamp.sigma,2)
    if currPy <= numParams
        tmpY = eval(['samps.mu_' paramNames{currPy}]);
    else
        tmpY = samps.sigma(:,currPy-numParams);
    end
    tmpY_d = tmpY(logical(samps.divergent__));
    tmpY = tmpY(~logical(samps.divergent__));
    for currPx = 1:numParams + size(samps.sigma,2)
        subplot(numParams+size(samps.sigma,2), numParams+size(samps.sigma,2),(currPy-1)*(numParams + size(samps.sigma,2)) + currPx); hold on;

        if currPy == currPx
            h = histogram(tmpY, 30, 'FaceColor', 'c', 'normalization', 'probability');
            histogram(tmpY_d, h.BinEdges, 'FaceColor', 'm', 'normalization', 'probability')
        else       
            if currPx <= numParams
                tmpX = eval(['samps.mu_' paramNames{currPx}]);
            else
                tmpX = samps.sigma(:,currPx-numParams);
            end
            tmpX_d = tmpX(logical(samps.divergent__));
            tmpX = tmpX(~logical(samps.divergent__));
            scatter(tmpX, tmpY, 1, 'c', 'filled')
            scatter(tmpX_d, tmpY_d, 1, 'm', 'filled')
        end

        if currPx == 1
            if currPy <= numParams
               ylabel(paramNames{currPy})
            else
               ylabel(['sigma' '(' num2str(currPy-numParams) ')'])
            end
        end
        if currPy == numParams + size(samps.sigma,2)
            if currPx <= numParams
               xlabel(paramNames{currPx})
            else
               xlabel(['sigma' '(' num2str(currPx-numParams) ')'])
            end
        end

    end
end
screen = get(0,'Screensize');
screen(4) = screen(4) - 100;
set(dFig, 'Position', screen)



if p.Results.saveFlag
    %save the full samples

    sampFile = [sheet category '_', p.Results.modelName];
    saveFile = [sampFile '.mat'];
    eval([sampFile,  ' = samples;']);
    saveFigurePDF(posterior,[savePath sheet category '_' p.Results.modelName '_posteriors'])
    saveFigurePDF(dFig,[savePath sheet category '_' p.Results.modelName '_diverg'])
    save([savePath saveFile], sampFile, 'dayList', 'summary');
        
end
    %save([savePath saveFile], sampFile, 'paramEsts', 'dayList', 'tbl');
    
    
