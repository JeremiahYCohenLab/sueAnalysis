function [paramEsts] = stan_qLearningFit_tranWeight(xlFile, animal, category, varargin)

p = inputParser;
% default parameters if none given
p.addParameter('revForFlag', 0)
p.addParameter('bernFlag', 1)
p.addParameter('fixedParams', 0)
p.addParameter('paramType', ['fiveParamOfourStart_preS'])
p.addParameter('paramInds', [])
p.addParameter('paramNames', {'aN', 'aP', 'aF', 'beta'})
p.addParameter('modelName', ['fourParam'])
p.addParameter('iter', 2000)
p.addParameter('warmup', [])
p.addParameter('saveFlag', 1)
p.addParameter('rwdProbs', [90 50 10])
p.addParameter('tranWin', 15)
p.addParameter('simFlag', 0)
p.addParameter('rSeed', 46578492)
p.addParameter('params', []);
p.parse(varargin{:});

pHigh = p.Results.rwdProbs(1);
pLow = p.Results.rwdProbs(2);
probDiff = p.Results.rwdProbs(1) - p.Results.rwdProbs(3);
tranWin = p.Results.tranWin;

if isempty(p.Results.paramInds)
    paramInds = [1:length(p.Results.paramNames)];
else
    paramInds = p.Results.paramInds;
end

if isempty(p.Results.warmup)
    warmup = max(floor(p.Results.iter/2),1);
else
    warmup = p.Results.warmup;
end

[root, sep] = currComputer();
if p.Results.bernFlag
    savePath = [root animal sep animal 'sorted' sep 'stan' sep 'bernoulli' sep p.Results.modelName sep];
else
    savePath = [root animal sep animal 'sorted' sep 'stan' sep p.Results.modelName sep];
end
if ~exist(savePath)
    mkdir(savePath);
end

%if fitting model to simulated data for model recovery, extract relevant params
if p.Results.simFlag
    if isempty(p.Results.params)
        if p.Results.bernFlag
            modelPath = [root animal sep animal 'sorted' sep 'stan' sep 'bernoulli' sep p.Results.modelName sep animal...
                    category '_' p.Results.modelName '.mat'];
        else
            modelPath = [root animal sep animal 'sorted' sep 'stan' sep p.Results.modelName sep animal...
                        category '_' p.Results.modelName '.mat'];
        end
        t = generateStanModelTerms_opMD(p.Results.modelName, modelPath, [], 0);
        params = t.params;
    else
        params = p.Results.params;
    end
end

%load session list
[~, dayList, ~] = xlsread(xlFile, animal);
[~,col] = find(~cellfun(@isempty,strfind(dayList, category)) == 1);
dayList = dayList(2:end,col);
endInd = find(cellfun(@isempty,dayList),1);
if ~isempty(endInd)
    dayList = dayList(1:endInd-1,:);
end

for i = 1:length(dayList)
    sessionName = dayList{i};
    filename = [sessionName '.asc'];
    [behSessionData, unCorrectedBlockSwitch, out] = loadBehavioralData(filename, p.Results.revForFlag);
    behavStruct = parseBehavioralData(behSessionData, unCorrectedBlockSwitch);

    if p.Results.simFlag
        [~, allRewards, allChoices, blockProbs, blockSwitch] = runSim_dF(p.Results.modelName, params, 1000,...
            p.Results.rSeed+i, p.Results.rwdProbs);
        rwdProb_L = nan(1, length(allChoices));
        rwdProb_R = nan(1, length(allChoices));
        for j = 2:length(blockSwitch)
            rwdProb_L(blockSwitch(j-1):blockSwitch(j)-1) = blockProbs(j-1, 1);
            rwdProb_R(blockSwitch(j-1):blockSwitch(j)-1) = blockProbs(j-1, 2);
        end
        rwdProb_L(blockSwitch(end):length(allChoices)) = blockProbs(end, 1);
        rwdProb_R(blockSwitch(end):length(allChoices)) = blockProbs(end, 2);
        transTmp{i} = zeros(1, 1000);
    else
        allChoices = behavStruct.allChoices;
        blockSwitch = behavStruct.blockSwitch;
        allRewards = behavStruct.allRewards;
        rwdProb_R = [behSessionData(behavStruct.responseInds).rewardProbR]; 
        rwdProb_L = [behSessionData(behavStruct.responseInds).rewardProbL];
        transTmp{i} = zeros(1, length(behavStruct.allChoices));
    end
        

    choiceTmp{i} = allChoices;
    if p.Results.bernFlag
        choiceTmp{i}(choiceTmp{i} == -1) = 0;
    else
        choiceTmp{i} = choiceTmp{i} + 1;
        choiceTmp{i}(choiceTmp{i} == 0) = 1;
    end
    outcomeTmp{i} = abs(allRewards); 
    Tsesh(i,1) = length(outcomeTmp{i});

    for j = 2:(length(blockSwitch) - 1)
        tmpInd = blockSwitch(j);
        if tmpInd+tranWin <= length(allChoices) & tmpInd-tranWin > 0
            if rwdProb_R(tmpInd-1) == pHigh & rwdProb_R(tmpInd) == 10 & any(diff(rwdProb_L(tmpInd-tranWin:tmpInd)) == probDiff)
                transTmp{i}(tmpInd:tmpInd+tranWin-1) = 1;
            elseif rwdProb_R(tmpInd-1) == pLow & rwdProb_R(tmpInd) == 10 & any(diff(rwdProb_L(tmpInd-tranWin:tmpInd)) == probDiff)
                transTmp{i}(tmpInd:tmpInd+tranWin-1) = 1;
            elseif rwdProb_L(tmpInd-1) == pHigh & rwdProb_L(tmpInd) == 10 & any(diff(rwdProb_R(tmpInd-tranWin:tmpInd)) == probDiff)
                transTmp{i}(tmpInd:tmpInd+tranWin-1) = 1;
            elseif rwdProb_L(tmpInd-1) == pLow & rwdProb_L(tmpInd) == 10 & any(diff(rwdProb_R(tmpInd-tranWin:tmpInd)) == probDiff)
                transTmp{i}(tmpInd:tmpInd+tranWin-1) = 1;
            end   
        end
    end

    
end

T = max(Tsesh);
N = length(dayList);
choice = zeros(N, T);
outcome = zeros(N, T);
trans = zeros(N, T);

for i = 1:N
    choice(i, 1:Tsesh(i)) = choiceTmp{i};
    outcome(i, 1:Tsesh(i)) = outcomeTmp{i};
    trans(i, 1:Tsesh(i)) = transTmp{i};
end

%create data structure to feed into stan model
if p.Results.fixedParams
    [allParams, modelNames, ~] = xlsread('stanParams.xlsx', animal);
    [~,col] = find(~cellfun(@isempty,strfind(modelNames, p.Results.paramType)) == 1);
    params = allParams(paramInds,col);
    session_dat = struct('N',N,'T',T, 'Tsesh', Tsesh, 'choice', choice, 'outcome', outcome, 'params', params, 'trans', trans);
else
    session_dat = struct('N',N,'T',T, 'Tsesh', Tsesh, 'choice', choice, 'outcome', outcome, 'trans', trans);
end

%run the stan model
if p.Results.bernFlag
    filePath = ['C:\Users\cooper\Documents\gitHubRepositories\cooperAnalysis\matlabCode\operantMatching\learningModels\stan\bernoulli\'];
else
    filePath = ['C:\Users\cooper\Documents\githubRepositories\cooperAnalysis\matlabCode\operantMatching\learningModels\stan\'];
end
switch p.Results.modelName
    case 'sixParam_absPePeAN_bi_tW'
        fit = stan('file',[filePath 'stan_qLearning_6params_absPePeAN_bi_tW.stan'],'data',session_dat,'verbose', true,...
            'iter', p.Results.iter, 'warmup', warmup, 'working_dir', savePath);
    case 'sevenParam_absPePeAN_bi_bias_tW'
        fit = stan('file',[filePath 'stan_qLearning_7params_absPePeAN_bi_bias_tW.stan'],'data',session_dat,'verbose', true,...
            'iter', p.Results.iter, 'warmup', warmup, 'working_dir', savePath);
    case 'sixParam_absPePeAN_exp_tW'
        fit = stan('file',[filePath 'stan_qLearning_6params_absPePeAN_exp_tW.stan'],'data',session_dat,'verbose', true,...
            'iter', p.Results.iter, 'warmup', warmup, 'working_dir', savePath);
    case 'sevenParam_absPePeAN_exp_bias_tW'
        fit = stan('file',[filePath 'stan_qLearning_7params_absPePeAN_exp_bias_tW.stan'],'data',session_dat,'verbose', true,...
            'iter', p.Results.iter, 'warmup', warmup, 'working_dir', savePath);
    otherwise
        error([p.Results.modelName ' does not exist'])
end

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
        if length(tmpCount) == length(baseCount) + 4
            doneFlag = 1;
        end
    end
end
delete([savePath 'diaryTmp.txt'])

%extract samples from the stan fit object
samples = [];
while isempty(samples)
    samples = fit.extract('permuted',true);
end
if isfield(samples, 'y_pred')
    samples = rmfield(samples, 'y_pred');
end
%[~, tbl] = fit.print();

%generate best estimates of parameters
paramEsts = [];
if p.Results.fixedParams
    tmp = eval(['samples.mu_' p.Results.paramNames{setdiff([1:5], paramInds)}]);
    [n,e] = histcounts(tmp, 50);
    [~, maxInd] = max(n);
    paramEsts = median(tmp(tmp > e(maxInd) & tmp < e(maxInd+1)));
else
    for i = 1:length(paramInds)
        tmp = eval(['samples.mu_' p.Results.paramNames{i}]);
        [n,e] = histcounts(tmp, 50);
        [~, maxInd] = max(n);
        paramEsts(i) = median(tmp(tmp > e(maxInd) & tmp < e(maxInd+1)));
    end
end

%plot the distributions of the mouse-level parameters
figure; 
if p.Results.fixedParams
    histogram(eval(['samples.mu_' p.Results.paramNames{setdiff([1:5], paramInds)}]), 100,...
            'Normalization', 'Probability', 'FaceColor', 'k')
        set(gca,'tickdir', 'out')
        xlabel( p.Results.paramNames{setdiff([1:5], paramInds)})
else
    numParams = length(paramInds);
    blue = [0 1 1];
    purp = [0.7 0 1];
    colors = [linspace(blue(1),purp(1),numParams)', linspace(blue(2),purp(2),numParams)', linspace(blue(3),purp(3),numParams)'];
    for i = 1:numParams
        subplot(1,numParams,i); hold on;
        histogram(eval(['samples.mu_' p.Results.paramNames{paramInds(i)}]) , 100,...
            'Normalization', 'Probability', 'FaceColor', colors(i,:), 'EdgeColor', 'none')
        set(gca,'tickdir', 'out') 
        title(p.Results.paramNames{paramInds(i)})
    end
end
titleTxt = strrep([animal ' - ' p.Results.modelName], '_', ' ');
suptitle(titleTxt);
set(gcf,'Renderer', 'Painters')

if p.Results.saveFlag
    %save the full samples
    if p.Results.fixedParams
        sampFile = [animal category, '_', p.Results.modelName, '_', p.Results.paramNames{setdiff([1:5], paramInds)}];
        saveFile = [sampFile '.mat'];
        eval([sampFile,  ' = samples;']);
    else
        sampFile = [animal category '_', p.Results.modelName];
        saveFile = [sampFile '.mat'];
        eval([sampFile,  ' = samples;']);
    end
    %save([savePath saveFile], sampFile, 'paramEsts', 'dayList', 'tbl');
    save([savePath saveFile], sampFile, 'paramEsts', 'dayList');
    
    saveFigurePDF(gcf,[savePath animal category '_' p.Results.modelName  '_posteriors'])
end
    
