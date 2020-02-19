function [tbl, samples] = stan_qLearningFit(xlFile, animal, category, varargin)

p = inputParser;
% default parameters if none given
p.addParameter('revForFlag', 0);
p.addParameter('bernFlag', 1);
p.addParameter('fixedParams', 0);
p.addParameter('params', ['fiveParamOfourStart_preS']);
p.addParameter('paramInds', []);
p.addParameter('paramNames', {'aN', 'aP', 'aF', 'beta', 'v'});
p.addParameter('modelName', ['fourParam']);
p.addParameter('iter', 2000);
p.addParameter('warmup', []);
p.addParameter('saveFlag', 1)
p.parse(varargin{:});

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

    choiceTmp{i} = behavStruct.allChoices;
    if p.Results.bernFlag
        choiceTmp{i}(choiceTmp{i} == -1) = 0;
    else
        choiceTmp{i} = choiceTmp{i} + 1;
        choiceTmp{i}(choiceTmp{i} == 0) = 1;
    end
    outcomeTmp{i} = abs(behavStruct.allRewards); 
    Tsesh(i,1) = length(outcomeTmp{i});
end

T = max(Tsesh);
N = length(dayList);
choice = zeros(N, T);
outcome = zeros(N, T);

for i = 1:N
    choice(i, 1:Tsesh(i)) = choiceTmp{i};
    outcome(i, 1:Tsesh(i)) = outcomeTmp{i};
end

%create data structure to feed into stan model
if p.Results.fixedParams
    [allParams, modelNames, ~] = xlsread('stanParams.xlsx', animal);
    [~,col] = find(~cellfun(@isempty,strfind(modelNames, p.Results.params)) == 1);
    params = allParams(paramInds,col);
    session_dat = struct('N',N,'T',T, 'Tsesh', Tsesh, 'choice', choice, 'outcome', outcome, 'params', params);
else
    session_dat = struct('N',N,'T',T, 'Tsesh', Tsesh, 'choice', choice, 'outcome', outcome);
end

%run the stan model
if p.Results.bernFlag
    filePath = ['C:\Users\cooper\Documents\gitHubRepositories\cooperAnalysis\matlabCode\operantMatching\learningModels\stan\bernoulli\'];
else
    filePath = ['C:\Users\cooper\Documents\githubRepositories\cooperAnalysis\matlabCode\operantMatching\learningModels\stan\'];
end
switch p.Results.modelName
    case 'twoParam'
        fit = stan('file',[filePath 'stan_qLearning_2params.stan'],'data',session_dat,'verbose',true,...
            'iter', p.Results.iter, 'warmup', warmup, 'working_dir', savePath);
    case 'threeParam'
        fit = stan('file',[filePath 'stan_qLearning_3params.stan'],'data',session_dat,'verbose',true,...
            'iter', p.Results.iter, 'warmup', warmup, 'working_dir', savePath);
    case 'fourParam'
        fit = stan('file',[filePath 'stan_qLearning_4params.stan'],'data',session_dat,'verbose',true,...
            'iter', p.Results.iter, 'warmup', warmup, 'working_dir', savePath);
   case 'fourParam_manSoft'
        fit = stan('file',[filePath 'stan_qLearning_4params_manSoft.stan'],'data',session_dat,'verbose',true,...
            'iter', p.Results.iter, 'warmup', warmup, 'working_dir', savePath);
    case 'fiveParam_k'
        fit = stan('file',[filePath 'stan_qLearning_5params_k.stan'],'data',session_dat,'verbose', true,...
            'iter', p.Results.iter, 'warmup', warmup, 'working_dir', savePath);
    case 'fourParamFixed'
        scriptName = ['stan_qLearning_4params_fixedParams_' ...
            p.Results.paramNames{setdiff([1:4], paramInds)} '.stan'];
        fit = stan('file',[filePath scriptName],'data',session_dat,'verbose', true,...
            'iter', p.Results.iter, 'warmup', warmup, 'working_dir', savePath);
    case 'fourParamO'
        fit = stan('file',[filePath 'stan_qLearning_4params_opponency.stan'],'data',session_dat,'verbose', true,...
            'iter', p.Results.iter, 'warmup', warmup, 'working_dir', savePath);
    case 'fiveParamO'
        fit = stan('file',[filePath 'stan_qLearning_5params_opponency.stan'],'data',session_dat,'verbose', true,...
            'iter', p.Results.iter, 'warmup', warmup, 'working_dir', savePath);
    case 'sixParamO_rBarStart'
        fit = stan('file',[filePath 'stan_qLearning_6params_opponency_rBarStart.stan'],'data',session_dat,'verbose', true,...
            'iter', p.Results.iter, 'warmup', warmup, 'working_dir', savePath);
    case 'fiveParamO_peUpdate'
        fit = stan('file',[filePath 'stan_qLearning_5params_opponency_peUpdate.stan'],'data',session_dat,'verbose', true,...
            'iter', p.Results.iter, 'warmup', warmup, 'working_dir', savePath);    
    case 'sixParamO_rBarStart_peUpdate'
        fit = stan('file',[filePath 'stan_qLearning_6params_opponency_rBarStart_peUpdate.stan'],'data',session_dat,'verbose', true,...
            'iter', p.Results.iter, 'warmup', warmup, 'working_dir', savePath); 
    case 'fiveParam_peBeta'
        fit = stan('file',[filePath 'stan_qLearning_5params_peBeta.stan'],'data',session_dat,'verbose', true,...
            'iter', p.Results.iter, 'warmup', warmup, 'working_dir', savePath);
    case 'fiveParam_peBeta_avg'
        fit = stan('file',[filePath 'stan_qLearning_5params_peBeta_avg.stan'],'data',session_dat,'verbose', true,...
            'iter', p.Results.iter, 'warmup', warmup, 'working_dir', savePath);
    case 'fiveParam_peBeta_fixedMax'
        fit = stan('file',[filePath 'stan_qLearning_5params_peBeta_fixedMax.stan'],'data',session_dat,'verbose', true,...
            'iter', p.Results.iter, 'warmup', warmup, 'working_dir', savePath);
    case 'sixParam_peBeta'
        fit = stan('file',[filePath 'stan_qLearning_6params_peBeta.stan'],'data',session_dat,'verbose', true,...
            'iter', p.Results.iter, 'warmup', warmup, 'working_dir', savePath);
    case 'sixParam_pePeBeta'
        fit = stan('file',[filePath 'stan_qLearning_6params_pePeBeta.stan'],'data',session_dat,'verbose', true,...
            'iter', p.Results.iter, 'warmup', warmup, 'working_dir', savePath);
    case 'sixParam_peBeta_diff'
        fit = stan('file',[filePath 'stan_qLearning_6params_peBeta_diff.stan'],'data',session_dat,'verbose', true,...
            'iter', p.Results.iter, 'warmup', warmup, 'working_dir', savePath);
    case 'sixParam_absPePeBeta'
        fit = stan('file',[filePath 'stan_qLearning_6params_absPePeBeta.stan'],'data',session_dat,'verbose', true,...
            'iter', p.Results.iter, 'warmup', warmup, 'working_dir', savePath);
    case 'sevenParam_peBeta_k'
        fit = stan('file',[filePath 'stan_qLearning_7params_peBeta_k.stan'],'data',session_dat,'verbose', true,...
            'iter', p.Results.iter, 'warmup', warmup, 'working_dir', savePath);
    case 'fourParam_rBeta_scale'
        fit = stan('file',[filePath 'stan_qLearning_4params_rBeta_scale.stan'],'data',session_dat,'verbose', true,...
            'iter', p.Results.iter, 'warmup', warmup, 'working_dir', savePath);
    case 'fourParam_rBeta_confQ'
        fit = stan('file',[filePath 'stan_qLearning_4params_rBeta_confQ.stan'],'data',session_dat,'verbose', true,...
            'iter', p.Results.iter, 'warmup', warmup, 'working_dir', savePath);
    case 'fiveParam_rBeta_confQ'
        fit = stan('file',[filePath 'stan_qLearning_5params_rBeta_confQ.stan'],'data',session_dat,'verbose', true,...
            'iter', p.Results.iter, 'warmup', warmup, 'working_dir', savePath);
    case 'fiveParam_rBeta_scale'
        fit = stan('file',[filePath 'stan_qLearning_5params_rBeta_scale.stan'],'data',session_dat,'verbose', true,...
            'iter', p.Results.iter, 'warmup', warmup, 'working_dir', savePath);
    case 'sixParam_rBeta_scale'
        fit = stan('file',[filePath 'stan_qLearning_6params_rBeta_scale.stan'],'data',session_dat,'verbose', true,...
            'iter', p.Results.iter, 'warmup', warmup, 'working_dir', savePath );
    case 'sixParam_rBeta_scale_min'
        fit = stan('file',[filePath 'stan_qLearning_6params_rBeta_scale_min.stan'],'data',session_dat,'verbose', true,...
            'iter', p.Results.iter, 'warmup', warmup, 'working_dir', savePath );
    case 'sixParam_rBeta_scale_max'
        fit = stan('file',[filePath 'stan_qLearning_6params_rBeta_scale_max.stan'],'data',session_dat,'verbose', true,...
            'iter', p.Results.iter, 'warmup', warmup, 'working_dir', savePath );
    case 'sixParam_rBeta_scale_initR'
        fit = stan('file',[filePath 'stan_qLearning_6params_rBeta_scale_initR.stan'],'data',session_dat,'verbose', true,...
            'iter', p.Results.iter, 'warmup', warmup, 'working_dir', savePath );
    case 'sixParam_rBeta_scale_kappa'
        fit = stan('file',[filePath 'stan_qLearning_6params_rBeta_scale_kappa.stan'],'data',session_dat,'verbose', true,...
            'iter', p.Results.iter, 'warmup', warmup, 'working_dir', savePath );
    case 'sixParam_rBeta_rRPE'
        fit = stan('file',[filePath 'stan_qLearning_6params_rBeta_rRPE.stan'],'data',session_dat,'verbose', true,...
            'iter', p.Results.iter, 'warmup', warmup, 'working_dir', savePath);
    case 'sixParam_rBeta_rRPE_noMin'
        fit = stan('file',[filePath 'stan_qLearning_6params_rBeta_rRPE_noMin.stan'],'data',session_dat,'verbose', true,...
            'iter', p.Results.iter, 'warmup', warmup, 'working_dir', savePath);
    case 'sixParam_rBeta_rV'
        fit = stan('file',[filePath 'stan_qLearning_6params_rBeta_rV.stan'],'data',session_dat,'verbose', true,...
            'iter', p.Results.iter, 'warmup', warmup, 'working_dir', savePath);
    case 'fourParam_pePeAN_noScale'
        fit = stan('file',[filePath 'stan_qLearning_4params_pePeAN_noScale.stan'],'data',session_dat,'verbose', true,...
            'iter', p.Results.iter, 'warmup', warmup, 'working_dir', savePath); 
    case 'fiveParam_pePeAN_noScale'
        fit = stan('file',[filePath 'stan_qLearning_5params_pePeAN_noScale.stan'],'data',session_dat,'verbose', true,...
            'iter', p.Results.iter, 'warmup', warmup, 'working_dir', savePath); 
    case 'sixParam_pePeAN'
        fit = stan('file',[filePath 'stan_qLearning_6params_pePeAN.stan'],'data',session_dat,'verbose', true,...
            'iter', p.Results.iter, 'warmup', warmup, 'working_dir', savePath);
    case 'fiveParam_pePeAN_lag'
        fit = stan('file',[filePath 'stan_qLearning_5params_pePeAN_lag.stan'],'data',session_dat,'verbose', true,...
            'iter', p.Results.iter, 'warmup', warmup, 'working_dir', savePath);
    case 'sixParam_pePeAN_lag'
        fit = stan('file',[filePath 'stan_qLearning_6params_pePeAN_lag.stan'],'data',session_dat,'verbose', true,...
            'iter', p.Results.iter, 'warmup', warmup, 'working_dir', savePath);
    case 'fiveParam_absPePeAN'
        fit = stan('file',[filePath 'stan_qLearning_5params_absPePeAN.stan'],'data',session_dat,'verbose', true,...
            'iter', p.Results.iter, 'warmup', warmup, 'working_dir', savePath); 
    case 'sixParam_absPePeAN'
        fit = stan('file',[filePath 'stan_qLearning_6params_absPePeAN.stan'],'data',session_dat,'verbose', true,...
            'iter', p.Results.iter, 'warmup', warmup, 'working_dir', savePath);
    case 'sixParam_absPePeAN_bi'
        fit = stan('file',[filePath 'stan_qLearning_6params_absPePeAN_bi.stan'],'data',session_dat,'verbose', true,...
            'iter', p.Results.iter, 'warmup', warmup, 'working_dir', savePath);
    case 'sixParam_absPePeAN_biSep'
        fit = stan('file',[filePath 'stan_qLearning_6params_absPePeAN_biSep.stan'],'data',session_dat,'verbose', true,...
            'iter', p.Results.iter, 'warmup', warmup, 'working_dir', savePath);
    case 'sevenParam_absPePeAN_biSep_f'
        fit = stan('file',[filePath 'stan_qLearning_7params_absPePeAN_biSep_f.stan'],'data',session_dat,'verbose', true,...
            'iter', p.Results.iter, 'warmup', warmup, 'working_dir', savePath);
    case 'sevenParam_absPePeAN_bi_k'
        fit = stan('file',[filePath 'stan_qLearning_7params_absPePeAN_bi_k.stan'],'data',session_dat,'verbose', true,...
            'iter', p.Results.iter, 'warmup', warmup, 'working_dir', savePath);
    case 'sevenParam_absPePeLR'
        fit = stan('file',[filePath 'stan_qLearning_7params_absPePeLR.stan'],'data',session_dat,'verbose', true,...
            'iter', p.Results.iter, 'warmup', warmup, 'working_dir', savePath);
    case 'sevenParam_pePeAN_k'
        fit = stan('file',[filePath 'stan_qLearning_7params_pePeAN_k.stan'],'data',session_dat,'verbose', true,...
            'iter', p.Results.iter, 'warmup', warmup, 'working_dir', savePath);
    case 'sevenParam_peLR'
        fit = stan('file',[filePath 'stan_qLearning_7params_peLR.stan'],'data',session_dat,'verbose', true,...
            'iter', p.Results.iter, 'warmup', warmup, 'working_dir', savePath);
    case 'sevenParam_pePeLR'
        fit = stan('file',[filePath 'stan_qLearning_7params_pePeLR.stan'],'data',session_dat,'verbose', true,...
            'iter', p.Results.iter, 'warmup', warmup, 'working_dir', savePath);
    case 'eightParam_absPePeLR_k'
        fit = stan('file',[filePath 'stan_qLearning_8params_absPePeLR_k.stan'],'data',session_dat,'verbose', true,...
            'iter', p.Results.iter, 'warmup', warmup, 'working_dir', savePath);
    case 'sixParam_rAN'
        fit = stan('file',[filePath 'stan_qLearning_6params_rAN.stan'],'data',session_dat,'verbose', true,...
            'iter', p.Results.iter, 'warmup', warmup, 'working_dir', savePath);
    case 'sixParam_rBeta_oppo'
        fit = stan('file',[filePath 'stan_qLearning_6params_rBeta_oppo.stan'],'data',session_dat,'verbose', true,...
            'iter', p.Results.iter, 'warmup', warmup, 'working_dir', savePath);
    case 'sevenParam_rBeta_oppo_rStart'
        fit = stan('file',[filePath 'stan_qLearning_7params_rBeta_oppo_rStart.stan'],'data',session_dat,'verbose', true,...
            'iter', p.Results.iter, 'warmup', warmup, 'working_dir', savePath);
    case 'sixParam_rBeta_oppo_rCont'
        fit = stan('file',[filePath 'stan_qLearning_6params_rBeta_oppo_rCont.stan'],'data',session_dat,'verbose', true,...
            'iter', p.Results.iter, 'warmup', warmup, 'working_dir', savePath);
    case 'eightParam_rBeta_pePeAN'
        fit = stan('file',[filePath 'stan_qLearning_8params_rBeta_pePeAN.stan'],'data',session_dat,'verbose', true,...
            'iter', p.Results.iter, 'warmup', warmup, 'working_dir', savePath);
    case 'eightParam_rBeta_absPePeAN'
        fit = stan('file',[filePath 'stan_qLearning_8params_rBeta_absPePeAN.stan'],'data',session_dat,'verbose', true,...
            'iter', p.Results.iter, 'warmup', warmup, 'working_dir', savePath);
    case 'tenParam_rBeta_pePeAN_kb'
        fit = stan('file',[filePath 'stan_qLearning_10params_rBeta_pePeAN_kb.stan'],'data',session_dat,'verbose', true,...
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
    paramEsts = median(eval(['samples.mu_' p.Results.paramNames{setdiff([1:5], paramInds)}]));
else
    for i = 1:length(paramInds)
        paramEsts(i)  = median(eval(['samples.mu_' p.Results.paramNames{i}]));
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
    
