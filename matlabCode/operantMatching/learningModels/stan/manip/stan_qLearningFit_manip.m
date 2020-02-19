function [fit, samples, tbl] = stan_qLearningFit_manip(xlFile, animal, pre, post, varargin)

p = inputParser;
% default parameters if none given
p.addParameter('revForFlag', 0);
p.addParameter('bernFlag', 1);
p.addParameter('changeFlag', 0);
p.addParameter('saveFlag', 1);
p.addParameter('paramNames', {'aN', 'aP', 'aF', 'beta', 'v'});
p.addParameter('modelName', ['fiveParamO']);
p.addParameter('iter', 2000);
p.addParameter('warmup', []);
p.parse(varargin{:});

if isempty(p.Results.warmup)
    warmup = max(floor(p.Results.iter/2),1);
else
    warmup = p.Results.warmup;
end

if contains(p.Results.modelName, 'hyper')
    hyperFlag = 1;
else
    hyperFlag = 0;
end

[root, sep] = currComputer();

[~, dayList, ~] = xlsread(xlFile, animal);
[~,colPre] = find(~cellfun(@isempty,strfind(dayList, pre)) == 1);
dayListPre = dayList(2:end,colPre);
endIndPre = find(cellfun(@isempty,dayListPre),1);
if ~isempty(endIndPre)
    dayListPre = dayListPre(1:endIndPre-1,:);
end
[~,colPost] = find(~cellfun(@isempty,strfind(dayList, post)) == 1);
dayListPost = dayList(2:end,colPost);
endIndPost = find(cellfun(@isempty,dayListPost),1);
if ~isempty(endIndPost)
    dayListPost = dayListPost(1:endIndPost-1,:);
end
dayList = [dayListPre dayListPost];

if p.Results.bernFlag
    savePath = [root animal sep animal 'sorted' sep 'stan' sep 'manip' sep 'bernoulli' sep p.Results.modelName '_manip' sep];
else
    savePath = [root animal sep animal 'sorted' sep 'stan' sep 'manip' sep p.Results.modelName '_manip' sep];
end
if ~exist(savePath)
    mkdir(savePath);
end

for i = 1:length(dayListPre)
    sessionName = dayListPre{i};
    filename = [sessionName '.asc'];
    [behSessionData, unCorrectedBlockSwitch, out] = loadBehavioralData(filename, p.Results.revForFlag);
    behavStruct = parseBehavioralData(behSessionData, unCorrectedBlockSwitch);

    choiceTmpPre{i} = behavStruct.allChoices;
    if p.Results.bernFlag
        choiceTmpPre{i}(choiceTmpPre{i} == -1) = 0;
    else
        choiceTmpPre{i} = choiceTmpPre{i} + 1;
        choiceTmpPre{i}(choiceTmpPre{i} == 0) = 1;
    end
    outcomeTmpPre{i} = abs(behavStruct.allRewards); 
    TseshPre(i,1) = length(outcomeTmpPre{i});
end

for i = 1:length(dayListPost)
    sessionName = dayListPost{i};
    filename = [sessionName '.asc'];
    [behSessionData, unCorrectedBlockSwitch, out] = loadBehavioralData(filename, p.Results.revForFlag);
    behavStruct = parseBehavioralData(behSessionData, unCorrectedBlockSwitch);

    choiceTmpPost{i} = behavStruct.allChoices;
    if p.Results.bernFlag
        choiceTmpPost{i}(choiceTmpPost{i} == -1) = 0;
    else
        choiceTmpPost{i} = choiceTmpPost{i} + 1;
        choiceTmpPost{i}(choiceTmpPost{i} == 0) = 1;
    end
    outcomeTmpPost{i} = abs(behavStruct.allRewards); 
    TseshPost(i,1) = length(outcomeTmpPost{i});
end

T = max([TseshPre; TseshPost]);
N = length(dayListPre);
M = length(dayListPost);
choicePre = zeros(N,T);
outcomePre = zeros(N,T);
choicePost = zeros(M,T);
outcomePost = zeros(M,T);

for i = 1:N
    choicePre(i, 1:TseshPre(i)) = choiceTmpPre{i};
    outcomePre(i, 1:TseshPre(i)) = outcomeTmpPre{i};
end
for i = 1:M
    choicePost(i, 1:TseshPost(i)) = choiceTmpPost{i};
    outcomePost(i, 1:TseshPost(i)) = outcomeTmpPost{i};
end

session_dat = struct('N',N,'M',M,'T',T,'Tsesh', TseshPre, 'choice', choicePre, 'outcome', outcomePre,...
    'TseshM', TseshPost, 'choiceM', choicePost, 'outcomeM', outcomePost);
if p.Results.bernFlag
    filePath = ['C:\Users\cooper\Documents\githubRepositories\cooperAnalysis\matlabCode\operantMatching\learningModels\stan\manip\bernoulli\'];
else
    filePath = ['C:\Users\cooper\Documents\githubRepositories\cooperAnalysis\matlabCode\operantMatching\learningModels\stan\manip\'];
end
switch p.Results.modelName
    case 'twoParam'
        fit = stan('file',[filePath 'stan_qLearning_2params_manip.stan'],'data',session_dat,'verbose',true,...
            'iter', p.Results.iter, 'warmup', warmup, 'working_dir', savePath);
    case 'threeParam'
        fit = stan('file',[filePath 'stan_qLearning_3params_manip.stan'],'data',session_dat,'verbose',true,...
            'iter', p.Results.iter, 'warmup', warmup, 'working_dir', savePath);
    case 'fourParam'
        fit = stan('file',[filePath 'stan_qLearning_4params_manip.stan'],'data',session_dat,'verbose',true,...
            'iter', p.Results.iter, 'warmup', warmup, 'working_dir', savePath);
   case 'fourParam_manSoft'
        fit = stan('file',[filePath 'stan_qLearning_4params_manSoft_manip.stan'],'data',session_dat,'verbose',true,...
            'iter', p.Results.iter, 'warmup', warmup, 'working_dir', savePath);
    case 'fiveParam_k'
        fit = stan('file',[filePath 'stan_qLearning_5params_k_manip.stan'],'data',session_dat,'verbose', true,...
            'iter', p.Results.iter, 'warmup', warmup, 'working_dir', savePath);
    case 'fourParamFixed'
        scriptName = ['stan_qLearning_4params_fixedParams_' ...
            p.Results.paramNames{setdiff([1:4], paramInds)} '_manip.stan'];
        fit = stan('file',[filePath scriptName],'data',session_dat,'verbose', true,...
            'iter', p.Results.iter, 'warmup', warmup, 'working_dir', savePath);
    case 'fourParamO'
        fit = stan('file',[filePath 'stan_qLearning_4params_opponency_manip.stan'],'data',session_dat,'verbose', true,...
            'iter', p.Results.iter, 'warmup', warmup, 'working_dir', savePath);
    case 'fiveParamO'
        fit = stan('file',[filePath 'stan_qLearning_5params_opponency_manip.stan'],'data',session_dat,'verbose', true,...
            'iter', p.Results.iter, 'warmup', warmup, 'working_dir', savePath);
    case 'sixParamO_rBarStart'
        fit = stan('file',[filePath 'stan_qLearning_6params_opponency_rBarStart_manip.stan'],'data',session_dat,'verbose', true,...
            'iter', p.Results.iter, 'warmup', warmup, 'working_dir', savePath);
    case 'fiveParamO_peUpdate'
        fit = stan('file',[filePath 'stan_qLearning_5params_opponency_peUpdate_manip.stan'],'data',session_dat,'verbose', true,...
            'iter', p.Results.iter, 'warmup', warmup, 'working_dir', savePath);    
    case 'sixParamO_rBarStart_peUpdate'
        fit = stan('file',[filePath 'stan_qLearning_6params_opponency_rBarStart_peUpdate_manip.stan'],'data',session_dat,'verbose', true,...
            'iter', p.Results.iter, 'warmup', warmup, 'working_dir', savePath); 
    case 'fiveParam_peBeta'
        fit = stan('file',[filePath 'stan_qLearning_5params_peBeta_manip.stan'],'data',session_dat,'verbose', true,...
            'iter', p.Results.iter, 'warmup', warmup, 'working_dir', savePath);
    case 'fiveParam_peBeta_avg'
        fit = stan('file',[filePath 'stan_qLearning_5params_peBeta_avg_manip.stan'],'data',session_dat,'verbose', true,...
            'iter', p.Results.iter, 'warmup', warmup, 'working_dir', savePath);
    case 'fiveParam_peBeta_fixedMax'
        fit = stan('file',[filePath 'stan_qLearning_5params_peBeta_fixedMax_manip.stan'],'data',session_dat,'verbose', true,...
            'iter', p.Results.iter, 'warmup', warmup, 'working_dir', savePath);
    case 'sixParam_peBeta'
        fit = stan('file',[filePath 'stan_qLearning_6params_peBeta_manip.stan'],'data',session_dat,'verbose', true,...
            'iter', p.Results.iter, 'warmup', warmup, 'working_dir', savePath);
    case 'sixParam_pePeBeta'
        fit = stan('file',[filePath 'stan_qLearning_6params_pePeBeta_manip.stan'],'data',session_dat,'verbose', true,...
            'iter', p.Results.iter, 'warmup', warmup, 'working_dir', savePath);
    case 'sixParam_peBeta_diff'
        fit = stan('file',[filePath 'stan_qLearning_6params_peBeta_diff_manip.stan'],'data',session_dat,'verbose', true,...
            'iter', p.Results.iter, 'warmup', warmup, 'working_dir', savePath);
    case 'sixParam_absPePeBeta'
        fit = stan('file',[filePath 'stan_qLearning_6params_absPePeBeta_manip.stan'],'data',session_dat,'verbose', true,...
            'iter', p.Results.iter, 'warmup', warmup, 'working_dir', savePath);
    case 'sevenParam_peBeta_k'
        fit = stan('file',[filePath 'stan_qLearning_7params_peBeta_k_manip.stan'],'data',session_dat,'verbose', true,...
            'iter', p.Results.iter, 'warmup', warmup, 'working_dir', savePath);
    case 'fourParam_rBeta_scale'
        fit = stan('file',[filePath 'stan_qLearning_4params_rBeta_scale_manip.stan'],'data',session_dat,'verbose', true,...
            'iter', p.Results.iter, 'warmup', warmup, 'working_dir', savePath);
    case 'fourParam_rBeta_confQ'
        fit = stan('file',[filePath 'stan_qLearning_4params_rBeta_confQ_manip.stan'],'data',session_dat,'verbose', true,...
            'iter', p.Results.iter, 'warmup', warmup, 'working_dir', savePath);
    case 'fiveParam_rBeta_confQ'
        fit = stan('file',[filePath 'stan_qLearning_5params_rBeta_confQ_manip.stan'],'data',session_dat,'verbose', true,...
            'iter', p.Results.iter, 'warmup', warmup, 'working_dir', savePath);
    case 'fiveParam_rBeta_scale'
        fit = stan('file',[filePath 'stan_qLearning_5params_rBeta_scale_manip.stan'],'data',session_dat,'verbose', true,...
            'iter', p.Results.iter, 'warmup', warmup, 'working_dir', savePath);
    case 'sixParam_rBeta_scale'
        fit = stan('file',[filePath 'stan_qLearning_6params_rBeta_scale_manip.stan'],'data',session_dat,'verbose', true,...
            'iter', p.Results.iter, 'warmup', warmup, 'working_dir', savePath );
    case 'sixParam_rBeta_scale_min'
        fit = stan('file',[filePath 'stan_qLearning_6params_rBeta_scale_min_manip.stan'],'data',session_dat,'verbose', true,...
            'iter', p.Results.iter, 'warmup', warmup, 'working_dir', savePath );
    case 'sixParam_rBeta_scale_max'
        fit = stan('file',[filePath 'stan_qLearning_6params_rBeta_scale_max_manip.stan'],'data',session_dat,'verbose', true,...
            'iter', p.Results.iter, 'warmup', warmup, 'working_dir', savePath );
    case 'sixParam_rBeta_scale_initR'
        fit = stan('file',[filePath 'stan_qLearning_6params_rBeta_scale_initR_manip.stan'],'data',session_dat,'verbose', true,...
            'iter', p.Results.iter, 'warmup', warmup, 'working_dir', savePath );
    case 'sixParam_rBeta_scale_kappa'
        fit = stan('file',[filePath 'stan_qLearning_6params_rBeta_scale_kappa_manip.stan'],'data',session_dat,'verbose', true,...
            'iter', p.Results.iter, 'warmup', warmup, 'working_dir', savePath );
    case 'sixParam_rBeta_rRPE'
        fit = stan('file',[filePath 'stan_qLearning_6params_rBeta_rRPE_manip.stan'],'data',session_dat,'verbose', true,...
            'iter', p.Results.iter, 'warmup', warmup, 'working_dir', savePath);
    case 'sixParam_rBeta_rRPE_noMin'
        fit = stan('file',[filePath 'stan_qLearning_6params_rBeta_rRPE_noMin_manip.stan'],'data',session_dat,'verbose', true,...
            'iter', p.Results.iter, 'warmup', warmup, 'working_dir', savePath);
    case 'sixParam_rBeta_rV'
        fit = stan('file',[filePath 'stan_qLearning_6params_rBeta_rV_manip.stan'],'data',session_dat,'verbose', true,...
            'iter', p.Results.iter, 'warmup', warmup, 'working_dir', savePath);
    case 'fourParam_pePeAN_noScale'
        fit = stan('file',[filePath 'stan_qLearning_4params_pePeAN_noScale_manip.stan'],'data',session_dat,'verbose', true,...
            'iter', p.Results.iter, 'warmup', warmup, 'working_dir', savePath); 
    case 'fiveParam_pePeAN_noScale'
        fit = stan('file',[filePath 'stan_qLearning_5params_pePeAN_noScale_manip.stan'],'data',session_dat,'verbose', true,...
            'iter', p.Results.iter, 'warmup', warmup, 'working_dir', savePath); 
    case 'sixParam_pePeAN'
        fit = stan('file',[filePath 'stan_qLearning_6params_pePeAN_manip.stan'],'data',session_dat,'verbose', true,...
            'iter', p.Results.iter, 'warmup', warmup, 'working_dir', savePath);
    case 'fiveParam_pePeAN_lag'
        fit = stan('file',[filePath 'stan_qLearning_5params_pePeAN_lag_manip.stan'],'data',session_dat,'verbose', true,...
            'iter', p.Results.iter, 'warmup', warmup, 'working_dir', savePath);
    case 'sixParam_pePeAN_lag'
        fit = stan('file',[filePath 'stan_qLearning_6params_pePeAN_lag_manip.stan'],'data',session_dat,'verbose', true,...
            'iter', p.Results.iter, 'warmup', warmup, 'working_dir', savePath);
    case 'fiveParam_absPePeAN'
        fit = stan('file',[filePath 'stan_qLearning_5params_absPePeAN_manip.stan'],'data',session_dat,'verbose', true,...
            'iter', p.Results.iter, 'warmup', warmup, 'working_dir', savePath); 
    case 'sixParam_absPePeAN'
        fit = stan('file',[filePath 'stan_qLearning_6params_absPePeAN_manip.stan'],'data',session_dat,'verbose', true,...
            'iter', p.Results.iter, 'warmup', warmup, 'working_dir', savePath);
    case 'sixParam_absPePeAN_bi'
        fit = stan('file',[filePath 'stan_qLearning_6params_absPePeAN_bi_manip.stan'],'data',session_dat,'verbose', true,...
            'iter', p.Results.iter, 'warmup', warmup, 'working_dir', savePath);
    case 'sevenParam_absPePeAN_bi_k'
        fit = stan('file',[filePath 'stan_qLearning_7params_absPePeAN_bi_k_manip.stan'],'data',session_dat,'verbose', true,...
            'iter', p.Results.iter, 'warmup', warmup, 'working_dir', savePath);
    case 'sevenParam_absPePeLR'
        fit = stan('file',[filePath 'stan_qLearning_7params_absPePeLR_manip.stan'],'data',session_dat,'verbose', true,...
            'iter', p.Results.iter, 'warmup', warmup, 'working_dir', savePath);
    case 'sevenParam_pePeAN_k'
        fit = stan('file',[filePath 'stan_qLearning_7params_pePeAN_k_manip.stan'],'data',session_dat,'verbose', true,...
            'iter', p.Results.iter, 'warmup', warmup, 'working_dir', savePath);
    case 'sevenParam_peLR'
        fit = stan('file',[filePath 'stan_qLearning_7params_peLR_manip.stan'],'data',session_dat,'verbose', true,...
            'iter', p.Results.iter, 'warmup', warmup, 'working_dir', savePath);
    case 'sevenParam_pePeLR'
        fit = stan('file',[filePath 'stan_qLearning_7params_pePeLR_manip.stan'],'data',session_dat,'verbose', true,...
            'iter', p.Results.iter, 'warmup', warmup, 'working_dir', savePath);
    case 'eightParam_absPePeLR_k'
        fit = stan('file',[filePath 'stan_qLearning_8params_absPePeLR_k_manip.stan'],'data',session_dat,'verbose', true,...
            'iter', p.Results.iter, 'warmup', warmup, 'working_dir', savePath);
    case 'sixParam_rAN'
        fit = stan('file',[filePath 'stan_qLearning_6params_rAN_manip.stan'],'data',session_dat,'verbose', true,...
            'iter', p.Results.iter, 'warmup', warmup, 'working_dir', savePath);
    case 'sixParam_rBeta_oppo'
        fit = stan('file',[filePath 'stan_qLearning_6params_rBeta_oppo_manip.stan'],'data',session_dat,'verbose', true,...
            'iter', p.Results.iter, 'warmup', warmup, 'working_dir', savePath);
    case 'sevenParam_rBeta_oppo_rStart'
        fit = stan('file',[filePath 'stan_qLearning_7params_rBeta_oppo_rStart_manip.stan'],'data',session_dat,'verbose', true,...
            'iter', p.Results.iter, 'warmup', warmup, 'working_dir', savePath);
    case 'sixParam_rBeta_oppo_rCont'
        fit = stan('file',[filePath 'stan_qLearning_6params_rBeta_oppo_rCont_manip.stan'],'data',session_dat,'verbose', true,...
            'iter', p.Results.iter, 'warmup', warmup, 'working_dir', savePath);
    case 'eightParam_rBeta_pePeAN'
        fit = stan('file',[filePath 'stan_qLearning_8params_rBeta_pePeAN_manip.stan'],'data',session_dat,'verbose', true,...
            'iter', p.Results.iter, 'warmup', warmup, 'working_dir', savePath);
    case 'eightParam_rBeta_absPePeAN'
        fit = stan('file',[filePath 'stan_qLearning_8params_rBeta_absPePeAN_manip.stan'],'data',session_dat,'verbose', true,...
            'iter', p.Results.iter, 'warmup', warmup, 'working_dir', savePath);
    case 'tenParam_rBeta_pePeAN_kb'
        fit = stan('file',[filePath 'stan_qLearning_10params_rBeta_pePeAN_kb_manip.stan'],'data',session_dat,'verbose', true,...
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
%[~, tbl] = fit.print();


%generate best estimates of parameters
numParams = length(p.Results.paramNames);

%plot the distributions of the mouse-level parameters
figure; 
blue = [0 1 1];
purp = [0.7 0 1];
colors = [linspace(blue(1),purp(1),numParams*2)', linspace(blue(2),purp(2),numParams*2)', linspace(blue(3),purp(3),numParams*2)'];

if p.Results.changeFlag
    for i = 1:numParams
        subplot(2,numParams,i); hold on;
        histogram(eval(['samples.mu_' p.Results.paramNames{i}]) , 100,...
            'Normalization', 'Probability', 'FaceColor', colors(i,:), 'EdgeColor', 'none')
        set(gca,'tickdir', 'out') 
        title(p.Results.paramNames{i})

        subplot(2,numParams,i+numParams); hold on;
        histogram(eval(['samples.d_mu_' p.Results.paramNames{i}]) , 100,...
            'Normalization', 'Probability', 'FaceColor', colors(i+numParams,:), 'EdgeColor', 'none')
        set(gca,'tickdir', 'out') 
        title(['d_' p.Results.paramNames{i}], 'interpreter', 'none') 
    end
else
    for i = 1:numParams
        subplot(1,numParams,i); hold on;
        histogram(eval(['samples.mu_' p.Results.paramNames{i}]) , 100,...
            'Normalization', 'Probability', 'FaceColor', blue, 'EdgeColor', 'none')

        histogram(eval(['samples.d_mu_' p.Results.paramNames{i}]) , 100,...
            'Normalization', 'Probability', 'FaceColor', purp, 'EdgeColor', 'none')
        set(gca,'tickdir', 'out') 
        title(p.Results.paramNames{i})
    end
    legend('pre', 'post')
end
titleTxt = strrep([animal ' - ' p.Results.modelName], '_', ' ');
suptitle(titleTxt);
set(gcf,'Renderer', 'Painters')

if p.Results.saveFlag
    paramEsts = nan(2,length(numParams));
    for i = 1:numParams
        paramEsts(1,i)  = median(eval(['samples.mu_' p.Results.paramNames{i}]));
        paramEsts(2,i)  = median(eval(['samples.d_mu_' p.Results.paramNames{i}]));
    end

    %save the full samples
    if p.Results.saveFlag
        sampFile = [animal pre post '_', p.Results.modelName '_manip'];
        saveFile = [sampFile '.mat'];
        eval([sampFile,  ' = samples;']);
        %save([savePath saveFile], sampFile, 'paramEsts', 'dayList', 'tbl')
        save([savePath saveFile], sampFile, 'paramEsts', 'dayList')
    end
    
    saveFigurePDF(gcf,[savePath animal pre post '_' p.Results.modelName  '_posteriors'])
end