function [paramEsts] = stan_qLearningFit_dynaVol(xlFile, sheet, category, varargin)

p = inputParser;
% default parameters if none given
p.addParameter('revForFlag', 0);
p.addParameter('bernFlag', 1);
p.addParameter('nonfixedParams', 0);
p.addParameter('fixedParams', []);
% p.addParameter('paramNames',{'aNmin', 'aP', 'aF', 'aPE', 'v', 'beta'}); % animal level
% p.addParameter('modelName', '7params_absPePeAN_scale_int_bias_ord');
p.addParameter('modelName', '5params');
p.addParameter('cmpParam', 1);
% p.addParameter('modelName', 'vkf_fixV_kappa');
p.addParameter('iter', 10000);
p.addParameter('warmup', []);
p.addParameter('saveFlag', 1);
p.addParameter('maxTrial', 1000);
p.addParameter('numChains', 8);
p.addParameter('simFlag', 0);
p.addParameter('control', struct('delta', 0.85))
p.parse(varargin{:});

paramNames = getParamNames_dF(p.Results.modelName,0);
if ~p.Results.nonfixedParams
    paramInds = 1:length(paramNames);
    fullName = ['stan_qLearning_' p.Results.modelName 'CmpVol.stan'];
else
    paramInds = find(~contains(paramNames, p.Results.nonfixedParams));
    fullName = ['stan_qLearning_' p.Results.modelName '_' p.Results.nonfixedParams 'CmpVol.stan'];
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

if p.Results.simFlag
    if p.Results.bernFlag
        savePath = [root 'sim' sep 'stan' sep 'bernoulli' sep p.Results.modelName sep];
    else
        savePath = [root 'sim' sep 'stan' sep p.Results.modelName sep];
    end
    
    if ~exist(savePath)
        mkdir(savePath);
    end
    iteration = 50;
    % generate random paramters
    params = zeros(iteration, 4);
    params(:,1) = 0.9*betarnd(3, 5, iteration,1)+0.1; % lambda [0 1]
    params(:,2) = 10*betarnd(5, 2, iteration,1); % v0 [0 5]
%     params(:,2) = 10*5/7 * ones(iteration,1); % v0 [0 5]
    params(:,3) = 10*(0.9*betarnd(4, 2, iteration,1)+0.1); % omega [0 5]
    params(:,4) = 0.7*betarnd(6, 3, iteration,1)+0.3; % beta [0 1]
    params(:,5) = 0.9*betarnd(6, 3, iteration,1)+0.1; % aF/kappa [0 1]
    % simulation
    choice = zeros(iteration, p.Results.maxTrial);
    outcome = zeros(iteration, p.Results.maxTrial);
    for sim = 1:iteration
        %simulation
        if contains(p.Results.modelName, 'aF')
            [~, outcomeSim, choiceSim] = vkfSim_aF('params', params(sim,:),'randomSeed', sim,'maxTrials', p.Results.maxTrial, 'plotFlag', 0);
        else
            if contains(p.Results.modelName, 'kappa')
               [~, outcomeSim, choiceSim] = vkfSim_kappa('params', params(sim,:),'randomSeed', sim,'maxTrials', p.Results.maxTrial, 'plotFlag', 0); 
            else
                [~, outcomeSim, choiceSim] = vkfSim('params', params(sim,:),'randomSeed', sim,'maxTrials', p.Results.maxTrial, 'plotFlag', 0);
            end
        end
        choiceSim(choiceSim<0) = 0;
        outcomeSim = abs(outcomeSim);
        choice(sim,:) = choiceSim;
        outcome(sim,:) = outcomeSim;
    end
    
    T = p.Results.maxTrial;
    N = iteration;
    Tsesh = p.Results.maxTrial*ones(iteration,1);
    
    session_dat = struct('N',N,'T',T, 'Tsesh', Tsesh, 'choice', choice, 'outcome', outcome);
else
    if p.Results.bernFlag
        savePath = [root sheet sep sheet 'sorted' sep 'stan' sep 'bernoulli' sep 'compare' sep p.Results.modelName sep paramNames{p.Results.cmpParam} sep category sep];
    else
        savePath = [root sheet sep sheet 'sorted' sep 'stan' sep 'compare' sep p.Results.modelName sep paramNames{p.Results.cmpParam} sep category sep];
    end
    if ~exist(savePath)
        mkdir(savePath);
    end
    dayList = getDayList(xlFile, sheet, category);

    for i = 1:length(dayList)
        sessionName = dayList{i};
%         fprintf([sessionName '\n']);
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
        volTemp{i} = behavStruct.vol;
    end

    T = max(Tsesh);
    N = length(dayList);
    choice = zeros(N, T);
    outcome = zeros(N, T);
    vol = zeros(N, T);
    ITI = zeros(N,T);

    for i = 1:N
        choice(i, 1:Tsesh(i)) = choiceTmp{i};
        outcome(i, 1:Tsesh(i)) = outcomeTmp{i};
        vol(i, 1:Tsesh(i)) = volTemp{i};
        ITI(i, 1:Tsesh(i)) = ITItemp{i}/1000; % convert from ms to s
    end
    choice = choice(:,1:min([T p.Results.maxTrial]));
    outcome = outcome(:,1:min([T p.Results.maxTrial]));
    ITI = ITI(:,1:min([T p.Results.maxTrial]));
    Tsesh(Tsesh>p.Results.maxTrial) = p.Results.maxTrial;
    T = min([T p.Results.maxTrial]);
    
    session_dat = struct('N',N,'T',T, 'Tsesh', Tsesh, 'choice', choice, 'outcome', outcome, 'ITI', ITI, 'vol', vol, 'cmpParam', p.Results.cmpParam);
end
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
figure;
colors = [0, 0.8, 0.8;
          1, 0.3, 1];
tmpSamp = samples;
numParams = length(paramNames);
for currP = 1:length(paramNames)
    subplot(2,numParams, currP); hold on;
    mu_L = tmpSamp.(['mu_' paramNames{currP} '_L']);
    mu_H = tmpSamp.(['mu_' paramNames{currP} '_H']);
    if currP == p.Results.cmpParam
        edges = linspace(min([mu_L; mu_H]), max([mu_L; mu_H]), 50);
        histogram(mu_L, edges, 'FaceColor', colors(1,:), 'EdgeColor', 'none', 'Normalization', 'probability');
        histogram(mu_H, edges, 'FaceColor', colors(2,:), 'EdgeColor', 'none', 'Normalization', 'probability');
        
        legend({'Low', 'High'})
    else
       edges = linspace(min(mu_L), max(mu_L), 50);
        histogram(mu_L, edges, 'FaceColor', colors(1,:), 'EdgeColor', 'none', 'Normalization', 'probability')
        legend({'mutual'})
    end
    title(paramNames{currP});
end

muDiff = [tmpSamp.(['mu_' paramNames{p.Results.cmpParam} '_H'])] - [tmpSamp.(['mu_' paramNames{p.Results.cmpParam} '_L'])];
diff = tmpSamp.mu_diff;
subplot(2, 3, 4); hold on;
edges = linspace(min(muDiff), max(muDiff), 50);
histogram(muDiff, edges, 'FaceColor', colors(1,:), 'EdgeColor', 'none', 'Normalization', 'probability');
line([0 0], [0 0.1], 'LineStyle', '--', 'Color', [0.5 0.5 0.5],'LineWidth',2)
title([paramNames{p.Results.cmpParam} ' mu_H-mu_L'], 'Interpreter', 'none');

subplot(2, 3, 5); hold on;
edges = linspace(min(diff), max(diff), 50);
histogram(diff, edges, 'FaceColor', colors(1,:), 'EdgeColor', 'none', 'Normalization', 'probability');
p5__ = summary{{'mu_diff'}, {'p5_'}};
p95__ = summary{{'mu_diff'}, {'p95_'}};
line([p5__ p5__], [0 0.1], 'LineStyle', '--', 'Color', [1 0.5 0.5],'LineWidth',2)
line([p95__ p95__], [0 0.1], 'LineStyle', '--', 'Color', [0.5 0.5 1],'LineWidth',2)
title([paramNames{p.Results.cmpParam} ' muDiff']);

subplot(2,3,6); hold on;
loglikMean = mean(tmpSamp.log_likMean,2);
edges = linspace(min(loglikMean), max(loglikMean), 50);
histogram(loglikMean, edges, 'FaceColor', colors(1,:), 'EdgeColor', 'none', 'Normalization', 'probability');
title('logLL/trial')
screen = get(0,'Screensize');
screen(4) = screen(4) - 100;
set(gcf, 'Position', screen)


if p.Results.saveFlag
    %save the full samples

    sampFile = [sheet category '_', p.Results.modelName '_' paramNames{p.Results.cmpParam}];
    saveFile = [sampFile '.mat'];
    eval([sampFile,  ' = samples;']);
    saveFigurePDF(gcf,[savePath sheet category '_' p.Results.modelName '_' paramNames{p.Results.cmpParam} '_posteriors'])
    save([savePath saveFile], sampFile, 'dayList', 'summary');
        
end
    %save([savePath saveFile], sampFile, 'paramEsts', 'dayList', 'tbl');
    
    
