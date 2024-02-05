function [paramEsts] = stan_qLearningFit(xlFile, sheet, category, varargin)

p = inputParser;
% default parameters if none given
p.addParameter('revForFlag', 0);
p.addParameter('bernFlag', 1);
p.addParameter('nonfixedParams', 0);
p.addParameter('fixedParams', []);
% p.addParameter('paramNames',{'aNmin', 'aP', 'aF', 'aPE', 'v', 'beta'}); % animal level
% p.addParameter('modelName', '7params_absPePeAN_scale_int_bias_ord');
p.addParameter('modelName', '5params');
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
        savePath = [root sheet sep sheet 'sorted' sep 'stan' sep 'bernoulli' sep p.Results.modelName sep category sep];
    else
        savePath = [root sheet sep sheet 'sorted' sep 'stan' sep p.Results.modelName sep category sep];
    end
    if ~exist(savePath)
        mkdir(savePath);
    end
    [~, dayList, ~] = xlsread([root xlFile], sheet);
    col = cell(1,size(dayList,2));
    col(:) = {category};
    col = cellfun(@strcmp, dayList(1,:), col)>0;
    % [~,col] = find(~cellfun(@isempty,strfind(dayList, category)) == 1);
    dayList = dayList(2:end,col);
    endInd = find(cellfun(@isempty,dayList),1);
    if ~isempty(endInd)
        dayList = dayList(1:endInd-1,:);
    end

    for i = 1:length(dayList)
        sessionName = dayList{i};
        filename = [sessionName '.asc'];
%         fprintf([sessionName '\n']);
        behSessionData = loadBehavioralData(filename, p.Results.revForFlag);
        behavStruct = parseBehavioralData(behSessionData, p.Results.maxTrial);

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
    end
    transGuess = [0.4 0.3 0.3;
                  0.2 0.8 0;
                  0.2 0 0.8];
    emisGuess = [0.5 0.5;
                 1 0;
                 0 1];    
             
    

    T = max(Tsesh);
    N = length(dayList);
    choice = zeros(N, T);
    outcome = zeros(N, T);
    ITI = zeros(N,T);

    for i = 1:N
        choice(i, 1:Tsesh(i)) = choiceTmp{i};
        outcome(i, 1:Tsesh(i)) = outcomeTmp{i};
        ITI(i, 1:Tsesh(i)) = ITItemp{i}/1000; % convert from ms to s
    end
    choice = choice(:,1:min([T p.Results.maxTrial]));
    outcome = outcome(:,1:min([T p.Results.maxTrial]));
    ITI = ITI(:,1:min([T p.Results.maxTrial]));
    Tsesh(Tsesh>p.Results.maxTrial) = p.Results.maxTrial;
    T = min([T p.Results.maxTrial]);
    
    session_dat = struct('N',N,'T',T, 'Tsesh', Tsesh, 'choice', choice, 'outcome', outcome, 'ITI', ITI);
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
        if length(unique(tmp))>1
            if inds(i) < 50
                paramEsts(i) = median(tmp(tmp >= edgeTmp(inds(i)) & tmp < edgeTmp(inds(i)+1)));
            else
                paramEsts(i) = edgeTmp(inds(i));
            end
        else
            paramEsts(i) = unique(tmp);
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
        subplot(1,numParams+1,i); hold on;
        histogram(eval(['samples.mu_' paramNames{paramInds(i)}]) , 100,...
            'Normalization', 'Probability', 'FaceColor', colors(i,:), 'EdgeColor', 'none')
%         line([paramEsts(i) paramEsts(i)], [0 0.05], 'color', [0 0 0]);
        set(gca,'tickdir', 'out') 
        title(paramNames{paramInds(i)})
    end
    subplot(1,numParams+1,numParams+1); hold on;
    sumLL = mean(samples.log_lik, 2);
    histogram(sumLL , 100,...
        'Normalization', 'Probability', 'FaceColor', [0.5 0.5 0.5], 'EdgeColor', 'none')
%         line([paramEsts(i) paramEsts(i)], [0 0.05], 'color', [0 0 0]);
    set(gca,'tickdir', 'out') 
    title('logLL')
    
end
titleTxt = strrep([sheet ' - ' p.Results.modelName], '_', ' ');
sgtitle(titleTxt);
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
sgtitle(titleTxt);
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
            if ~isnan(str2double(sheet))
                sampFile = ['m' sheet category, '_', p.Results.modelName, '_', p.Results.nonfixedParams];
            else
                sampFile = [sheet category, '_', p.Results.modelName, '_', p.Results.nonfixedParams];
            end

            saveFile = [sampFile '.mat'];
            eval([sampFile,  ' = samples;']);
        else
            if ~isnan(str2double(sheet))
                sampFile = ['m' sheet category '_', p.Results.modelName];
            else
                sampFile = [sheet category '_', p.Results.modelName];
            end
            
            saveFile = [sampFile '.mat'];
            eval([sampFile,  ' = samples;']);
        end
        saveFigurePDF(pFig,[savePath sheet category '_' p.Results.modelName  '_posteriors'])
        saveFigurePDF(dFig,[savePath sheet category '_' p.Results.modelName  '_divergence'])
        save([savePath saveFile], sampFile, 'paramEsts', 'dayList', 'summary');
        
    end
    % save([savePath saveFile], sampFile, 'paramEsts', 'dayList', 'tbl');
    
end
    
