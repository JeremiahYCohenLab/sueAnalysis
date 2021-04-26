function [paramEsts] = stan_qLearningFit(xlFile, sheet, category, varargin)

p = inputParser;
% default parameters if none given
p.addParameter('revForFlag', 0);
p.addParameter('bernFlag', 1);
p.addParameter('nonfixedParams', 0);
p.addParameter('fixedParams', []);
p.addParameter('paramNames', {'aN', 'aP', 'aF', 'beta','bias'});
p.addParameter('modelName', '5params');
p.addParameter('iter', 20000);
p.addParameter('warmup', []);
p.addParameter('saveFlag', 1);
p.addParameter('maxTrial', 500);
p.parse(varargin{:});

if ~p.Results.nonfixedParams
    paramInds = 1:length(p.Results.paramNames);
    fullName = ['stan_qLearning_' p.Results.modelName '.stan'];
else
    paramInds = find(~contains(p.Results.paramNames, p.Results.nonfixedParams));
    fullName = ['stan_qLearning_' p.Results.modelName '_' p.Results.nonfixedParams '.stan'];
end

if isempty(p.Results.warmup)
    warmup = max(floor(p.Results.iter/2),1);
else
    warmup = p.Results.warmup;
end

[root, sep] = currComputer();
if p.Results.bernFlag
    savePath = [root sheet sep sheet 'sorted' sep 'stan' sep 'bernoulli' sep p.Results.modelName sep];
else
    savePath = [root sheet sep sheet 'sorted' sep 'stan' sep p.Results.modelName sep];
end
if ~exist(savePath)
    mkdir(savePath);
end

[~, dayList, ~] = xlsread([root xlFile], sheet);
[~,col] = find(contains(dayList, category) == 1);
% [~,col] = find(~cellfun(@isempty,strfind(dayList, category)) == 1);
dayList = dayList(2:end,col);
endInd = find(cellfun(@isempty,dayList),1);
if ~isempty(endInd)
    dayList = dayList(1:endInd-1,:);
end

for i = 1:length(dayList)
    sessionName = dayList{i};
    filename = [sessionName '.asc'];
    fprintf([sessionName '\n']);
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

T = max(Tsesh);
N = length(dayList);
choice = zeros(N, T);
outcome = zeros(N, T);
ITI = zeros(N,T);

for i = 1:N
    choice(i, 1:Tsesh(i)) = choiceTmp{i};
    outcome(i, 1:Tsesh(i)) = outcomeTmp{i};
    ITI(i, 1:Tsesh(i)) = ITItemp{i};
end
choice = choice(:,1:min([T p.Results.maxTrial]));
outcome = outcome(:,1:min([T p.Results.maxTrial]));
ITI = ITI(:,1:min([T p.Results.maxTrial]));
Tsesh(Tsesh>p.Results.maxTrial) = p.Results.maxTrial;
T = min([T p.Results.maxTrial]);
%create data structure to feed into stan model
session_dat = struct('N',N,'T',T, 'Tsesh', Tsesh, 'choice', choice, 'outcome', outcome, 'ITI', ITI);
if p.Results.nonfixedParams
    for j = 1:length(paramInds)
        session_dat.(p.Results.paramNames{paramInds(j)}) = p.Results.fixedParams(:,paramInds(j));
    end
end

%run the stan model
if p.Results.bernFlag
    filePath = 'C:\Users\zhixi\Documents\gitRepositories\sueAnalysis\matlabCode\operantMatching\learningModels\stan\bernoulli\';
else
    filePath = 'C:\Users\zhixi\Documents\gitRepositories\sueAnalysis\matlabCode\operantMatching\learningModels\stan\';
end
fit = stan('file',[filePath fullName],'data',session_dat,'verbose',true,...
            'iter', p.Results.iter, 'warmup', warmup, 'working_dir', savePath);
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
%[~, tbl] = fit.print();

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
        tmp = eval(['samples.mu_' p.Results.paramNames{i}]);
        allSamples = [allSamples tmp];
        edges{i} = linspace(min(tmp), max(tmp),50);
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
 figure2; 
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
        histogram(eval(['samples.mu_' p.Results.paramNames{paramInds(i)}]) , 100,...
            'Normalization', 'Probability', 'FaceColor', colors(i,:), 'EdgeColor', 'none')
        set(gca,'tickdir', 'out') 
        title(p.Results.paramNames{paramInds(i)})
    end
end
titleTxt = strrep([sheet ' - ' p.Results.modelName], '_', ' ');
suptitle(titleTxt);
set(gcf,'Renderer', 'Painters')

if p.Results.saveFlag
    %save the full samples
    if p.Results.nonfixedParams
        sampFile = [sheet category, '_', p.Results.modelName, '_', p.Results.nonfixedParams];
        saveFile = [sampFile '.mat'];
        eval([sampFile,  ' = samples;']);
    else
        sampFile = [sheet category '_', p.Results.modelName];
        saveFile = [sampFile '.mat'];
        eval([sampFile,  ' = samples;']);
    end
    %save([savePath saveFile], sampFile, 'paramEsts', 'dayList', 'tbl');
    save([savePath saveFile], sampFile, 'paramEsts', 'dayList');
    
    saveFigurePDF(gcf,[savePath sheet category '_' p.Results.modelName  '_posteriors'])
end
    
