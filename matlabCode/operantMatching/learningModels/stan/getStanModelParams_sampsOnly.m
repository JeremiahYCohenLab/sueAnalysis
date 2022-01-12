function [params, modelName, ll, noSession] = getStanModelParams_sampsOnly(animalName, category, modelName, numSamps, varargin)


p = inputParser;
% default parameters if none given
p.addParameter('biasFlag',1)
p.addParameter('sessionParamsFlag', 0)
p.addParameter('sessionName', [])
p.addParameter('plotFlag', 0)
p.parse(varargin{:});

[root, sep] = currComputer();
paramNames = getParamNames_dF(modelName, p.Results.biasFlag);
sampFile = [animalName category '_', modelName];
path = [root animalName sep animalName 'sorted' sep 'stan' sep 'bernoulli' sep modelName sep category sep];
modelPath = [path sampFile '.mat'];
%get model params
noSession = false;
params = nan(numSamps, length(paramNames));
ll = NaN;
load(modelPath);
if p.Results.sessionParamsFlag
    sessionInd = find(~cellfun(@isempty,strfind(dayList,p.Results.sessionName)));
    if isempty(sessionInd)
        noSession = true;
    return
    end
end

tmp = whos;
samples = eval(tmp(1).name);

inds = find(samples.divergent__<1);
inds = inds(randperm(sum(samples.divergent__<1), numSamps));

for currP = 1:length(paramNames)
    if p.Results.sessionParamsFlag
        tmp = eval(['samples.' paramNames{currP}]);
        tmp = tmp(:,sessionInd);
    else
        if ~strcmp(paramNames{currP}, 'bias')
            tmp = eval(['samples.mu_' paramNames{currP}]);
        else
            tmp = zeros(size(tmp));
        end
    end
    params(:,currP) = tmp(inds);
end
if p.Results.sessionParamsFlag
    ll = samples.log_lik(inds,sessionInd);
end

if p.Results.plotFlag
    figure; 
    colors = cool(length(paramNames));
    for i = 1:length(paramNames)
        subplot(1,length(paramNames),i);
        edges = linspace(min(params(:,i)),max(params(:,i)), 25);
        histogram(params(:,i), edges, 'FaceColor', colors(i,:),'Normalization', 'Probability');
        title(paramNames{i})
    end
end
end 
    
    
    
    
    

