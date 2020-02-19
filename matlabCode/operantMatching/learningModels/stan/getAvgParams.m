function [pE] = getAvgParams(varargin)

p = inputParser;
% default parameters if none given
p.addParameter('aList', [{'CG33'}, {'CG34'}, {'CG37'}, {'CG39'}, {'CG40'}]);
p.addParameter('beh', 'preS');
p.addParameter('modelName', ['fourParam']);
p.addParameter('paramNames', {'aN', 'aP','aF', 'beta'});
p.addParameter('bernFlag', 0)
p.addParameter('avgType', 'median');

[root, sep] = currComputer;

p.parse(varargin{:});
    
pE = nan(length(p.Results.aList), length(p.Results.paramNames));
for i = 1:length(p.Results.aList)
    animalName = p.Results.aList{i};
        if p.Results.bernFlag
            modelPath = [root animalName sep animalName 'sorted' sep 'stan' sep 'bernoulli' sep p.Results.modelName sep animalName...
            p.Results.beh '_' p.Results.modelName '.mat'];
        else
            modelPath = [root animalName sep animalName 'sorted' sep 'stan' sep p.Results.modelName sep animalName...
            p.Results.beh '_' p.Results.modelName '.mat'];
        end
    load(modelPath);
    for j = 1:length(p.Results.paramNames)
        samples = eval([animalName p.Results.beh '_' p.Results.modelName]);
        pE(i,j) = median(eval(['samples.mu_' p.Results.paramNames{j}]));
    end
end

switch p.Results.avgType
    case 'median'
        pE = median(pE);
    case 'mean'
        pE= mean(pE);
end

end

