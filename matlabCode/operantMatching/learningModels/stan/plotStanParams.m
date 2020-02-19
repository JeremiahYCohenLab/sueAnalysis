function plotStanParams(animals, beh, varargin)

p = inputParser;
% default parameters if none given
p.addParameter('revForFlag', 0);
p.addParameter('bernFlag', 0);
p.addParameter('paramNames', {'aN', 'aP', 'aF', 'v', 'betaScale'});
p.addParameter('modelName', ['fiveParam_rBeta_scale']);
p.parse(varargin{:});

[root, sep] = currComputer();
numParams = length(p.Results.paramNames);
numAnimals = length(animals);

for aI = 1:length(animals)
    if p.Results.bernFlag
        filePath = [root animals{aI} sep  animals{aI} 'sorted' sep 'stan' sep 'bernoulli' sep p.Results.modelName sep];
    else        
        filePath = [root animals{aI} sep  animals{aI} 'sorted' sep 'stan' sep p.Results.modelName sep];
    end
    mdlFolder = dir(filePath);
    behInd = (~cellfun(@isempty,strfind({mdlFolder.name}, beh)));
    fileInd = (~cellfun(@isempty,strfind({mdlFolder.name}, '.mat')));
    tmpInd = behInd & fileInd;
    load([filePath mdlFolder(tmpInd).name]);
    samps = eval([animals{aI} beh '_' p.Results.modelName]);

    for ind = 1:numParams
        avgParams{ind}(aI, 1) = eval(['mean(samps.mu_' p.Results.paramNames{ind} ')']);
        avgParams{ind}(aI, 2) = eval(['std(samps.mu_' p.Results.paramNames{ind} ')']);
    end
    
end

blue = [0 1 1];
purp = [0.7 0 1];
colors = [linspace(blue(1),purp(1),numAnimals)', linspace(blue(2),purp(2),numAnimals)', linspace(blue(3),purp(3),numAnimals)'];

titles = generateParamTitles(p.Results.paramNames);
x = linspace(0.5, 1.5, numAnimals);

figure;
for i = 1:numParams
    subplot(1,numParams,i); hold on;
    for j = 1:length(animals)
        errorbar(x(j), [avgParams{i}(j,1)], [avgParams{i}(j,2)],...
            'Color', colors(j,:), 'linewidth', 1.5);
    end
    xticks([1])
    xticklabels({' '})
    xlim([0.25 1.75])
    set(gca,'tickdir', 'out')
    if isempty(strfind(p.Results.paramNames{i}, 'beta'))
        ylim([0 1])
    end
    title(titles{i})
end
legend(animals)
titleTxt = strrep([p.Results.modelName], '_', ' ');
suptitle(titleTxt);
set(gcf,'Renderer', 'Painters')
