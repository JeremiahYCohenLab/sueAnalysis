function plotStanParams(animals, beh, varargin)

p = inputParser;
% default parameters if none given
p.addParameter('revForFlag', 0);
p.addParameter('bernFlag', 1);
p.addParameter('paramNames', {'aNscale', 'aNmin', 'aP', 'aF', 'aPE', 'beta'});
p.addParameter('modelName', ['sixParam_absPePeAN_bi']);
p.parse(varargin{:});

[root, sep] = currComputer();
numParams = length(p.Results.paramNames);
numAnimals = length(animals);

aiS = 0;
sInds = [];
for aI = 1:numAnimals
    if p.Results.bernFlag
        filePath = [root animals{aI} sep  animals{aI} 'sorted' sep 'stan' sep 'bernoulli' sep p.Results.modelName sep];
    else        
        filePath = [root animals{aI} sep  animals{aI} 'sorted' sep 'stan' sep p.Results.modelName sep];
    end
    mdlFolder = dir(filePath);
    behInd = (~cellfun(@isempty,strfind({mdlFolder.name}, beh)));
    fileInd = (~cellfun(@isempty,strfind({mdlFolder.name}, '.mat')));
    tmpInd = behInd & fileInd;
    if sum(tmpInd) == 0
        aiS = aiS + 1;
        sInds = [sInds aI];
    else
        load([filePath mdlFolder(tmpInd).name]);
        samps = eval([animals{aI} beh '_' p.Results.modelName]);

        for ind = 1:numParams
            tmp = eval(['samps.mu_' p.Results.paramNames{ind}]);
            [n,e] = histcounts(tmp, 50);
            [~, maxInd] = max(n);
            avgParams{ind}(aI-aiS, 1) = median(tmp(tmp > e(maxInd) & tmp < e(maxInd+1)));
            avgParams{ind}(aI-aiS, 2) = eval(['std(samps.mu_' p.Results.paramNames{ind} ')']);
             %avgParams{ind}(aI, 1) = median(tmp);   %for median params
        end
    end
end

numAnimals = numAnimals - aiS;
colorsA = cool(numAnimals);
colorsP = cool(numParams);
titles = generateParamTitles(p.Results.paramNames);
x = linspace(0.5, 1.5, numAnimals);

animals(sInds) = [];

figure;
for i = 1:numParams
    subplot(2,numParams,i); hold on;
    for j = 1:numAnimals
        errorbar(x(j), [avgParams{i}(j,1)], [avgParams{i}(j,2)],...
            'Color', colorsA(j,:), 'linewidth', 1.5);
    end
    xticks([])
    xticklabels({' '})
    xlim([0.25 1.75])
    set(gca,'tickdir', 'out')
    if isempty(strfind(p.Results.paramNames{i}, 'beta'))
        ylim([0 1])
    end
    title(titles{i})
end
legend(animals)
for i = 1:numParams
    subplot(2,numParams,numParams+i); hold on;
    histogram(avgParams{i}(:,1), 10, 'FaceColor', colorsP(i,:), 'Normalization', 'probability');
    set(gca, 'box', 'off', 'tickdir', 'out')
end
titleTxt = strrep([p.Results.modelName], '_', ' ');
suptitle(titleTxt);
set(gcf,'Renderer', 'Painters', 'position', [-1919 41 1920 963])
