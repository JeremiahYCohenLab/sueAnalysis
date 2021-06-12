function [paramMatx] = plotStanParams(animals, beh, varargin)

p = inputParser;
% default parameters if none given
p.addParameter('revForFlag', 0)
p.addParameter('bernFlag', 1)
p.addParameter('sessionParamsFlag', 0)
p.addParameter('plotFlag', 1)
p.addParameter('biasFlag', 0)
p.addParameter('modelName', ['sevenParam_absPePeAN_scale_int_bias_ord'])
p.parse(varargin{:});

[root, sep] = currComputer();
paramNames = getParamNames_dF(p.Results.modelName, 0);
numParams = length(paramNames);
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
    fileName = [animals{aI} beh '_' p.Results.modelName '.mat'];    
    fileInd = find(strcmp(fileName, {mdlFolder.name}));
    if sum(fileInd) == 0
        aiS = aiS + 1;
        sInds = [sInds aI];
    else
        mdl = load([filePath mdlFolder(fileInd).name]);
        samps = eval(['mdl.' animals{aI} beh '_' p.Results.modelName]);
        dayLists{aI} = mdl.dayList;

        for currP = 1:numParams
            tmp = eval(['samps.mu_' paramNames{currP}]);
            [n,e] = histcounts(tmp, 50);
            [~, maxInd] = max(n);
            avgParams{currP}(aI-aiS, 1) = median(tmp(tmp > e(maxInd) & tmp < e(maxInd+1)));
            avgParams{currP}(aI-aiS, 2) = std(tmp);
        end
        if p.Results.biasFlag
            inds = randperm(size(samps.bias, 2));
            tmp = samps.bias(:, inds(1));
            [n,e] = histcounts(tmp, 50);
            [~, maxInd] = max(n);
            avgParams{numParams+1}(aI-aiS, 1) = median(tmp(tmp > e(maxInd) & tmp < e(maxInd+1)));
            avgParams{numParams+1}(aI-aiS, 2) = std(tmp);
        end
        clear mdl;
    end
end

numAnimals = numAnimals - aiS;
colorsA = cool(numAnimals);
colorsP = cool(numParams);
titles = generateParamTitles(paramNames);
x = linspace(0.5, 1.5, numAnimals);

animals(sInds) = [];

paramMatx = horzcat(avgParams{:});
paramMatx = paramMatx(:, [1:2:size(paramMatx,2)]);

if p.Results.plotFlag
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
        if isempty(strfind(paramNames{i}, 'beta'))
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

    figure;
    [rho, pVal] = corr(paramMatx);
    imagesc(rho);
    colormap(cool(256));
    ax = gca;
    set(ax, 'tickdir', 'out', 'box', 'off')
    set(ax, 'XTick', [1:numParams], 'YTick', [1:numParams]);
    set(ax, 'XTickLabel', titles, 'YTickLabel', titles);
    cb = colorbar('Peer', ax);
    set(cb, 'tickdir', 'out');

    if any(any(pVal < 0.05))
        [sigInds_x, sigInds_y] = find(pVal < 0.05);
        for currS = 1:length(sigInds_x)
            text('Position', [sigInds_x(currS) sigInds_y(currS)  0],'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle',  'string', ...
            num2str(round(pVal(sigInds_x(currS), sigInds_y(currS)), 4)));
        end
    end

    suptitle(titleTxt);
    set(gcf,'Renderer', 'Painters', 'position', [-1520 182 955 765])
    
end
    
    
