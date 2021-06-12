function [paramMatx] = plotStanSessionParams(animals, varargin)

p = inputParser;
% default parameters if none given
p.addParameter('beh', 'good')
p.addParameter('bernFlag', 1)
p.addParameter('modelName','7params_absPePeAN_scale_int_bias_ord')
p.addParameter('plotFlag', 1)
p.addParameter('saveFigFlag', 0)
p.parse(varargin{:});

[root, sep] = currComputer();
paramNames = getParamNames_dF(p.Results.modelName, 1);
numParams = length(paramNames);
numAnimals = length(animals);

%get info for plotting
colors = cool(numParams);
titles = generateParamTitles(paramNames);

paramMatx = [];
sInds = [];
for aI = 1:numAnimals
    if p.Results.bernFlag
        filePath = [root animals{aI} sep  animals{aI} 'sorted' sep 'stan' sep 'bernoulli' sep...
            p.Results.modelName sep p.Results.beh sep animals{aI} p.Results.beh '_' p.Results.modelName '.mat'];
    else        
        filePath = [root animals{aI} sep  animals{aI} 'sorted' sep 'stan' sep...
            p.Results.modelName sep p.Results.beh sep animals{aI} p.Results.beh '_' p.Results.modelName '.mat'];
    end

    mdl = load(filePath);
    structName = [animals{aI} p.Results.beh '_' p.Results.modelName];
    samps = mdl.(structName);
    numSesh = length(mdl.dayList);

    avgParams = nan(numSesh, numParams);
    for currS = 1:numSesh
        for currP = 1:numParams
            tmp = samps.(paramNames{currP})(:,currS);
            if strcmp(paramNames{currP}, 'aPE') || strcmp(paramNames{currP}, 'v')
                tmp = log(tmp);
            end
            [n,e] = histcounts(tmp, 50);
            [~, maxInd] = max(n);
            avgParams(currS, currP) = median(tmp(tmp > e(maxInd) & tmp < e(maxInd+1)));
        end
    end
    paramMatx = [paramMatx; avgParams];
    clear mdl;
    
    if p.Results.plotFlag
        figure;
        for currP = 1:numParams
            subplot(1,numParams+3,currP)
            histogram(avgParams(:,currP), 5, 'facecolor', colors(currP,:))
            if strcmp(paramNames{currP}, 'aPE') || strcmp(paramNames{currP}, 'v')
                xlabel(['log(' titles{currP} ')'])
            else
                xlabel(titles{currP})
            end
            set(gca, 'tickdir', 'out', 'box', 'off')
            if currP == 1
                ylabel('count')
            end
        end

        subplot(1,numParams+3,[numParams+1:numParams+3])
        [rho, pVal] = corr(avgParams);
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
        xtickangle(20)
        ytickangle(45)

        titleTxt = [animals{aI} ' ' strrep([p.Results.modelName], '_', ' ')];
        suptitle(titleTxt);
        set(gcf,'Renderer', 'Painters', 'position', [-1928 278 1924 566])

        if p.Results.saveFigFlag
            saveFigurePDF(gcf,['C:\Users\cooper\Desktop\analysis\model param pdfs\' p.Results.modelName '_' animals{aI}])
        end
    end
    
end


%amend cell figures as one pdf
if p.Results.saveFigFlag
    dirTmp = dir('C:\Users\cooper\Desktop\analysis\model param pdfs\');
    for currFig = 3:length(dirTmp)
        append_pdfs(['C:\Users\cooper\Desktop\analysis\model param pdfs\' p.Results.modelName '_all.pdf'], ...
            [dirTmp(currFig).folder '\' dirTmp(currFig).name]);
    end
end 
