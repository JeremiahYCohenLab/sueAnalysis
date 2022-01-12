function [avgParams] = plotStanSessionParamsSim(modelName, subFolder, varargin)

p = inputParser;
% default parameters if none given
p.addParameter('bernFlag', 1)
p.addParameter('plotFlag', 1)
p.addParameter('saveFigFlag', 1)
p.parse(varargin{:});

bin1 = 1;
bin2 = 40;

[root, sep] = currComputer();
paramNames = getParamNames_dF(modelName, 1);
numParams = length(paramNames);

%get info for plotting
colors = cool(numParams);
titles = generateParamTitles(paramNames);

if p.Results.bernFlag
    filePath = [root 'sim' sep 'stan' sep 'bernoulli' sep...
        modelName sep subFolder sep 'sim_' modelName '.mat'];
    savePath = [root 'sim' sep 'stan' sep 'bernoulli' sep...
        modelName sep subFolder sep];
else        
    filePath = [root 'sim' sep  'stan' sep...
        modelName sep subFolder sep 'sim_' modelName '.mat'];
    savePath = [root 'sim' sep 'stan' sep ...
        modelName sep subFolder sep];
end

mdl = load(filePath);
structName = ['sim_' modelName];
samps = mdl.(structName);
numSesh = size(mdl.choice, 1);

avgParams = nan(numSesh, numParams);
for currS = 1:numSesh
    allSamples = [];
        % first bins
        edges = cell(1,length(paramNames));
        for i = 1:length(paramNames)
            if size(samps.(paramNames{i}),2) == 1
                tmp = samps.(paramNames{i});    
            else
                tmp = samps.(paramNames{i})(:,currS);
            end
            allSamples = [allSamples tmp];
            edges{i} = linspace(min(tmp), max(tmp)+0.0001,bin1+1);
        end
%         figure2;
%         subplot(1,2,1); histogram(allSamples(:,1), 30);
%         subplot(1,2,2); histogram(allSamples(:,2), 30);
        
        n = histcnd(allSamples,edges); %bin samples by multiple dimensions
        [~, inds] = myMaxAll(n); %find the bin with max num in bin
        % second bins
        edges2 = cell(1,length(paramNames));
        for i = 1:length(paramNames) %use previous best bin as newbin
            edgeTmp = edges{i};
            edges2{i} = linspace(edgeTmp(inds(i)),edgeTmp(inds(i)+1)+0.0001,bin2+1);
        end
        n = histcnd(allSamples,edges2); %bin samples by multiple dimensions
        [~, inds] = myMaxAll(n); %find the bin with max num in bin   
        for i = 1:length(paramNames) %use median in bin as best estimate
            tmp = allSamples(:,i);
            edgeTmp = edges2{i};
            if inds(i) < bin2
                avgParams(currS,i) = median(tmp(tmp >= edgeTmp(inds(i)) & tmp < edgeTmp(inds(i)+1)));
            else
                avgParams(currS,i) = median(tmp(tmp >= edgeTmp(inds(i)) & tmp <= edgeTmp(inds(i)+1)));
            end
        end
end

if p.Results.plotFlag
    
    figure;
    screenSize = get(0,'Screensize');
    screenSize(4) = screenSize(4) - 100;
    set(gcf, 'Position', screenSize)
    for currP = 1:numParams
        subplot(2,numParams+3,currP+numParams+3)
        histogram(avgParams(:,currP), 10, 'facecolor', colors(currP,:))
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
    
     for currP = 1:numParams
        subplot(2,numParams+3,currP)
        scatter(mdl.params(:,currP), avgParams(:,currP), 8, colors(currP,:), 'filled')
        lowB = min([mdl.params(:,currP); avgParams(:,currP)]);
        upB = max([mdl.params(:,currP); avgParams(:,currP)]);
        line([lowB upB], [lowB upB], 'color', [0.7 0.7 0.7]);
        xlim([lowB upB])
        ylim([lowB upB])
        [rho, pVal] = corr([mdl.params(:,currP), avgParams(:,currP)]);
        title(sprintf('%.2f p=%.2f', rho(1,2), pVal(1,2)));
        if strcmp(paramNames{currP}, 'aPE') || strcmp(paramNames{currP}, 'v')
            xlabel(['log(' titles{currP} ')'])
        else
            xlabel(titles{currP})
        end
        set(gca, 'tickdir', 'out', 'box', 'off')
        if currP == 1
            ylabel('estimation')
        end
    end
    
    

    subplot(2,numParams+3,2*numParams+4:2*numParams+6)
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
    title('Estimation')
    
    subplot(2,numParams+3,numParams+1:numParams+3)
    [rho, pVal] = corr(mdl.params);
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
    title('Simulation')
    
    titleTxt = ['sim ' strrep([modelName], '_', ' ') ' ' num2str(size(mdl.outcome, 2))];
    suptitle(titleTxt);
    %set(gcf,'Renderer', 'Painters', 'position', [-1928 278 1924 566])
    if p.Results.saveFigFlag
        saveFigurePDF(gcf,[savePath modelName '_' 'parameters'])
    end
    figure;
    scatterAll(mdl.params, paramNames, 8, 'm');
    suptitle('Simulation params');
    if p.Results.saveFigFlag
        saveFigurePDF(gcf,[savePath modelName '_' 'simulationParameters'])
    end
    figure;
    scatterAll(avgParams, paramNames, 8, 'c');
    suptitle('Estimated params')
    if p.Results.saveFigFlag
        saveFigurePDF(gcf,[savePath modelName '_' 'estimationParameters'])
    end    
    
    dFig = figure;
    paramNames = getParamNames_dF(modelName,0);
    numParams = length(paramNames);
    for currPy = 1:numParams*2
        if currPy <= numParams
            tmpY = eval(['samps.mu_' paramNames{currPy}]);
        else
            tmpY = samps.sigma(:,currPy-numParams);
        end
        tmpY_d = tmpY(logical(samps.divergent__));
        tmpY = tmpY(~logical(samps.divergent__));
        for currPx = 1:2*numParams
            subplot(2*numParams,2*numParams,[(currPy-1)*2*numParams + currPx]); hold on;

            if currPy == currPx
                h = histogram(tmpY, 30, 'FaceColor', 'c', 'normalization', 'probability');
                histogram(tmpY_d, h.BinEdges, 'FaceColor', 'm', 'normalization', 'probability')
            else       
                if currPx <= numParams
                    tmpX = eval(['samps.mu_' paramNames{currPx}]);
                else
                    tmpX = samps.sigma(:,currPx-numParams);
                end
                tmpX_d = tmpX(logical(samps.divergent__));
                tmpX = tmpX(~logical(samps.divergent__));
                scatter(tmpX, tmpY, 5, 'c', 'filled')
                scatter(tmpX_d, tmpY_d, 5, 'm', 'filled')
            end

            if currPx == 1
                if currPy <= numParams
                   ylabel(paramNames{currPy})
                else
                   ylabel(['sigma' '(' num2str(currPy-numParams) ')'])
                end
            end
            if currPy == 2*numParams
                if currPx <= numParams
                   xlabel(paramNames{currPx})
                else
                   xlabel(['sigma' '(' num2str(currPx-numParams) ')'])
                end
            end

        end
    end
    titleTxt = [titleTxt ' (divergence rate = ' num2str(sum(samps.divergent__)/length(samps.divergent__)) ')'];
    suptitle(titleTxt);
    set(dFig, 'renderer', 'painters', 'position', screenSize)
    if p.Results.saveFigFlag
        saveFigurePDF(gcf,[savePath modelName '_' 'divergence'])
    end
end
    



%amend cell figures as one pdf
% if p.Results.saveFigFlag
%     dirTmp = dir('C:\Users\cooper\Desktop\analysis\model param pdfs\');
%     for currFig = 3:length(dirTmp)
%         append_pdfs(['C:\Users\cooper\Desktop\analysis\model param pdfs\' p.Results.modelName '_all.pdf'], ...
%             [dirTmp(currFig).folder '\' dirTmp(currFig).name]);
%     end
% end 
