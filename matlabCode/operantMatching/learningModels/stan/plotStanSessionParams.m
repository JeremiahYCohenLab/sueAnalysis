function [paramMatx] = plotStanSessionParams(animals, varargin)

p = inputParser;
% default parameters if none given
p.addParameter('beh', 'good')
p.addParameter('bernFlag', 1)
% p.addParameter('modelName','7params_absPePeAN_scale_int_bias_ord')
p.addParameter('modelName', '5paramsHmm');
p.addParameter('plotFlag', 1)
p.addParameter('saveFigFlag', 1)
p.parse(varargin{:});

bin1 = 1;
bin2 = 35;

[root, sep] = currComputer();
paramNames = getParamNames_dF(p.Results.modelName, 0);
numParams = length(paramNames);
numAnimals = length(animals);

%get info for plotting
colors = cool(numParams);
titles = generateParamTitles(paramNames);

paramMatx = [];
sInds = [];
for aI = 1:numAnimals
    if p.Results.bernFlag
        
        % filePath = [root animals{aI} sep  animals{aI} 'sorted' sep 'stan' sep 'bernoulli' sep...
        %     p.Results.modelName sep p.Results.beh sep animals{aI} p.Results.beh '_' p.Results.modelName '.mat'];
        savePath = [root animals{aI} sep animals{aI} 'sorted' sep 'stan' sep 'bernoulli' sep p.Results.modelName sep p.Results.beh sep];
    
    
        if  ~isempty(regexp(animals{aI}, '^[A-Z]', 'once')) || ~isempty(regexp(animals{aI}, '^[a-z]', 'once'))
            sampFile = [animals{aI} p.Results.beh '_',  p.Results.modelName];
        else
            sampFile = ['m' animals{aI} p.Results.beh '_',  p.Results.modelName];
        end

        path = [root animals{aI} sep  animals{aI} 'sorted' sep 'stan' sep 'bernoulli' sep...
            p.Results.modelName sep p.Results.beh sep];
        filePath = [path sampFile '.mat'];
   
    else        
        filePath = [root animals{aI} sep  animals{aI} 'sorted' sep 'stan' sep...
            p.Results.modelName sep p.Results.beh sep animals{aI} p.Results.beh '_' p.Results.modelName '.mat'];
    end

    mdl = load(filePath);
    structName = sampFile;
    samps = mdl.(structName);
    numSesh = length(mdl.dayList);
    dates = zeros(size(mdl.dayList));
    for i = 1:length(dates)
        tmp = mdl.dayList{i};
        tmp = str2double(tmp(end-3:end));
        if isnan(tmp)
            tmp = mdl.dayList{i};
            tmp = str2double(tmp(end-4:end-1));
        end
        dates(i) = tmp;
    end

    [~, ind] = sort(dates);
    ind(ind) = 1:length(dates);

    avgParams = nan(numSesh, numParams);
    for currS = 1:numSesh
        allSamples = [];
            % first bins
            edges = cell(1,length(paramNames));
            for i = 1:length(paramNames)
                if size(samps.(paramNames{i}),2) == 1 % if there's only animal level
                    tmp = samps.(paramNames{i});    
                else
                    tmp = samps.(paramNames{i})(:,currS);
                end
                
                if strcmp(paramNames{i}, 'aPE') || strcmp(paramNames{i}, 'v') % log scale
                    tmp = log(tmp);
                end
                
                allSamples = [allSamples tmp];
                edges{i} = linspace(min(tmp), max(tmp)+0.0001,bin1+1);
            end
%             figure2;
%             subplot(1,2,1); histogram(allSamples(:,1), 30);
%             subplot(1,2,2); histogram(allSamples(:,3), 30);

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
    paramMatx = [paramMatx; avgParams];
    clear mdl;
    
    if p.Results.plotFlag
        figure;
        for currP = 1:numParams
            subplot(2,numParams+3,currP)
            histogram(avgParams(:,currP), 5, 'facecolor', colors(currP,:))
            if strcmp(paramNames{currP}, 'aPE') || strcmp(paramNames{currP}, 'v')
                title(['log(' titles{currP} ')'])
            else
                title(titles{currP})
            end
            set(gca, 'tickdir', 'out', 'box', 'off')
            if currP == 1
                ylabel('count')
            end
            
            subplot(2, numParams+3, currP+numParams+3)
            scatter(ind, avgParams(:, currP), 10, colors(currP,:), 'filled');
            [rho, pval] = corr(ind, avgParams(:, currP));
            xlabel('days')
            title(sprintf('%0.2f p:%0.2f', rho, pval));
            
        end
        
        subplot(2,numParams+3, [2*numParams+4]);
        if sum(contains(fieldnames(samps),'log_likMean'))>0
            scatter(ind, mean([samps.log_likMean]), 10, [0.3, 0.3, 0.3], 'filled');
            xlabel('days')
            title('trialLL');
        end
        
        subplot(2,numParams+3, [2*numParams+5]);
        sumLL = mean(samps.log_lik,2);
         histogram(sumLL , 100, 'Normalization', 'Probability', 'FaceColor', [0.5 0.5 0.5], 'EdgeColor', 'none')
        set(gca,'tickdir', 'out') 
        title('meanLL/session')
        
        subplot(2,numParams+3,[numParams+1:numParams+3])
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

        titleTxt = [animals{aI} ' ' p.Results.beh ' ' strrep([p.Results.modelName], '_', ' ')];
        sgtitle(titleTxt);
        set(gcf,'Renderer', 'Painters', 'position', [-1928 278 1924 1066])

        if p.Results.saveFigFlag
            saveFigurePDF(gcf,[savePath animals{aI} p.Results.beh '_' p.Results.modelName '_parameters.pdf']);
        end
        
        dFig = figure;
        for currPy = 1:numParams + size(samps.sigma,2)
            if currPy <= numParams
                tmpY = eval(['samps.mu_' paramNames{currPy}]);
            else
                tmpY = samps.sigma(:,currPy-numParams);
            end
            tmpY_d = tmpY(logical(samps.divergent__));
            tmpY = tmpY(~logical(samps.divergent__));
            for currPx = 1:numParams + size(samps.sigma,2)
                subplot(numParams+size(samps.sigma,2), numParams+size(samps.sigma,2),(currPy-1)*(numParams + size(samps.sigma,2)) + currPx); hold on;

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
                if currPy == numParams + size(samps.sigma,2)
                    if currPx <= numParams
                       xlabel(paramNames{currPx})
                    else
                       xlabel(['sigma' '(' num2str(currPx-numParams) ')'])
                    end
                end

            end
        end
        titleTxt = [titleTxt ' (divergence rate = ' num2str(sum(samps.divergent__)/length(samps.divergent__)) ')'];
        sgtitle(titleTxt);
        screenSize = get(0,'Screensize');
        screenSize(4) = screenSize(4) - 100;
        set(dFig, 'renderer', 'painters', 'position', screenSize)
        if p.Results.saveFigFlag
            saveFigurePDF(gcf,[savePath animals{aI} p.Results.beh '_' p.Results.modelName '_' 'divergence.pdf'])
        end
        cFig = figure;
        scatterAll(avgParams, paramNames, 7, 'm')
        sgtitle(titleTxt);
        if p.Results.saveFigFlag
            saveFigurePDF(cFig,[savePath animals{aI} p.Results.beh '_' p.Results.modelName '_' 'estimatedParams.pdf'])
        end      
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
