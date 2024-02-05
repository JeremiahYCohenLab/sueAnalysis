function photometryAnalysisSession(session, modelName)
% load data
[root, sep] = currComputer();
pd = parseSessionString_df(session, root, sep);
dataPath = [pd.sortedFolder session '_photometryCombinewithKH.mat'];
category = 'good';
numSamps = 2000;
numBins = 4;
psthBinNum = 4;
if ~exist(dataPath, 'file')
    fprintf([session 'no photometry data \n'])
    return
else
    load(dataPath);
end

myKernel = ones(2,1);
myKernel = myKernel/sum(myKernel);
LCmatChoiceFilteredG = [];
mPFCmatChoiceFilteredG = [];
LCNmatChoiceFilteredG = [];

for j = 1:size(mPFCmatChoiceG,1)
    temp = conv(LCmatChoiceG(j,:), myKernel);
    temp = temp((floor(0.5*length(myKernel))+1):(end-floor(0.5*length(myKernel))));
    LCmatChoiceFilteredG(j,:) = temp;

    temp = conv(LCNmatChoiceG(j,:), myKernel);
    temp = temp((floor(0.5*length(myKernel))+1):(end-floor(0.5*length(myKernel))));
    LCNmatChoiceFilteredG(j,:) = temp;
    
    temp = conv(mPFCmatChoiceG(j,:), myKernel);
    temp = temp((floor(0.5*length(myKernel))+1):(end-floor(0.5*length(myKernel))));
    mPFCmatChoiceFilteredG(j,:) = temp;        
end

% LCmatChoiceFilteredR = [];
% mPFCmatChoiceFilteredR = [];
% NACmatChoiceFilteredR = [];
% 
% for j = 1:size(LCmatChoiceR,1)
%     temp = conv(LCmatChoiceR(j,:), myKernel);
%     temp = temp((floor(0.5*length(myKernel))+1):(end-floor(0.5*length(myKernel))));
%     LCmatChoiceFilteredR(j,:) = temp;
% 
%     temp = conv(NACmatChoiceR(j,:), myKernel);
%     temp = temp((floor(0.5*length(myKernel))+1):(end-floor(0.5*length(myKernel))));
%     NACmatChoiceFilteredR(j,:) = temp;
% 
%     temp = conv(mPFCmatChoiceR(j,:), myKernel);
%     temp = temp((floor(0.5*length(myKernel))+1):(end-floor(0.5*length(myKernel))));
%     mPFCmatChoiceFilteredR(j,:) = temp;        
% end

if stepSize > 1
    midPointsFilter = midPoints + 0.5*length(myKernel)*stepSize/1000;
else
    midPointsFilter = midPoints + 0.5*length(myKernel)*stepSize;
end
midPointsFilter = midPointsFilter(1:length(temp));

os = behAnalysisNoPlot_opMD(session, 'simpleFlag', 1);

focusWins = {[-2 -1] [0.3 1.8]};
rwdWin = focusWins{2};
names = {'mPFC', 'LCN', 'LC'};
if ~strcmp(modelName, 'none')
    params = getStanModelParams_sampsOnly(pd.animalName, category, modelName, numSamps, 'sessionName', session);
    t = inferModelVar(session, params, modelName);
    pe = t.pe;
end
for c  = 1
    if c == 1
        currC = 'G';
    else
        if c == 2
            currC = 'R';
        end
    end

    for i = 1:length(names)
        currfig = figure;
        sgtitle([session names{i} currC])
        
        eval(['currMat =' names{i} 'mat' currC ';']);
        eval(['currMatChoice =' names{i} 'matChoice' currC ';'])
    
        subplot(4,6,1);
        m1 = currMat(os.responseInds(os.allRewards~=0),:);
        m2 = currMat(os.responseInds(os.allRewards==0),:);
        plotFilled(midPoints, m1, 'b');
        plotFilled(midPoints, m2, 'r');
        legend({'rwd', '', 'norwd', ''})
        
        subplot(4,6,2);
        m1 = currMat(os.responseInds(os.changeChoice_Inds),:);
        m2 = currMat(os.responseInds(os.stayChoice_Inds),:);
        plotFilled(midPoints, m1, 'b');
        plotFilled(midPoints, m2, 'r');
        legend({'change', '', 'stay', ''})
        
        subplot(4,6,3);
        m1 = currMat(os.CSplus,:);
        m2 = currMat(os.CSminus,:);
        plotFilled(midPoints, m1, 'b');
        plotFilled(midPoints, m2, 'r');
        legend({'plus', '', 'minus', ''})
        
        subplot(4,6,4);
        m1 = currMat(os.responseInds(os.allChoices==1),:);
        m2 = currMat(os.responseInds(os.allChoices~=1),:);
        plotFilled(midPoints, m1, 'b');
        plotFilled(midPoints, m2, 'r');
        legend({'right', '', 'left', ''})
        
        subplot(4,6,5);
        m1 = currMat(os.responseInds(intersect(find(os.allChoices==1), os.rwd_Inds)),:);
        m2 = currMat(os.responseInds(intersect(find(os.allChoices==1), os.nrwd_Inds)),:);
        plotFilled(midPoints, m1, 'b');
        plotFilled(midPoints, m2, 'r');
        legend({'Rrwd', '', 'Rnorwd', ''})
        
        subplot(4,6,6);
        m1 = currMat(os.responseInds(intersect(find(os.allChoices==-1), os.rwd_Inds)),:);
        m2 = currMat(os.responseInds(intersect(find(os.allChoices==-1), os.nrwd_Inds)),:);
        plotFilled(midPoints, m1, 'b');
        plotFilled(midPoints, m2, 'r');
        legend({'Rrwd', '', 'Rnorwd', ''})
        
        %% need model
        if strcmp(modelName, 'none')
            screen = get(0,'Screensize');
            screen(4) = screen(4) - 100;
            set(currfig, 'Position', screen)                             
        else
            eval(['currFiltered =' names{i} 'matChoiceFiltered' currC ';'])
            target = pe;
            subplot(4,3,[4,7,10]); hold on;
            [~, sortedInd] = sort(pe);
            imagesc(midPointsFilter, 1:size(currFiltered,1), currFiltered(sortedInd,:));
            hold on;
            plot([0.001*os.rwdDelay 0.001*os.rwdDelay], [0 size(currFiltered,1)],'Color', 'w', 'LineWidth', 2, 'LineStyle', '--');
            xlim([-1 2.5])
            myMap = [[linspace(0, 1, 50)', linspace(0, 1, 50)', linspace(1, 1, 50)'];
            [linspace(1, 1, 100)', linspace(1, 0, 100)', linspace(1, 0, 100)']];
            colormap(myMap)
            colorbar
            % maxScale = max(abs(waveformsSession), [], 'all');
            clim([min(currFiltered,[],"all"), max(currFiltered,[],"all")])
    
            
            currWin = focusWins{1};
            bl = mean(currMat(:, midPoints>=currWin(1)&midPoints<currWin(2)), 2);
            for j = 2:length(focusWins)
                currWin = focusWins{j};
                focusMean = mean(currMat(:, midPoints>=currWin(1)&midPoints<currWin(2)), 2) - bl;
                subplot(4,6,[11,17,23]); hold on;
                imagesc(currWin, [1 length(focusMean)], zscore(focusMean(sortedInd)));
                colormap(myMap);
                title('-bl')
                focusMean = mean(currMat(:, midPoints>=currWin(1)&midPoints<currWin(2)), 2);
                subplot(4,6,[12, 18, 24]); hold on;
                imagesc(currWin, [1 length(focusMean)], zscore(focusMean(sortedInd)));
                colormap(myMap);
            end
        
       
            
            edges = [linspace(min(target)-0.01, 0, 0.5*psthBinNum+1) linspace(0, max(target)+0.01, 0.5*psthBinNum+1)];
            edges = edges([1:0.5*psthBinNum 0.5*psthBinNum+2:end]);
            colorPSTH = [1 0 0;
                         1 0.4 0.4;
                         0.4 0.4 1;
                         0 0 1];
            subplot(4,3, [5 8]); hold on;
            for j = 1:psthBinNum
                plotFilled(midPointsFilter, currFiltered(target>=edges(j)&target<edges(j+1),:), colorPSTH(j,:));
            end
            
        
            meanSignal = zeros(1, numBins);
            semSignal = zeros(1, numBins);
            target = t.pe;
            edges = linspace(min(target)-0.01, max(target)+0.01, numBins+1);
            signal = zscore(mean(currMatChoice(:, midPoints>=rwdWin(1)&midPoints<rwdWin(2)), 2));
            signalBl = mean(currMatChoice(:, midPoints>=rwdWin(1)&midPoints<rwdWin(2)), 2);
            blWin = focusWins{1};
            bl = mean(currMatChoice(:, midPoints>=blWin(1)&midPoints<blWin(2)), 2);
            signalBl = zscore(signalBl - bl);
            for j = 1:numBins
                meanPe(j) = mean(target(target>=edges(j) & target<edges(j+1)));
                meanSignal(j) = mean(signal(target>=edges(j) & target<edges(j+1)));
                semSignal(j) = sem(signal(target>=edges(j) & target<edges(j+1)));
                meanSignalBl(j) = mean(signalBl(target>=edges(j) & target<edges(j+1)));
                semSignalBl(j) = sem(signalBl(target>=edges(j) & target<edges(j+1)));
            end
    
            subplot(4,6,21); hold on;
            plot(meanPe, meanSignal, 'LineWidth', 2, 'Color', 'k');
            patch([meanPe, flip(meanPe)], [meanSignal-semSignal flip(meanSignal+semSignal)], 'k', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
            xlabel('RPE');
            ylabel('dF/F zscored')
            set(gca, 'TickDir', 'out');
            set(gca, 'Box', 'off');
            set(gca, 'XTick', -1:0.5:1);
            set(gca, 'YTick', -0.5:0.5:1.0);
            title('raw');
        
        
            subplot(4,6,22); hold on;
            plot(meanPe, meanSignalBl, 'LineWidth', 2, 'Color', 'k');
            patch([meanPe, flip(meanPe)], [meanSignalBl-semSignalBl  flip(meanSignalBl+semSignalBl)], 'k', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
            xlabel('RPE');
            ylabel('dF/F zscored')
            set(gca, 'TickDir', 'out');
            set(gca, 'Box', 'off');
            set(gca, 'XTick', -1:0.5:1);
            set(gca, 'YTick', -0.5:0.5:1.0);
            title('-bl');      
    
        end
        screen = get(0,'Screensize');
        screen(4) = screen(4) - 100;
        set(currfig, 'Position', screen)
        savePath = pd.saveFigFolder;
        saveFigurePDF(currfig,[savePath session '_' names{i} '_' currC '.pdf']);
    end
end
end