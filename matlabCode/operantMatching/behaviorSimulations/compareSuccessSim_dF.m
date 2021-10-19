function compareSuccessSim_dF(xlFile, sheet, category, varargin)

p = inputParser;
% default parameters if none given
p.addParameter('rwdProbs', [90 50 10])
p.addParameter('maxTrials', 350)
p.addParameter('randomSeed', 89465)
p.addParameter('modelNames', [{'fiveParam_bias'}, {'sevenParam_absPePeAN_scale_int_bias_ord'}])
p.addParameter('mdlBeh', 'clean')
p.addParameter('bernFlag', [])
p.addParameter('samps', 1000)
p.addParameter('runsPerSamp', 1)
p.addParameter('sessionParamsFlag', 0)
p.addParameter('compareParam', []) %model ind param ind
p.addParameter('lesionFlag', 0)
p.addParameter('varName', {'peBar', 'aN'})
p.addParameter('revForFlag', 0)
p.parse(varargin{:})

%get info about computer
[root, sep] = currComputer();

%extract parameters from given input
rSeed = p.Results.randomSeed;
numVars = length(p.Results.varName);
numMdls = length(p.Results.modelNames);
if p.Results.sessionParamsFlag
    biasFlag = 1;
else
    biasFlag = 0;
end

%get session list
[~, dayList, ~] = xlsread(xlFile, sheet);
[~,col] = find(strcmp(dayList, category));
dayList = dayList(2:end, col);
endInd = find(cellfun(@isempty,dayList),1);
if ~isempty(endInd)
    dayList = dayList(1:endInd-1,:);
end

%find number of animals
prevAnimal = [];
numA = 0;
for currS = 1:length(dayList)
    sessionName = dayList{currS};
    [animal, ~] = strtok(sessionName, 'd'); 
    animal = animal(2:end);
    if strcmp(animal, prevAnimal) == 0
        numA = numA + 1;
        prevAnimal = animal;
    end
end

%if simulating a lesion, add the lesion model to the list
if p.Results.lesionFlag
    numMdls = numMdls + 1;
    allVarsL = cell(numA, numVars);
    modelNames = [p.Results.modelNames p.Results.modelNames{p.Results.compareParam(1)}];
else
    modelNames = p.Results.modelNames;
end

%generate flags for bernoulli models
if isempty(p.Results.bernFlag)
    bernFlag = ones(1, numMdls);
else
    bernFlag = p.Results.bernFlag;
end

%initalize reward, correct, and choice counts
sumRwds = zeros(numA, 1);
sumChoices = zeros(numA, 1);
sumCorrect = zeros(numA, 1);
sumRwds_sim = zeros(numA, numMdls);
sumChoices_sim = zeros(numA, numMdls);
sumCorrect_sim = zeros(numA, numMdls);
sumRwds_c = zeros(numA, 1);
sumChoices_c = zeros(numA, 1);
sumCorrect_c = zeros(numA, 1);
sumRwds_r = zeros(numA, 1);
sumChoices_r = zeros(numA, 1);
sumCorrect_r = zeros(numA, 1);
sumRwds_v = zeros(numA, 1);
sumChoices_v = zeros(numA, 1);
sumCorrect_v = zeros(numA, 1);
allVars = cell(numA, numVars);


aInd = 0;
prevAnimal = [];
for currS = 1:length(dayList)
    sessionName = dayList{currS};
    [animal, ~] = strtok(sessionName, 'd'); 
    animal = animal(2:end);
    if strcmp(animal, prevAnimal) == 0
        fprintf('Simulating animal %s \n', animal);
        aInd = aInd+1;
    end

    [behSessionData, unCorrectedBlockSwitch, out] = loadBehavioralData([sessionName '.asc'], p.Results.revForFlag);
    s = parseBehavioralData(behSessionData, unCorrectedBlockSwitch);

    sumRwds(aInd) = sumRwds(aInd) + sum(abs(s.allRewards)); 
    sumChoices(aInd) = sumChoices(aInd) + length(s.allRewards);

    rwdProbL = [behSessionData(s.responseInds).rewardProbL]; rwdProbR = [behSessionData(s.responseInds).rewardProbR];
    sumCorrect(aInd) = sumCorrect(aInd) + sum(rwdProbL(logical(s.allChoice_L)) >= rwdProbR(logical(s.allChoice_L))) + ...
                        sum(rwdProbR(logical(s.allChoice_R)) >= rwdProbL(logical(s.allChoice_R)));
    
    if strcmp(animal, prevAnimal) == 0 || p.Results.sessionParamsFlag
        for currM = 1:numMdls
            if bernFlag(currM)
                modelPath = [root animal sep animal 'sorted' sep 'stan' sep 'bernoulli' sep modelNames{currM}...
                    sep animal p.Results.mdlBeh '_' modelNames{currM} '.mat'];
            else
                modelPath = [root animal sep animal 'sorted' sep 'stan' sep modelNames{currM}...
                    sep animal p.Results.mdlBeh '_' modelNames{currM} '.mat'];
            end
            t = getStanModelParams_samps(modelNames{currM}, modelPath, p.Results.samps, 'sessionParamsFlag', p.Results.sessionParamsFlag,...
                'sessionName', sessionName, 'biasFlag', biasFlag);
            params{currM} = t.params;
            if p.Results.lesionFlag & currM == numMdls
                params{currM}(:, p.Results.compareParam(2)) = 0;
            end
        end
        if ~isempty(p.Results.compareParam)
            paramC(aInd) = median(params{p.Results.compareParam(1)}(:,p.Results.compareParam(2)));
        end

        for currSamp = 1:p.Results.samps
            for i = 1:p.Results.runsPerSamp
                rSeed = rSeed + 1;
                for currM = 1:numMdls
                    [t, allRewards, allChoices, blockProbs, blockSwitch] = runSim_dF(modelNames{currM}, params{currM}(currSamp,:), ...
                                                                                        p.Results.maxTrials, rSeed, p.Results.rwdProbs);
                    trialProbs = nan(length(allChoices), 2);
                    for j = 2:length(blockSwitch)
                        trialProbs(blockSwitch(j-1):blockSwitch(j)-1, 1) = blockProbs(j-1, 1);
                        trialProbs(blockSwitch(j-1):blockSwitch(j)-1, 2) = blockProbs(j-1, 2);
                    end
                    trialProbs(blockSwitch(end):length(allChoices), 1) = blockProbs(end, 1);
                    trialProbs(blockSwitch(end):length(allChoices), 2) = blockProbs(end, 2);                                                               
                    sumCorrect_sim(aInd,currM) = sumCorrect_sim(aInd,currM) + sum(trialProbs(allChoices == -1, 1) >= trialProbs(allChoices == -1, 2)) + ...
                        sum(trialProbs(allChoices == 1, 2) >= trialProbs(allChoices == 1, 1));                                                                 
                    sumRwds_sim(aInd,currM) = sumRwds_sim(aInd,currM) + sum(abs(allRewards));
                    sumChoices_sim(aInd,currM) = sumChoices_sim(aInd,currM) + length(allChoices);
                    
                    for currV = 1:numVars
                        if isfield(t, p.Results.varName{currV})
                            if p.Results.lesionFlag & currM == numMdls
                                allVarsL{aInd, currV} = [allVars{aInd, currV} NaN t.(p.Results.varName{currV})'];
                            else
                                allVars{aInd, currV} = [allVars{aInd, currV} NaN t.(p.Results.varName{currV})'];
                            end
                        end
                    end
                end
            end
        end
        if strcmp(animal, prevAnimal) == 0
            prevAnimal = animal;
            for currR = 1:100
                %find reward rate of clairvoyant mouse
                [allRewards, allChoices, blockProbs, blockSwitch, correctArray] = clairvoyant_simNoPlot('maxTrials', p.Results.maxTrials,...
                'randomSeed', rSeed, 'rwdProbs', p.Results.rwdProbs);
                sumCorrect_c(aInd) = sumCorrect_c(aInd) + sum(correctArray); 
                sumRwds_c(aInd) = sumRwds_c(aInd) + sum(abs(allRewards));
                sumChoices_c(aInd) = sumChoices_c(aInd) + length(allChoices);

                %find reward rate of random mouse
                [allRewards, allChoices, ~, ~, correctArray] = random_simNoPlot('maxTrials', p.Results.maxTrials,...
                'randomSeed', rSeed, 'rwdProbs', p.Results.rwdProbs);
                sumCorrect_r(aInd) = sumCorrect_r(aInd) + sum(correctArray); 
                sumRwds_r(aInd) = sumRwds_r(aInd) + sum(abs(allRewards));
                sumChoices_r(aInd) = sumChoices_r(aInd) + length(allChoices);
                
                %find reward rate of optimized VKF
                [allRewards, allChoices, ~, ~, correctArray] = optimizedVKF_simNoPlot('maxTrials', p.Results.maxTrials,...
                'randomSeed', rSeed, 'rwdProbs', p.Results.rwdProbs);
                sumCorrect_v(aInd) = sumCorrect_v(aInd) + sum(correctArray); 
                sumRwds_v(aInd) = sumRwds_v(aInd) + sum(abs(allRewards));
                sumChoices_v(aInd) = sumChoices_v(aInd) + length(allChoices);
            end
        end
    end
end

%% find averages

correctRate = sumCorrect./sumChoices;
rewardRate = sumRwds./sumChoices;
successRate = correctRate + rewardRate;

correctRate_c = sumCorrect_c./sumChoices_c;
rewardRate_c = sumRwds_c./sumChoices_c ;
successRate_c = correctRate_c + rewardRate_c;

correctRate_r = sumCorrect_r./sumChoices_r;
rewardRate_r = sumRwds_r./sumChoices_r ;
successRate_r = correctRate_r + rewardRate_r;

correctRate_v = sumCorrect_v./sumChoices_v;
rewardRate_v = sumRwds_v./sumChoices_v ;
successRate_v = correctRate_v + rewardRate_v;

correctRate_sim = sumCorrect_sim./sumChoices_sim;
rewardRate_sim = sumRwds_sim./sumChoices_sim ;
successRate_sim = correctRate_sim + rewardRate_sim;

%% plot rates
colors = cool(numMdls + 4);
figure; 
set(gcf, 'renderer', 'painters', 'position', [-1919 41 1920 963])

subplot(1,3,1); hold on;
if length(correctRate) == 1
    scatter(0, correctRate, 2000, colors(1,:),  'filled')
    text(0, correctRate, num2str(correctRate), 'horizontalalignment', 'center');
    scatter(1, correctRate_c, 2000, colors(2,:),  'filled')
    text(1, correctRate_c, num2str(correctRate_c), 'horizontalalignment', 'center');
    scatter(2, correctRate_r, 2000, colors(3,:),  'filled')
    text(2, correctRate_r, num2str(correctRate_r), 'horizontalalignment', 'center');
    scatter(3, correctRate_v, 2000, colors(4,:),  'filled')
    text(3, correctRate_v, num2str(correctRate_v), 'horizontalalignment', 'center');
else
    errorbar(0, mean(correctRate), std(correctRate), 'Color', colors(1,:), 'linewidth', 2);
    text(0, mean(correctRate), num2str(mean(correctRate)), 'horizontalalignment', 'center');
    errorbar(1, mean(correctRate_c), std(correctRate_c), 'Color', colors(2,:), 'linewidth', 2);
    text(1, mean(correctRate_c), num2str(mean(correctRate_c)), 'horizontalalignment', 'center');
    errorbar(2, mean(correctRate_r), std(correctRate_r), 'Color', colors(3,:), 'linewidth', 2);
    text(2, mean(correctRate_r), num2str(mean(correctRate_r)), 'horizontalalignment', 'center');
    errorbar(3, mean(correctRate_v), std(correctRate_v), 'Color', colors(4,:), 'linewidth', 2);
    text(3, mean(correctRate_v), num2str(mean(correctRate_v)), 'horizontalalignment', 'center');
    
    [~, pVal, ~, st] = ttest2(correctRate, correctRate_r);
    plot([0 2], [mean(correctRate) + std(correctRate) + 0.05 mean(correctRate) + std(correctRate) + 0.05], '-k', 'linewidth', 2)
    text('position', [1 mean(correctRate) + std(correctRate) + 0.05 0], 'verticalalignment', 'bottom', 'horizontalalignment', 'center',...
        'string', sprintf(['t stat = ' num2str(st.tstat)  '\np = ' num2str(pVal)]))
end


tickLbls = {'actual', 'clairvoyant', 'random', 'vkf'};
for currM = 1:numMdls
    if length(correctRate) == 1
        scatter(currM+3, correctRate_sim(currM), 2000, colors(currM+4,:), 'filled')
        text(currM+3, correctRate_sim(currM), num2str(correctRate_sim(currM)), 'horizontalalignment', 'center');
    else
        errorbar(currM+3, mean(correctRate_sim(:,currM)), std(correctRate_sim(:,currM)), 'Color', colors(currM+4,:), 'linewidth', 2); 
        text(currM+3, mean(correctRate_sim(:,currM)), num2str(mean(correctRate_sim(:,currM))), 'horizontalalignment', 'center');
    end
    tickLbls = [tickLbls {modelNames{currM}}];
end
tickLbls = strrep(tickLbls, '_', ' ');
xlim([-0.5 numMdls+3.5])
ylabel('correct rate')
xticks([0:numMdls+3])
xticklabels(tickLbls)
xtickangle(20)
set(gca, 'tickdir', 'out')

subplot(1,3,2); hold on;
if length(rewardRate) == 1
    scatter(0, rewardRate, 2000, colors(1,:),  'filled')
    text(0, rewardRate, num2str(rewardRate), 'horizontalalignment', 'center');
    scatter(1, rewardRate_c, 2000, colors(2,:),  'filled')
    text(1, rewardRate_c, num2str(rewardRate_c), 'horizontalalignment', 'center');
    scatter(2, rewardRate_r, 2000, colors(3,:),  'filled')
    text(2, rewardRate_r, num2str(rewardRate_r), 'horizontalalignment', 'center');
    scatter(3, rewardRate_v, 2000, colors(4,:),  'filled')
    text(3, rewardRate_v, num2str(rewardRate_v), 'horizontalalignment', 'center');
else
    errorbar(0, mean(rewardRate), std(rewardRate), 'Color', colors(1,:), 'linewidth', 2);
    text(0, mean(rewardRate), num2str(mean(rewardRate)), 'horizontalalignment', 'center');
    errorbar(1, mean(rewardRate_c), std(rewardRate_c), 'Color', colors(2,:), 'linewidth', 2);
    text(1, mean(rewardRate_c), num2str(mean(rewardRate_c)), 'horizontalalignment', 'center');
    errorbar(2, mean(rewardRate_r), std(rewardRate_r), 'Color', colors(3,:), 'linewidth', 2);
    text(2, mean(rewardRate_r), num2str(mean(rewardRate_r)), 'horizontalalignment', 'center');
    errorbar(3, mean(rewardRate_v), std(rewardRate_v), 'Color', colors(4,:), 'linewidth', 2);
    text(3, mean(rewardRate_v), num2str(mean(rewardRate_v)), 'horizontalalignment', 'center');
    
    [~, pVal, ~, st] = ttest2(rewardRate, rewardRate_r);
    plot([0 2], [mean(rewardRate) + std(rewardRate) + 0.05 mean(rewardRate) + std(rewardRate) + 0.05], '-k', 'linewidth', 2)
    text('position', [1 mean(rewardRate) + std(rewardRate) + 0.05 0], 'verticalalignment', 'bottom', 'horizontalalignment', 'center',...
        'string', sprintf(['t stat = ' num2str(st.tstat)  '\np = ' num2str(pVal)]))
end

for currM = 1:numMdls
    if length(rewardRate) == 1
        scatter(currM+3, rewardRate_sim(currM), 2000, colors(currM+4,:), 'filled')
        text(currM+3, rewardRate_sim(currM), num2str(rewardRate_sim(currM)), 'horizontalalignment', 'center');
    else
        errorbar(currM+3, mean(rewardRate_sim(:,currM)), std(rewardRate_sim(:,currM)), 'Color', colors(currM+4,:), 'linewidth', 2); 
        text(currM+3, mean(rewardRate_sim(:,currM)), num2str(mean(rewardRate_sim(:,currM))), 'horizontalalignment', 'center');
    end
end
tickLbls = strrep(tickLbls, '_', ' ');
xlim([-0.5 numMdls+3.5])
ylabel('reward rate')
xticks([0:numMdls+3])
xticklabels(tickLbls)
xtickangle(20)
set(gca, 'tickdir', 'out')

subplot(1,3,3); hold on;
if length(successRate) == 1
    scatter(0, successRate, 2000, colors(1,:),  'filled')
    text(0, successRate, num2str(successRate), 'horizontalalignment', 'center');
    scatter(1, successRate_c, 2000, colors(2,:),  'filled')
    text(1, successRate_c, num2str(successRate_c), 'horizontalalignment', 'center');
    scatter(2, successRate_r, 2000, colors(3,:),  'filled')
    text(2, successRate_r, num2str(successRate_r), 'horizontalalignment', 'center');
    scatter(3, successRate_v, 2000, colors(4,:),  'filled')
    text(3, successRate_v, num2str(successRate_v), 'horizontalalignment', 'center');
else
    errorbar(0, mean(successRate), std(successRate), 'Color', colors(1,:), 'linewidth', 2);
    text(0, mean(successRate), num2str(mean(successRate)), 'horizontalalignment', 'center');
    errorbar(1, mean(successRate_c), std(successRate_c), 'Color', colors(2,:), 'linewidth', 2);
    text(1, mean(successRate_c), num2str(mean(successRate_c)), 'horizontalalignment', 'center');
    errorbar(2, mean(successRate_r), std(successRate_r), 'Color', colors(3,:), 'linewidth', 2);
    text(2, mean(successRate_r), num2str(mean(successRate_r)), 'horizontalalignment', 'center');
    errorbar(3, mean(successRate_v), std(successRate_v), 'Color', colors(4,:), 'linewidth', 2);
    text(3, mean(successRate_v), num2str(mean(successRate_v)), 'horizontalalignment', 'center');
    
    [~, pVal, ~, st] = ttest2(successRate, successRate_r);
    plot([0 2], [mean(successRate) + std(successRate) + 0.05 mean(successRate) + std(successRate) + 0.05], '-k', 'linewidth', 2)
    text('position', [1 mean(successRate) + std(successRate) + 0.05 0], 'verticalalignment', 'bottom', 'horizontalalignment', 'center',...
        'string', sprintf(['t stat = ' num2str(st.tstat)  '\np = ' num2str(pVal)]))
end

for currM = 1:numMdls
    if length(successRate) == 1
        scatter(currM+3, successRate_sim(currM), 2000, colors(currM+4,:), 'filled')
        text(currM+3, successRate_sim(currM), num2str(successRate_sim(currM)), 'horizontalalignment', 'center');
    else
        errorbar(currM+3, mean(successRate_sim(:,currM)), std(successRate_sim(:,currM)), 'Color', colors(currM+4,:), 'linewidth', 2); 
        text(currM+3, mean(successRate_sim(:,currM)), num2str(mean(successRate_sim(:,currM))), 'horizontalalignment', 'center');
    end
end
tickLbls = strrep(tickLbls, '_', ' ');
xlim([-0.5 numMdls+3.5])
ylabel('rwd rate + correct rate')
xticks([0:numMdls+3])
xticklabels(tickLbls)
xtickangle(20)
set(gca, 'tickdir', 'out')




%% for comparing param

if ~isempty(p.Results.compareParam)
    paramName = getParamNames_dF(modelNames{p.Results.compareParam(1)}, 0);
    paramName = paramName{p.Results.compareParam(2)};
    for currV = 1:numVars
        figure; 
        set(gcf, 'renderer', 'painters', 'position', [-1625 144 1279 808]);
        subplot(2,3,1); hold on;
        colors = cool(aInd);
        scatter(paramC, successRate, [], colors, 'filled')
        mdl = fitlm(paramC', successRate);
        x = linspace(min(paramC), max(paramC), 100);
        y = mdl.Coefficients.Estimate(1) + mdl.Coefficients.Estimate(2) * x;
        plot(x, y, '-k', 'linewidth', 1.5);
        xlabel([paramName ' estimate'])
        ylabel('actual rwd rate + correct rate')
        xl = xlim; yl = ylim;
        text('Position', [xl(2) yl(2)  0], 'HorizontalAlignment', 'right', 'VerticalAlignment', 'top', 'string', ...
            sprintf([' R^2 = ' num2str(mdl.Rsquared.Ordinary)]));
        text('Position', [xl(2) yl(2) - 0.02  0], 'HorizontalAlignment', 'right', 'VerticalAlignment', 'top', 'string', ...
            sprintf([' p = ' num2str(mdl.Coefficients.pValue(2))]));


        subplot(2,3,2); hold on;
        for currA = 1:aInd
            changeRate(currA) = nanmean(abs(diff(allVars{currA, currV})));
        end
        scatter(changeRate, successRate, [], colors, 'filled')
        mdl = fitlm(changeRate', successRate);
        x = linspace(min(changeRate), max(changeRate), 100);
        y = mdl.Coefficients.Estimate(1) + mdl.Coefficients.Estimate(2) * x;
        plot(x, y, '-k', 'linewidth', 1.5);
        xlabel([p.Results.varName{currV} ' mean change per trial'])
        ylabel('actual rwd rate + correct rate')
        xl = xlim; yl = ylim;
        text('Position', [xl(2) yl(2)  0], 'HorizontalAlignment', 'right', 'VerticalAlignment', 'top', 'string', ...
            sprintf([' R^2 = ' num2str(mdl.Rsquared.Ordinary)]));
        text('Position', [xl(2) yl(2) - 0.02  0], 'HorizontalAlignment', 'right', 'VerticalAlignment', 'top', 'string', ...
            sprintf([' p = ' num2str(mdl.Coefficients.pValue(2))]));
        set(gca, 'tickdir', 'out')

        subplot(2,3,3); hold on;
        scatter(changeRate, paramC, [], colors, 'filled')
        mdl = fitlm(changeRate', paramC);
        y = mdl.Coefficients.Estimate(1) + mdl.Coefficients.Estimate(2) * x;
        plot(x, y, '-k', 'linewidth', 1.5);
        xlabel([p.Results.varName{currV} ' mean change per trial'])
        ylabel([paramName ' estimate'])
        xl = xlim; yl = ylim;
        text('Position', [xl(1) yl(2)  0], 'HorizontalAlignment', 'left', 'VerticalAlignment', 'top', 'string', ...
            sprintf([' R^2 = ' num2str(mdl.Rsquared.Ordinary)]));
        text('Position', [xl(1) yl(2) - 0.02  0], 'HorizontalAlignment', 'left', 'VerticalAlignment', 'top', 'string', ...
            sprintf([' p = ' num2str(mdl.Coefficients.pValue(2))]));
        set(gca, 'tickdir', 'out')   


        subplot(2,3,4); hold on;
        colors = cool(aInd);
        scatter(paramC, successRate_sim(:,p.Results.compareParam(1)), [], colors, 'filled')
        mdl = fitlm(paramC', successRate_sim(:,p.Results.compareParam(1)));
        x = linspace(min(paramC), max(paramC), 100);
        y = mdl.Coefficients.Estimate(1) + mdl.Coefficients.Estimate(2) * x;
        plot(x, y, '-k', 'linewidth', 1.5);
        xlabel([paramName ' estimate'])
        ylabel('simulated rwd rate + correct rate')
        xl = xlim; yl = ylim;
        text('Position', [xl(2) yl(2)  0], 'HorizontalAlignment', 'right', 'VerticalAlignment', 'top', 'string', ...
            sprintf([' R^2 = ' num2str(mdl.Rsquared.Ordinary)]));
        text('Position', [xl(2) yl(2) - 0.02  0], 'HorizontalAlignment', 'right', 'VerticalAlignment', 'top', 'string', ...
            sprintf([' p = ' num2str(mdl.Coefficients.pValue(2))]));


        subplot(2,3,5); hold on;
        scatter(changeRate, successRate_sim(:,p.Results.compareParam(1)), [], colors, 'filled')
        mdl = fitlm(changeRate', successRate_sim(:,p.Results.compareParam(1)));
        x = linspace(min(changeRate), max(changeRate), 100);
        y = mdl.Coefficients.Estimate(1) + mdl.Coefficients.Estimate(2) * x;
        plot(x, y, '-k', 'linewidth', 1.5);
        xlabel([p.Results.varName{currV} ' mean change per trial'])
        ylabel('simulated rwd rate + correct rate')
        xl = xlim; yl = ylim;
        text('Position', [xl(2) yl(2)  0], 'HorizontalAlignment', 'right', 'VerticalAlignment', 'top', 'string', ...
            sprintf([' R^2 = ' num2str(mdl.Rsquared.Ordinary)]));
        text('Position', [xl(2) yl(2) - 0.02  0], 'HorizontalAlignment', 'right', 'VerticalAlignment', 'top', 'string', ...
            sprintf([' p = ' num2str(mdl.Coefficients.pValue(2))]));
        set(gca, 'tickdir', 'out')

        titleTxt = [strrep(modelNames{p.Results.compareParam(1)}, '_', ' ') ' - ' paramName ' & ' p.Results.varName{currV}];
        suptitle(titleTxt);
    end
    
end

%% plotting effect of simulated lesion

if p.Results.lesionFlag
    figure; hold on;
    set(gcf, 'renderer', 'painters', 'position', [-1625 499 559 453]);
    for currA = 1:numA
        plot([1 2], [successRate_sim(currA, p.Results.compareParam(1)) successRate_sim(currA, end)], '-', 'color', colors(currA,:))
    end
    xticks([1 2])
    xticklabels({paramName, 'lesion'})
    ylabel('simulated rwd rate + correct rate')
    xlim([0.5 2.5])
    set(gca, 'tickdir', 'out')
    title(strrep(modelNames{p.Results.compareParam(1)}, '_', ' '))
end




