function compareTransitionSimVar_opMD(xlFile, sheet, category, varargin)

p = inputParser;
% default parameters if none given
p.addParameter('rwdProbs', [90 50 10]);
p.addParameter('params', []);
p.addParameter('maxTrials', 350);
p.addParameter('runs', []);
p.addParameter('randomSeed', 982);
p.addParameter('modelName', 'fiveParam_rBeta');
p.addParameter('figFlag', 0);
p.addParameter('bernFlag', 0);
p.parse(varargin{:});

pHigh = p.Results.rwdProbs(1);
pLow = p.Results.rwdProbs(2);

if isempty(p.Results.params)
    [root, sep] = currComputer();
    if p.Results.bernFlag
        modelPath = [root sheet sep sheet 'sorted' sep 'stan' sep 'bernoulli' sep p.Results.modelName sep sheet...
                category '_' p.Results.modelName '.mat'];
    else
        modelPath = [root sheet sep sheet 'sorted' sep 'stan' sep p.Results.modelName sep sheet...
                    category '_' p.Results.modelName '.mat'];
    end
    t = generateStanModelTerms_opMD(p.Results.modelName, modelPath, [], 0);
    params = t.params;
else
    params = p.Results.params;
end


[wsls, transLow, transHigh] = transitionWsLs_opMD(xlFile, sheet, category, p.Results.rwdProbs);
close;

if isempty(p.Results.runs)
    [~, dayList, ~] = xlsread(xlFile, sheet);
    [~,col] = find(~cellfun(@isempty,strfind(dayList, category)) == 1);
    runs = length(dayList(2:end,col));
else
    runs = p.Results.runs;
end

transLowSim = [];
transHighSim = [];
changeChoiceMatxLowSim = [];
changeChoiceMatxHighSim = [];
prevRwdMatxLowSim = [];
prevRwdMatxHighSim = [];
varLowSim = [];
varHighSim = [];
range = 20;
rSeed = p.Results.randomSeed;
errorCount = 0;

for i = 1:runs
    if rem(i,20) == 0
        fprintf('Running simulation %d of %d \n', i, runs);
    end
    rSeed = rSeed + 1;
    switch p.Results.modelName
        case 'fourParam'
            [allRewards, allChoices, blockProbs, blockSwitch] = qLearningModel_simNoPlot('params', params,...
                'maxTrials', p.Results.maxTrials, 'randomSeed', rSeed, 'rwdProbs', p.Results.rwdProbs);
        case 'fiveParam_k'
            [allRewards, allChoices, blockProbs, blockSwitch] = qLearningModel_k_simNoPlot('params', params,...
                'maxTrials', p.Results.maxTrials, 'randomSeed', rSeed, 'rwdProbs', p.Results.rwdProbs);
        case 'fiveParamO'
            [~, allRewards, allChoices, blockProbs, blockSwitch] = qLearningModel_oppo_simNoPlot('params', params,...
                'maxTrials', p.Results.maxTrials, 'randomSeed', rSeed, 'rwdProbs', p.Results.rwdProbs);
        case 'fiveParam_rBeta_scale'
            [~, allRewards, allChoices, blockProbs, blockSwitch] = qLearningModel_rBeta_simNoPlot('params', params,...
                'maxTrials', p.Results.maxTrials, 'randomSeed', rSeed, 'rwdProbs', p.Results.rwdProbs);
        case 'fiveParam_rBeta_kappa'
            [~, allRewards, allChoices, blockProbs, blockSwitch] = qLearningModel_rBetaKappa_simNoPlot('params', params,...
                'maxTrials', p.Results.maxTrials, 'randomSeed', rSeed, 'rwdProbs', p.Results.rwdProbs);
        case 'sixParamO_rBarStart'
            [~, allRewards, allChoices, blockProbs, blockSwitch] = qLearningModel_oppo_rBarStart_simNoPlot('params', params,...
                'maxTrials', p.Results.maxTrials, 'randomSeed', rSeed, 'rwdProbs', p.Results.rwdProbs);
        case 'sixParam_rBeta_oppo'
            [~, allRewards, allChoices, blockProbs, blockSwitch] = qLearningModel_rBeta_oppo_simNoPlot('params', params,...
                'maxTrials', p.Results.maxTrials, 'randomSeed', rSeed, 'rwdProbs', p.Results.rwdProbs);
        case 'sixParam_peBeta'
            [~, allRewards, allChoices, blockProbs, blockSwitch] = qLearningModel_peBeta_simNoPlot('params', params,...
                'maxTrials', p.Results.maxTrials, 'randomSeed', rSeed, 'rwdProbs', p.Results.rwdProbs);
        case 'sevenParam_peBeta_k'
            [~, allRewards, allChoices, blockProbs, blockSwitch] = qLearningModel_peBeta_k_simNoPlot('params', params,...
                'maxTrials', p.Results.maxTrials, 'randomSeed', rSeed, 'rwdProbs', p.Results.rwdProbs);
        case 'sixParam_pePeBeta'
            [~, allRewards, allChoices, blockProbs, blockSwitch] = qLearningModel_pePeBeta_simNoPlot('params', params,...
                'maxTrials', p.Results.maxTrials, 'randomSeed', rSeed, 'rwdProbs', p.Results.rwdProbs);
        case 'sixParam_rAN'
            [~, allRewards, allChoices, blockProbs, blockSwitch] = qLearningModel_rAN_simNoPlot('params', params,...
                'maxTrials', p.Results.maxTrials, 'randomSeed', rSeed, 'rwdProbs', p.Results.rwdProbs);
        case 'sixParam_pePeAN'
            [~, allRewards, allChoices, blockProbs, blockSwitch] = qLearningModel_pePeAN_simNoPlot('params', params,...
                'maxTrials', p.Results.maxTrials, 'randomSeed', rSeed, 'rwdProbs', p.Results.rwdProbs);
        case 'sixParam_pePeAN_lag'
            [~, allRewards, allChoices, blockProbs, blockSwitch] = qLearningModel_pePeAN_lag_simNoPlot('params', params,...
                'maxTrials', p.Results.maxTrials, 'randomSeed', rSeed, 'rwdProbs', p.Results.rwdProbs);
        case 'sixParam_absPePeAN'
            [~, var, allRewards, allChoices, blockProbs, blockSwitch] = qLearningModel_absPePeAN_simNoPlot('params', params,...
                'maxTrials', p.Results.maxTrials, 'randomSeed', rSeed, 'rwdProbs', p.Results.rwdProbs);
        case 'sevenParam_absPePeLR'
            [~, var, allRewards, allChoices, blockProbs, blockSwitch] = qLearningModel_absPePeLR_simNoPlot('params', params,...
                'maxTrials', p.Results.maxTrials, 'randomSeed', rSeed, 'rwdProbs', p.Results.rwdProbs); 
        case 'fiveParam_absPePeAN'
            [~, allRewards, allChoices, blockProbs, blockSwitch] = qLearningModel_absPePeAN_noMin_simNoPlot('params', params,...
                'maxTrials', p.Results.maxTrials, 'randomSeed', rSeed, 'rwdProbs', p.Results.rwdProbs);
        case 'sevenParam_pePeAN_k'
            [~, allRewards, allChoices, blockProbs, blockSwitch] = qLearningModel_pePeAN_k_simNoPlot('params', params,...
                'maxTrials', p.Results.maxTrials, 'randomSeed', rSeed, 'rwdProbs', p.Results.rwdProbs);
        case 'sevenParam_peLR'
            [~, allRewards, allChoices, blockProbs, blockSwitch] = qLearningModel_peLR_simNoPlot('params', params,...
                'maxTrials', p.Results.maxTrials, 'randomSeed', rSeed, 'rwdProbs', p.Results.rwdProbs);
        case 'eightParam_rBeta_pePeAN'
            [~, allRewards, allChoices, blockProbs, blockSwitch] = qLearningModel_rBeta_pePeAN_simNoPlot('params', params,...
                'maxTrials', p.Results.maxTrials, 'randomSeed', rSeed, 'rwdProbs', p.Results.rwdProbs);
        case 'eightParam_rBeta_absPePeAN'
            [~, allRewards, allChoices, blockProbs, blockSwitch] = qLearningModel_rBeta_absPePeAN_simNoPlot('params', params,...
                'maxTrials', p.Results.maxTrials, 'randomSeed', rSeed, 'rwdProbs', p.Results.rwdProbs);
    end
    
    var = var';
    if length(blockProbs) == length(blockSwitch)
        trialProbs = nan(length(allChoices), 2);
        for j = 2:length(blockSwitch)
            trialProbs(blockSwitch(j-1):blockSwitch(j)-1, 1) = blockProbs(j-1, 1);
            trialProbs(blockSwitch(j-1):blockSwitch(j)-1, 2) = blockProbs(j-1, 2);
        end
        trialProbs(blockSwitch(end):length(allChoices), 1) = blockProbs(end, 1);
        trialProbs(blockSwitch(end):length(allChoices), 2) = blockProbs(end, 2);
        prevRewardsBin = [0 abs(allRewards(1:end-1))];
        changeChoice = [0 abs(diff(allChoices)) > 0];

        for j = 2:(length(blockSwitch) - 1)
            tmpInd = blockSwitch(j);
            if trialProbs(tmpInd-1, 2) == pHigh & trialProbs(tmpInd-1, 1) == 10 & trialProbs(tmpInd, 2) == 10 & trialProbs(tmpInd, 1) == pHigh
                    if (tmpInd - range - 1) > 0 & length(allChoices) > (tmpInd + range)
                        transHighSim = [transHighSim; allChoices((tmpInd-range+1):(tmpInd+range))];
                        changeChoiceMatxHighSim = [changeChoiceMatxHighSim; changeChoice((tmpInd-range+1):(tmpInd+range))];
                        prevRwdMatxHighSim = [prevRwdMatxHighSim; prevRewardsBin((tmpInd-range+1):(tmpInd+range))];
                        varHighSim = [varHighSim; var((tmpInd-range+1):(tmpInd+range))];
                    end
            elseif trialProbs(tmpInd-1,2) == pLow & trialProbs(tmpInd-1, 1) == 10 & trialProbs(tmpInd, 2) == 10 & trialProbs(tmpInd, 1) == pHigh
                    if (tmpInd - range - 1) > 0 & length(allChoices) > (tmpInd + range)
                        transLowSim = [transLowSim; allChoices((tmpInd-range+1):(tmpInd+range))];
                        changeChoiceMatxLowSim = [changeChoiceMatxLowSim; changeChoice((tmpInd-range+1):(tmpInd+range))];
                        prevRwdMatxLowSim = [prevRwdMatxLowSim; prevRewardsBin((tmpInd-range+1):(tmpInd+range))]; 
                        varLowSim = [varLowSim; var((tmpInd-range+1):(tmpInd+range))];
                    end
            elseif trialProbs(tmpInd-1, 1) == pHigh & trialProbs(tmpInd-1, 2) == 10 & trialProbs(tmpInd, 1) == 10 & trialProbs(tmpInd, 2) == pHigh
                    if (tmpInd - range - 1) > 0 & length(allChoices) > (tmpInd + range)
                        transHighSim = [transHighSim; (allChoices((tmpInd-range+1):(tmpInd+range))*-1)];
                        changeChoiceMatxHighSim = [changeChoiceMatxHighSim; changeChoice((tmpInd-range+1):(tmpInd+range))];
                        prevRwdMatxHighSim = [prevRwdMatxHighSim; prevRewardsBin((tmpInd-range+1):(tmpInd+range))];
                        varHighSim = [varHighSim; var((tmpInd-range+1):(tmpInd+range))];
                    end 
            elseif trialProbs(tmpInd-1, 1) == pLow & trialProbs(tmpInd-1, 2) == 10 & trialProbs(tmpInd, 1) == 10 & trialProbs(tmpInd, 2) == pHigh
                    if (tmpInd - range - 1) > 0 & length(allChoices) > (tmpInd + range)
                        transLowSim = [transLowSim; (allChoices((tmpInd-range+1):(tmpInd+range))*-1)];
                        changeChoiceMatxLowSim = [changeChoiceMatxLowSim; changeChoice((tmpInd-range+1):(tmpInd+range))];
                        prevRwdMatxLowSim = [prevRwdMatxLowSim; prevRewardsBin((tmpInd-range+1):(tmpInd+range))];
                        varLowSim = [varLowSim; var((tmpInd-range+1):(tmpInd+range))];
                    end
            end
        end
    else
        errorCount = errorCount + 1;
    end
        
end

lowAvg = mean(transLow,1);
highAvg = mean(transHigh,1);
lowAvgSim = mean(transLowSim,1);
highAvgSim = mean(transHighSim,1);


% blue = [0 1 1];
% purp = [0.7 0 1];
% colors = [linspace(blue(1),purp(1),4)', linspace(blue(2),purp(2),4)', linspace(blue(3),purp(3),4)'];
colors = [0 1 1; 0.7 0 1];


figure; 
subplot(2,4,1:2); hold on
x = [-range+1:range];
plot(x,lowAvg, '-', 'Color', colors(1,:))
plot(x,highAvg, '-', 'Color', colors(2,:))

% xx = [1:range+1];
% mdlFitLow = singleExpFitInt(xx,lowAvg(range:end));
% mdlFitHigh = singleExpFitInt(xx,highAvg(range:end));
% expConvLow = mdlFitLow.a*exp(-(mdlFitLow.b)*(1:range+1)) + mdlFitLow.c;
% expConvHigh = mdlFitHigh.a*exp(-(mdlFitHigh.b)*(1:range+1)) + mdlFitHigh.c;
% plot([0:range], expConvLow, '--', 'Color', colors(1,:))
% plot([0:range], expConvHigh,'--', 'Color', colors(2,:))

ylim([-1 1])
linetype = 'k';
vline(0, linetype)
set(gca, 'tickdir', 'out')
xlabel('Trials from switch')
ylabel('Choice average')
legend('medium -> low', 'high -> low')
title('actual')


subplot(2,4,3:4); hold on
plot(x,lowAvgSim, '-', 'Color', colors(1,:))
plot(x,highAvgSim, '-', 'Color', colors(2,:))

% mdlFitLowSim = singleExpFitInt(xx,lowAvgSim(range:end));
% mdlFitHighSim = singleExpFitInt(xx,highAvgSim(range:end));
% expConvLowSim = mdlFitLowSim.a*exp(-(mdlFitLowSim.b)*(1:range+1)) + mdlFitLowSim.c;
% expConvHighSim = mdlFitHighSim.a*exp(-(mdlFitHighSim.b)*(1:range+1)) + mdlFitHighSim.c;
% plot([0:range], expConvLowSim, '--', 'Color', colors(1,:))
% plot([0:range], expConvHighSim,'--', 'Color', colors(2,:))

ylim([-1 1])
vline(0, linetype)
set(gca, 'tickdir', 'out')
xlabel('Trials from switch')
legend('sim: medium -> low', 'sim: high -> low')
title('simulation')

% subplot(2,2,2); hold on
% taus = [1/mdlFitLow.b 1/mdlFitHigh.b 1/mdlFitLowSim.b 1/mdlFitHighSim.b];
% scatter([1:4], taus, 'k', 'filled')
% xticks([1:4])
% xticklabels({'med->low', 'high->low', 'sim: med->low', 'sim: high->low'})
% ylabel('\tau')
% xlim([0 5])
% set(gca, 'tickdir', 'out')

subplot(2,4,6:7); hold on
plot(x,[highAvg - lowAvg], '-k')
plot(x,[highAvgSim - lowAvgSim], '-', 'Color', [0.7 0.7 0.7])
plot([-range range], [0 0], ':k');
title([ 'Difference in choice probabilities'])
ylabel('Choice average difference')
xlabel('Trials from switch')
legend('high - medium', 'sim: high - medium')
set(gca, 'tickdir', 'out')
y = get(gca, 'ylim');
plot([0 0], y, ':k')

titleTxt = strrep(['Choice at block transitions: ' sheet ' ' category ' - ' p.Results.modelName], '_', ' ');
suptitle(titleTxt); 
set(gcf, 'renderer', 'painters', 'position', [-1919 1 1920 1004])

%look at win-stay lose-shift around transitions
for tInd = 1:range*2
    lS_low(tInd) = sum(changeChoiceMatxLowSim(find(prevRwdMatxLowSim(:,tInd)==0), tInd))/sum(prevRwdMatxLowSim(:,tInd)==0);
    sem_lS_low(tInd) = sem_bernoulli(sum(changeChoiceMatxLowSim(find(prevRwdMatxLowSim(:,tInd)==0), tInd)), sum(prevRwdMatxLowSim(:,tInd)==0));
    lS_high(tInd) = sum(changeChoiceMatxHighSim(find(prevRwdMatxHighSim(:,tInd)==0), tInd))/sum(prevRwdMatxHighSim(:,tInd)==0);
    sem_lS_high(tInd) = sem_bernoulli(sum(changeChoiceMatxHighSim(find(prevRwdMatxHighSim(:,tInd)==0), tInd)), sum(prevRwdMatxHighSim(:,tInd)==0));

    wS_low(tInd) = 1 - ((sum(changeChoiceMatxLowSim(find(prevRwdMatxLowSim(:,tInd)==1), tInd)))/sum(prevRwdMatxLowSim(:,tInd)==1));
    sem_wS_low(tInd) = sem_bernoulli(sum(~changeChoiceMatxLowSim(find(prevRwdMatxLowSim(:,tInd)==1), tInd)), sum(prevRwdMatxLowSim(:,tInd)==1));
    wS_high(tInd) = 1 - ((sum(changeChoiceMatxHighSim(find(prevRwdMatxHighSim(:,tInd)==1), tInd)))/sum(prevRwdMatxHighSim(:,tInd)==1));
    sem_wS_high(tInd) = sem_bernoulli(sum(~changeChoiceMatxHighSim(find(prevRwdMatxHighSim(:,tInd)==1), tInd)), sum(prevRwdMatxHighSim(:,tInd)==1));
end

figure;
subplot(2,2,1); hold on
errorbar(x, wsls.wS_low, wsls.sem_wS_low, 'Color', colors(1,:), 'linewidth', 1.3)
errorbar(x, wsls.wS_high, wsls.sem_wS_high, 'Color', colors(2,:), 'linewidth', 1.3)
ylabel('Win-stay')
xlabel('Trials from switch')
yws = ylim;
plot([0 0], yws, '--k')
set(gca, 'tickdir', 'out')
title('actual')

subplot(2,2,2); hold on
errorbar(x, wsls.lS_low, wsls.sem_lS_low, 'Color', colors(1,:), 'linewidth', 1.3)
errorbar(x, wsls.lS_high, wsls.sem_lS_low, 'Color', colors(2,:), 'linewidth', 1.3)
ylabel('Lose-shift')
xlabel('Trials from switch')
legend('medium -> low', 'high -> low')
yls = ylim;
plot([0 0], yls, '--k')
set(gca, 'tickdir', 'out')
title('actual')

subplot(2,2,3); hold on
errorbar(x, wS_low, sem_wS_low, 'Color', colors(1,:), 'linewidth', 1.3)
errorbar(x, wS_high, sem_wS_high, 'Color', colors(2,:), 'linewidth', 1.3)
ylabel('Win-stay')
xlabel('Trials from switch')
ylim([yws])
plot([0 0], yws, '--k')
set(gca, 'tickdir', 'out')
title('simulated')

subplot(2,2,4); hold on
errorbar(x, lS_low, sem_lS_low, 'Color', colors(1,:), 'linewidth', 1.3)
errorbar(x, lS_high, sem_lS_low, 'Color', colors(2,:), 'linewidth', 1.3)
ylabel('Lose-shift')
xlabel('Trials from switch')
legend('sim: medium -> low', 'sim: high -> low')
ylim([yls])
plot([0 0], yls, '--k')
set(gca, 'tickdir', 'out')
title('simulated')

titleTxt = strrep(['Win-Stay Lose-Shift: ' sheet ' ' category ' - ' p.Results.modelName], '_', ' ');
suptitle(titleTxt); 
set(gcf, 'renderer', 'painters', 'position', [-1919 1 1920 1004])


figure; hold on
plot(x,mean(varLowSim), '-', 'Color', colors(1,:))
plot(x,mean(varHighSim), '-', 'Color', colors(2,:))

ylim([-1 1])
vline(0, linetype)
set(gca, 'tickdir', 'out')
xlabel('Trials from switch')
legend('sim: medium -> low', 'sim: high -> low')
title('simulation variable')

end

