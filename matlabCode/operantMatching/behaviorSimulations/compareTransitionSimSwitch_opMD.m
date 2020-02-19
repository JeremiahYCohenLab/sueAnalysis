function compareTransitionSimSwitch_opMD(xlFile, sheet, category, varargin)

p = inputParser;
% default parameters if none given
p.addParameter('params', []);
p.addParameter('bernFlag', 0);
p.addParameter('rwdProbs', [90 50 10]);
p.addParameter('maxTrials', 350);
p.addParameter('runs', []);
p.addParameter('randomSeed', 77582);
p.addParameter('modelName', 'sixParam_pePeAN_lag');
p.addParameter('figFlag', 0);
p.parse(varargin{:});


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

[transMatx, mdlFit] = transitionAnalysisAll_opM(xlFile, sheet, category, p.Results.rwdProbs(1));
close;
[transMatx_s, mdlFit_s] = transitionAnalysisAll_opM(xlFile, sheet, 'switch', p.Results.rwdProbs(2));
close;


if isempty(p.Results.runs)
    [~, dayList, ~] = xlsread(xlFile, sheet);
    [~,col] = find(~cellfun(@isempty,strfind(dayList, category)) == 1);
    runs = length(dayList(2:end,col));
else
    runs = p.Results.runs;
end

transMatx_sim = [];
transMatx_sim_s = [];
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
                'maxTrials', p.Results.maxTrials, 'randomSeed', rSeed, 'rwdProbs', p.Results.rwdProbs, 'taskType', 'switch');
        case 'fiveParam_k'
            [allRewards, allChoices, blockProbs, blockSwitch] = qLearningModel_k_simNoPlot('params', params,...
                'maxTrials', p.Results.maxTrials, 'randomSeed', rSeed, 'rwdProbs', p.Results.rwdProbs, 'taskType', 'switch');
        case 'fiveParamO'
            [~, allRewards, allChoices, blockProbs, blockSwitch] = qLearningModel_oppo_simNoPlot('params', params,...
                'maxTrials', p.Results.maxTrials, 'randomSeed', rSeed, 'rwdProbs', p.Results.rwdProbs, 'taskType', 'switch');
        case 'fiveParam_rBeta_scale'
            [~, allRewards, allChoices, blockProbs, blockSwitch] = qLearningModel_rBeta_simNoPlot('params', params,...
                'maxTrials', p.Results.maxTrials, 'randomSeed', rSeed, 'rwdProbs', p.Results.rwdProbs, 'taskType', 'switch');
        case 'fiveParam_rBeta_kappa'
            [~, allRewards, allChoices, blockProbs, blockSwitch] = qLearningModel_rBetaKappa_simNoPlot('params', params,...
                'maxTrials', p.Results.maxTrials, 'randomSeed', rSeed, 'rwdProbs', p.Results.rwdProbs, 'taskType', 'switch');
        case 'sixParam_rBeta_oppo'
            [~, allRewards, allChoices, blockProbs, blockSwitch] = qLearningModel_rBeta_oppo_simNoPlot('params', params,...
                'maxTrials', p.Results.maxTrials, 'randomSeed', rSeed, 'rwdProbs', p.Results.rwdProbs, 'taskType', 'switch');
        case 'sixParam_peBeta'
            [~, allRewards, allChoices, blockProbs, blockSwitch] = qLearningModel_peBeta_simNoPlot('params', params,...
                'maxTrials', p.Results.maxTrials, 'randomSeed', rSeed, 'rwdProbs', p.Results.rwdProbs, 'taskType', 'switch');
        case 'sevenParam_peBeta_k'
            [~, allRewards, allChoices, blockProbs, blockSwitch] = qLearningModel_peBeta_k_simNoPlot('params', params,...
                'maxTrials', p.Results.maxTrials, 'randomSeed', rSeed, 'rwdProbs', p.Results.rwdProbs, 'taskType', 'switch');
        case 'sixParam_pePeBeta'
            [~, allRewards, allChoices, blockProbs, blockSwitch] = qLearningModel_pePeBeta_simNoPlot('params', params,...
                'maxTrials', p.Results.maxTrials, 'randomSeed', rSeed, 'rwdProbs', p.Results.rwdProbs, 'taskType', 'switch');
        case 'sixParam_rAN'
            [~, allRewards, allChoices, blockProbs, blockSwitch] = qLearningModel_rAN_simNoPlot('params', params,...
                'maxTrials', p.Results.maxTrials, 'randomSeed', rSeed, 'rwdProbs', p.Results.rwdProbs, 'taskType', 'switch');
        case 'sixParam_peAN'
            [~, allRewards, allChoices, blockProbs, blockSwitch] = qLearningModel_peAN_simNoPlot('params', params,...
                'maxTrials', p.Results.maxTrials, 'randomSeed', rSeed, 'rwdProbs', p.Results.rwdProbs, 'taskType', 'switch');
        case 'sixParam_pePeAN'
            [~, allRewards, allChoices, blockProbs, blockSwitch] = qLearningModel_pePeAN_simNoPlot('params', params,...
                'maxTrials', p.Results.maxTrials, 'randomSeed', rSeed, 'rwdProbs', p.Results.rwdProbs, 'taskType', 'switch');
        case 'sixParam_pePeAN_lag'
            [~, allRewards, allChoices, blockProbs, blockSwitch] = qLearningModel_pePeAN_lag_simNoPlot('params', params,...
                'maxTrials', p.Results.maxTrials, 'randomSeed', rSeed, 'rwdProbs', p.Results.rwdProbs, 'taskType', 'switch');
        case 'sixParam_absPePeAN'
            [~, allRewards, allChoices, blockProbs, blockSwitch] = qLearningModel_absPePeAN_simNoPlot('params', params,...
                'maxTrials', p.Results.maxTrials, 'randomSeed', rSeed, 'rwdProbs', p.Results.rwdProbs, 'taskType', 'switch');
        case 'sevenParam_peAN_k'
            [~, allRewards, allChoices, blockProbs, blockSwitch] = qLearningModel_peAN_k_simNoPlot('params', params,...
                'maxTrials', p.Results.maxTrials, 'randomSeed', rSeed, 'rwdProbs', p.Results.rwdProbs, 'taskType', 'switch');
        case 'sevenParam_pePeAN_k'
            [~, allRewards, allChoices, blockProbs, blockSwitch] = qLearningModel_pePeAN_k_simNoPlot('params', params,...
                'maxTrials', p.Results.maxTrials, 'randomSeed', rSeed, 'rwdProbs', p.Results.rwdProbs, 'taskType', 'switch');
        case 'sevenParam_peLR'
            [~, allRewards, allChoices, blockProbs, blockSwitch] = qLearningModel_peLR_simNoPlot('params', params,...
                'maxTrials', p.Results.maxTrials, 'randomSeed', rSeed, 'rwdProbs', p.Results.rwdProbs, 'taskType', 'switch');
        case 'eightParam_rBeta_peAN'
            [~, allRewards, allChoices, blockProbs, blockSwitch] = qLearningModel_rBeta_peAN_simNoPlot('params', params,...
                'maxTrials', p.Results.maxTrials, 'randomSeed', rSeed, 'rwdProbs', p.Results.rwdProbs, 'taskType', 'switch');
        case 'eightParam_rBeta_pePeAN'
            [~, allRewards, allChoices, blockProbs, blockSwitch] = qLearningModel_rBeta_pePeAN_simNoPlot('params', params,...
                'maxTrials', p.Results.maxTrials, 'randomSeed', rSeed, 'rwdProbs', p.Results.rwdProbs, 'taskType', 'switch');
    end
    
    if length(blockProbs) == length(blockSwitch)
        trialProbs = nan(length(allChoices), 2);
        for j = 2:length(blockSwitch)
            trialProbs(blockSwitch(j-1):blockSwitch(j)-1, 1) = blockProbs(j-1, 1);
            trialProbs(blockSwitch(j-1):blockSwitch(j)-1, 2) = blockProbs(j-1, 2);
        end
        trialProbs(blockSwitch(end):length(allChoices), 1) = blockProbs(end, 1);
        trialProbs(blockSwitch(end):length(allChoices), 2) = blockProbs(end, 2);    
        
        for j = 2:(length(blockSwitch) - 1)
            tmpInd = blockSwitch(j);
            if blockProbs(j,1) == p.Results.rwdProbs(1)
                    if (tmpInd - range - 1) > 0 & length(allChoices) > (tmpInd + range)
                        transMatx_sim = [transMatx_sim; allChoices((tmpInd-range+1):(tmpInd+range))];
                    end
            elseif blockProbs(j,2) == p.Results.rwdProbs(1)
                    if (tmpInd - range - 1) > 0 & length(allChoices) > (tmpInd + range)
                        transMatx_sim = [transMatx_sim; (allChoices((tmpInd-range+1):(tmpInd+range)))*-1];
                    end
            elseif blockProbs(j,1) == p.Results.rwdProbs(2)
                    if (tmpInd - range - 1) > 0 & length(allChoices) > (tmpInd + range)
                        transMatx_sim_s = [transMatx_sim_s; allChoices((tmpInd-range+1):(tmpInd+range))];
                    end
            elseif blockProbs(j,2) == p.Results.rwdProbs(2)
                    if (tmpInd - range - 1) > 0 & length(allChoices) > (tmpInd + range)
                        transMatx_sim_s = [transMatx_sim_s; (allChoices((tmpInd-range+1):(tmpInd+range)))*-1];
                    end
            end
        end
    else
        errorCount = errorCount + 1
    end
    
        
end

choiceAvg = mean(transMatx,1);
choiceAvg_s = mean(transMatx_s,1);
choiceAvg_sim = mean(transMatx_sim,1);
choiceAvg_sim_s = mean(transMatx_sim_s,1);

blue = [0 1 1];
purp = [0.7 0 1];

figure; 
subplot(1,3,1); hold on
x = [-range+1:range];
plot(x,choiceAvg, '-', 'Color', blue)
plot(x,choiceAvg_s, '-', 'Color', purp)
title([ 'Choice at block transitions'])
ylabel('Choice average')
xlabel('Trials from switch')
legend([{'high->low'} {'med->low'}]); 
set(gca, 'tickdir', 'out')
ylim([-1 1])
linetype = 'k';
vline(0, linetype)


subplot(1,3,2); hold on
plot(x,choiceAvg_sim, '-', 'Color', blue)
plot(x,choiceAvg_sim_s, '-', 'Color', purp)
title([ 'Choice at block transitions -- sim'])
ylabel('Choice average')
xlabel('Trials from switch')
legend([{'high->low'} {'med->low'}]); 
set(gca, 'tickdir', 'out')
ylim([-1 1])
linetype = 'k';
vline(0, linetype)

subplot(1,3,3); hold on
plot(x,[choiceAvg - choiceAvg_s], '-k')
plot(x,[choiceAvg_sim - choiceAvg_sim_s], '-', 'Color', [0.5 0.5 0.5])
ylabel('Difference in choice average')
xlabel('Trials from switch')
legend([{'actual'} {'simulated'}]); 
set(gca, 'tickdir', 'out')
ylim([-1 1])
linetype = 'k';
vline(0, linetype)
plot(x, zeros(1,length(x)), ':k');

set(gcf, 'renderer', 'painters', 'position', [-1812 476 1654 390])
end

