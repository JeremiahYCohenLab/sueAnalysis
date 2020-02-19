function compareTransitionSim_opM(xlFile, sheet, category, varargin)

p = inputParser;
% default parameters if none given
p.addParameter('probs', [90 50]);
p.addParameter('bestParams', []);
p.addParameter('maxTrials', 350);
p.addParameter('runs', [30]);
p.addParameter('randomSeed', 77582);
p.addParameter('modelType', 'rBeta');
p.addParameter('figFlag', 0);
p.parse(varargin{:});

pHigh = p.Results.probs(1);
pLow = p.Results.probs(2);

[transLow, transHigh, mdlFitLow, mdlFitHigh] = transitionAnalysisAll_opMD(xlFile, sheet, category, p.Results.probs);
close;


transLowSim = [];
transHighSim = [];
range = 20;
rSeed = p.Results.randomSeed;

for i = 1:p.Results.runs
    fprintf('Running simulation %d of %d \n', i, p.Results.runs);
    rSeed = rSeed + 1;
    switch p.Results.modelType
        case 'rBeta'
            [~, allRewards, allChoices, blockProbs, blockSwitch] = qLearningModel_rBeta_sim('bestParams', p.Results.bestParams,...
                'maxTrials', p.Results.maxTrials, 'randomSeed', rSeed, 'rwdProb', [pHigh 10], 'taskType', 'coupled');
            close;
    end
    
    trialProbs = nan(length(allChoices), 2);
    for j = 2:length(blockSwitch)
        trialProbs(blockSwitch(j-1):blockSwitch(j)-1, 1) = blockProbs(j-1, 1);
        trialProbs(blockSwitch(j-1):blockSwitch(j)-1, 2) = blockProbs(j-1, 2);
    end
    trialProbs(blockSwitch(end):length(allChoices), 1) = blockProbs(end, 1);
    trialProbs(blockSwitch(end):length(allChoices), 2) = blockProbs(end, 2);
    
    for j = 2:(length(blockSwitch) - 1)
        tmpInd = blockSwitch(j);
        if trialProbs(tmpInd-1, 2) == pHigh
                if (tmpInd - range - 1) > 0 & length(allChoices) > (tmpInd + range)
                    transHighSim = [transHighSim; allChoices((tmpInd-range+1):(tmpInd+range))];
                end
        elseif trialProbs(tmpInd-1, 1) == pHigh
                if (tmpInd - range - 1) > 0 & length(allChoices) > (tmpInd + range)
                    transHighSim = [transHighSim; (allChoices((tmpInd-range+1):(tmpInd+range))*-1)];
                end 
        end
    end
    
   switch p.Results.modelType
        case 'rBeta'
            [~, allRewards, allChoices, blockProbs, blockSwitch] = qLearningModel_rBeta_sim('bestParams', p.Results.bestParams,...
                'maxTrials', p.Results.maxTrials, 'randomSeed', rSeed, 'taskType', 'switch');
            close;
    end
    
    trialProbs = nan(length(allChoices), 2);
    for j = 2:length(blockSwitch)
        trialProbs(blockSwitch(j-1):blockSwitch(j)-1, 1) = blockProbs(j-1, 1);
        trialProbs(blockSwitch(j-1):blockSwitch(j)-1, 2) = blockProbs(j-1, 2);
    end
    trialProbs(blockSwitch(end):length(allChoices), 1) = blockProbs(end, 1);
    trialProbs(blockSwitch(end):length(allChoices), 2) = blockProbs(end, 2);
    
    for j = 2:(length(blockSwitch) - 1)
        tmpInd = blockSwitch(j);
        if trialProbs(tmpInd-1, 2) == pLow
                if (tmpInd - range - 1) > 0 & length(allChoices) > (tmpInd + range)
                    transLowSim = [transLowSim; allChoices((tmpInd-range+1):(tmpInd+range))];
                end
        elseif trialProbs(tmpInd-1, 1) == pLow
                if (tmpInd - range - 1) > 0 & length(allChoices) > (tmpInd + range)
                    transLowSim = [transLowSim; (allChoices((tmpInd-range+1):(tmpInd+range))*-1)];
                end 
        end
    end 
end

lowAvg = mean(transLow,1);
highAvg = mean(transHigh,1);
lowAvgSim = mean(transLowSim,1);
highAvgSim = mean(transHighSim,1);


figure; hold on
x = [-range+1:range];
plot(x,lowAvg, '-r')
plot(x,lowAvgSim, '--r')
plot(x,highAvg, '-b')
plot(x,highAvgSim, '--b')
title([ 'Choice at block transitions'])
ylabel('Choice average')
xlabel('Trials from switch')
legend('medium -> low, actual', 'medium -> low, sim', 'high -> low, actual', 'high -> low, sim')
linetype = 'k';
vline(0, linetype)
set(gca, 'tickdir', 'out')

end

