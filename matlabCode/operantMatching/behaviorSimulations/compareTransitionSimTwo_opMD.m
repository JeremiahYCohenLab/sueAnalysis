function compareTransitionSimTwo_opMD(xlFile, sheet, category, varargin)

p = inputParser;
% default parameters if none given
p.addParameter('rwdProbs', [50 10]);
p.addParameter('bestParams', []);
p.addParameter('maxTrials', 350);
p.addParameter('runs', []);
p.addParameter('randomSeed', 77582);
p.addParameter('modelType', [{'rBeta'} {'fourParam'}]);
p.addParameter('figFlag', 0);
p.parse(varargin{:});

pHigh = p.Results.rwdProbs(1);
pLow = p.Results.rwdProbs(2);

[transLow, transHigh, mdlFitLow, mdlFitHigh] = transitionAnalysisAll_opMD(xlFile, sheet, category, p.Results.rwdProbs);
close;

if isempty(p.Results.runs)
    [~, dayList, ~] = xlsread(xlFile, sheet);
    [~,col] = find(~cellfun(@isempty,strfind(dayList, category)) == 1);
    runs = length(dayList(2:end,col));
else
    runs = p.Results.runs;
end

transLowSim_one = [];
transLowSim_two = [];
transHighSim_one = [];
transHighSim_two = [];
range = 20;
rSeed = p.Results.randomSeed;
errorCount = 0;

for i = 1:runs
    fprintf('Running simulation %d of %d \n', i, runs);
    rSeed = rSeed + 1;
    switch p.Results.modelType{1}
        case 'fourParam'
            [allRewards, allChoices, blockProbs, blockSwitch] = qLearningModel_simNoPlot('bestParams', p.Results.bestParams{1},...
                'maxTrials', p.Results.maxTrials, 'randomSeed', rSeed, 'rwdProbs', p.Results.rwdProbs);
        case 'rBeta'
            [~, allRewards, allChoices, blockProbs, blockSwitch] = qLearningModel_rBeta_simNoPlot('bestParams', p.Results.bestParams{1},...
                'maxTrials', p.Results.maxTrials, 'randomSeed', rSeed, 'rwdProbs', p.Results.rwdProbs);
        case 'rBeta_kappa'
            [~, allRewards, allChoices, blockProbs, blockSwitch] = qLearningModel_rBetaKappa_simNoPlot('bestParams', p.Results.bestParams{1},...
                'maxTrials', p.Results.maxTrials, 'randomSeed', rSeed, 'rwdProbs', p.Results.rwdProbs);
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
            if trialProbs(tmpInd-1, 2) == pHigh & trialProbs(tmpInd-1, 1) == 10 & trialProbs(tmpInd, 2) == 10 & trialProbs(tmpInd, 1) == pHigh
                    if (tmpInd - range - 1) > 0 & length(allChoices) > (tmpInd + range)
                        transHighSim_one = [transHighSim_one; allChoices((tmpInd-range+1):(tmpInd+range))];
                    end
            elseif trialProbs(tmpInd-1,2) == pLow & trialProbs(tmpInd-1, 1) == 10 & trialProbs(tmpInd, 2) == 10 & trialProbs(tmpInd, 1) == pHigh
                    if (tmpInd - range - 1) > 0 & length(allChoices) > (tmpInd + range)
                        transLowSim_one = [transLowSim_one; allChoices((tmpInd-range+1):(tmpInd+range))];
                    end
            elseif trialProbs(tmpInd-1, 1) == pHigh & trialProbs(tmpInd-1, 2) == 10 & trialProbs(tmpInd, 1) == 10 & trialProbs(tmpInd, 2) == pHigh
                    if (tmpInd - range - 1) > 0 & length(allChoices) > (tmpInd + range)
                        transHighSim_one = [transHighSim_one; (allChoices((tmpInd-range+1):(tmpInd+range))*-1)];
                    end 
            elseif trialProbs(tmpInd-1, 1) == pLow & trialProbs(tmpInd-1, 2) == 10 & trialProbs(tmpInd, 1) == 10 & trialProbs(tmpInd, 2) == pHigh
                    if (tmpInd - range - 1) > 0 & length(allChoices) > (tmpInd + range)
                        transLowSim_one = [transLowSim_one; (allChoices((tmpInd-range+1):(tmpInd+range))*-1)];
                    end
            end
        end
    else
        errorCount = errorCount + 1;
    end
    
    switch p.Results.modelType{2}
        case 'fourParam'
            [allRewards, allChoices, blockProbs, blockSwitch] = qLearningModel_simNoPlot('bestParams', p.Results.bestParams{2},...
                'maxTrials', p.Results.maxTrials, 'randomSeed', rSeed, 'rwdProbs', p.Results.rwdProbs);
        case 'rBeta'
            [~, allRewards, allChoices, blockProbs, blockSwitch] = qLearningModel_rBeta_simNoPlot('bestParams', p.Results.bestParams{2},...
                'maxTrials', p.Results.maxTrials, 'randomSeed', rSeed, 'rwdProbs', p.Results.rwdProbs);
        case 'rBeta_kappa'
            [~, allRewards, allChoices, blockProbs, blockSwitch] = qLearningModel_rBetaKappa_simNoPlot('bestParams', p.Results.bestParams{2},...
                'maxTrials', p.Results.maxTrials, 'randomSeed', rSeed, 'rwdProbs', p.Results.rwdProbs);
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
            if trialProbs(tmpInd-1, 2) == pHigh & trialProbs(tmpInd-1, 1) == 10 & trialProbs(tmpInd, 2) == 10 & trialProbs(tmpInd, 1) == pHigh
                    if (tmpInd - range - 1) > 0 & length(allChoices) > (tmpInd + range)
                        transHighSim_two = [transHighSim_two; allChoices((tmpInd-range+1):(tmpInd+range))];
                    end
            elseif trialProbs(tmpInd-1,2) == pLow & trialProbs(tmpInd-1, 1) == 10 & trialProbs(tmpInd, 2) == 10 & trialProbs(tmpInd, 1) == pHigh
                    if (tmpInd - range - 1) > 0 & length(allChoices) > (tmpInd + range)
                        transLowSim_two = [transLowSim_two; allChoices((tmpInd-range+1):(tmpInd+range))];
                    end
            elseif trialProbs(tmpInd-1, 1) == pHigh & trialProbs(tmpInd-1, 2) == 10 & trialProbs(tmpInd, 1) == 10 & trialProbs(tmpInd, 2) == pHigh
                    if (tmpInd - range - 1) > 0 & length(allChoices) > (tmpInd + range)
                        transHighSim_two = [transHighSim_two; (allChoices((tmpInd-range+1):(tmpInd+range))*-1)];
                    end 
            elseif trialProbs(tmpInd-1, 1) == pLow & trialProbs(tmpInd-1, 2) == 10 & trialProbs(tmpInd, 1) == 10 & trialProbs(tmpInd, 2) == pHigh
                    if (tmpInd - range - 1) > 0 & length(allChoices) > (tmpInd + range)
                        transLowSim_two = [transLowSim_two; (allChoices((tmpInd-range+1):(tmpInd+range))*-1)];
                    end
            end
        end
    else
        errorCount = errorCount + 1;
    end
    
        
end

lowAvg = mean(transLow,1);
highAvg = mean(transHigh,1);
lowAvgSim_one = mean(transLowSim_one,1);
highAvgSim_one = mean(transHighSim_one,1);
lowAvgSim_two = mean(transLowSim_one,1);
highAvgSim_two = mean(transHighSim_one,1);


blue = [0 1 1];
purp = [0.7 0 1];
colors = [linspace(blue(1),purp(1),6)', linspace(blue(2),purp(2),6)', linspace(blue(3),purp(3),6)'];

figure; 
subplot(1,2,1); hold on
x = [-range+1:range];
plot(x,lowAvg, '-', 'Color', colors(1,:))
plot(x,lowAvgSim_one, '-', 'Color', colors(2,:))
plot(x,lowAvgSim_two, '-', 'Color', colors(3,:))
plot(x,highAvg, '-', 'Color', colors(4,:))
plot(x,highAvgSim_one, '-', 'Color', colors(5,:))
plot(x,highAvgSim_two, '-', 'Color', colors(6,:))
title([ 'Choice at block transitions'])
ylabel('Choice average')
xlabel('Trials from switch')
legend('medium -> low', 'sim1: m -> l', 'sim2: m -> l', 'high -> low', 'sim1: h -> l', 'sim2: h -> l')


xx = [1:range+1];
mdlFitLowSim_one = singleExpFitInt(xx,lowAvgSim_one(range:end));
mdlFitLowSim_two = singleExpFitInt(xx,lowAvgSim_two(range:end));
mdlFitHighSim_one = singleExpFitInt(xx,highAvgSim_one(range:end));
mdlFitHighSim_two = singleExpFitInt(xx,highAvgSim_two(range:end));
expConvLow = mdlFitLow.a*exp(-(mdlFitLow.b)*(1:range+1)) + mdlFitLow.c;
expConvLowSim_one = mdlFitLowSim_one.a*exp(-(mdlFitLowSim_one.b)*(1:range+1)) + mdlFitLowSim_one.c;
expConvLowSim_two = mdlFitLowSim_two.a*exp(-(mdlFitLowSim_two.b)*(1:range+1)) + mdlFitLowSim_two.c;
expConvHigh = mdlFitHigh.a*exp(-(mdlFitHigh.b)*(1:range+1)) + mdlFitHigh.c;
expConvHighSim_one = mdlFitHighSim_one.a*exp(-(mdlFitHighSim_one.b)*(1:range+1)) + mdlFitHighSim_one.c;
expConvHighSim_two = mdlFitHighSim_two.a*exp(-(mdlFitHighSim_two.b)*(1:range+1)) + mdlFitHighSim_two.c;
plot([0:range], expConvLow, '--', 'Color', colors(1,:))
plot([0:range], expConvLowSim_one, '--', 'Color', colors(2,:))
plot([0:range], expConvLowSim_two, '--', 'Color', colors(3,:))
plot([0:range], expConvHigh,'--', 'Color', colors(4,:))
plot([0:range], expConvHighSim_one,'--', 'Color', colors(5,:))
plot([0:range], expConvHighSim_two,'--', 'Color', colors(6,:))

linetype = 'k';
vline(0, linetype)
ylim([-1 1])
set(gca, 'tickdir', 'out')

subplot(1,2,2); hold on
taus = [1/mdlFitLow.b 1/mdlFitHigh.b 1/mdlFitLowSim_one.b 1/mdlFitHighSim_one.b 1/mdlFitLowSim_two.b 1/mdlFitHighSim_two.b];
scatter([1:6], taus, 'k', 'filled')
xticks([1:6])
xticklabels({'med->low', 'high->low', 'sim1:\n med->low', 'sim1:\n high->low', 'sim2:\n med->low', 'sim2:\n high->low'})
ylabel('\tau')
xlim([0 7])
set(gca, 'tickdir', 'out')

end

