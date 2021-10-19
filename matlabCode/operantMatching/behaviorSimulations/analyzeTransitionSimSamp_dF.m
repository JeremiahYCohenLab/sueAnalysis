function [t] = analyzeTransitionSimSamp_dF(xlFile, sheet, category, varargin)

p = inputParser;
% default parameters if none given
p.addParameter('rwdProbs', [90 50 10])
p.addParameter('tranWin', 10)
p.addParameter('maxTrials', 1000)
p.addParameter('runs', 1)
p.addParameter('samps', 2000)
p.addParameter('randomSeed', 98773)
p.addParameter('modelName', 'sevenParam_absPePeAN_scale_int_bias_ord')
p.addParameter('modelBeh', 'clean')
p.addParameter('bernFlag', 1)
p.addParameter('sessionParamsFlag', 0)
p.addParameter('lesionInd', [])
p.addParameter('lesionVal', 0)
p.parse(varargin{:});

[root, sep] = currComputer();

%set task conditions
pHigh = p.Results.rwdProbs(1);
pLow = p.Results.rwdProbs(2);
probDiffH = pHigh - p.Results.rwdProbs(3);
tranWin = p.Results.tranWin;

%extract session list from excel sheet
[~, dayList, ~] = xlsread(xlFile, sheet);
[~,col] = find(strcmp(dayList, category));
dayList = dayList(2:end, col);
endInd = find(cellfun(@isempty,dayList),1);
if ~isempty(endInd)
    dayList = dayList(1:endInd-1,:);
end

[~, transMed, transHigh, ~, ~, rwdMed, rwdHigh] = transitionAnalysis_opMD(xlFile, sheet, category, p.Results.rwdProbs, tranWin);

%initialize values for simulation
rSeed = p.Results.randomSeed;
errorCount = 0;
runs = p.Results.runs;

%set params for transition analysis
range = 15;
postTranWin = 15;
thresh = 0.4;

%initialize matrices for looking at features of transition behavior
medSlopeSim = [];
highSlopeSim = [];
zeroCrossSim = [];
threshCrossHighSim = [];
threshCrossLowSim = [];
transMedSim = [];
transHighSim = [];
transMedAvgSim = [];
transHighAvgSim = [];

%set params for exponential fit to choice prob curves
tauSim = [];
rwdHxTauSim = [];
ft = fittype('a*exp((-1/b)*x)');
sp = [0.5 7];

animals = [];
prevAnimal = [];
for currSesh = 1:length(dayList)
    sessionName = dayList{currSesh};
    [animal, ~] = strtok(sessionName, 'd'); 
    animal = animal(2:end);
    if strcmp(animal, prevAnimal) == 0 || p.Results.sessionParamsFlag
        animals = [animals {animal}];
        if strcmp(animal, prevAnimal) == 0
            fprintf('Simulating animal %s \n', animal);
        end

        if p.Results.bernFlag
            modelPath = [root animal sep animal 'sorted' sep 'stan' sep 'bernoulli' sep p.Results.modelName sep animal...
                    p.Results.modelBeh '_' p.Results.modelName '.mat'];
        else
            modelPath = [root animal sep animal 'sorted' sep 'stan' sep p.Results.modelName sep animal...
                        p.Results.modelBeh '_' p.Results.modelName '.mat'];
        end
        if p.Results.sessionParamsFlag
            [t, ~] = getStanModelParams_samps(p.Results.modelName, modelPath, p.Results.samps, 'sessionParamsFlag', 1,...
                'sessionName', sessionName, 'biasFlag', 1);
        else
            [t, ~] = getStanModelParams_samps(p.Results.modelName, modelPath, p.Results.samps);
        end
        if isempty(p.Results.lesionInd)
            params = t.params;
        else
            params = t.params;
            params(:, p.Results.lesionInd) = p.Results.lesionVal;
        end
        
        tmpHigh = []; tmpMed = []; tmpRwdHigh = []; tmpRwdMed = []; tmpChoiceHigh = []; tmpChoiceMed = [];
        for currSamp = 1:p.Results.samps
            for i = 1:runs
                if rem(i,20) == 0
                    fprintf('Running simulation %d of %d (%d total) \n', i, runs, p.Results.runs);
                end
                rSeed = rSeed + 1;
                [~,allRewards, allChoices, blockProbs, blockSwitch] = runSim_dF(p.Results.modelName, params(currSamp,:),...
                            p.Results.maxTrials, rSeed, p.Results.rwdProbs);

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
                        if tmpInd-tranWin > 0 & tmpInd+tranWin <= length(allChoices)
                            if trialProbs(tmpInd-1, 2) == pHigh & trialProbs(tmpInd, 2) == 10 & any(diff(trialProbs(tmpInd-tranWin:tmpInd, 1)) == probDiffH)
                                if (tmpInd - range - 1) > 0 & length(allChoices) >= (tmpInd + range)
                                    tmpHigh = [tmpHigh; allChoices(tmpInd:(tmpInd+postTranWin))];
                                    tmpRwdHigh = [tmpRwdHigh; allRewards(tmpInd-9:tmpInd)];
                                    tmpChoiceHigh = [tmpChoiceHigh; allChoices(tmpInd-9:tmpInd)];
                                end
                            elseif trialProbs(tmpInd-1, 2) == pLow & trialProbs(tmpInd, 2) == 10 & any(diff(trialProbs(tmpInd-tranWin:tmpInd, 1)) == probDiffH)
                                if (tmpInd - range - 1) > 0 & length(allChoices) >= (tmpInd + range)
                                    tmpMed = [tmpMed; allChoices(tmpInd:(tmpInd+postTranWin))];
                                    tmpRwdMed = [tmpRwdMed; allRewards(tmpInd-9:tmpInd)];
                                    tmpChoiceMed = [tmpChoiceMed; allChoices(tmpInd-9:tmpInd)];
                                end
                            elseif trialProbs(tmpInd-1, 1) == pHigh & trialProbs(tmpInd, 1) == 10 & any(diff(trialProbs(tmpInd-tranWin:tmpInd, 2)) == probDiffH)
                                if (tmpInd - range - 1) > 0 & length(allChoices) >= (tmpInd + range)
                                    tmpHigh = [tmpHigh; allChoices(tmpInd:(tmpInd+postTranWin))*-1];
                                    tmpRwdHigh = [tmpRwdHigh; allRewards(tmpInd-9:tmpInd)*-1];
                                    tmpChoiceHigh = [tmpChoiceHigh; allChoices(tmpInd-9:tmpInd)*-1];
                                end 
                            elseif trialProbs(tmpInd-1, 1) == pLow & trialProbs(tmpInd, 1) == 10 & any(diff(trialProbs(tmpInd-tranWin:tmpInd, 2)) == probDiffH)
                                if (tmpInd - range - 1) > 0 & length(allChoices) >= (tmpInd + range)
                                    tmpMed = [tmpMed; allChoices(tmpInd:(tmpInd+postTranWin))*-1];
                                    tmpRwdMed = [tmpRwdMed; allRewards(tmpInd-9:tmpInd)*-1];
                                    tmpChoiceMed = [tmpChoiceMed; allChoices(tmpInd-9:tmpInd)*-1];
                                end
                            end
                        end
                    end
                else
                    errorCount = errorCount + 1
                end
            end
        end
        
        if ~isempty(tmpMed) & ~isempty(tmpHigh)
            %min slope (max steepness)
            tmpMed(tmpMed==-1) = 0;
            tmpHigh(tmpHigh==-1) = 0;
            
            transMedSim = [transMedSim; tmpMed];
            transHighSim = [transHighSim; tmpHigh];
            
            %fit exponentials to transitions by rwd hx
            choiceX = sum([tmpChoiceMed; tmpChoiceHigh], 2);
            xInds = find(choiceX == 10);
            transAll = [tmpMed; tmpHigh];
            transAll = transAll(xInds, :);
            rwdHx = sum([tmpRwdMed; tmpRwdHigh], 2);
            rwdHx = rwdHx(xInds);
            sortInds{1} = find(rwdHx <= 5); sortInds{2} = find(rwdHx > 5);
            x = [0:postTranWin];
            rwdHxTauSim = [rwdHxTauSim; nan(1,2)];
            for currR = 1:2
                mdl = fit(x', mean(transAll(sortInds{currR},:))', ft, 'start', sp);
                rwdHxTauSim(end,currR) = mdl.b;
            end
            
            tmpMed = mean(tmpMed);
            tmpHigh = mean(tmpHigh);
            medSlopeSim = [medSlopeSim abs(min(diff(tmpMed)))];
            highSlopeSim = [highSlopeSim abs(min(diff(tmpHigh)))];
            
            %zero crossing point
            tmpDiff = [tmpHigh - tmpMed];
            zeroInd = find(tmpDiff(2:end) < 0, 1);                    %find point where trace crosses zero
            if isempty(zeroInd) || zeroInd > postTranWin - 1
                zeroCrossSim = [zeroCrossSim NaN];
            else
                zeroInd = zeroInd + 1;
                x = [zeroInd-1 zeroInd];                    %find line between points on sides of zero
                y = tmpDiff(x);
                c = [[1; 1]  x(:)]\y(:);                    %calculate parameter vector
                zeroCrossSim = [zeroCrossSim (0 - c(1)) / c(2)];  %x = (y-b)/ m
            end
            
            
            %thresh crossing point
            threshInd = find(tmpHigh < thresh, 1);
            if isempty(threshInd) || threshInd > postTranWin || threshInd == 1
                threshCrossHighSim = [threshCrossHighSim NaN];
            else
                x = [threshInd-1 threshInd];                                 %find line between points on sides of zero
                y = tmpHigh(x);
                c = [[1; 1]  x(:)]\y(:);                                     %calculate parameter vector
                threshCrossHighSim = [threshCrossHighSim (thresh - c(1)) / c(2)];  %x = (y-b)/ m
            end
            threshInd = find(tmpMed < thresh, 1);
            if isempty(threshInd) || threshInd > postTranWin || threshInd == 1
                threshCrossLowSim = [threshCrossLowSim NaN];
            else
                x = [threshInd-1 threshInd];                               %find line between points on sides of zero
                y = tmpMed(x);
                c = [[1; 1]  x(:)]\y(:);                                   %calculate parameter vector
                threshCrossLowSim = [threshCrossLowSim (thresh - c(1)) / c(2)];  %x = (y-b)/ m
            end
            
            %fit exponential curve
            x = [0:postTranWin];
            tauSim = [tauSim; nan(1,2)];
            medMdl = fit(x', tmpMed', ft, 'start', sp);
            tauSim(end,1) = medMdl.b;
            highMdl = fit(x', tmpHigh', ft, 'start', sp);
            tauSim(end,2) = highMdl.b;
            
            transMedAvgSim = [transMedAvgSim; tmpMed];
            transHighAvgSim = [transHighAvgSim; tmpHigh];
     
        end
    end
    prevAnimal = animal;
end

%average actual data within transition window
medAvg = mean(transMed(:,range:range+postTranWin));
highAvg = mean(transHigh(:,range:range+postTranWin));

%% plot results
figure; 
set(gcf, 'renderer', 'painters', 'position', [-1811 389 1641 469])
subplot(1,4,1); hold on;
histogram(zeroCrossSim, 15, 'FaceColor', 'c', 'Normalization', 'probability')
ylabel('probability')
xlabel('trial at 0 cross')
numNoCross = sum(isnan(zeroCrossSim))/length(zeroCrossSim);
legend([num2str(numNoCross*100) '% never cross'])
badAnimals = animals(isnan(zeroCrossSim));

tmpDiff = [highAvg - medAvg];
zeroInd = find(tmpDiff < 0, 1);
if zeroInd == 1
    zeroInd = find(tmpDiff < 0, 2);
    zeroInd = zeroInd(2);
end
x = [zeroInd-1 zeroInd];                 
y = tmpDiff(x);
c = [[1; 1]  x(:)]\y(:);
zeroCross = (0 - c(1)) / c(2);

yl = ylim;
plot([zeroCross zeroCross], [0 yl(2)], '--k', 'linewidth', 2)
set(gca, 'tickdir', 'out', 'box', 'off')


subplot(1,4,2); hold on;
colors = cool(length(medSlopeSim));
scatter(medSlopeSim, highSlopeSim, [], colors, 'filled')
scatter(abs(min(diff(medAvg))), abs(min(diff(highAvg))), 100, 'k', 'filled');
xl = xlim; yl = ylim;
plot([0 max([xl yl])], [0 max([xl yl])], '--k') 
ylim([0 max([xl yl])]);
xlim([0 max([xl yl])]);
xlabel('steepest medium slope')
ylabel('steepest high slope')
set(gca, 'tickdir', 'out')

subplot(1,4,3); hold on;
threshInd = find(highAvg < thresh, 1);
if isempty(threshInd)
    threshCrossHigh = NaN;
else
    x = [threshInd-1 threshInd];                        %find line between points on sides of zero
    y = highAvg(x);
    c = [[1; 1]  x(:)]\y(:);                            %calculate parameter vector
    threshCrossHigh = (thresh - c(1)) / c(2);    %x = (y-b)/ m
end
threshInd = find(medAvg < thresh, 1);
if isempty(threshInd)
    threshCrossLow = NaN;
else
    x = [threshInd-1 threshInd];                 	   %find line between points on sides of zero
    y = medAvg(x);
    c = [[1; 1]  x(:)]\y(:);                           %calculate parameter vector
    threshCrossLow = (thresh - c(1)) / c(2);    %x = (y-b)/ m
end

colors = cool(length(threshCrossLowSim));
scatter(threshCrossLowSim, threshCrossHighSim, [], colors, 'filled')
scatter(threshCrossLow, threshCrossHigh, 100, 'k', 'filled');
xl = xlim; yl = ylim;
plot([0 max([xl yl])], [0 max([xl yl])], '--k') 
ylim([0 max([xl yl])]);
xlim([0 max([xl yl])]);
xlabel('thresh cross medium')
ylabel('thresh cross high')
set(gca, 'tickdir', 'out')

numNoCross = sum(isnan(threshCrossLowSim) | isnan(threshCrossHighSim)) / length(threshCrossLowSim);
legTxt = [{['thresh = ' num2str(thresh)]}, {[num2str(numNoCross*100) '% never cross']}];
legend(legTxt)

subplot(1,4,4); hold on;
colors = cool(length(tauSim));
scatter(tauSim(:,1), tauSim(:,2), [], colors, 'filled')
xl = xlim; yl = ylim;
plot([0 max([xl yl])], [0 max([xl yl])], '--k') 
ylim([0 max([xl yl])]);
xlim([0 max([xl yl])]);
xlabel('medium \tauSim')
ylabel('high \tauSim')

%get tau for choice probabilities across animals' actual behavior
x = [0:postTranWin];
tau = nan(1,2); tauCI = nan(1,2);
medMdl = fit(x', medAvg', ft, 'start', sp);
tau(1) = medMdl.b;
tmp = confint(medMdl); 
tauCI(1) = medMdl.b - tmp(1,2);
highMdl = fit(x', highAvg', ft, 'start', sp);
tau(2) = highMdl.b;
tmp = confint(highMdl); 
tauCI(2) = highMdl.b - tmp(1,2);

%plot results from actual behavior
errorbar(tau(1), tau(2), tauCI(2), tauCI(2), tauCI(1), tauCI(1), 'o',...
    'color','k', 'markerfacecolor', 'k', 'markersize', 5)
set(gca, 'tickdir', 'out', 'box', 'off')

%find probability of actual beh belonging to simulated multivariate distribution
mahalDist = mahal(tau, tauSim);
tauProb = 1 - chi2cdf(mahalDist, 2);
legend(['prob actual = ' num2str(tauProb)])


%% plot taus from transitions sorted by rwd hx
figure;
set(gcf, 'renderer', 'painters', 'position', [-1440 320 885 527])
subplot(1,2,1); hold on;
for currA = 1:length(rwdHxTauSim)
    plot([1:2], rwdHxTauSim(currA, :), 'color', colors(currA,:))
end

transAll = [transMed; transHigh];
choiceX = sum(transAll(:,range-9:range), 2);
xInds = find(choiceX == 10);
transAll = transAll(xInds, :);
rwdHx = [rwdMed; rwdHigh];
rwdHx = sum(rwdHx(:,range-9:range), 2);
rwdHx = rwdHx(xInds);
sortInds{1} = find(rwdHx <= 5); sortInds{2} = find(rwdHx > 5);
x = [0:postTranWin];
for currR = 1:2
    mdl = fit(x', mean(transAll(sortInds{currR},range:range+postTranWin))', ft, 'start', sp);
    rwdHxTau(currR) = mdl.b;
    tmp = confint(mdl); 
    rwdHxTauCI(currR) = mdl.b - tmp(1,2);
end
errorbar([1:2], rwdHxTau, rwdHxTauCI, '-o', 'color', 'k', 'markerfacecolor', 'k', 'markersize', 5)

xticks([1:2])
xticklabels({'low', 'high'})
xlim([0.5 2.5])
xlabel('transition by rwd hx')
ylabel('\tau')
set(gca, 'tickdir', 'out', 'box', 'off')

%find probability of actual beh belonging to simulated multivariate distribution
mahalDist = mahal(rwdHxTau, rwdHxTauSim);
tauProb = 1 - chi2cdf(mahalDist, 2);
legend(['prob actual = ' num2str(tauProb)])

subplot(1,2,2); hold on;
scatter(rwdHxTauSim(:,1), rwdHxTauSim(:,2), [], colors, 'filled')
errorbar(rwdHxTau(1), rwdHxTau(2), rwdHxTauCI(2), rwdHxTauCI(2), rwdHxTauCI(1), rwdHxTauCI(1), 'o',...
    'color','k', 'markerfacecolor', 'k', 'markersize', 5)
xl = xlim; yl = ylim;
plot([0 max([xl yl])], [0 max([xl yl])], '--k')
ylim([0 max([xl yl])]);
xlim([0 max([xl yl])]);
xlabel('low \tau')
ylabel('high \tau')
set(gca, 'tickdir', 'out', 'box', 'off')

titleTxt = strrep([sheet ' ' category ' - ' p.Results.modelName], '_', ' ');
suptitle(titleTxt);


%% analyze decay constants

%get tau for choice probabilities average across animals' simulations
tauSimAvg = nan(1,2); tauSimAvgCI = nan(1,2);
medAvgSim = mean(transMedSim);
highAvgSim = mean(transHighSim);
medMdl = fit(x', medAvgSim', ft, 'start', sp);
tauSimAvg(1) = medMdl.b;
tmp = confint(medMdl); 
tauSimAvgCI(1) = medMdl.b - tmp(1,2);
highMdl = fit(x', highAvgSim', ft, 'start', sp);
tauSimAvg(2) = highMdl.b;
tmp = confint(highMdl); 
tauSimAvgCI(2) = highMdl.b - tmp(1,2);

%save relevant data in output structure
t.zeroCrossSim = zeroCrossSim; 
t.zeroCross = zeroCross;
t.threshCrossLowSim = threshCrossLowSim;
t.threshCrossHighSim = threshCrossHighSim;
t.threshCrossLow = threshCrossLow;
t.threshCrossHigh = threshCrossHigh;
t.tauSim = tauSim;
t.tau = tau;
t.tauCI = tauCI;
t.tauSimAvg = tauSimAvg;
t.tauSimAvgCI = tauSimAvgCI;
t.medAvg = medAvg;
t.highAvg = highAvg;
t.medAvgSim = transMedAvgSim;
t.highAvgSim = transHighAvgSim;
