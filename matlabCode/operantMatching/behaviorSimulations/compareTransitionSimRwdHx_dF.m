function compareTransitionSimRwdHx_dF(xlFile, sheet, category, varargin)

p = inputParser;
% default parameters if none given
p.addParameter('rwdProbs', [90 50 10])
p.addParameter('tranWin', 5)
p.addParameter('postTranWin', 15)
p.addParameter('tMax', 10)
p.addParameter('numBins', 5)
p.addParameter('maxTrials', 350)
p.addParameter('runs', 100) %max runs per samp
p.addParameter('samps', 1000)
p.addParameter('randomSeed', 7474124)
p.addParameter('modelName', 'sevenParam_absPePeAN_scale_int_bias_ord')
p.addParameter('mdlBeh', 'clean')
p.addParameter('bernFlag', 1)
p.addParameter('sessionParamsFlag', 0)
p.addParameter('biasFlag', 0)
p.parse(varargin{:});

[root, sep] = currComputer();

pHigh = p.Results.rwdProbs(1);
pLow = p.Results.rwdProbs(2);
probDiffH = pHigh - p.Results.rwdProbs(3);
tranWin = p.Results.tranWin;
consecChoice = p.Results.tMax;

%extract session list from excel sheet
[~, dayList, ~] = xlsread(xlFile, sheet);
[~,col] = find(strcmp(dayList, category));
dayList = dayList(2:end, col);
endInd = find(cellfun(@isempty,dayList),1);
if ~isempty(endInd)
    dayList = dayList(1:endInd-1,:);
end

[transChoiceMatx, rwdMatx, rwdKern, binE, rwdHxInds, tau, tauCI, sList] = ...
    transitionAnalysisRwdHx_opMD(xlFile, sheet, category, 'tranWin', p.Results.tranWin, 'postTranWin', p.Results.postTranWin,...
    'numBins', p.Results.numBins', 'tMax', p.Results.tMax, 'plotFlag', 0);

if p.Results.sessionParamsFlag
    numC = length(dayList);
else
    aList = cellfun(@(x) x(2:5), sList, 'uniformoutput', 0);
    numC = length(unique(aList));
end

transChoiceMatxSim = cell(1,numC);
rwdMatxSim = cell(1,numC);
changeChoiceMatxSim = cell(1,numC);
prevRwdMatxSim = cell(1,numC);
range = 15;
rSeed = p.Results.randomSeed;
errorCount = 0;

animals = [];
prevAnimal = [];
currSim = 0;

for currSesh = 1:length(dayList)
    sessionName = dayList{currSesh};
    [animal, ~] = strtok(sessionName, 'd'); 
    animal = animal(2:end);
    if strcmp(animal, prevAnimal) == 0 || p.Results.sessionParamsFlag
        animals = [animals {animal}];
        if strcmp(animal, prevAnimal) == 0 
            fprintf('Simulating animal %s \n', animal);
        end
        
        if p.Results.sessionParamsFlag
            fracSesh = sum(~cellfun(@isempty, regexp(sList, sessionName))) / length(sList);
            runs = ceil(fracSesh * p.Results.runs);
        else
            fracSesh = sum(~cellfun(@isempty, regexp(sList, animal))) / length(sList);
            runs = ceil(fracSesh * p.Results.runs);
        end
        
        if runs > 0
            currSim = currSim + 1;
            %get matrix of parameter samples
            if p.Results.bernFlag
                modelPath = [root animal sep animal 'sorted' sep 'stan' sep 'bernoulli' sep p.Results.modelName sep animal...
                        p.Results.mdlBeh '_' p.Results.modelName '.mat'];
            else
                modelPath = [root animal sep animal 'sorted' sep 'stan' sep p.Results.modelName sep animal...
                            p.Results.mdlBeh '_' p.Results.modelName '.mat'];
            end
            if p.Results.sessionParamsFlag
                [t, ~] = getStanModelParams_samps(p.Results.modelName, modelPath, p.Results.samps, 'sessionParamsFlag', 1,...
                    'sessionName', sessionName, 'biasFlag', p.Results.biasFlag);
            else
                [t, ~] = getStanModelParams_samps(p.Results.modelName, modelPath, p.Results.samps);
            end

            for currSamp = 1:p.Results.samps
                for i = 1:runs
                    rSeed = rSeed + 1;
                    [~,allRewards, allChoices, blockProbs, blockSwitch] = runSim_dF(p.Results.modelName, t.params(currSamp,:),...
                                350, rSeed, p.Results.rwdProbs);

                    if length(blockProbs) == length(blockSwitch)
                        trialProbs = nan(length(allChoices), 2);
                        for j = 2:length(blockSwitch)
                            trialProbs(blockSwitch(j-1):blockSwitch(j)-1, 1) = blockProbs(j-1, 1);
                            trialProbs(blockSwitch(j-1):blockSwitch(j)-1, 2) = blockProbs(j-1, 2);
                        end
                        trialProbs(blockSwitch(end):length(allChoices), 1) = blockProbs(end, 1);
                        trialProbs(blockSwitch(end):length(allChoices), 2) = blockProbs(end, 2);
                        allRewardsBin = abs(allRewards);
                        prevRewardsBin = [0 allRewardsBin(1:end-1)];
                        changeChoice = [0 abs(diff(allChoices)) > 0];

                        for j = 2:(length(blockSwitch) - 1)
                            tmpInd = blockSwitch(j);
                            if tmpInd-tranWin > 0 & tmpInd-consecChoice > 0 & tmpInd+tranWin <= length(allChoices)
                                if sum(allChoices(tmpInd-consecChoice+1:tmpInd)) == consecChoice & any(diff(trialProbs(tmpInd-tranWin:tmpInd, 1)) == probDiffH)
                                    if (tmpInd - range - 1) > 0 & length(allChoices) >= (tmpInd + range)

                                    end
                                elseif sum(allChoices(tmpInd-consecChoice+1:tmpInd)) == -consecChoice & any(diff(trialProbs(tmpInd-tranWin:tmpInd, 2)) == probDiffH)
                                    if (tmpInd - range - 1) > 0 & length(allChoices) >= (tmpInd + range)
                                        transChoiceMatxSim{currSim} = [transChoiceMatxSim{currSim}; allChoices((tmpInd-range+1):(tmpInd+range))*-1];
                                        changeChoiceMatxSim{currSim} = [changeChoiceMatxSim{currSim}; changeChoice((tmpInd-range+1):(tmpInd+range))];
                                        prevRwdMatxSim{currSim} = [prevRwdMatxSim{currSim}; abs(prevRewardsBin((tmpInd-range+1):(tmpInd+range)))];
                                        rwdMatxSim{currSim} = [rwdMatxSim{currSim}; prevRewardsBin((tmpInd-range+1):(tmpInd+range))];
                                    end
                                end
                            end
                        end
                    else
                        errorCount = errorCount + 1
                    end
                end
            end
        end
    end
    prevAnimal = animal;
end

%convert from cell to matrix (saves time)
transChoiceMatxSim = cell2mat(transChoiceMatxSim');
rwdMatxSim = cell2mat(rwdMatxSim');

transChoiceMatxSim(transChoiceMatxSim==-1) = 0;

%get exponentially weight rwd hist prior to transition, find indeces
rwdHxSim = sum(rwdMatxSim(:, range-consecChoice+1:range) .* repmat(rwdKern, size(rwdMatxSim, 1), 1), 2);
[~, ~, rwdHxIndsSim] = histcounts(rwdHxSim, binE); 

%fit exponentials to the choice average curves for simulated data
ft = fittype('a*exp((-1/b)*x)');
sp = [0.5 7];
x = [0:p.Results.postTranWin];
tauSim = [];
for currR = 1:p.Results.numBins
    tmp = nanmean(transChoiceMatxSim(rwdHxIndsSim==currR, range:range+p.Results.postTranWin));
    mdl = fit(x', tmp', ft, 'start', sp);
    tauSim(currR) = mdl.b;
    tmp = confint(mdl); 
    tauCIsim(currR) = mdl.b - tmp(1,2);
end

x = [-range+1:range];
colors = cool(p.Results.numBins);
figure; 
subplot(2,2,1); hold on;
for currR = 1:p.Results.numBins
    if any(rwdHxInds == currR)
        plotFilledBern(x, transChoiceMatx(rwdHxInds==currR,:), colors(currR, :));
    end
end
plot([0 0], [0 1], ':k')
plot([x(1) x(end)], [0.5 0.5], ':k')
ylim([0 1])
xlabel('trials from switch')
ylabel('Choice probability')
title('actual')
set(gca, 'tickdir', 'out')

subplot(2,2,2); hold on;
for currR = 1:p.Results.numBins
    if any(rwdHxInds == currR)
        plotFilledBern(x, transChoiceMatxSim(rwdHxIndsSim==currR,:), colors(currR, :));
    end
end
plot([0 0], [0 1], ':k')
plot([x(1) x(end)], [0.5 0.5], ':k')
ylim([0 1])
xlabel('trials from switch')
ylabel('Choice probability')
title('simulated')
set(gca, 'tickdir', 'out')


subplot(2,2,3); hold on;
x = binE(1:end-1) + diff(binE)/2;
errorbar(x, tau, tauCI, '-k', 'linewidth', 2);
xlim([binE(1) binE(end)]);
xlabel('rwd hx')
ylabel('\tau')
set(gca, 'tickdir', 'out', 'box' , 'off')

subplot(2,2,4); hold on;
errorbar(x, tauSim, tauCIsim, '-k', 'linewidth', 2);
xlim([binE(1) binE(end)]);
xlabel('rwd hx')
ylabel('\tau')
set(gca, 'tickdir', 'out', 'box' , 'off')

suptitle([sheet ' ' category ' ' strrep(p.Results.modelName, '_', ' ')])
set(gcf, 'renderer', 'painters', 'position', [-1465 97 989  868])
