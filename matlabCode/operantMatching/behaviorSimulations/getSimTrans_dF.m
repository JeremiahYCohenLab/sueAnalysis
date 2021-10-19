function [transMed, transHigh] = getSimTrans_dF(xlFile, sheet, category, varargin)

p = inputParser;
% default parameters if none given
p.addParameter('rwdProbs', [90 50 10])
p.addParameter('tranWin', 5)
p.addParameter('range', 15)
p.addParameter('maxTrials', 350)
p.addParameter('runs', 1) %max runs per samp
p.addParameter('samps', 1000)
p.addParameter('randomSeed', 98773)
p.addParameter('modelName', 'sevenParam_absPePeAN_scale_int_bias_ord')
p.addParameter('modelBeh', 'clean')
p.addParameter('bernFlag', 1)
p.addParameter('sessionParamsFlag', 0)
p.addParameter('biasFlag', 0)
p.parse(varargin{:});

[root, sep] = currComputer();

pHigh = p.Results.rwdProbs(1);
pLow = p.Results.rwdProbs(2);
probDiffH = pHigh - p.Results.rwdProbs(3);
tranWin = p.Results.tranWin;
rSeed = p.Results.randomSeed;
range = p.Results.range;

[wsls, transMed, transHigh, sList, numTrials, rwdMed, rwdHigh] = ...
    transitionAnalysis_opMD(xlFile, sheet, category, p.Results.rwdProbs, tranWin);

%extract session list from excel sheet
[~, dayList, ~] = xlsread(xlFile, sheet);
[~,col] = find(strcmp(dayList, category));
dayList = dayList(2:end, col);
endInd = find(cellfun(@isempty,dayList),1);
if ~isempty(endInd)
    dayList = dayList(1:endInd-1,:);
end

transMed = [];
transHigh = [];

prevAnimal = [];
for currSesh = 1:length(dayList)
    sessionName = dayList{currSesh};
    [animal, ~] = strtok(sessionName, 'd'); 
    animal = animal(2:end);
    if strcmp(animal, prevAnimal) == 0 || p.Results.sessionParamsFlag
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
            %get matrix of parameter samples
            if p.Results.bernFlag
                modelPath = [root animal sep animal 'sorted' sep 'stan' sep 'bernoulli' sep p.Results.modelName sep animal...
                        p.Results.modelBeh '_' p.Results.modelName '.mat'];
            else
                modelPath = [root animal sep animal 'sorted' sep 'stan' sep p.Results.modelName sep animal...
                            p.Results.modelBeh '_' p.Results.modelName '.mat'];
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
                                p.Results.maxTrials, rSeed, p.Results.rwdProbs);

                    if length(blockProbs) == length(blockSwitch)
                        trialProbs = nan(length(allChoices), 2);
                        for j = 2:length(blockSwitch)
                            trialProbs(blockSwitch(j-1):blockSwitch(j)-1, 1) = blockProbs(j-1, 1);
                            trialProbs(blockSwitch(j-1):blockSwitch(j)-1, 2) = blockProbs(j-1, 2);
                        end
                        trialProbs(blockSwitch(end):length(allChoices), 1) = blockProbs(end, 1);
                        trialProbs(blockSwitch(end):length(allChoices), 2) = blockProbs(end, 2);
                        prevRewards = [0 allRewards(1:end-1)];
                        changeChoice = [0 abs(diff(allChoices)) > 0];

                        for j = 2:(length(blockSwitch) - 1)
                            tmpInd = blockSwitch(j);
                            if tmpInd-tranWin > 0 & tmpInd+tranWin <= length(allChoices)
                                if trialProbs(tmpInd-1, 2) == pHigh & trialProbs(tmpInd, 2) == 10 & any(diff(trialProbs(tmpInd-tranWin:tmpInd, 1)) == probDiffH)
                                    if (tmpInd - range - 1) > 0 & length(allChoices) >= (tmpInd + range)
                                        transHigh = [transHigh; allChoices((tmpInd-range+1):(tmpInd+range))];
                                    end
                                elseif trialProbs(tmpInd-1, 2) == pLow & trialProbs(tmpInd, 2) == 10 & any(diff(trialProbs(tmpInd-tranWin:tmpInd, 1)) == probDiffH)
                                    if (tmpInd - range - 1) > 0 & length(allChoices) >= (tmpInd + range)
                                        transMed = [transMed; allChoices((tmpInd-range+1):(tmpInd+range))];
                                    end
                                elseif trialProbs(tmpInd-1, 1) == pHigh & trialProbs(tmpInd, 1) == 10 & any(diff(trialProbs(tmpInd-tranWin:tmpInd, 2)) == probDiffH)
                                    if (tmpInd - range - 1) > 0 & length(allChoices) >= (tmpInd + range)
                                        transHigh = [transHigh; (allChoices((tmpInd-range+1):(tmpInd+range))*-1)];
                                    end 
                                elseif trialProbs(tmpInd-1, 1) == pLow & trialProbs(tmpInd, 1) == 10 & any(diff(trialProbs(tmpInd-tranWin:tmpInd, 2)) == probDiffH)
                                    if (tmpInd - range - 1) > 0 & length(allChoices) >= (tmpInd + range)
                                        transMed = [transMed; (allChoices((tmpInd-range+1):(tmpInd+range))*-1)];
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    end
    prevAnimal = animal;
end

transMed(transMed == -1) = 0;
transHigh(transHigh == -1) = 0;