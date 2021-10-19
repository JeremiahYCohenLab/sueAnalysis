function [zeroCross, zeroCross_actual, threshCrossLow, threshCrossHigh, threshCrossLow_actual, threshCrossHigh_actual] = analyzeTransitionSim_dF(xlFile, sheet, category, beh, varargin)

p = inputParser;
% default parameters if none given
p.addParameter('rwdProbs', [90 50 10])
p.addParameter('tranWin', 5)
p.addParameter('params', [])
p.addParameter('maxTrials', 350)
p.addParameter('runs', 400)
p.addParameter('randomSeed', 98773)
p.addParameter('modelName', 'sevenParam_absPePeAN_scale_int_bias_ord')
p.addParameter('bernFlag', 1)
p.addParameter('sessionFlag', 0)
p.parse(varargin{:});

[root, sep] = currComputer();

if ~isempty(p.Results.params)
    params = p.Results.params;
end

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

[~, transLow, transHigh, ~, ~] = transitionAnalysis_opMD(xlFile, sheet, category, p.Results.rwdProbs, tranWin);

%initialize values for simulation
rSeed = p.Results.randomSeed;
errorCount = 0;
runs = p.Results.runs;

%set params for transition analysis
range = 20;
postTranWin = 20;
thresh = 0.2;

%initialize matrices for looking at features of transition behavior
lowSlope = [];
highSlope = [];
zeroCross = [];
threshCrossHigh = [];
threshCrossLow = [];

animals = [];
prevAnimal = [];
for currSesh = 1:length(dayList)
    sessionName = dayList{currSesh};
    [animal, ~] = strtok(sessionName, 'd'); 
    animal = animal(2:end);
    if strcmp(animal, prevAnimal) == 0 || p.Results.sessionFlag
        animals = [animals {animal}];
        if strcmp(animal, prevAnimal) == 0
            fprintf('Simulating animal %s \n', animal);
        end
        
        if isempty(p.Results.params) & runs > 0
            if p.Results.bernFlag
                modelPath = [root animal sep animal 'sorted' sep 'stan' sep 'bernoulli' sep p.Results.modelName sep animal...
                        beh '_' p.Results.modelName '.mat'];
            else
                modelPath = [root animal sep animal 'sorted' sep 'stan' sep p.Results.modelName sep animal...
                            beh '_' p.Results.modelName '.mat'];
            end
            t = getStanModelParams_mode(p.Results.modelName, modelPath, sessionName, p.Results.sessionFlag);
            params = t.params;
        end

        tmpHigh = []; tmpLow = [];
        for i = 1:runs
            if rem(i,20) == 0
                fprintf('Running simulation %d of %d (%d total) \n', i, runs, p.Results.runs);
            end
            rSeed = rSeed + 1;
            [~,allRewards, allChoices, blockProbs, blockSwitch] = runSim_dF(p.Results.modelName, params,...
                        1000, rSeed, p.Results.rwdProbs);
            
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
                    if tmpInd-tranWin > 0 & tmpInd+tranWin <= length(allChoices)
                        if trialProbs(tmpInd-1, 2) == pHigh & trialProbs(tmpInd, 2) == 10 & any(diff(trialProbs(tmpInd-tranWin:tmpInd, 1)) == probDiffH)
                            if (tmpInd - range - 1) > 0 & length(allChoices) >= (tmpInd + range)
                                tmpHigh = [tmpHigh; allChoices(tmpInd:(tmpInd+postTranWin))];
                            end
                        elseif trialProbs(tmpInd-1, 2) == pLow & trialProbs(tmpInd, 2) == 10 & any(diff(trialProbs(tmpInd-tranWin:tmpInd, 1)) == probDiffH)
                            if (tmpInd - range - 1) > 0 & length(allChoices) >= (tmpInd + range)
                                tmpLow = [tmpLow; allChoices(tmpInd:(tmpInd+postTranWin))];
                            end
                        elseif trialProbs(tmpInd-1, 1) == pHigh & trialProbs(tmpInd, 1) == 10 & any(diff(trialProbs(tmpInd-tranWin:tmpInd, 2)) == probDiffH)
                            if (tmpInd - range - 1) > 0 & length(allChoices) >= (tmpInd + range)
                                tmpHigh = [tmpHigh; (allChoices(tmpInd:(tmpInd+postTranWin))*-1)];
                            end 
                        elseif trialProbs(tmpInd-1, 1) == pLow & trialProbs(tmpInd, 1) == 10 & any(diff(trialProbs(tmpInd-tranWin:tmpInd, 2)) == probDiffH)
                            if (tmpInd - range - 1) > 0 & length(allChoices) >= (tmpInd + range)
                                tmpLow = [tmpLow; (allChoices(tmpInd:(tmpInd+postTranWin))*-1)];
                            end
                        end
                    end
                end
            else
                errorCount = errorCount + 1
            end
        end
        
        if ~isempty(tmpLow) & ~isempty(tmpHigh)
            %min slope (max steepness)
            tmpLow(tmpLow==-1) = 0;
            tmpHigh(tmpHigh==-1) = 0;
            tmpLow = mean(tmpLow);
            tmpHigh = mean(tmpHigh);
            lowSlope = [lowSlope abs(min(diff(tmpLow)))];
            highSlope = [highSlope abs(min(diff(tmpHigh)))];
            
            %zero crossing point
            tmpDiff = [tmpHigh - tmpLow];
            zeroInd = find(tmpDiff(2:end) < 0, 1);                    %find point where trace crosses zero
            if isempty(zeroInd) || zeroInd > postTranWin - 1
                zeroCross = [zeroCross NaN];
            else
                zeroInd = zeroInd + 1;
                x = [zeroInd-1 zeroInd];                    %find line between points on sides of zero
                y = tmpDiff(x);
                c = [[1; 1]  x(:)]\y(:);                    %calculate parameter vector
                zeroCross = [zeroCross (0 - c(1)) / c(2)];  %x = (y-b)/ m
            end
            
            
            %thresh crossing point
            threshInd = find(tmpHigh < thresh, 1);
            if isempty(threshInd) || threshInd > postTranWin || threshInd == 1
                threshCrossHigh = [threshCrossHigh NaN];
            else
                x = [threshInd-1 threshInd];                                 %find line between points on sides of zero
                y = tmpHigh(x);
                c = [[1; 1]  x(:)]\y(:);                                     %calculate parameter vector
                threshCrossHigh = [threshCrossHigh (thresh - c(1)) / c(2)];  %x = (y-b)/ m
            end
            threshInd = find(tmpLow < thresh, 1);
            if isempty(threshInd) || threshInd > postTranWin || threshInd == 1
                threshCrossLow = [threshCrossLow NaN];
            else
                x = [threshInd-1 threshInd];                               %find line between points on sides of zero
                y = tmpLow(x);
                c = [[1; 1]  x(:)]\y(:);                                   %calculate parameter vector
                threshCrossLow = [threshCrossLow (thresh - c(1)) / c(2)];  %x = (y-b)/ m
            end
        end
    end
    prevAnimal = animal;
end

%average actual data within transition window
lowAvg = mean(transLow(:,range:range+postTranWin));
highAvg = mean(transHigh(:,range:range+postTranWin));

figure; 
subplot(1,3,1); hold on;
histogram(zeroCross, 15, 'FaceColor', 'c', 'Normalization', 'probability')
ylabel('probability')
xlabel('trial at 0 cross')
numNoCross = sum(isnan(zeroCross))/length(zeroCross);
legend([num2str(numNoCross*100) '% never cross'])
badAnimals = animals(isnan(zeroCross));

tmpDiff = [highAvg - lowAvg];
zeroInd = find(tmpDiff < 0, 1);
if zeroInd == 1
    zeroInd = find(tmpDiff < 0, 2);
    zeroInd = zeroInd(2);
end
x = [zeroInd-1 zeroInd];                 
y = tmpDiff(x);
c = [[1; 1]  x(:)]\y(:);
zeroCross_actual = (0 - c(1)) / c(2);

yl = ylim;
plot([zeroCross_actual zeroCross_actual], [0 yl(2)], '--k', 'linewidth', 2)
set(gca, 'tickdir', 'out', 'box', 'off')


subplot(1,3,2); hold on;
colors = cool(length(lowSlope));
scatter(lowSlope, highSlope, [], colors, 'filled')
scatter(abs(min(diff(lowAvg))), abs(min(diff(highAvg))), 100, 'k', 'filled');
xl = xlim; yl = ylim;
plot([0 max([xl yl])], [0 max([xl yl])], '--k') 
ylim([0 max([xl yl])]);
xlim([0 max([xl yl])]);
xlabel('steepest medium slope')
ylabel('steepest high slope')
set(gca, 'tickdir', 'out')

subplot(1,3,3); hold on;
threshInd = find(highAvg < thresh, 1);
if isempty(threshInd)
    threshCrossHigh_actual = NaN;
else
    x = [threshInd-1 threshInd];                        %find line between points on sides of zero
    y = highAvg(x);
    c = [[1; 1]  x(:)]\y(:);                            %calculate parameter vector
    threshCrossHigh_actual = (thresh - c(1)) / c(2);    %x = (y-b)/ m
end
threshInd = find(lowAvg < thresh, 1);
if isempty(threshInd)
    threshCrossLow_actual = NaN;
else
    x = [threshInd-1 threshInd];                 	   %find line between points on sides of zero
    y = lowAvg(x);
    c = [[1; 1]  x(:)]\y(:);                           %calculate parameter vector
    threshCrossLow_actual = (thresh - c(1)) / c(2);    %x = (y-b)/ m
end

colors = cool(length(threshCrossLow));
scatter(threshCrossLow, threshCrossHigh, [], colors, 'filled')
scatter(threshCrossLow_actual, threshCrossHigh_actual, 100, 'k', 'filled');
xl = xlim; yl = ylim;
plot([0 max([xl yl])], [0 max([xl yl])], '--k') 
ylim([0 max([xl yl])]);
xlim([0 max([xl yl])]);
xlabel('thresh cross medium')
ylabel('thresh cross high')
set(gca, 'tickdir', 'out')

numNoCross = sum(isnan(threshCrossLow) | isnan(threshCrossHigh)) / length(threshCrossLow);
legTxt = [{['thresh = ' num2str(thresh)]}, {[num2str(numNoCross*100) '% never cross']}];
legend(legTxt)

set(gcf, 'renderer', 'painters', 'position', [-1796 421 1558 463])
titleTxt = strrep([sheet ' ' category ' - ' p.Results.modelName], '_', ' ');
suptitle(titleTxt);

