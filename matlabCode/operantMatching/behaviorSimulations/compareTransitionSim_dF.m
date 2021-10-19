function [zeroCross, zeroCross_actual] = compareTransitionSim_dF(xlFile, sheet, category, beh, varargin)

p = inputParser;
% default parameters if none given
p.addParameter('rwdProbs', [90 50 10])
p.addParameter('tranWin', 5)
p.addParameter('params', [])
p.addParameter('maxTrials', 350)
p.addParameter('runs', 1000)
p.addParameter('randomSeed', 98773)
p.addParameter('modelName', 'sevenParam_absPePeAN_scale_int_bias_ord')
p.addParameter('bernFlag', 1)
p.addParameter('sessionParamsFlag', 0)
p.parse(varargin{:});

[root, sep] = currComputer();

if ~isempty(regexp(p.Results.modelName, 'dbm')) | ~isempty(regexp(p.Results.modelName, 'fbm'))
    dbmFlag = 1;
else
    dbmFlag = 0;
end

if ~isempty(p.Results.params)
    params = p.Results.params;
end

pHigh = p.Results.rwdProbs(1);
pMed = p.Results.rwdProbs(2);
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

[wsls, transMed, transHigh, sList, numTrials, rwdMed, rwdHigh] = transitionAnalysis_opMD(xlFile, sheet, category, p.Results.rwdProbs, tranWin);
    
if p.Results.sessionParamsFlag
    numC = length(dayList);
else
    aList = cellfun(@(x) x(2:5), sList, 'uniformoutput', 0);
    numC = length(unique(aList));
end
transMedSim = cell(1,numC);
transHighSim = cell(1,numC);
changeChoiceMatxMedSim = cell(1,numC);
changeChoiceMatxHighSim = cell(1,numC);
prevRwdMatxMedSim = cell(1,numC);
prevRwdMatxHighSim = cell(1,numC);
rwdMatxMedSim = cell(1,numC);
rwdMatxHighSim = cell(1,numC);
range = 15;
rSeed = p.Results.randomSeed;
errorCount = 0;

%initialize matrices for looking at features of transition behavior
zeroCross = [];
medSlope = [];
highSlope = [];

%set params for transition analysis
postTranWin = 15;

animals = [];
prevAnimal = [];
cInd = 0;
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
            cInd = cInd + 1; %only increase index if session or animal has instances of real transitions
            
            if isempty(p.Results.params)
                if dbmFlag
                   modelPath = [root animal sep animal 'sorted' sep 'bayesianModels' sep sessionName...
                                   sep p.Results.modelName sep 'mdlStruct.mat'];
                    load(modelPath);
                    params = mdlStruct.bestParams;
                else
                    if p.Results.bernFlag
                        modelPath = [root animal sep animal 'sorted' sep 'stan' sep 'bernoulli' sep p.Results.modelName sep animal...
                                beh '_' p.Results.modelName '.mat'];
                    else
                        modelPath = [root animal sep animal 'sorted' sep 'stan' sep p.Results.modelName sep animal...
                                    beh '_' p.Results.modelName '.mat'];
                    end
                    t = getStanModelParams_mode(p.Results.modelName, modelPath, sessionName, p.Results.sessionParamsFlag);
                    params = t.params;
                end
            end
        end

        for currSim = 1:runs
            rSeed = rSeed + 1;
            if dbmFlag
                [~, allRewards, allChoices, blockProbs, blockSwitch] = bayesianModelSim(p.Results.modelName,...
                    params, rSeed, 'maxTrials', p.Results.maxTrials, 'rwdProbs', ...
                    p.Results.rwdProbs);
            else
                [~, allRewards, allChoices, blockProbs, blockSwitch] = runSim_dF(p.Results.modelName, params,...
                            1000, rSeed, p.Results.rwdProbs);
            end
            
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
                                transHighSim{cInd} = [transHighSim{cInd}; allChoices((tmpInd-range+1):(tmpInd+range))];
                                changeChoiceMatxHighSim{cInd} = [changeChoiceMatxHighSim{cInd}; changeChoice((tmpInd-range+1):(tmpInd+range))];
                                prevRwdMatxHighSim{cInd} = [prevRwdMatxHighSim{cInd}; abs(prevRewards((tmpInd-range+1):(tmpInd+range)))];
                                rwdMatxHighSim{cInd} = [rwdMatxHighSim{cInd}; prevRewards((tmpInd-range+1):(tmpInd+range))];
                            end
                        elseif trialProbs(tmpInd-1, 2) == pMed & trialProbs(tmpInd, 2) == 10 & any(diff(trialProbs(tmpInd-tranWin:tmpInd, 1)) == probDiffH)
                            if (tmpInd - range - 1) > 0 & length(allChoices) >= (tmpInd + range)
                                transMedSim{cInd} = [transMedSim{cInd}; allChoices((tmpInd-range+1):(tmpInd+range))];
                                changeChoiceMatxMedSim{cInd} = [changeChoiceMatxMedSim{cInd}; changeChoice((tmpInd-range+1):(tmpInd+range))];
                                prevRwdMatxMedSim{cInd} = [prevRwdMatxMedSim{cInd}; abs(prevRewards((tmpInd-range+1):(tmpInd+range)))];
                                rwdMatxMedSim{cInd} = [rwdMatxMedSim{cInd}; prevRewards((tmpInd-range+1):(tmpInd+range))];
                            end
                        elseif trialProbs(tmpInd-1, 1) == pHigh & trialProbs(tmpInd, 1) == 10 & any(diff(trialProbs(tmpInd-tranWin:tmpInd, 2)) == probDiffH)
                            if (tmpInd - range - 1) > 0 & length(allChoices) >= (tmpInd + range)
                                transHighSim{cInd} = [transHighSim{cInd}; (allChoices((tmpInd-range+1):(tmpInd+range))*-1)];
                                changeChoiceMatxHighSim{cInd} = [changeChoiceMatxHighSim{cInd}; changeChoice((tmpInd-range+1):(tmpInd+range))];
                                prevRwdMatxHighSim{cInd} = [prevRwdMatxHighSim{cInd}; abs(prevRewards((tmpInd-range+1):(tmpInd+range)))];
                                rwdMatxHighSim{cInd} = [rwdMatxHighSim{cInd}; prevRewards((tmpInd-range+1):(tmpInd+range))*-1];
                            end 
                        elseif trialProbs(tmpInd-1, 1) == pMed & trialProbs(tmpInd, 1) == 10 & any(diff(trialProbs(tmpInd-tranWin:tmpInd, 2)) == probDiffH)
                            if (tmpInd - range - 1) > 0 & length(allChoices) >= (tmpInd + range)
                                transMedSim{cInd} = [transMedSim{cInd}; (allChoices((tmpInd-range+1):(tmpInd+range))*-1)];
                                changeChoiceMatxMedSim{cInd} = [changeChoiceMatxMedSim{cInd}; changeChoice((tmpInd-range+1):(tmpInd+range))];
                                prevRwdMatxMedSim{cInd} = [prevRwdMatxMedSim{cInd}; abs(prevRewards((tmpInd-range+1):(tmpInd+range)))];
                                rwdMatxMedSim{cInd} = [rwdMatxMedSim{cInd}; prevRewards((tmpInd-range+1):(tmpInd+range))*-1];
                            end
                        end
                    end
                end
            else
                errorCount = errorCount + 1
            end
        end
    end
    prevAnimal = animal;
end

%convert from cell to matrix (saves time)
transMedSim = cell2mat(transMedSim');
transHighSim = cell2mat(transHighSim');
changeChoiceMatxMedSim = cell2mat(changeChoiceMatxMedSim');
changeChoiceMatxHighSim = cell2mat(changeChoiceMatxHighSim');
prevRwdMatxMedSim = cell2mat(prevRwdMatxMedSim');
prevRwdMatxHighSim = cell2mat(prevRwdMatxHighSim');
rwdMatxMedSim = cell2mat(rwdMatxMedSim');
rwdMatxHighSim = cell2mat(rwdMatxHighSim');

%get choice probabilities for actual data
medAvg = mean(transMed,1);
highAvg = mean(transHigh,1);

transMedSim(transMedSim==-1) = 0;
transHighSim(transHighSim==-1) = 0;
medAvgSim = mean(transMedSim,1);
highAvgSim = mean(transHighSim,1);


% rwdMatxMedSim(rwdMatxMedSim==-1) = 0;
% rwdMatxHighSim(rwdMatxHighSim==-1) = 0;

colors = cool(4);

figure; 
subplot(4,4,1:2); hold on
x = [-range+1:range];
plotFilledBern(x, transMed, colors(1,:));
plotFilledBern(x, transHigh, colors(3,:));
plot([0 0], [0 1], ':k')
plot([x(1) x(end)], [0.5 0.5], ':k')
ylim([0 1])
set(gca, 'tickdir', 'out')
xlabel('Trials from switch')
ylabel('Choice probability')
legend('medium -> med', '', 'high -> med', '')
title('actual')


subplot(4,4,3:4); hold on
plot(x, medAvgSim, '-', 'Color', colors(2,:), 'linewidth', 1.5)
plot(x, highAvgSim, '-', 'Color', colors(4,:), 'linewidth', 1.5)
plot([0 0], [0 1], ':k')
plot([x(1) x(end)], [0.5 0.5], ':k')
ylim([0 1])
set(gca, 'tickdir', 'out')
xlabel('Trials from switch')
legend('sim: medium -> med', 'sim: high -> med')
title('simulation')


subplot(4,4,[5:6]); hold on;
tmp = [rwdMed; rwdHigh];
rwdHx = sum(tmp(:, range-9:range), 2);
medInds = find(rwdHx < 3 & rwdHx > 0);
highInds = find(rwdHx > 6);
tmp = [transMed; transHigh];
plotFilledBern(x, tmp(medInds,:), colors(1,:));
plotFilledBern(x, tmp(highInds,:), colors(3,:));
plot([0 0], [0 1], ':k')
plot([x(1) x(end)], [0.5 0.5], ':k')
set(gca, 'tickdir', 'out', 'box', 'off')
xlabel('Trials from switch')
legend('med rwd hist', '', 'high rwd hist', '')

subplot(4,4,[7:8]); hold on;
tmp = [rwdMatxMedSim; rwdMatxHighSim];
rwdHx = sum(tmp(:, range-9:range), 2);
medInds = find(rwdHx < 3 & rwdHx > 0);
highInds = find(rwdHx > 6);
tmp = [transMedSim; transHighSim];
plot(x, mean(tmp(medInds,:)), 'Color', colors(2,:), 'linewidth', 1.5)
plot(x, mean(tmp(highInds,:)), 'Color', colors(3,:), 'linewidth', 1.5)
plot([0 0], [0 1], ':k')
plot([x(1) x(end)], [0.5 0.5], ':k')
set(gca, 'tickdir', 'out', 'box', 'off')
xlabel('Trials from switch')
legend('sim: med rwd hist', 'sim: high rwd hist')


subplot(4,4,[9:10]); hold on;
colors2 = cool(3);
transMedInds = zeros(1,size(transMed, 1));
for currT = 1:size(transMed, 1)
    if sum(transMed(currT, range-9:range)) == 10
        transMedInds(currT) = 1;
    end
end
transMedX = transMed(logical(transMedInds), :);
rwdHx = rwdMed(logical(transMedInds), :);
rwdHx = sum(rwdHx(:, range-9:range), 2);
sortIndsI = find(rwdHx < 4);
sortIndsII = find(rwdHx == 5);
sortIndsIII = find(rwdHx > 6);

plotFilled(x, transMedX(sortIndsI, :), colors2(1,:));
plotFilled(x, transMedX(sortIndsII, :), colors2(2,:));
plotFilled(x, transMedX(sortIndsIII, :), colors2(3,:));
plot([0 0], [0 1], ':k')
plot([x(1) x(end)], [0.5 0.5], ':k')
ylim([0 1])
ylabel('choice average')
title('medium transitions by rwd hist')
set(gca, 'tickdir', 'out')


subplot(4,4,[11:12]); hold on;
transMedInds = zeros(1,size(transMedSim, 1));
for currT = 1:size(transMedSim, 1)
    if sum(transMedSim(currT, range-9:range)) == 10
        transMedInds(currT) = 1;
    end
end
transMedX = transMedSim(logical(transMedInds), :);
rwdHx = rwdMatxMedSim(logical(transMedInds), :);
rwdHx = sum(rwdHx(:, range-9:range), 2);
sortIndsI = find(rwdHx < 4);
sortIndsII = find(rwdHx == 5);
sortIndsIII = find(rwdHx > 6);
plotFilled(x, transMedX(sortIndsI, :), colors2(1,:));
plotFilled(x, transMedX(sortIndsII, :), colors2(2,:));
plotFilled(x, transMedX(sortIndsIII, :), colors2(3,:));
plot([0 0], [0 1], ':k')
plot([x(1) x(end)], [0.5 0.5], ':k')
ylim([0 1])
title('medium transitions by rwd hist')
set(gca, 'tickdir', 'out')


subplot(4,4,13:14); hold on
plot(x, medAvg, '-', 'Color', colors(1,:), 'linewidth', 1.5)
plot(x, highAvg, '-', 'Color',  colors(3,:), 'linewidth', 1.5)
plot(x, medAvgSim, '-', 'Color', colors(2,:), 'linewidth', 1.5)
plot(x, highAvgSim, '-', 'Color', colors(4,:), 'linewidth', 1.5)
plot([0 0], [0 1], ':k')
plot([x(1) x(end)], [0.5 0.5], ':k')
ylim([0 1])
set(gca, 'tickdir', 'out')
mdlL = fitlm([medAvgSim], medAvg');
mdlH = fitlm([highAvgSim], highAvg');
text('Position', [min(x) min(ylim)  0], 'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom', 'string', ...
        sprintf(['med R^2 = ' num2str(mdlL.Rsquared.Ordinary) '\nhigh R^2 = ' num2str(mdlH.Rsquared.Ordinary)]));
    
  
subplot(4,4,15:16); hold on
plot(x,[highAvg - medAvg], '-k')
plot(x,[highAvgSim - medAvgSim], '-', 'Color', [0.7 0.7 0.7])
plot([-range range], [0 0], ':k');
ylabel('Difference in choice probabilities')
xlabel('Trials from switch')
legend('high - medium', 'sim: high - medium')
set(gca, 'tickdir', 'out')
y = get(gca, 'ylim');
plot([0 0], y, ':k')
mdl = fitlm([highAvgSim - medAvgSim], [highAvg - medAvg]');
text('Position', [min(x) min(ylim)  0], 'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom', 'string', ...
        sprintf(['R^2 = ' num2str(mdl.Rsquared.Ordinary)]));

titleTxt = strrep(['Choice at block transitions: ' sheet ' ' category ' - ' p.Results.modelName], '_', ' ');
suptitle(titleTxt); 
set(gcf, 'renderer', 'painters', 'position', [-1492 42 1006 954])

%% look at win-stay lose-shift around transitions

%reduce choices to those that occur relative to the higher side
for rInd = 1:size(transHighSim,1)
    tmp = [];
    for tInd = 1:range*2-1
        if transHighSim(rInd, tInd:tInd+1) == [1 1] | transHighSim(rInd, tInd:tInd+1) == [1 0]
            tmp = [tmp tInd tInd+1];
        end
    end
    tmp = setdiff([1:range*2], tmp); 
    prevRwdMatxHighSim(rInd, tmp) = NaN;
    changeChoiceMatxHighSim(rInd, tmp) = NaN;
end
for rInd = 1:size(transMedSim,1)
    tmp = [];
    for tInd = 1:range*2-1
        if transMedSim(rInd, tInd:tInd+1) == [1 1] | transMedSim(rInd, tInd:tInd+1) == [1 0]
            tmp = [tmp tInd tInd+1];
        end
    end
    tmp = setdiff([1:range*2], tmp); 
    prevRwdMatxMedSim(rInd, tmp) = NaN;
    changeChoiceMatxMedSim(rInd, tmp) = NaN;
end


for tInd = 1:range*2
    lS_med(tInd) = sum(changeChoiceMatxMedSim(find(prevRwdMatxMedSim(:,tInd)==0), tInd))/sum(prevRwdMatxMedSim(:,tInd)==0);
    lS_high(tInd) = sum(changeChoiceMatxHighSim(find(prevRwdMatxHighSim(:,tInd)==0), tInd))/sum(prevRwdMatxHighSim(:,tInd)==0);

    wS_med(tInd) = 1 - ((sum(changeChoiceMatxMedSim(find(prevRwdMatxMedSim(:,tInd)==1), tInd)))/sum(prevRwdMatxMedSim(:,tInd)==1));
    wS_high(tInd) = 1 - ((sum(changeChoiceMatxHighSim(find(prevRwdMatxHighSim(:,tInd)==1), tInd)))/sum(prevRwdMatxHighSim(:,tInd)==1));
end

figure;
subplot(2,2,1); hold on
errorbar(x, wsls.wS_med, wsls.sem_wS_med, 'Color', colors(1,:), 'linewidth', 1.3)
errorbar(x, wsls.wS_high, wsls.sem_wS_high, 'Color', colors(3,:), 'linewidth', 1.3)
ylabel('Win-stay')
xlabel('Trials from switch')
yws = ylim;
plot([0 0], yws, '--k')
set(gca, 'tickdir', 'out')
title('actual')

subplot(2,2,2); hold on
errorbar(x, wsls.lS_med, wsls.sem_lS_med, 'Color', colors(1,:), 'linewidth', 1.3)
errorbar(x, wsls.lS_high, wsls.sem_lS_med, 'Color', colors(3,:), 'linewidth', 1.3)
ylabel('Lose-shift')
xlabel('Trials from switch')
legend('medium -> med', 'high -> med')
yls = ylim;
plot([0 0], yls, '--k')
set(gca, 'tickdir', 'out')
title('actual')

subplot(2,2,3); hold on
plot(x, wS_med, 'Color', colors(2,:), 'linewidth', 1.3)
plot(x, wS_high, 'Color', colors(4,:), 'linewidth', 1.3)
ylabel('Win-stay')
xlabel('Trials from switch')
ylim([yws])
plot([0 0], yws, '--k')
set(gca, 'tickdir', 'out')
title('simulated')

subplot(2,2,4); hold on
plot(x, lS_med, 'Color', colors(2,:), 'linewidth', 1.3)
plot(x, lS_high, 'Color', colors(4,:), 'linewidth', 1.3)
ylabel('Lose-shift')
xlabel('Trials from switch')
legend('sim: medium -> med', 'sim: high -> med')
ylim([yls])
plot([0 0], yls, '--k')
set(gca, 'tickdir', 'out')
title('simulated')

titleTxt = strrep(['Win-Stay Lose-Shift: ' sheet ' ' category ' - ' p.Results.modelName], '_', ' ');
suptitle(titleTxt); 
set(gcf, 'renderer', 'painters', 'position', [-1492 42 1006 954])

end

