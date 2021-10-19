function [zeroCross, zeroCross_actual] = compareTransitionSimSamp_dF(xlFile, sheet, category, varargin)

p = inputParser;
% default parameters if none given
p.addParameter('rwdProbs', [90 50 10])
p.addParameter('tranWin', 10)
p.addParameter('maxTrials', 350)
p.addParameter('runs', 100) %max runs per samp
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
changeChoiceMatxLowSim = cell(1,numC);
changeChoiceMatxHighSim = cell(1,numC);
prevRwdMatxLowSim = cell(1,numC);
prevRwdMatxHighSim = cell(1,numC);
rwdMatxMedSim = cell(1,numC);
rwdMatxHighSim = cell(1,numC);
range = 15;
rSeed = p.Results.randomSeed;
errorCount = 0;

%initialize matrices for looking at features of transition behavior
zeroCross = [];
lowSlope = [];
highSlope = [];

%set params for transition analysis
postTranWin = 15;
trialNum = repmat([1:range+1], 1, numC*2);
transType = repmat([ones(1,range+1) ones(1,range+1)*2], 1, numC);
mouse = [];
for currA = 1:numC
    mouse = [mouse ones(1, (range+1)*2)*currA];
end
allTrans = nan(1, length(mouse));

animals = [];
prevAnimal = [];
currSim = 0;
aInd = 0;
for currSesh = 1:length(dayList)
    sessionName = dayList{currSesh};
    [animal, ~] = strtok(sessionName, 'd'); 
    animal = animal(2:end);
    if strcmp(animal, prevAnimal) == 0 || p.Results.sessionParamsFlag
        animals = [animals {animal}];
        if strcmp(animal, prevAnimal) == 0 
            fprintf('Simulating animal %s \n', animal);
        end
        if aInd > 0
            tmpInds = [find(isnan(allTrans), 1) : find(isnan(allTrans), 1) + (range+1)*2 - 1];
            allTrans(tmpInds) = [mean(transMedSim{aInd}(:,range:range+range)) mean(transHighSim{aInd}(:,range:range+range))];
        end
        
        if p.Results.sessionParamsFlag
            fracSesh = sum(~cellfun(@isempty, regexp(sList, sessionName))) / length(sList);
            runs = ceil(fracSesh * p.Results.runs);
        else
            fracSesh = sum(~cellfun(@isempty, regexp(sList, animal))) / length(sList);
            runs = ceil(fracSesh * p.Results.runs);
        end
        
        if runs > 0
            aInd = aInd + 1;
            currSim = currSim + 1;
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
                    [~, allRewards, allChoices, blockProbs, blockSwitch] = runSim_dF(p.Results.modelName, t.params(currSamp,:),...
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
                                        transHighSim{currSim} = [transHighSim{currSim}; allChoices((tmpInd-range+1):(tmpInd+range))];
                                        changeChoiceMatxHighSim{currSim} = [changeChoiceMatxHighSim{currSim}; changeChoice((tmpInd-range+1):(tmpInd+range))];
                                        prevRwdMatxHighSim{currSim} = [prevRwdMatxHighSim{currSim}; abs(prevRewards((tmpInd-range+1):(tmpInd+range)))];
                                        rwdMatxHighSim{currSim} = [rwdMatxHighSim{currSim}; prevRewards((tmpInd-range+1):(tmpInd+range))];
                                    end
                                elseif trialProbs(tmpInd-1, 2) == pLow & trialProbs(tmpInd, 2) == 10 & any(diff(trialProbs(tmpInd-tranWin:tmpInd, 1)) == probDiffH)
                                    if (tmpInd - range - 1) > 0 & length(allChoices) >= (tmpInd + range)
                                        transMedSim{currSim} = [transMedSim{currSim}; allChoices((tmpInd-range+1):(tmpInd+range))];
                                        changeChoiceMatxLowSim{currSim} = [changeChoiceMatxLowSim{currSim}; changeChoice((tmpInd-range+1):(tmpInd+range))];
                                        prevRwdMatxLowSim{currSim} = [prevRwdMatxLowSim{currSim}; abs(prevRewards((tmpInd-range+1):(tmpInd+range)))];
                                        rwdMatxMedSim{currSim} = [rwdMatxMedSim{currSim}; prevRewards((tmpInd-range+1):(tmpInd+range))];
                                    end
                                elseif trialProbs(tmpInd-1, 1) == pHigh & trialProbs(tmpInd, 1) == 10 & any(diff(trialProbs(tmpInd-tranWin:tmpInd, 2)) == probDiffH)
                                    if (tmpInd - range - 1) > 0 & length(allChoices) >= (tmpInd + range)
                                        transHighSim{currSim} = [transHighSim{currSim}; (allChoices((tmpInd-range+1):(tmpInd+range))*-1)];
                                        changeChoiceMatxHighSim{currSim} = [changeChoiceMatxHighSim{currSim}; changeChoice((tmpInd-range+1):(tmpInd+range))];
                                        prevRwdMatxHighSim{currSim} = [prevRwdMatxHighSim{currSim}; abs(prevRewards((tmpInd-range+1):(tmpInd+range)))];
                                        rwdMatxHighSim{currSim} = [rwdMatxHighSim{currSim}; prevRewards((tmpInd-range+1):(tmpInd+range))*-1];
                                    end 
                                elseif trialProbs(tmpInd-1, 1) == pLow & trialProbs(tmpInd, 1) == 10 & any(diff(trialProbs(tmpInd-tranWin:tmpInd, 2)) == probDiffH)
                                    if (tmpInd - range - 1) > 0 & length(allChoices) >= (tmpInd + range)
                                        transMedSim{currSim} = [transMedSim{currSim}; (allChoices((tmpInd-range+1):(tmpInd+range))*-1)];
                                        changeChoiceMatxLowSim{currSim} = [changeChoiceMatxLowSim{currSim}; changeChoice((tmpInd-range+1):(tmpInd+range))];
                                        prevRwdMatxLowSim{currSim} = [prevRwdMatxLowSim{currSim}; abs(prevRewards((tmpInd-range+1):(tmpInd+range)))];
                                        rwdMatxMedSim{currSim} = [rwdMatxMedSim{currSim}; prevRewards((tmpInd-range+1):(tmpInd+range))*-1];
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

%run lme
% tmpInds = [find(isnan(allTrans), 1) : find(isnan(allTrans), 1) + (range+1)*2 - 1];
% allTrans(tmpInds) = [mean(transMedSim{aInd}(:,range:range+range)) mean(transHighSim{aInd}(:,range:range+range))];
% tt = table(zscore(allTrans)', trialNum', transType', mouse', ...
%         'VariableNames', {'choiceProbs', 'trial', 'trans', 'mouse'});
% tt.trans = nominal(tt.trans);
% tt.mouse = nominal(tt.mouse);
% 
% mdl = fitlme(tt, 'choiceProbs~trial*trans + (trial*trans|mouse)')
% stats = anova(mdl)


%convert from cell to matrix (saves time)
transMedSim = cell2mat(transMedSim');
transHighSim = cell2mat(transHighSim');
changeChoiceMatxLowSim = cell2mat(changeChoiceMatxLowSim');
changeChoiceMatxHighSim = cell2mat(changeChoiceMatxHighSim');
prevRwdMatxLowSim = cell2mat(prevRwdMatxLowSim');
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

%%
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
ylabel('Choice probability')
legend('medium -> low', '', 'high -> low', '')
title('actual')


subplot(4,4,3:4); hold on
plot(x, medAvgSim, '-', 'Color', colors(2,:), 'linewidth', 1.5)
plot(x, highAvgSim, '-', 'Color', colors(4,:), 'linewidth', 1.5)
plot([0 0], [0 1], ':k')
plot([x(1) x(end)], [0.5 0.5], ':k')
ylim([0 1])
set(gca, 'tickdir', 'out')
legend('sim: medium -> low', 'sim: high -> low')
title('simulation')


subplot(4,4,[5:6]); hold on;
tmp = [rwdMed; rwdHigh];
rwdHx = sum(tmp(:, range-9:range), 2);
lowInds = find(rwdHx < 3 & rwdHx > 0);
highInds = find(rwdHx > 6);
tmp = [transMed; transHigh];
legTxt = [];
if ~isempty(lowInds)
    plotFilledBern(x, tmp(lowInds,:), colors(1,:));
    legTxt = [{'low rwd hist', ''}];
end
if ~isempty(highInds)
    plotFilledBern(x, tmp(highInds,:), colors(3,:));
    legTxt = [legTxt {'high rwd hist', ''}];
end
plot([0 0], [0 1], ':k')
plot([x(1) x(end)], [0.5 0.5], ':k')
set(gca, 'tickdir', 'out', 'box', 'off')
legend(legTxt)

subplot(4,4,[7:8]); hold on;
tmp = [rwdMatxMedSim; rwdMatxHighSim];
rwdHx = sum(tmp(:, range-9:range), 2);
lowInds = find(rwdHx < 3 & rwdHx > 0);
highInds = find(rwdHx > 6);
tmp = [transMedSim; transHighSim];
legTxt = [];
if ~isempty(lowInds)
    plot(x, mean(tmp(lowInds,:)), 'Color', colors(2,:), 'linewidth', 1.5)
    legTxt = [{'sim: low rwd hist'}];
end
if ~isempty(highInds)
    plot(x, mean(tmp(highInds,:)), 'Color', colors(3,:), 'linewidth', 1.5)
    legTxt = [legTxt {'sim: high rwd hist'}];
end
plot([0 0], [0 1], ':k')
plot([x(1) x(end)], [0.5 0.5], ':k')
set(gca, 'tickdir', 'out', 'box', 'off')
legend(legTxt)


%look at transitions by reward history
subplot(4,4,[9:10]); hold on;
colors2 = cool(2);
transAll = [transMed; transHigh];
choiceX = sum(transAll(:,range-9:range), 2);
xInds = find(choiceX == 10);
transAll = transAll(xInds, :);
rwdHx = [rwdMed; rwdHigh];
rwdHx = sum(rwdHx(:,range-9:range), 2);
rwdHx = rwdHx(xInds);
sortIndsI = find(rwdHx <= 5);
sortIndsII = find(rwdHx > 5);
plotFilled(x, transAll(sortIndsI, :), colors2(1,:));
plotFilled(x, transAll(sortIndsII, :), colors2(2,:));
plot([0 0], [0 1], ':k')
plot([x(1) x(end)], [0.5 0.5], ':k')
ylim([0 1])
ylabel('choice average')
title('medium transitions by rwd hist')
set(gca, 'tickdir', 'out')


subplot(4,4,[11:12]); hold on;
transAll = [transMedSim; transHighSim];
choiceX = sum(transAll(:,range-9:range), 2);
xInds = find(choiceX == 10);
transAll = transAll(xInds, :);
rwdHx = [rwdMatxMedSim; rwdMatxHighSim];
rwdHx = sum(rwdHx(:,range-9:range), 2);
rwdHx = rwdHx(xInds);
sortIndsI = find(rwdHx <= 5);
sortIndsII = find(rwdHx > 5);
plotFilled(x, transAll(sortIndsI, :), colors2(1,:));
plotFilled(x, transAll(sortIndsII, :), colors2(2,:));
plot([0 0], [0 1], ':k')
plot([x(1) x(end)], [0.5 0.5], ':k')
ylim([0 1])
title('medium transitions by rwd hist')
set(gca, 'tickdir', 'out')


subplot(4,4,[13:14]); hold on
plot(x, medAvg, '-', 'Color', colors(1,:), 'linewidth', 1.5)
plot(x, highAvg, '-', 'Color',  colors(3,:), 'linewidth', 1.5)
plot(x, medAvgSim, '-', 'Color', colors(2,:), 'linewidth', 1.5)
plot(x, highAvgSim, '-', 'Color', colors(4,:), 'linewidth', 1.5)
plot([0 0], [0 1], ':k')
plot([x(1) x(end)], [0.5 0.5], ':k')
ylim([0 1])
set(gca, 'tickdir', 'out')
medMdl = fitlm([medAvgSim], medAvg');
highMdl = fitlm([highAvgSim], highAvg');
xlabel('Trials from switch')
text('Position', [min(x) min(ylim)  0], 'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom', 'string', ...
        sprintf(['low R^2 = ' num2str(medMdl.Rsquared.Ordinary) '\nhigh R^2 = ' num2str(highMdl.Rsquared.Ordinary)]));

subplot(4,4,[15:16]); hold on
actualDiff = [highAvg - medAvg];
plot(x, actualDiff, '-k')
simDiff = [highAvgSim - medAvgSim];
plot(x, simDiff, '-', 'Color', [0.7 0.7 0.7])
plot([-range range], [0 0], ':k');
ylabel('Difference in choice probabilities')
xlabel('Trials from switch')
legend('high - medium', 'sim: high - medium')
set(gca, 'tickdir', 'out')
y = get(gca, 'ylim');
plot([0 0], y, ':k')
diffMdl = fitlm([highAvgSim - medAvgSim], [highAvg - medAvg]');
text('Position', [min(x) min(ylim)  0], 'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom', 'string', ...
        sprintf(['R^2 = ' num2str(diffMdl.Rsquared.Ordinary)]));

titleTxt = strrep(['Choice at block transitions: ' sheet ' ' category ' - ' p.Results.modelName], '_', ' ');
suptitle(titleTxt); 
set(gcf, 'renderer', 'painters', 'position', [-1439 42 1007 954])

%%
%make struct of linear models and choice probability curves
t.diffMdl = diffMdl; t.medMdl = medMdl; t.highMdl = highMdl;
t.actualDiff = actualDiff; t.simDiff = simDiff;
t.highAvg = highAvg; t.medAvg = medAvg;
t.highAvgSim = highAvgSim; t.medAvgSim = medAvgSim;
t.transMed = transMed; t.transHigh = transHigh

%save models and choice probability curves
savePath = [root 'transitionData' sep p.Results.modelName sep sheet sep 'transitionMdl' sep];
if ~exist(savePath)
    mkdir(savePath);
end
save([savePath sheet category '_transitionSim.mat'], 't');

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
    prevRwdMatxLowSim(rInd, tmp) = NaN;
    changeChoiceMatxLowSim(rInd, tmp) = NaN;
end


for tInd = 1:range*2
    lS_med(tInd) = sum(changeChoiceMatxLowSim(find(prevRwdMatxLowSim(:,tInd)==0), tInd))/sum(prevRwdMatxLowSim(:,tInd)==0);
    lS_high(tInd) = sum(changeChoiceMatxHighSim(find(prevRwdMatxHighSim(:,tInd)==0), tInd))/sum(prevRwdMatxHighSim(:,tInd)==0);

    wS_med(tInd) = 1 - ((sum(changeChoiceMatxLowSim(find(prevRwdMatxLowSim(:,tInd)==1), tInd)))/sum(prevRwdMatxLowSim(:,tInd)==1));
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
legend('medium -> low', 'high -> low')
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
legend('sim: medium -> low', 'sim: high -> low')
ylim([yls])
plot([0 0], yls, '--k')
set(gca, 'tickdir', 'out')
title('simulated')

titleTxt = strrep(['Win-Stay Lose-Shift: ' sheet ' ' category ' - ' p.Results.modelName], '_', ' ');
suptitle(titleTxt); 
set(gcf, 'renderer', 'painters', 'position', [-1439 42 1007 954])


