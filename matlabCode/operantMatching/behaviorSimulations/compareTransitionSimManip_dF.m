function compareTransitionSimManip_dF(xlFile, sheet, pre, post, varargin)

p = inputParser;
% default parameters if none given
p.addParameter('rwdProbs', [90 50 10])
p.addParameter('tranWin', 3)
p.addParameter('maxTrials', 350)
p.addParameter('runs', 3000)
p.addParameter('randomSeed', 98207)
p.addParameter('modelName', 'eightParam_absPePeAN_scale_int_bias')
p.addParameter('params', [])
p.addParameter('figFlag', 0)
p.addParameter('bernFlag', 1)
p.addParameter('manipMdlFlag', 1)
p.addParameter('lesionInd', [])
p.addParameter('lesionVal', 0)
p.parse(varargin{:});

[root, sep] = currComputer();

pHigh = p.Results.rwdProbs(1);
pLow = p.Results.rwdProbs(2);
probDiffH = pHigh - 10;
tranWin = p.Results.tranWin;
beh = [{pre} {post}];

[wsls{1}, transLow_actual{1}, transHigh_actual{1}, aNames{1}] = transitionAnalysis_opMD(xlFile, sheet, pre, p.Results.rwdProbs, p.Results.tranWin);
[wsls{2}, transLow_actual{2}, transHigh_actual{2}, aNames{2}] = transitionAnalysis_opMD(xlFile, sheet, post, p.Results.rwdProbs, p.Results.tranWin);

transLow = cell(1,2);
transHigh = cell(1,2);
changeChoiceMatxLow= cell(1,2);
changeChoiceMatxHigh = cell(1,2);
prevRwdMatxLow = cell(1,2);
prevRwdMatxHigh = cell(1,2);

range = 20;
rSeed = p.Results.randomSeed;

for currP = 1:2
        
    [~, seshList, ~] = xlsread(xlFile, sheet);
    [~,col] = find(strcmp(seshList, beh{currP}));
    seshList = seshList(2:end, col);
    endInd = find(cellfun(@isempty,seshList),1);
    if ~isempty(endInd)
        seshList = seshList(1:endInd-1,:);
    end
    
    prevAnimals = cell(1);
    for currSesh = 1:length(seshList)
        sessionName = seshList{currSesh};
        [animal, ~] = strtok(sessionName, 'd'); 
        animal = animal(2:end);
        if sum(~cellfun(@isempty, regexp(prevAnimals, animal))) == 0
            if currP == 1
             fprintf('Simulating animal %s pre \n', animal);
            else
                fprintf('Simulating animal %s post \n', animal);
            end
            fracTrans =  sum(~cellfun(@isempty, regexp(aNames{currP}, animal))) / length(aNames{currP});
            runs = ceil(fracTrans * p.Results.runs);
            %fracSesh = sum(~cellfun(@isempty, regexp(seshList, animal))) / length(seshList);
            %runs = ceil(fracSesh * p.Results.runs);
            
            if isempty(p.Results.params)
                if p.Results.manipMdlFlag
                    if p.Results.bernFlag
                        modelPath = [root animal sep animal 'sorted' sep 'stan' sep 'manip' sep 'bernoulli' sep p.Results.modelName '_manip'...
                            sep animal pre post '_' p.Results.modelName '_manip.mat'];
                    else
                        modelPath = [root animal sep animal 'sorted' sep 'stan' sep 'manip' sep p.Results.modelName '_manip'...
                            sep animal pre post '_' p.Results.modelName '_manip.mat'];
                    end
                    tmp = load(modelPath);
                    params = tmp.paramEsts;
                    clear tmp;
                    
                    if ~isempty(p.Results.lesionInd)
                        params(2,:) = params(1,:);
                        params(2, p.Results.lesionInd) = p.Results.lesionVal;
                    end
                else
                    if p.Results.bernFlag
                        modelPath = [root animal sep animal 'sorted' sep 'stan' sep 'bernoulli' sep p.Results.modelName...
                            sep animal pre '_' p.Results.modelName '.mat'];
                    else
                        modelPath = [root animal sep animal 'sorted' sep 'stan' sep 'manip' sep p.Results.modelName...
                            sep animal pre '_' p.Results.modelName '.mat'];
                    end
                    tmp = load(modelPath);
                    params = tmp.paramEsts;
                    params(2,:) = params(1,:);
                    params(2, p.Results.lesionInd) = p.Results.lesionVal;
                    clear tmp;
                end
            else
                params = p.Results.params;
            end

            for run = 1:runs
                if rem(run,100) == 0
                    fprintf('Running simulation %d of %d \n', run, runs);
                end
                rSeed = rSeed + 1;
                [~,allRewards, allChoices, blockProbs, blockSwitch] = runSim_dF(p.Results.modelName, params(currP,:),...
                    1000, rSeed, p.Results.rwdProbs);

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
                    if tmpInd-tranWin > 0
                        if trialProbs(tmpInd-1, 2) == pHigh & trialProbs(tmpInd, 2) == 10 & any(diff(trialProbs(tmpInd-tranWin:tmpInd, 1)) == probDiffH)
                                if (tmpInd - range - 1) > 0 & length(allChoices) > (tmpInd + range)
                                    transHigh{currP} = [transHigh{currP}; allChoices((tmpInd-range+1):(tmpInd+range))];
                                    changeChoiceMatxHigh{currP} = [changeChoiceMatxHigh{currP}; changeChoice((tmpInd-range+1):(tmpInd+range))];
                                    prevRwdMatxHigh{currP} = [prevRwdMatxHigh{currP}; prevRewardsBin((tmpInd-range+1):(tmpInd+range))];
                                end
                        elseif trialProbs(tmpInd-1, 2) == pLow & trialProbs(tmpInd, 2) == 10 & any(diff(trialProbs(tmpInd-tranWin:tmpInd, 1)) == probDiffH)
                                if (tmpInd - range - 1) > 0 & length(allChoices) > (tmpInd + range)
                                    transLow{currP} = [transLow{currP}; allChoices((tmpInd-range+1):(tmpInd+range))];
                                    changeChoiceMatxLow{currP} = [changeChoiceMatxLow{currP}; changeChoice((tmpInd-range+1):(tmpInd+range))];
                                    prevRwdMatxLow{currP} = [prevRwdMatxLow{currP}; prevRewardsBin((tmpInd-range+1):(tmpInd+range))]; 
                                end
                        elseif trialProbs(tmpInd-1, 1) == pHigh & trialProbs(tmpInd, 1) == 10 & any(diff(trialProbs(tmpInd-tranWin:tmpInd, 2)) == probDiffH)
                                if (tmpInd - range - 1) > 0 & length(allChoices) > (tmpInd + range)
                                    transHigh{currP} = [transHigh{currP}; (allChoices((tmpInd-range+1):(tmpInd+range))*-1)];
                                    changeChoiceMatxHigh{currP} = [changeChoiceMatxHigh{currP}; changeChoice((tmpInd-range+1):(tmpInd+range))];
                                    prevRwdMatxHigh{currP} = [prevRwdMatxHigh{currP}; prevRewardsBin((tmpInd-range+1):(tmpInd+range))];
                                end 
                        elseif trialProbs(tmpInd-1, 1) == pLow & trialProbs(tmpInd, 1) == 10 & any(diff(trialProbs(tmpInd-tranWin:tmpInd, 2)) == probDiffH)
                                if (tmpInd - range - 1) > 0 & length(allChoices) > (tmpInd + range)
                                    transLow{currP} = [transLow{currP}; (allChoices((tmpInd-range+1):(tmpInd+range))*-1)];
                                    changeChoiceMatxLow{currP} = [changeChoiceMatxLow{currP}; changeChoice((tmpInd-range+1):(tmpInd+range))];
                                    prevRwdMatxLow{currP} = [prevRwdMatxLow{currP}; prevRewardsBin((tmpInd-range+1):(tmpInd+range))]; 
                                end
                        end
                    end
                end
            end
            prevAnimals = [prevAnimals {animal}];
        end
    end
end

%find averages and sem of actual transition behavior
lowAvg_actual{1} = mean(transLow_actual{1},1);
highAvg_actual{1} = mean(transHigh_actual{1},1);
lowAvg_actual{2} = mean(transLow_actual{2},1);
highAvg_actual{2} = mean(transHigh_actual{2},1);


%find averages of simulated behavior
transLow{1}(transLow{1} == -1) = 0;
transHigh{1}(transHigh{1} == -1) = 0;
transLow{2}(transLow{2} == -1) = 0;
transHigh{2}(transHigh{2} == -1) = 0;

lowAvg{1} = mean(transLow{1},1);
highAvg{1} = mean(transHigh{1},1);
lowAvg{2} = mean(transLow{2},1);
highAvg{2} = mean(transHigh{2},1);

colors = cool(4);
legTxt = {'pre', 'post'};
figure;
x = [-range+1:range];
for currP = 1:2
    %plot actual transition behavior for pre/post
    subplot(3,2,currP); hold on
    plotFilledBern(x, transLow_actual{currP}, colors(currP,:));
    plotFilledBern(x, transHigh_actual{currP}, colors(currP+2,:));
    plot([-range range], [0.5 0.5], ':k')
    plot([0 0], [0 1], ':k')
    ylim([0 1])
    set(gca, 'tickdir', 'out')
    ylabel('choice probability')
    legend(' medium -> low', 'high -> low')
    title(['actual ' legTxt{currP}])

    %plot transition behavior for simualted pre/post
    subplot(3,2,currP+2); hold on
    plot(x, lowAvg{currP}, '-', 'Color', colors(currP,:), 'linewidth', 2)
    plot(x, highAvg{currP}, '-', 'Color', colors(currP+2,:), 'linewidth', 2)
    plot([-range range], [0.5 0.5], ':k')
    plot([0 0], [0 1], ':k')
    ylim([0 1])
    ylabel('choice probability')
    set(gca, 'tickdir', 'out')
    if isempty(p.Results.lesionInd)
        title(['simulated ' legTxt{currP}])
    else
        title('simulated lesion')
    end

    %plot differences in choice probs between conditions for actual and sim pre/post
    subplot(3,2,currP+4); hold on
    plot(x,[highAvg_actual{currP} - lowAvg_actual{currP}], '-k', 'linewidth', 2)
    plot(x,[highAvg{currP} - lowAvg{currP}], '-', 'Color', [0.7 0.7 0.7], 'linewidth', 2)
    plot([-range range], [0 0], ':k')
    xlabel('Trials from switch')
    ylabel('Choice probability difference')
    legend('high - medium', 'sim: high - medium')
    set(gca, 'tickdir', 'out')
    y = get(gca, 'ylim');
    plot([0 0], y, ':k')
    mdl = fitlm([highAvg{currP} - lowAvg{currP}], [highAvg_actual{currP} - lowAvg_actual{currP}]');
    text('Position', [min(x) min(ylim)  0], 'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom', 'string', ...
            sprintf(['R^2 = ' num2str(mdl.Rsquared.Ordinary)]));
end

set(gcf, 'renderer', 'painters', 'position', [-1467 42 1024 954])
if isempty(p.Results.lesionInd)
    titleTxt = strrep([sheet ' ' p.Results.modelName ' transition behavior'], '_', ' ');
else
    titleTxt = strrep([sheet ' ' p.Results.modelName ' transition behavior - simulated lesion'], '_', ' ');
end
suptitle(titleTxt)

%% plot decay constants from exp fits
%figure; hold on;
%x = [0:15];
%ft = fittype('a*exp((-1/b)*x)+c');
%sp = [0.5 7 0.1];

%mdl = fit(x',lowAvg_actual{1}(range:range+15)',ft, 'start', sp);
%lowTau_actual(1) = mdl.b;
%tmp = confint(mdl); 
%lowTau_actual_CI(1) = mdl.b - tmp(1,2);
%mdl = fit(x',lowAvg_actual{2}(range:range+15)',ft, 'start', sp);
%lowTau_actual(2) = mdl.b;
%tmp = confint(mdl); 
%lowTau_actual_CI(2) = mdl.b - tmp(1,2);
%errorbar(lowTau_actual, lowTau_actual_CI, '-', 'Color', [0.7 0.7 0.7], 'linewidth', 1.5)


%% look at win-stay lose-shift around transitions
for currP = 1:2
    for tInd = 1:range*2
        lS_low{currP}(tInd) = sum(changeChoiceMatxLow{currP}(find(prevRwdMatxLow{currP}(:,tInd)==0), tInd))/sum(prevRwdMatxLow{currP}(:,tInd)==0);
        sem_lS_low{currP}(tInd) = sem_bernoulli(sum(changeChoiceMatxLow{currP}(find(prevRwdMatxLow{currP}(:,tInd)==0), tInd)), sum(prevRwdMatxLow{currP}(:,tInd)==0));
        lS_high{currP}(tInd) = sum(changeChoiceMatxHigh{currP}(find(prevRwdMatxHigh{currP}(:,tInd)==0), tInd))/sum(prevRwdMatxHigh{currP}(:,tInd)==0);
        sem_lS_high{currP}(tInd) = sem_bernoulli(sum(changeChoiceMatxHigh{currP}(find(prevRwdMatxHigh{currP}(:,tInd)==0), tInd)), sum(prevRwdMatxHigh{currP}(:,tInd)==0));

        wS_low{currP}(tInd) = 1 - ((sum(changeChoiceMatxLow{currP}(find(prevRwdMatxLow{currP}(:,tInd)==1), tInd)))/sum(prevRwdMatxLow{currP}(:,tInd)==1));
        sem_wS_low{currP}(tInd) = sem_bernoulli(sum(~changeChoiceMatxLow{currP}(find(prevRwdMatxLow{currP}(:,tInd)==1), tInd)), sum(prevRwdMatxLow{currP}(:,tInd)==1));
        wS_high{currP}(tInd) = 1 - ((sum(changeChoiceMatxHigh{currP}(find(prevRwdMatxHigh{currP}(:,tInd)==1), tInd)))/sum(prevRwdMatxHigh{currP}(:,tInd)==1));
        sem_wS_high{currP}(tInd) = sem_bernoulli(sum(~changeChoiceMatxHigh{currP}(find(prevRwdMatxHigh{currP}(:,tInd)==1), tInd)), sum(prevRwdMatxHigh{currP}(:,tInd)==1));
    end
end

figure;
for currP = 1:2
    %actual win-stay pre/post
    subplot(2,4,(currP-1)*2+1); hold on
    errorbar(x, wsls{currP}.wS_low, wsls{currP}.sem_wS_low, 'Color', colors(currP,:), 'linewidth', 1.3)
    errorbar(x, wsls{currP}.wS_high, wsls{currP}.sem_wS_high, 'Color', colors(currP+2,:), 'linewidth', 1.3)
    ylabel('win-stay')
    yws = ylim;
    plot([0 0], yws, '--k')
    set(gca, 'tickdir', 'out')
    title('actual')

    %actual lose-shift pre/post
    subplot(2,4,(currP-1)*2+2); hold on
    errorbar(x, wsls{currP}.lS_low, wsls{currP}.sem_lS_low, 'Color', colors(currP,:), 'linewidth', 1.3)
    errorbar(x, wsls{currP}.lS_high, wsls{currP}.sem_lS_high, 'Color', colors(currP+2,:), 'linewidth', 1.3)
    ylabel('lose-shift')
    yls = ylim;
    plot([0 0], yws, '--k')
    set(gca, 'tickdir', 'out')
    title('actual')

    %simulated win-stay pre/post
    subplot(2,4,(currP-1)*2+5); hold on
    errorbar(x, wS_low{currP}, sem_wS_low{currP}, 'Color', colors(currP,:), 'linewidth', 1.3)
    errorbar(x, wS_high{currP}, sem_wS_high{currP}, 'Color', colors(currP+2,:), 'linewidth', 1.3)
    xlabel('Trials from switch')
    ylabel('win-stay')
    ylim([yws])
    plot([0 0], yws, '--k')
    set(gca, 'tickdir', 'out')
    title('simulated')

    %simulated lose-shift pre
    subplot(2,4,(currP-1)*2+6); hold on
    errorbar(x, lS_low{currP}, sem_lS_low{currP}, 'Color', colors(currP,:), 'linewidth', 1.3)
    errorbar(x, lS_high{currP}, sem_lS_low{currP}, 'Color', colors(currP+2,:), 'linewidth', 1.3)
    xlabel('Trials from switch')
    ylabel('lose-shift')
    ylim([yls])
    plot([0 0], yls, '--k')
    set(gca, 'tickdir', 'out')
    legend('medium->low', 'high->low')
    title('simulated')
end
set(gcf, 'renderer', 'painters', 'position', [-1467 42 1024 954])
suptitle(titleTxt)

end

