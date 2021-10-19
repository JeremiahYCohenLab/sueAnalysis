function compareTransitionSimManipSamp_dF(xlFile, sheet, pre, post, varargin)

p = inputParser;
% default parameters if none given
p.addParameter('rwdProbs', [90 50 10])
p.addParameter('tranWin', 5)
p.addParameter('maxTrials', 350)
p.addParameter('runs', 10) %max runs per sample per condition
p.addParameter('randomSeed', 698512)
p.addParameter('modelName', 'sevenParam_absPePeAN_scale_int_bias_ord')
p.addParameter('params', [])
p.addParameter('sessionParamsFlag', 0)
p.addParameter('biasFlag', 0)
p.addParameter('bernFlag', 1)
p.addParameter('manipMdlFlag', 0)
p.addParameter('samps', 1000)
p.addParameter('lesionInd', [])
p.addParameter('lesionVal', 0)
p.parse(varargin{:});

[root, sep] = currComputer();

pHigh = p.Results.rwdProbs(1);
pMed = p.Results.rwdProbs(2);
probDiffH = pHigh - 10;
tranWin = p.Results.tranWin;

if isempty(p.Results.lesionInd)
    beh = [{pre} {post}];
else
    beh = [{pre} {pre}];
end


[wsls{1}, transMed{1}, transHigh{1}, aNames{1}] = transitionAnalysis_opMD(xlFile, sheet, pre, p.Results.rwdProbs, p.Results.tranWin);
[wsls{2}, transMed{2}, transHigh{2}, aNames{2}] = transitionAnalysis_opMD(xlFile, sheet, post, p.Results.rwdProbs, p.Results.tranWin);

transMedSim = cell(1,2);
transHighSim = cell(1,2);
changeChoiceMatxMed = cell(1,2);
changeChoiceMatxHigh = cell(1,2);
prevRwdMatxMed = cell(1,2);
prevRwdMatxHigh = cell(1,2);

range = 15;
rSeed = p.Results.randomSeed;
runTot = 0;
for currP = 1:2
    
    if currP == 1 || isempty(p.Results.lesionInd)
        [~, seshList, ~] = xlsread(xlFile, sheet);
        [~,col] = find(strcmp(seshList, beh{currP}));
        seshList = seshList(2:end, col);
        endInd = find(cellfun(@isempty,seshList),1);
        if ~isempty(endInd)
            seshList = seshList(1:endInd-1,:);
        end
    end
    
    prevAnimals = cell(1);
    for currSesh = 1:length(seshList)
        sessionName = seshList{currSesh};
        [animal, ~] = strtok(sessionName, 'd'); 
        animal = animal(2:end);
        if sum(~cellfun(@isempty, regexp(prevAnimals, animal))) == 0
            fprintf('Simulating animal %s %s \n', animal, beh{currP});
            fracTrans =  sum(~cellfun(@isempty, regexp(aNames{currP}, animal))) / length(aNames{currP});
            runs = ceil(fracTrans * p.Results.runs);
            runTot = runTot + runs*p.Results.samps;
            
            if isempty(p.Results.params)
                if p.Results.manipMdlFlag
                    if p.Results.bernFlag
                        modelPath = [root animal sep animal 'sorted' sep 'stan' sep 'manip' sep 'bernoulli' sep p.Results.modelName '_manip'...
                            sep animal pre post '_' p.Results.modelName '_manip.mat'];
                    else
                        modelPath = [root animal sep animal 'sorted' sep 'stan' sep 'manip' sep p.Results.modelName '_manip'...
                            sep animal pre post '_' p.Results.modelName '_manip.mat'];
                    end
                    if p.Results.sessionParamsFlag
                        t = getStanModelParamsManip_samps(p.Results.modelName, modelPath, p.Results.samps, 'sessionParamsFlag', 1,...
                            'sessionName', sessionName, 'biasFlag', p.Results.biasFlag);
                    else
                        t = getStanModelParamsManip_samps(p.Results.modelName, modelPath, p.Results.samps);
                    end
                    params{1} = t.paramsPre;
                    if ~isempty(p.Results.lesionInd)
                        params{2} = t.paramsPre;
                        params{2}(:, p.Results.lesionInd) = p.Results.lesionVal;
                    else
                        params{2} = t.paramsPost;
                    end
                else
                    if p.Results.bernFlag
                        modelPath = [root animal sep animal 'sorted' sep 'stan' sep 'bernoulli' sep p.Results.modelName...
                            sep animal beh{currP} '_' p.Results.modelName '.mat'];
                    else
                        modelPath = [root animal sep animal 'sorted' sep 'stan' sep 'manip' sep p.Results.modelName...
                            sep animal beh{currP} '_' p.Results.modelName '.mat'];
                    end
                    if p.Results.sessionParamsFlag
                        t = getStanModelParams_samps(p.Results.modelName, modelPath, p.Results.samps, 'sessionParamsFlag', 1,...
                            'sessionName', sessionName);
                    else
                        t = getStanModelParams_samps(p.Results.modelName, modelPath, p.Results.samps);
                    end
                    params{currP} = t.params;
                    if ~isempty(p.Results.lesionInd) & currP == 2
                        params{currP}(:,p.Results.lesionInd) = p.Results.lesionVal;
                    end
                end
            else
                params = p.Results.params;
            end
            for currSamp = 1:p.Results.samps
                for run = 1:runs
                    if rem(run,100) == 0
                        fprintf('Running simulation %d of %d \n', run, runs);
                    end
                    rSeed = rSeed + 1;
                    [~,allRewards, allChoices, blockProbs, blockSwitch] = runSim_dF(p.Results.modelName, params{currP}(run,:),...
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
                                        transHighSim{currP} = [transHighSim{currP}; allChoices((tmpInd-range+1):(tmpInd+range))];
                                        changeChoiceMatxHigh{currP} = [changeChoiceMatxHigh{currP}; changeChoice((tmpInd-range+1):(tmpInd+range))];
                                        prevRwdMatxHigh{currP} = [prevRwdMatxHigh{currP}; prevRewardsBin((tmpInd-range+1):(tmpInd+range))];
                                    end
                            elseif trialProbs(tmpInd-1, 2) == pMed & trialProbs(tmpInd, 2) == 10 & any(diff(trialProbs(tmpInd-tranWin:tmpInd, 1)) == probDiffH)
                                    if (tmpInd - range - 1) > 0 & length(allChoices) > (tmpInd + range)
                                        transMedSim{currP} = [transMedSim{currP}; allChoices((tmpInd-range+1):(tmpInd+range))];
                                        changeChoiceMatxMed{currP} = [changeChoiceMatxMed{currP}; changeChoice((tmpInd-range+1):(tmpInd+range))];
                                        prevRwdMatxMed{currP} = [prevRwdMatxMed{currP}; prevRewardsBin((tmpInd-range+1):(tmpInd+range))]; 
                                    end
                            elseif trialProbs(tmpInd-1, 1) == pHigh & trialProbs(tmpInd, 1) == 10 & any(diff(trialProbs(tmpInd-tranWin:tmpInd, 2)) == probDiffH)
                                    if (tmpInd - range - 1) > 0 & length(allChoices) > (tmpInd + range)
                                        transHighSim{currP} = [transHighSim{currP}; (allChoices((tmpInd-range+1):(tmpInd+range))*-1)];
                                        changeChoiceMatxHigh{currP} = [changeChoiceMatxHigh{currP}; changeChoice((tmpInd-range+1):(tmpInd+range))];
                                        prevRwdMatxHigh{currP} = [prevRwdMatxHigh{currP}; prevRewardsBin((tmpInd-range+1):(tmpInd+range))];
                                    end 
                            elseif trialProbs(tmpInd-1, 1) == pMed & trialProbs(tmpInd, 1) == 10 & any(diff(trialProbs(tmpInd-tranWin:tmpInd, 2)) == probDiffH)
                                    if (tmpInd - range - 1) > 0 & length(allChoices) > (tmpInd + range)
                                        transMedSim{currP} = [transMedSim{currP}; (allChoices((tmpInd-range+1):(tmpInd+range))*-1)];
                                        changeChoiceMatxMed{currP} = [changeChoiceMatxMed{currP}; changeChoice((tmpInd-range+1):(tmpInd+range))];
                                        prevRwdMatxMed{currP} = [prevRwdMatxMed{currP}; prevRewardsBin((tmpInd-range+1):(tmpInd+range))]; 
                                    end
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
medAvg{1} = mean(transMed{1},1);
highAvg{1} = mean(transHigh{1},1);
medAvg{2} = mean(transMed{2},1);
highAvg{2} = mean(transHigh{2},1);


%find averages of simulated behavior
transMedSim{1}(transMedSim{1} == -1) = 0;
transHighSim{1}(transHighSim{1} == -1) = 0;
transMedSim{2}(transMedSim{2} == -1) = 0;
transHighSim{2}(transHighSim{2} == -1) = 0;

medAvgSim{1} = mean(transMedSim{1},1);
highAvgSim{1} = mean(transHighSim{1},1);
medAvgSim{2} = mean(transMedSim{2},1);
highAvgSim{2} = mean(transHighSim{2},1);

colors = cool(4);
figure;
x = [-range+1:range];
for currP = 1:2
    %plot actual transition behavior for pre/post
    subplot(4,2,currP); hold on
    plotFilledBern(x, transMed{currP}, colors(currP,:));
    plotFilledBern(x, transHigh{currP}, colors(currP+2,:));
    plot([-range range], [0.5 0.5], ':k')
    plot([0 0], [0 1], ':k')
    ylim([0 1])
    set(gca, 'tickdir', 'out')
    ylabel('choice probability')
    legend(' medium -> low', '', 'high -> low', '')
    title(['actual ' beh{currP}])

    %plot transition behavior for simualted pre/post
    subplot(4,2,currP+2); hold on
    plot(x, medAvgSim{currP}, '-', 'Color', colors(currP,:), 'linewidth', 2)
    plot(x, highAvgSim{currP}, '-', 'Color', colors(currP+2,:), 'linewidth', 2)
    plot([-range range], [0.5 0.5], ':k')
    plot([0 0], [0 1], ':k')
    ylim([0 1])
    ylabel('choice probability')
    set(gca, 'tickdir', 'out')
    if isempty(p.Results.lesionInd) || currP == 1
        title(['simulated ' beh{currP}])
    else
        title('simulated lesion')
    end

    %plot differences in choice probs between conditions for actual and sim pre/post
    subplot(4,2,currP+4); hold on
    actualDiff{currP} = [highAvg{currP} - medAvg{currP}];
    plot(x, actualDiff{currP}, '-k', 'linewidth', 2)
    simDiff{currP} = [highAvgSim{currP} - medAvgSim{currP}];
    plot(x, simDiff{currP}, '-', 'Color', [0.7 0.7 0.7], 'linewidth', 2)
    plot([-range range], [0 0], ':k')
    xlabel('Trials from switch')
    ylabel('Choice probability difference')
    legend('high - medium', 'sim: high - medium')
    set(gca, 'tickdir', 'out')
    y = get(gca, 'ylim');
    plot([0 0], y, ':k')
    diffMdl{currP} = fitlm([highAvgSim{currP} - medAvgSim{currP}], [highAvg{currP} - medAvg{currP}]');
    text('Position', [min(x) min(ylim)  0], 'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom', 'string', ...
            sprintf(['R^2 = ' num2str(diffMdl{currP}.Rsquared.Ordinary)]));
        
    subplot(4,2,currP+6); hold on
    plot(x, medAvg{currP}, '-', 'Color',  colors(currP,:));
    plot(x, highAvg{currP}, '-', 'Color',  colors(currP+2,:));
    plot(x, medAvgSim{currP}, '-', 'Color', colors(currP,:), 'linewidth', 2)
    plot(x, highAvgSim{currP}, '-', 'Color', colors(currP+2,:), 'linewidth', 2)
    plot([-range range], [0.5 0.5], ':k')
    plot([0 0], [0 1], ':k')

    medMdl{currP} = fitlm(medAvgSim{currP}, medAvg{currP}');
    text('Position', [min(x) min(ylim) + 0.1  0], 'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom', 'string', ...
            sprintf(['low R^2 = ' num2str(medMdl{currP}.Rsquared.Ordinary)]));
    highMdl{currP} = fitlm(highAvgSim{currP}, highAvg{currP}');
    text('Position', [min(x) min(ylim)  0], 'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom', 'string', ...
            sprintf(['high R^2 = ' num2str(highMdl{currP}.Rsquared.Ordinary)]));
        
    ylim([0 1])
    set(gca, 'tickdir', 'out')
    ylabel('choice probability')
    title(['actual v sim' beh{currP}])
end

set(gcf, 'renderer', 'painters', 'position', [-1467 42 1024 954])
if isempty(p.Results.lesionInd)
    titleTxt = strrep([sheet ' ' p.Results.modelName ' transition behavior'], '_', ' ');
else
    titleTxt = strrep([sheet ' ' p.Results.modelName ' transition behavior - simulated lesion'], '_', ' ');
end
suptitle(titleTxt)

t.diffMdl = diffMdl; t.medMdl = medMdl; t.highMdl = highMdl;
t.actualDiff = actualDiff; t.simDiff = simDiff;
t.highAvg = highAvg; t.medAvg = medAvg;
t.highAvgSim = highAvgSim; t.medAvgSim = medAvgSim;

%save models and choice probability curves
if p.Results.manipMdlFlag
    savePath = [root 'transitionData' sep p.Results.modelName sep sheet sep 'transitionManipMdl' sep];
else
    savePath = [root 'transitionData' sep p.Results.modelName sep sheet sep 'transitionManip' sep];
end
if ~exist(savePath)
    mkdir(savePath);
end
if isempty(p.Results.lesionInd)
    save([savePath beh{1} beh{2} '_transitionSim.mat'], 't');
else
    save([savePath beh{1} beh{2} '_transitionSimLesion.mat'], 't');
end

%% plot decay constants from exp fits
%figure; hold on;
%x = [0:15];
%ft = fittype('a*exp((-1/b)*x)+c');
%sp = [0.5 7 0.1];

%mdl = fit(x',medAvg{1}(range:range+15)',ft, 'start', sp);
%lowTau_actual(1) = mdl.b;
%tmp = confint(mdl); 
%lowTau_actual_CI(1) = mdl.b - tmp(1,2);
%mdl = fit(x',medAvg{2}(range:range+15)',ft, 'start', sp);
%lowTau_actual(2) = mdl.b;
%tmp = confint(mdl); 
%lowTau_actual_CI(2) = mdl.b - tmp(1,2);
%errorbar(lowTau_actual, lowTau_actual_CI, '-', 'Color', [0.7 0.7 0.7], 'linewidth', 1.5)


%% look at win-stay lose-shift around transitions
for currP = 1:2
    for tInd = 1:range*2
        lS_low{currP}(tInd) = sum(changeChoiceMatxMed{currP}(find(prevRwdMatxMed{currP}(:,tInd)==0), tInd))/sum(prevRwdMatxMed{currP}(:,tInd)==0);
        sem_lS_low{currP}(tInd) = sem_bernoulli(sum(changeChoiceMatxMed{currP}(find(prevRwdMatxMed{currP}(:,tInd)==0), tInd)), sum(prevRwdMatxMed{currP}(:,tInd)==0));
        lS_high{currP}(tInd) = sum(changeChoiceMatxHigh{currP}(find(prevRwdMatxHigh{currP}(:,tInd)==0), tInd))/sum(prevRwdMatxHigh{currP}(:,tInd)==0);
        sem_lS_high{currP}(tInd) = sem_bernoulli(sum(changeChoiceMatxHigh{currP}(find(prevRwdMatxHigh{currP}(:,tInd)==0), tInd)), sum(prevRwdMatxHigh{currP}(:,tInd)==0));

        wS_low{currP}(tInd) = 1 - ((sum(changeChoiceMatxMed{currP}(find(prevRwdMatxMed{currP}(:,tInd)==1), tInd)))/sum(prevRwdMatxMed{currP}(:,tInd)==1));
        sem_wS_low{currP}(tInd) = sem_bernoulli(sum(~changeChoiceMatxMed{currP}(find(prevRwdMatxMed{currP}(:,tInd)==1), tInd)), sum(prevRwdMatxMed{currP}(:,tInd)==1));
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

fprintf('Run total: %d \n', runTot);

end

