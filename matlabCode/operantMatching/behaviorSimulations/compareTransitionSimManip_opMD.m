function compareTransitionSimManip_opMD(xlFile, sheet, pre, post, varargin)

p = inputParser;
% default parameters if none given
p.addParameter('rwdProbs', [90 50 10]);
p.addParameter('tranWin', 3);
p.addParameter('maxTrials', 350);
p.addParameter('runs', 2000);
p.addParameter('randomSeed', 98207);
p.addParameter('modelName', 'sixParam_absPePeAN_bi');
p.addParameter('figFlag', 0);
p.addParameter('bernFlag', 1);
p.parse(varargin{:});

[root, sep] = currComputer();

pHigh = p.Results.rwdProbs(1);
pLow = p.Results.rwdProbs(2);
probDiffH = pHigh - 10;
probDiffM = pLow - 10;
tranWin = p.Results.tranWin;
beh = [{pre} {post}];

[wslsPre, transLowPre_actual, transHighPre_actual] = transitionAnalysis_opMD(xlFile, sheet, pre, p.Results.rwdProbs, p.Results.tranWin);
[wslsPost, transLowPost_actual, transHighPost_actual] = transitionAnalysis_opMD(xlFile, sheet, post, p.Results.rwdProbs, p.Results.tranWin);

transLowPre = [];
transHighPre = [];
changeChoiceMatxLowPre= [];
changeChoiceMatxHighPre = [];
prevRwdMatxLowPre = [];
prevRwdMatxHighPre = [];

transLowPost = [];
transHighPost = [];
changeChoiceMatxLowPost= [];
changeChoiceMatxHighPost = [];
prevRwdMatxLowPost = [];
prevRwdMatxHighPost = [];


range = 20;
rSeed = p.Results.randomSeed;

for currP = 1:2
        
    [~, seshList, ~] = xlsread(xlFile, sheet);
    [~,col] = find(~cellfun(@isempty,strfind(seshList, beh{currP})) == 1);
    seshList = seshList(2:end, col);
    endInd = find(cellfun(@isempty,seshList),1);
    if ~isempty(endInd)
        seshList = seshList(1:endInd-1,:);
    end
    
    prevAnimal = [];
    for currSesh = 1:length(seshList)
        sessionName = seshList{currSesh};
        [animal, ~] = strtok(sessionName, 'd'); 
        animal = animal(2:end);
        if strcmp(animal, prevAnimal) == 0
            fprintf('Simulating animal %s \n', animal);
            fracSesh = sum(~cellfun(@isempty, regexp(seshList, animal))) / length(seshList);
            runs = ceil(fracSesh * p.Results.runs);
            
            if p.Results.bernFlag
                modelPath = [root animal sep animal 'sorted' sep 'stan' sep 'manip' sep 'bernoulli' sep p.Results.modelName '_manip'...
                    sep animal pre post '_' p.Results.modelName '_manip.mat'];
            else
                modelPath = [root animal sep animal 'sorted' sep 'stan' sep 'manip' sep p.Results.modelName '_manip'...
                    sep animal pre post '_' p.Results.modelName '_manip.mat'];
            end
            load(modelPath);
            params = paramEsts;

            for run = 1:runs
                if rem(run,100) == 0
                    fprintf('Running simulation %d of %d \n', run, runs);
                end
                rSeed = rSeed + 1;
                switch p.Results.modelName
                    case 'fourParam'
                        [allRewards, allChoices, blockProbs, blockSwitch] = qLearningModel_simNoPlot('params', params(currP,:),...
                            'maxTrials', p.Results.maxTrials, 'randomSeed', rSeed, 'rwdProbs', p.Results.rwdProbs);
                    case 'fiveParam_k'
                        [allRewards, allChoices, blockProbs, blockSwitch] = qLearningModel_k_simNoPlot('params', params(currP,:),...
                            'maxTrials', p.Results.maxTrials, 'randomSeed', rSeed, 'rwdProbs', p.Results.rwdProbs);
                    case 'fiveParamO'
                        [~, allRewards, allChoices, blockProbs, blockSwitch] = qLearningModel_oppo_simNoPlot('params', params(currP,:),...
                            'maxTrials', p.Results.maxTrials, 'randomSeed', rSeed, 'rwdProbs', p.Results.rwdProbs);
                    case 'fiveParam_rBeta_scale'
                        [~, allRewards, allChoices, blockProbs, blockSwitch] = qLearningModel_rBeta_simNoPlot('params', params(currP,:),...
                            'maxTrials', p.Results.maxTrials, 'randomSeed', rSeed, 'rwdProbs', p.Results.rwdProbs);
                    case 'fiveParam_rBeta_kappa'
                        [~, allRewards, allChoices, blockProbs, blockSwitch] = qLearningModel_rBetaKappa_simNoPlot('params', params(currP,:),...
                            'maxTrials', p.Results.maxTrials, 'randomSeed', rSeed, 'rwdProbs', p.Results.rwdProbs);
                    case 'sixParamO_rBarStart'
                        [~, allRewards, allChoices, blockProbs, blockSwitch] = qLearningModel_oppo_rBarStart_simNoPlot('params', params(currP,:),...
                            'maxTrials', p.Results.maxTrials, 'randomSeed', rSeed, 'rwdProbs', p.Results.rwdProbs);
                    case 'sixParam_rBeta_oppo'
                        [~, allRewards, allChoices, blockProbs, blockSwitch] = qLearningModel_rBeta_oppo_simNoPlot('params', params(currP,:),...
                            'maxTrials', p.Results.maxTrials, 'randomSeed', rSeed, 'rwdProbs', p.Results.rwdProbs);
                    case 'sixParam_peBeta'
                        [~, allRewards, allChoices, blockProbs, blockSwitch] = qLearningModel_peBeta_simNoPlot('params', params(currP,:),...
                            'maxTrials', p.Results.maxTrials, 'randomSeed', rSeed, 'rwdProbs', p.Results.rwdProbs);
                    case 'sevenParam_peBeta_k'
                        [~, allRewards, allChoices, blockProbs, blockSwitch] = qLearningModel_peBeta_k_simNoPlot('params', params(currP,:),...
                            'maxTrials', p.Results.maxTrials, 'randomSeed', rSeed, 'rwdProbs', p.Results.rwdProbs);
                    case 'sixParam_pePeBeta'
                        [~, allRewards, allChoices, blockProbs, blockSwitch] = qLearningModel_pePeBeta_simNoPlot('params', params(currP,:),...
                            'maxTrials', p.Results.maxTrials, 'randomSeed', rSeed, 'rwdProbs', p.Results.rwdProbs);
                    case 'sixParam_rAN'
                        [~, allRewards, allChoices, blockProbs, blockSwitch] = qLearningModel_rAN_simNoPlot('params', params(currP,:),...
                            'maxTrials', p.Results.maxTrials, 'randomSeed', rSeed, 'rwdProbs', p.Results.rwdProbs);
                    case 'sixParam_pePeAN'
                        [~, allRewards, allChoices, blockProbs, blockSwitch] = qLearningModel_pePeAN_simNoPlot('params', params(currP,:),...
                            'maxTrials', p.Results.maxTrials, 'randomSeed', rSeed, 'rwdProbs', p.Results.rwdProbs);
                    case 'sixParam_pePeAN_lag'
                        [~, allRewards, allChoices, blockProbs, blockSwitch] = qLearningModel_pePeAN_lag_simNoPlot('params', params(currP,:),...
                            'maxTrials', p.Results.maxTrials, 'randomSeed', rSeed, 'rwdProbs', p.Results.rwdProbs);
                    case 'sixParam_absPePeAN'
                        [~, allRewards, allChoices, blockProbs, blockSwitch] = qLearningModel_absPePeAN_simNoPlot('params', params(currP,:),...
                            'maxTrials', p.Results.maxTrials, 'randomSeed', rSeed, 'rwdProbs', p.Results.rwdProbs);
                    case 'sixParam_absPePeAN_bi'
                        [~, allRewards, allChoices, blockProbs, blockSwitch] = qLearningModel_absPePeAN_bi_simNoPlot('params', params(currP,:),...
                            'maxTrials', p.Results.maxTrials, 'randomSeed', rSeed, 'rwdProbs', p.Results.rwdProbs);
                    case 'sixParam_absPePeAN_biSep'
                        [~, allRewards, allChoices, blockProbs, blockSwitch] = qLearningModel_absPePeAN_biSep_simNoPlot('params', params(currP,:),...
                            'maxTrials', p.Results.maxTrials, 'randomSeed', rSeed, 'rwdProbs', p.Results.rwdProbs);
                    case 'sevenParam_absPePeAN_bi_k'
                        [~, allRewards, allChoices, blockProbs, blockSwitch] = qLearningModel_absPePeAN_bi_k_simNoPlot('params', params(currP,:),...
                            'maxTrials', p.Results.maxTrials, 'randomSeed', rSeed, 'rwdProbs', p.Results.rwdProbs);
                    case 'fiveParam_absPePeAN'
                        [~, allRewards, allChoices, blockProbs, blockSwitch] = qLearningModel_absPePeAN_noMin_simNoPlot('params', params(currP,:),...
                            'maxTrials', p.Results.maxTrials, 'randomSeed', rSeed, 'rwdProbs', p.Results.rwdProbs);
                    case 'sevenParam_absPePeLR'
                        [~, allRewards, allChoices, blockProbs, blockSwitch] = qLearningModel_absPePeLR_simNoPlot('params', params(currP,:),...
                            'maxTrials', p.Results.maxTrials, 'randomSeed', rSeed, 'rwdProbs', p.Results.rwdProbs);
                    case 'eightParam_absPePeLR_k'
                        [~, allRewards, allChoices, blockProbs, blockSwitch] = qLearningModel_absPePeLR_k_simNoPlot('params', params(currP,:),...
                            'maxTrials', p.Results.maxTrials, 'randomSeed', rSeed, 'rwdProbs', p.Results.rwdProbs);
                    case 'sevenParam_pePeAN_k'
                        [~, allRewards, allChoices, blockProbs, blockSwitch] = qLearningModel_pePeAN_k_simNoPlot('params', params(currP,:),...
                            'maxTrials', p.Results.maxTrials, 'randomSeed', rSeed, 'rwdProbs', p.Results.rwdProbs);
                    case 'sevenParam_peLR'
                        [~, allRewards, allChoices, blockProbs, blockSwitch] = qLearningModel_peLR_simNoPlot('params', params(currP,:),...
                            'maxTrials', p.Results.maxTrials, 'randomSeed', rSeed, 'rwdProbs', p.Results.rwdProbs);
                    case 'eightParam_rBeta_pePeAN'
                        [~, allRewards, allChoices, blockProbs, blockSwitch] = qLearningModel_rBeta_pePeAN_simNoPlot('params', params(currP,:),...
                            'maxTrials', p.Results.maxTrials, 'randomSeed', rSeed, 'rwdProbs', p.Results.rwdProbs);
                    case 'eightParam_rBeta_absPePeAN'
                        [~, allRewards, allChoices, blockProbs, blockSwitch] = qLearningModel_rBeta_absPePeAN_simNoPlot('params', params(currP,:),...
                            'maxTrials', p.Results.maxTrials, 'randomSeed', rSeed, 'rwdProbs', p.Results.rwdProbs);
                end

                trialProbs = nan(length(allChoices), 2);
                for j = 2:length(blockSwitch)
                    trialProbs(blockSwitch(j-1):blockSwitch(j)-1, 1) = blockProbs(j-1, 1);
                    trialProbs(blockSwitch(j-1):blockSwitch(j)-1, 2) = blockProbs(j-1, 2);
                end
                trialProbs(blockSwitch(end):length(allChoices), 1) = blockProbs(end, 1);
                trialProbs(blockSwitch(end):length(allChoices), 2) = blockProbs(end, 2);
                prevRewardsBin = [0 abs(allRewards(1:end-1))];
                changeChoice = [0 abs(diff(allChoices)) > 0];

                if currP == 1
                    for j = 2:(length(blockSwitch) - 1)
                        tmpInd = blockSwitch(j);
                        if tmpInd-tranWin > 0 & tmpInd+tranWin <= length(allChoices)
                         %   if trialProbs(tmpInd-1, 2) == pHigh & trialProbs(tmpInd-1, 1) == 10 & trialProbs(tmpInd, 2) == 10 & trialProbs(tmpInd, 1) == pHigh
                            if any(diff(trialProbs(tmpInd-tranWin:tmpInd+tranWin, 2)) == -probDiffH) & any(diff(trialProbs(tmpInd-tranWin:tmpInd+tranWin, 1)) == probDiffH)
                                    if (tmpInd - range - 1) > 0 & length(allChoices) > (tmpInd + range)
                                        transHighPre = [transHighPre; allChoices((tmpInd-range+1):(tmpInd+range))];
                                        changeChoiceMatxHighPre = [changeChoiceMatxHighPre; changeChoice((tmpInd-range+1):(tmpInd+range))];
                                        prevRwdMatxHighPre = [prevRwdMatxHighPre; prevRewardsBin((tmpInd-range+1):(tmpInd+range))];
                                    end
                         %   elseif trialProbs(tmpInd-1,2) == pLow & trialProbs(tmpInd-1, 1) == 10 & trialProbs(tmpInd, 2) == 10 & trialProbs(tmpInd, 1) == pHigh
                            elseif any(diff(trialProbs(tmpInd-tranWin:tmpInd+tranWin, 2)) == -probDiffM) & any(diff(trialProbs(tmpInd-tranWin:tmpInd+tranWin, 1)) == probDiffH)
                                    if (tmpInd - range - 1) > 0 & length(allChoices) > (tmpInd + range)
                                        transLowPre = [transLowPre; allChoices((tmpInd-range+1):(tmpInd+range))];
                                        changeChoiceMatxLowPre = [changeChoiceMatxLowPre; changeChoice((tmpInd-range+1):(tmpInd+range))];
                                        prevRwdMatxLowPre = [prevRwdMatxLowPre; prevRewardsBin((tmpInd-range+1):(tmpInd+range))]; 
                                    end
                         %   elseif trialProbs(tmpInd-1, 1) == pHigh & trialProbs(tmpInd-1, 2) == 10 & trialProbs(tmpInd, 1) == 10 & trialProbs(tmpInd, 2) == pHigh
                            elseif any(diff(trialProbs(tmpInd-tranWin:tmpInd+tranWin, 1)) == -probDiffH) & any(diff(trialProbs(tmpInd-tranWin:tmpInd+tranWin, 2)) == probDiffH)
                                    if (tmpInd - range - 1) > 0 & length(allChoices) > (tmpInd + range)
                                        transHighPre = [transHighPre; (allChoices((tmpInd-range+1):(tmpInd+range))*-1)];
                                        changeChoiceMatxHighPre = [changeChoiceMatxHighPre; changeChoice((tmpInd-range+1):(tmpInd+range))];
                                        prevRwdMatxHighPre = [prevRwdMatxHighPre; prevRewardsBin((tmpInd-range+1):(tmpInd+range))];
                                    end 
                         %   elseif trialProbs(tmpInd-1, 1) == pLow & trialProbs(tmpInd-1, 2) == 10 & trialProbs(tmpInd, 1) == 10 & trialProbs(tmpInd, 2) == pHigh
                            elseif any(diff(trialProbs(tmpInd-tranWin:tmpInd+tranWin, 1)) == -probDiffM) & any(diff(trialProbs(tmpInd-tranWin:tmpInd+tranWin, 2)) == probDiffH)
                                    if (tmpInd - range - 1) > 0 & length(allChoices) > (tmpInd + range)
                                        transLowPre = [transLowPre; (allChoices((tmpInd-range+1):(tmpInd+range))*-1)];
                                        changeChoiceMatxLowPre = [changeChoiceMatxLowPre; changeChoice((tmpInd-range+1):(tmpInd+range))];
                                        prevRwdMatxLowPre = [prevRwdMatxLowPre; prevRewardsBin((tmpInd-range+1):(tmpInd+range))]; 
                                    end
                            end
                        end
                    end
                else
                    for j = 2:(length(blockSwitch) - 1)
                        tmpInd = blockSwitch(j);
                        if tmpInd-tranWin > 0 & tmpInd+tranWin <= length(allChoices)
                         %   if trialProbs(tmpInd-1, 2) == pHigh & trialProbs(tmpInd-1, 1) == 10 & trialProbs(tmpInd, 2) == 10 & trialProbs(tmpInd, 1) == pHigh
                            if any(diff(trialProbs(tmpInd-tranWin:tmpInd+tranWin, 2)) == -probDiffH) & any(diff(trialProbs(tmpInd-tranWin:tmpInd+tranWin, 1)) == probDiffH)
                                    if (tmpInd - range - 1) > 0 & length(allChoices) > (tmpInd + range)
                                        transHighPost = [transHighPost; allChoices((tmpInd-range+1):(tmpInd+range))];
                                        changeChoiceMatxHighPost = [changeChoiceMatxHighPost; changeChoice((tmpInd-range+1):(tmpInd+range))];
                                        prevRwdMatxHighPost = [prevRwdMatxHighPost; prevRewardsBin((tmpInd-range+1):(tmpInd+range))];
                                    end
                         %   elseif trialProbs(tmpInd-1,2) == pLow & trialProbs(tmpInd-1, 1) == 10 & trialProbs(tmpInd, 2) == 10 & trialProbs(tmpInd, 1) == pHigh
                            elseif any(diff(trialProbs(tmpInd-tranWin:tmpInd+tranWin, 2)) == -probDiffM) & any(diff(trialProbs(tmpInd-tranWin:tmpInd+tranWin, 1)) == probDiffH)
                                    if (tmpInd - range - 1) > 0 & length(allChoices) > (tmpInd + range)
                                        transLowPost = [transLowPost; allChoices((tmpInd-range+1):(tmpInd+range))];
                                        changeChoiceMatxLowPost = [changeChoiceMatxLowPost; changeChoice((tmpInd-range+1):(tmpInd+range))];
                                        prevRwdMatxLowPost = [prevRwdMatxLowPost; prevRewardsBin((tmpInd-range+1):(tmpInd+range))]; 
                                    end
                         %   elseif trialProbs(tmpInd-1, 1) == pHigh & trialProbs(tmpInd-1, 2) == 10 & trialProbs(tmpInd, 1) == 10 & trialProbs(tmpInd, 2) == pHigh
                            elseif any(diff(trialProbs(tmpInd-tranWin:tmpInd+tranWin, 1)) == -probDiffH) & any(diff(trialProbs(tmpInd-tranWin:tmpInd+tranWin, 2)) == probDiffH)
                                    if (tmpInd - range - 1) > 0 & length(allChoices) > (tmpInd + range)
                                        transHighPost = [transHighPost; (allChoices((tmpInd-range+1):(tmpInd+range))*-1)];
                                        changeChoiceMatxHighPost = [changeChoiceMatxHighPost; changeChoice((tmpInd-range+1):(tmpInd+range))];
                                        prevRwdMatxHighPost = [prevRwdMatxHighPost; prevRewardsBin((tmpInd-range+1):(tmpInd+range))];
                                    end 
                         %   elseif trialProbs(tmpInd-1, 1) == pLow & trialProbs(tmpInd-1, 2) == 10 & trialProbs(tmpInd, 1) == 10 & trialProbs(tmpInd, 2) == pHigh
                            elseif any(diff(trialProbs(tmpInd-tranWin:tmpInd+tranWin, 1)) == -probDiffM) & any(diff(trialProbs(tmpInd-tranWin:tmpInd+tranWin, 2)) == probDiffH)
                                    if (tmpInd - range - 1) > 0 & length(allChoices) > (tmpInd + range)
                                        transLowPost = [transLowPost; (allChoices((tmpInd-range+1):(tmpInd+range))*-1)];
                                        changeChoiceMatxLowPost = [changeChoiceMatxLowPost; changeChoice((tmpInd-range+1):(tmpInd+range))];
                                        prevRwdMatxLowPost = [prevRwdMatxLowPost; prevRewardsBin((tmpInd-range+1):(tmpInd+range))]; 
                                    end
                            end
                        end
                    end
                end
            end
        end
        prevAnimal = animal;
    end
end

%find averages and sem of actual transition behavior
lowAvgPre_actual = mean(transLowPre_actual,1);
lowSemPre_actual = sem_bernoulli(sum(transLowPre_actual==1), repmat(size(transLowPre_actual,1), 1, range*2));
highAvgPre_actual = mean(transHighPre_actual,1);
highSemPre_actual = sem_bernoulli(sum(transHighPre_actual==1), repmat(size(transHighPre_actual,1), 1, range*2));
lowAvgPost_actual = mean(transLowPost_actual,1);
lowSemPost_actual = sem_bernoulli(sum(transLowPost_actual==1), repmat(size(transLowPost_actual,1), 1, range*2));
highAvgPost_actual = mean(transHighPost_actual,1);
highSemPost_actual = sem_bernoulli(sum(transHighPost_actual==1), repmat(size(transHighPost_actual,1), 1, range*2));


%find averages of simulated behavior
lowAvgPre = mean(transLowPre,1);
highAvgPre = mean(transHighPre,1);
lowAvgPost = mean(transLowPost,1);
highAvgPost = mean(transHighPost,1);

colors = cool(4);

figure;
%plot transition behavior for pre
x = [-range+1:range];
subplot(3,2,1); hold on
plot(x, lowAvgPre_actual, '-', 'Color', colors(1,:), 'linewidth', 2)
plot(x, highAvgPre_actual, '-', 'Color', colors(3,:), 'linewidth', 2)
plot([-range range], [0 0], ':k')
ylim([-1 1])
linetype = 'k';
vline(0, linetype)
set(gca, 'tickdir', 'out')
legend(' medium -> low', 'high -> low')
title('actual pre')

%plot transition behavior for post
subplot(3,2,2); hold on
plot(x, lowAvgPost_actual, '-', 'Color', colors(2,:), 'linewidth', 2)
plot(x, highAvgPost_actual, '-', 'Color', colors(4,:), 'linewidth', 2)
plot([-range range], [0 0], ':k')
ylim([-1 1])
vline(0, linetype)
set(gca, 'tickdir', 'out')
legend('medium -> low', 'high -> low')
title('actual post')

%plot transition behavior for simualted pre
subplot(3,2,3); hold on
plot(x, lowAvgPre, '-', 'Color', colors(1,:), 'linewidth', 2)
plot(x, highAvgPre, '-', 'Color', colors(3,:), 'linewidth', 2)
plot([-range range], [0 0], ':k')
ylim([-1 1])
vline(0, linetype)
set(gca, 'tickdir', 'out')
title('simulated pre')

%plot transition behavior for simualted post
subplot(3,2,4); hold on
plot(x, lowAvgPost, '-', 'Color', colors(2,:), 'linewidth', 2)
plot(x, highAvgPost, '-', 'Color', colors(4,:), 'linewidth', 2)
plot([-range range], [0 0], ':k')
ylim([-1 1])
vline(0, linetype)
set(gca, 'tickdir', 'out')
title('simulated post')

%plot differences in choice probs between conditions for actual and sim pre
subplot(3,2,5); hold on
plot(x,[highAvgPre_actual - lowAvgPre_actual], '-k', 'linewidth', 2)
plot(x,[highAvgPre - lowAvgPre], '-', 'Color', [0.7 0.7 0.7], 'linewidth', 2)
plot([-range range], [0 0], ':k')
xlabel('Trials from switch')
ylabel('Choice average difference')
legend('high - medium', 'sim: high - medium')
set(gca, 'tickdir', 'out')
y = get(gca, 'ylim');
plot([0 0], y, ':k')

%plot differences in choice probs between conditions for actual and sim post
subplot(3,2,6); hold on
plot(x,[highAvgPost_actual - lowAvgPost_actual], '-k', 'linewidth', 2)
plot(x,[highAvgPost - lowAvgPost], '-', 'Color', [0.7 0.7 0.7], 'linewidth', 2)
plot([-range range], [0 0], ':k')
xlabel('Trials from switch')
ylabel('Choice average difference')
set(gca, 'tickdir', 'out')
y = get(gca, 'ylim');
plot([0 0], y, ':k')


set(gcf, 'renderer', 'painters', 'position', [-1467 42 1024 954])
titleTxt = strrep([sheet ' ' p.Results.modelName ' transition behavior'], '_', ' ');
suptitle(titleTxt)


%look at win-stay lose-shift around transitions
for tInd = 1:range*2
    lS_low_pre(tInd) = sum(changeChoiceMatxLowPre(find(prevRwdMatxLowPre(:,tInd)==0), tInd))/sum(prevRwdMatxLowPre(:,tInd)==0);
    sem_lS_low_pre(tInd) = sem_bernoulli(sum(changeChoiceMatxLowPre(find(prevRwdMatxLowPre(:,tInd)==0), tInd)), sum(prevRwdMatxLowPre(:,tInd)==0));
    lS_high_pre(tInd) = sum(changeChoiceMatxHighPre(find(prevRwdMatxHighPre(:,tInd)==0), tInd))/sum(prevRwdMatxHighPre(:,tInd)==0);
    sem_lS_high_pre(tInd) = sem_bernoulli(sum(changeChoiceMatxHighPre(find(prevRwdMatxHighPre(:,tInd)==0), tInd)), sum(prevRwdMatxHighPre(:,tInd)==0));

    wS_low_pre(tInd) = 1 - ((sum(changeChoiceMatxLowPre(find(prevRwdMatxLowPre(:,tInd)==1), tInd)))/sum(prevRwdMatxLowPre(:,tInd)==1));
    sem_wS_low_pre(tInd) = sem_bernoulli(sum(~changeChoiceMatxLowPre(find(prevRwdMatxLowPre(:,tInd)==1), tInd)), sum(prevRwdMatxLowPre(:,tInd)==1));
    wS_high_pre(tInd) = 1 - ((sum(changeChoiceMatxHighPre(find(prevRwdMatxHighPre(:,tInd)==1), tInd)))/sum(prevRwdMatxHighPre(:,tInd)==1));
    sem_wS_high_pre(tInd) = sem_bernoulli(sum(~changeChoiceMatxHighPre(find(prevRwdMatxHighPre(:,tInd)==1), tInd)), sum(prevRwdMatxHighPre(:,tInd)==1));

    lS_low_post(tInd) = sum(changeChoiceMatxLowPost(find(prevRwdMatxLowPost(:,tInd)==0), tInd))/sum(prevRwdMatxLowPost(:,tInd)==0);
    sem_lS_low_post(tInd) = sem_bernoulli(sum(changeChoiceMatxLowPost(find(prevRwdMatxLowPost(:,tInd)==0), tInd)), sum(prevRwdMatxLowPost(:,tInd)==0));
    lS_high_post(tInd) = sum(changeChoiceMatxHighPost(find(prevRwdMatxHighPost(:,tInd)==0), tInd))/sum(prevRwdMatxHighPost(:,tInd)==0);
    sem_lS_high_post(tInd) = sem_bernoulli(sum(changeChoiceMatxHighPost(find(prevRwdMatxHighPost(:,tInd)==0), tInd)), sum(prevRwdMatxHighPost(:,tInd)==0));

    wS_low_post(tInd) = 1 - ((sum(changeChoiceMatxLowPost(find(prevRwdMatxLowPost(:,tInd)==1), tInd)))/sum(prevRwdMatxLowPost(:,tInd)==1));
    sem_wS_low_post(tInd) = sem_bernoulli(sum(~changeChoiceMatxLowPost(find(prevRwdMatxLowPost(:,tInd)==1), tInd)), sum(prevRwdMatxLowPost(:,tInd)==1));
    wS_high_post(tInd) = 1 - ((sum(changeChoiceMatxHighPost(find(prevRwdMatxHighPost(:,tInd)==1), tInd)))/sum(prevRwdMatxHighPost(:,tInd)==1));
    sem_wS_high_post(tInd) = sem_bernoulli(sum(~changeChoiceMatxHighPost(find(prevRwdMatxHighPost(:,tInd)==1), tInd)), sum(prevRwdMatxHighPost(:,tInd)==1));
end


figure;

%actual win-stay pre
subplot(2,4,1); hold on
errorbar(x, wslsPre.wS_low, wslsPre.sem_wS_low, 'Color', colors(1,:), 'linewidth', 1.3)
errorbar(x, wslsPre.wS_high, wslsPre.sem_wS_high, 'Color', colors(3,:), 'linewidth', 1.3)
ylabel('win-stay')
yws = ylim;
plot([0 0], yws, '--k')
set(gca, 'tickdir', 'out')

%actual lose-shift pre
subplot(2,4,2); hold on
errorbar(x, wslsPre.lS_low, wslsPre.sem_lS_low, 'Color', colors(1,:), 'linewidth', 1.3)
errorbar(x, wslsPre.lS_high, wslsPre.sem_lS_high, 'Color', colors(3,:), 'linewidth', 1.3)
ylabel('lose-shift')
yls = ylim;
plot([0 0], yws, '--k')
set(gca, 'tickdir', 'out')

%actual win-stay post
subplot(2,4,3); hold on
errorbar(x, wslsPost.wS_low, wslsPost.sem_wS_low, 'Color', colors(1,:), 'linewidth', 1.3)
errorbar(x, wslsPost.wS_high, wslsPost.sem_wS_high, 'Color', colors(3,:), 'linewidth', 1.3)
ylabel('win-stay')
ylim([yws])
plot([0 0], yws, '--k')
set(gca, 'tickdir', 'out')

%actual lose-shift post
subplot(2,4,4); hold on
errorbar(x, wslsPost.lS_low, wslsPost.sem_lS_low, 'Color', colors(1,:), 'linewidth', 1.3)
errorbar(x, wslsPost.lS_high, wslsPost.sem_lS_high, 'Color', colors(3,:), 'linewidth', 1.3)
ylabel('lose-shift')
ylim([yls])
plot([0 0], yls, '--k')
set(gca, 'tickdir', 'out')

%simulated win-stay pre
subplot(2,4,5); hold on
errorbar(x, wS_low_pre, sem_wS_low_pre, 'Color', colors(1,:), 'linewidth', 1.3)
errorbar(x, wS_high_pre, sem_wS_high_pre, 'Color', colors(3,:), 'linewidth', 1.3)
xlabel('Trials from switch')
ylabel('win-stay')
ylim([yws])
plot([0 0], yws, '--k')
set(gca, 'tickdir', 'out')

%simulated lose-shift pre
subplot(2,4,6); hold on
errorbar(x, lS_low_pre, sem_lS_low_pre, 'Color', colors(1,:), 'linewidth', 1.3)
errorbar(x, lS_high_pre, sem_lS_low_pre, 'Color', colors(3,:), 'linewidth', 1.3)
xlabel('Trials from switch')
ylabel('lose-shift')
ylim([yls])
plot([0 0], yls, '--k')
set(gca, 'tickdir', 'out')

%simulated win-stay post
subplot(2,4,7); hold on
errorbar(x, wS_low_post, sem_wS_low_post, 'Color', colors(2,:), 'linewidth', 1.3)
errorbar(x, wS_high_post, sem_wS_high_post, 'Color', colors(4,:), 'linewidth', 1.3)
ylabel('Win-stay')
ylim([yws])
plot([0 0], yws, '--k')
set(gca, 'tickdir', 'out')

%simulated lose-shift post
subplot(2,4,8); hold on
errorbar(x, lS_low_post, sem_lS_low_post, 'Color', colors(2,:), 'linewidth', 1.3)
errorbar(x, lS_high_post, sem_lS_low_post, 'Color', colors(4,:), 'linewidth', 1.3)
ylabel('Lose-shift')
xlabel('Trials from switch')
ylim([yls])
plot([0 0], yls, '--k')
set(gca, 'tickdir', 'out')


set(gcf, 'renderer', 'painters', 'position', [-1467 42 1024 954])
suptitle(titleTxt)

end

