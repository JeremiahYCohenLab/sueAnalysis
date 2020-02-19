function compareWsLsSim_opMD(xlFile, sheet, category, varargin)

p = inputParser;
% default parameters if none given
p.addParameter('rwdProbs', [90 50 10]);
p.addParameter('params', []);
p.addParameter('maxTrials', 350);
p.addParameter('runs', []);
p.addParameter('randomSeed', 77582);
p.addParameter('modelName', [{'fiveParam_rBeta_scale'} {'fourParam'} {'fiveParamO'}]);
p.addParameter('taskType', 'decoupled');
p.addParameter('revForFlag', 0);
p.addParameter('tMax', 10);
p.parse(varargin{:});

numMdls = length(p.Results.modelName);
tMax = p.Results.tMax;
pHigh = p.Results.rwdProbs(1);
pLow = p.Results.rwdProbs(2);

if isempty(p.Results.params)
    [root, sep] = currComputer();
    for i =1:numMdls
        modelPath = [root sheet sep sheet 'sorted' sep 'stan' sep p.Results.modelName{i} sep sheet...
                    category '_' p.Results.modelName{i} '.mat'];
        t = generateStanModelTerms_opMD(p.Results.modelName, modelPath, [], 0, p.Results.revForFlag);
        params{i} = t.params;
    end
else
    params = p.Results.params;
end

if regexp(p.Results.taskType, 'switch')
    [wS_rwdHx, lS_rwdHx, wS, lS] = wslsRwdCount_opMD(xlFile, sheet, 'switch', 'tMax', p.Results.tMax, 'revForFlag', p.Results.revForFlag);
else
    [wS_rwdHx, lS_rwdHx, wS, lS] = wslsRwdCount_opMD(xlFile, sheet, category, 'tMax', p.Results.tMax, 'revForFlag', p.Results.revForFlag);
end

if isempty(p.Results.runs)
    [~, dayList, ~] = xlsread(xlFile, sheet);
    [~,col] = find(~cellfun(@isempty,strfind(dayList, category)) == 1);
    runs = length(dayList(2:end,col));
else
    runs = p.Results.runs;
end

combinedRwdHx = cell(1,numMdls);
combinedRwds = cell(1,numMdls);
combinedChangeChoice = cell(1,numMdls);
rSeed = p.Results.randomSeed;

for i = 1:runs
    fprintf('Running simulation %d of %d \n', i, runs);
    rSeed = rSeed + 1;
    
    for j = 1:numMdls
        switch p.Results.modelName{j}
            case 'fourParam'
                [allRewards, allChoices, blockProbs, blockSwitch] = qLearningModel_simNoPlot('params', params{j},...
                    'maxTrials', p.Results.maxTrials, 'randomSeed', rSeed, 'rwdProbs', p.Results.rwdProbs, 'taskType', p.Results.taskType);
            case 'fiveParamO'
                [~, allRewards, allChoices, blockProbs, blockSwitch] = qLearningModel_oppo_simNoPlot('params', params{j},...
                    'maxTrials', p.Results.maxTrials, 'randomSeed', rSeed, 'rwdProbs', p.Results.rwdProbs, 'taskType', p.Results.taskType);
            case 'fiveParam_rBeta_scale'
                [~, allRewards, allChoices, blockProbs, blockSwitch] = qLearningModel_rBeta_simNoPlot('params', params{j},...
                    'maxTrials', p.Results.maxTrials, 'randomSeed', rSeed, 'rwdProbs', p.Results.rwdProbs, 'taskType', p.Results.taskType);
            case 'fiveParam_rBeta_kappa'
                [~, allRewards, allChoices, blockProbs, blockSwitch] = qLearningModel_rBetaKappa_simNoPlot('params', params{j},...
                    'maxTrials', p.Results.maxTrials, 'randomSeed', rSeed, 'rwdProbs', p.Results.rwdProbs, 'taskType', p.Results.taskType);
        end


        rwdMatx = [];
        for k=1:tMax
            rwdMatx(k,:) = [nan(1,k) allRewards(1:end-k)];
        end
        rwdHx = nansum(rwdMatx, 1);
        rwdHx(1) = 0;

        combinedRwdHx{j} = [combinedRwdHx{j} rwdHx(1:end-2)];
        combinedRwds{j} = [combinedRwds{j} allRewards(2:end-1)];
        changeChoice = [abs(diff(allChoices)) > 0];
        combinedChangeChoice{j} = [combinedChangeChoice{j} changeChoice(2:end)];
    end
    
        
end



for i = 1:numMdls
    wS_sim(i) = 1 - (sum(combinedChangeChoice{i}(combinedRwds{i}==1))/sum(combinedRwds{i}==1)); 
    lS_sim(i) = sum(combinedChangeChoice{i}(combinedRwds{i}==0))/sum(combinedRwds{i}==0);
    for j = 1:tMax+1
        tmpInds = logical(combinedRwdHx{i} == j-1);
        wS_rwdHx_sim(i,j) = 1 - (sum(combinedChangeChoice{i}(combinedRwds{i}==1 & tmpInds))/sum(combinedRwds{i}==1 & tmpInds));
        lS_rwdHx_sim(i,j) = sum(combinedChangeChoice{i}(combinedRwds{i}==0 & tmpInds))/sum(combinedRwds{i}==0 & tmpInds);
    end 
end


blue = [0 1 1];
purp = [0.7 0 1];
colors = [linspace(blue(1),purp(1),numMdls+1)', linspace(blue(2),purp(2),numMdls+1)',...
    linspace(blue(3),purp(3),numMdls+1)'];

figure;
subplot(2,2,1); hold on;
plot([0:tMax], wS_rwdHx, 'linewidth', 2, 'Color', colors(1,:))
for i = 1:numMdls
    plot([0:tMax], wS_rwdHx_sim(i,:), 'linewidth', 2, 'Color', colors(i+1,:))
end
xlim([-0.5 tMax+0.5])
ylabel('probability')
title('win-stay')
set(gca, 'tickdir', 'out')

subplot(2,2,2); hold on;
plot([0:tMax], lS_rwdHx, 'linewidth', 2, 'Color', colors(1,:))
legTmp = {'actual'};
for i = 1:numMdls
    plot([0:tMax], lS_rwdHx_sim(i,:), 'linewidth', 2, 'Color', colors(i+1,:))
    legTmp = [legTmp {p.Results.modelName{i}}];
end
xlim([-0.5 tMax+0.5])
ylabel('probability')
title('lose-shift')
set(gca, 'tickdir', 'out')
legend(legTmp)

subplot(2,2,3); hold on;
scatter(0, wS, 200, colors(1,:),  'filled')
tickLbls = {'actual'};
for i = 1:numMdls
    scatter(i, wS_sim(i), 200, colors(i+1,:),  'filled')
    tickLbls = [tickLbls {p.Results.modelName{i}}];
end
tickLbls = strrep(tickLbls, '_', ' ');
xlim([-0.5 numMdls+0.5])
ylabel('win-stay')
xticks([0:numMdls])
xticklabels(tickLbls)
set(gca, 'tickdir', 'out')

subplot(2,2,4); hold on;
scatter(0, lS, 200, colors(1,:),  'filled')
tickLbls = {'actual'};
for i = 1:numMdls
    scatter(i, lS_sim(i), 200, colors(i+1,:),  'filled')
    tickLbls = [tickLbls {p.Results.modelName{i}}];
end
tickLbls = strrep(tickLbls, '_', ' ');
xlim([-0.5 numMdls+0.5])
ylabel('lose-shift')
xticks([0:numMdls])
xticklabels(tickLbls)
set(gca, 'tickdir', 'out')


end

