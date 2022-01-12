%% settings
xlFile = 'inhibitionAll.xlsx';
sheet = 'allinhibition';
col = 'inhibition';
modelName = '5params';
[root, sep] = currComputer();
dayList = getDayList(xlFile, sheet, col);
%% without model
len = 2:5;
lickLats = cell(1,length(len));
lickRates = cell(1,length(len));
lickRatesRwd = cell(1,length(len));
lickRatesNoRwd = cell(1,length(len));
lcChoice = cell(1,length(len));
wsChoice = cell(1,length(len));
lickLatsAfterNoRwd = cell(1,length(len));
% 
for sess = 1:length(dayList)
    os = behAnalysisNoPlot_opMD(dayList{sess},'simpleFlag',1);
    switches = zeros(1, length(os.allChoices));
    switches(os.changeChoice_Inds) = 1;
    lickInds = cell(1,length(len));
    lickRwdInds = cell(1, length(len));
    lickNoRwdInds = cell(1, length(len));
    for i = 1:length(len)
        % generate choice and reward history
        choicesHis = conv(os.allChoices, ones(1,len(i)));
        choicesHis = choicesHis(1:end-(len(i)-1));
        rewardsHis = conv(os.allRewards, ones(1,len(i)-1));
        rewardsHis = rewardsHis(1:end-(len(i)-2));
        noRewardsHis = conv(os.allNoRewards, ones(1,len(i)-1));
        noRewardsHis = noRewardsHis(1:end-(len(i)-2));
        % detect len consecutive choices and len-1 consecutive rewards
        conChoicesInds = find(abs(choicesHis)>=len(i));
        conPreRewardInds = find(abs(rewardsHis)>=(len(i)-1))+1;
        lickInds{i} = intersect(conChoicesInds, conPreRewardInds); 
        lickNoRwdInds{i} = intersect(conChoicesInds, intersect(conPreRewardInds, os.nrwd_Inds)); 
        lickRwdInds{i} = intersect(conChoicesInds, intersect(conPreRewardInds, os.rwd_Inds)); 
    end
    for i = 1:(length(len)-1)
        lickInds{i} = setxor(lickInds{i},lickInds{i+1}); % take longer ones away from shorter groups
        lickRwdInds{i} = setxor(lickRwdInds{i},lickRwdInds{i+1});
        lickNoRwdInds{i} = setxor(lickNoRwdInds{i},lickNoRwdInds{i+1});
    end
    for i = 1:length(len)
        lickLats{i} = [lickLats{i} os.lickLatLogZ(lickInds{i})]; % take lickLatz of certain groups
        lickRates{i} = [lickRates{i} os.lickRateZ(lickInds{i})]; % take lick rate of certain groups
        lickRatesNoRwd{i} = [lickRatesNoRwd{i} os.lickRateRwdZ(lickNoRwdInds{i})]; % take cons lick Rate of certain groups
        lickRatesRwd{i} = [lickRatesRwd{i} os.lickRateRwdZ(lickRwdInds{i})]; % take cons lick Rate of certain groups
        tempInds = lickNoRwdInds{i};
        tempInds = tempInds(tempInds<length(os.allChoices)); % take away the last trial
        lcChoice{i} = [lcChoice{i} switches(tempInds+1)]; % take next choice
        lickLatsAfterNoRwd{i} = [lickLatsAfterNoRwd{i} os.lickLatLogZ(tempInds+1)];
        tempInds = lickRwdInds{i};
        tempInds = tempInds(tempInds<length(os.allChoices)); % take away the last trial
        wsChoice{i} = [wsChoice{i} 1-switches(tempInds+1)]; % take next choice
    end
end

meanLats = cellfun(@(x) mean(x, 'omitnan'), lickLats);
meanRates = cellfun(@(x) mean(x, 'omitnan'), lickRates);
meanRatesNoRwd = cellfun(@(x) mean(x, 'omitnan'), lickRatesNoRwd);
meanRatesRwd = cellfun(@(x) mean(x, 'omitnan'), lickRatesRwd);
meanLc = cellfun(@(x) mean(x, 'omitnan'), lcChoice);
meanWs = cellfun(@(x) mean(x, 'omitnan'), wsChoice);

semLats = cellfun(@(x) sem(x), lickLats);
semRates = cellfun(@(x) sem(x), lickRates);
semRatesNoRwd = cellfun(@(x) sem(x), lickRatesNoRwd);
semRatesRwd = cellfun(@(x) sem(x), lickRatesRwd);
semDiff = sqrt(semRatesNoRwd.^2 + semRatesRwd.^2);
semLc = cellfun(@(x) sem_bern(x), lcChoice);
semWs = cellfun(@(x) sem_bern(x), wsChoice);

meanLatsNoRwdStay = zeros(1,length(len));
semLatsNoRwdsStay = zeros(1,length(len));
meanLatsNoRwdSwitch = zeros(1,length(len));
semLatsNoRwdsSwitch = zeros(1,length(len));
for i = 1:length(len)
    meanLatsNoRwdStay(i) = mean(lickLatsAfterNoRwd{i}(lcChoice{i}<1), 'omitnan');
    semLatsNoRwdsStay(i) = sem(lickLatsAfterNoRwd{i}(lcChoice{i}<1));
    meanLatsNoRwdSwitch(i) = mean(lickLatsAfterNoRwd{i}(lcChoice{i}>0), 'omitnan');
    semLatsNoRwdsSwitch(i) = sem(lickLatsAfterNoRwd{i}(lcChoice{i}>0));
end

%%
figure;
suptitle([sheet '  ' col])
subplot(4,6,1); hold on;
plot(len, meanLats, 'c', 'lineWidth', 2);
patch([len flip(len)], [meanLats + semLats, flip(meanLats - semLats)], 'c', 'FaceAlpha', 0.4, 'EdgeColor', 'none');
title('lickLat')

subplot(4,6,2); hold on;
plot(len, meanRates, 'color', [0.5 0.5 1], 'lineWidth', 2);
patch([len flip(len)], [meanRates + semRates, flip(meanRates - semRates)], [0.5 0.5 1], 'FaceAlpha', 0.4, 'EdgeColor', 'none');
title('lickRate')

subplot(4,6,3); hold on;
plot(len, meanRatesRwd, 'm', 'lineWidth', 2);
patch([len flip(len)], [meanRatesRwd + semRatesRwd, flip(meanRatesRwd - semRatesRwd)], 'm', 'FaceAlpha', 0.4, 'EdgeColor', 'none');
title('lickRateRwd')

subplot(4,6,4); hold on;
plot(len, meanRatesNoRwd, 'b', 'lineWidth', 2);
patch([len flip(len)], [meanRatesNoRwd + semRatesNoRwd, flip(meanRatesNoRwd - semRatesNoRwd)], 'b', 'FaceAlpha', 0.4, 'EdgeColor', 'none');
title('lickRateNoRwd')

subplot(4,6,5); hold on;
plot(len, meanRatesRwd - meanRatesNoRwd, 'r', 'lineWidth', 2);
patch([len flip(len)], [meanRatesRwd - meanRatesNoRwd + semDiff, flip(meanRatesRwd - meanRatesNoRwd - semDiff)], 'r', 'FaceAlpha', 0.4, 'EdgeColor', 'none');
title('lickRateRwd-NoRwd')

subplot(4,6,[6 12]); hold on;
plot(len, meanLatsNoRwdStay, 'b', 'lineWidth', 2);
patch([len flip(len)], [meanLatsNoRwdStay + semLatsNoRwdsStay, flip(meanLatsNoRwdStay - semLatsNoRwdsStay)], 'b', 'FaceAlpha', 0.4, 'EdgeColor', 'none');
plot(len, meanLatsNoRwdSwitch, 'm', 'lineWidth', 2);
patch([len flip(len)], [meanLatsNoRwdSwitch + semLatsNoRwdsSwitch, flip(meanLatsNoRwdSwitch - semLatsNoRwdsSwitch)], 'm', 'FaceAlpha', 0.4, 'EdgeColor', 'none');

plot(len, meanLats, 'color', [0.5 0.5 0.5], 'lineWidth', 2);
patch([len flip(len)], [meanLats + semLats, flip(meanLats - semLats)], [0.5 0.5 0.5], 'FaceAlpha', 0.4, 'EdgeColor', 'none');
title('lickLatsAfterNoRwd')

subplot(4,6,7); hold on;
plot(len, meanLc, 'b', 'lineWidth', 2);
patch([len flip(len)], [meanLc + semLc, flip(meanLc - semLc)], 'b', 'FaceAlpha', 0.4, 'EdgeColor', 'none');
title('lc')

subplot(4,6,8); hold on;
plot(len, meanWs, 'r', 'lineWidth', 2);
patch([len flip(len)], [meanWs + semWs, flip(meanWs - semWs)], 'r', 'FaceAlpha', 0.4, 'EdgeColor', 'none');
title('ws')
%% with model
numSamps = 200;
combinePe = [];
combineQc = [];
combineConf = [];
combineLickLat = [];
combineLickRates = [];
combineLickRatesRwd = [];
combineSw = [];

for sess = 1:length(dayList)
    os = behAnalysisNoPlot_opMD(dayList{sess},'simpleFlag',1);
    switches = zeros(1, length(os.allChoices));
    switches(os.changeChoice_Inds) = 1;
    [animalName, date] = strtok(dayList{sess}, 'd'); 
    animalName = animalName(2:end);
    date = date(1:9);
    [params, modelName, ~, noSession] = getStanModelParams_sampsOnly(animalName, col, modelName, numSamps, 'sessionParamsFlag', 1, 'sessionName', dayList{sess});
    if noSession
        fprintf([dayList{sess} ' no good behavior \n'])
        continue
    end
    t = inferModelVar(dayList{sess}, params, modelName);
    combinePe = [combinePe; t.pe];
    combineConf = [combineConf; 2*t.probChoice-1];
    qChosen = zeros(length(os.allChoices),1);
    qChosen(os.allChoices>0) = t.Q(os.allChoices>0,2);
    qChosen(os.allChoices<0) = t.Q(os.allChoices<0,1);
    combineQc = [combineQc; zscore(qChosen)];
    combineLickLat = [combineLickLat, os.lickLatZ];
    combineLickRates = [combineLickRates, os.lickRateZ];
    combineLickRatesRwd = [combineLickRatesRwd, os.lickRateRwdZ];
    combineSw = [combineSw, switches];
end
%% bin pe, plot lickRate after stay
numBins = 20;
combinePeStay = combinePe(combineSw<1);
combineLickRatesRwdStay = combineLickRatesRwd(combineSw<1);
edges = unique([linspace(min(combinePeStay), 0, 0.5*numBins+1), linspace(0, max(combinePeStay), 0.5*numBins+1)]);
rateMeans = zeros(numBins,1);
rateSems = zeros(numBins,1);
peMeans = zeros(numBins,1);
%     spikeMeansLate = zeros(numBins,1);
%     spikeSemsLate = zeros(numBins,1);
for k = 1:numBins
    if k < numBins
        rateTemp = combineLickRatesRwdStay(combinePeStay >= edges(k) & combinePeStay < edges(k+1));
        peMeans(k) = mean(combinePeStay(combinePeStay >= edges(k) & combinePeStay < edges(k+1)));
    else
        rateTemp = combineLickRatesRwdStay(combinePeStay >= edges(k) & combinePeStay <= edges(k+1));
        peMeans(k) = mean(combinePeStay(combinePeStay >= edges(k) & combinePeStay <= edges(k+1)));
    end
    rateMeans(k) = mean(rateTemp, 'omitnan');
    rateSems(k) = sem(rateTemp);
end
%%
subplot(4,6,13); hold on;
plot(peMeans(peMeans<0), rateMeans(peMeans<0), 'b', 'lineWidth', 2);
patch([peMeans(peMeans<0); flip(peMeans(peMeans<0))], [rateMeans(peMeans<0) + rateSems(peMeans<0); flip(rateMeans(peMeans<0) - rateSems(peMeans<0))], 'b', 'FaceAlpha', 0.4, 'EdgeColor', 'none');
title('lickRateNoRwd-pe')

subplot(4,6,14); hold on;
plot(peMeans(peMeans>0), rateMeans(peMeans>0), 'm', 'lineWidth', 2);
patch([peMeans(peMeans>0); flip(peMeans(peMeans>0))], [rateMeans(peMeans>0) + rateSems(peMeans>0); flip(rateMeans(peMeans>0) - rateSems(peMeans>0))], 'm', 'FaceAlpha', 0.4, 'EdgeColor', 'none');
title('lickRateRwd-pe')

%% bin pe, plot lickRate after switch
numBins = 12;
combinePeSwitch = combinePe(combineSw>0);
combineLickRatesRwdSwitch = combineLickRatesRwd(combineSw>0);
edges = unique([linspace(min(combinePeSwitch), 0, 0.5*numBins+1), linspace(0, max(combinePeSwitch), 0.5*numBins+1)]);
rateMeans = zeros(numBins,1);
rateSems = zeros(numBins,1);
peMeans = zeros(numBins,1);
%     spikeMeansLate = zeros(numBins,1);
%     spikeSemsLate = zeros(numBins,1);
for k = 1:numBins
    if k < numBins
        rateTemp = combineLickRatesRwdSwitch(combinePeSwitch >= edges(k) & combinePeSwitch < edges(k+1));
        peMeans(k) = mean(combinePeSwitch(combinePeSwitch >= edges(k) & combinePeSwitch < edges(k+1)));
    else
        rateTemp = combineLickRatesRwdSwitch(combinePeSwitch >= edges(k) & combinePeSwitch <= edges(k+1));
        peMeans(k) = mean(combinePeSwitch(combinePeSwitch >= edges(k) & combinePeSwitch <= edges(k+1)));
    end
    rateMeans(k) = mean(rateTemp, 'omitnan');
    rateSems(k) = sem(rateTemp);
end
%%
subplot(4,6,15); hold on;
plot(peMeans(peMeans<0), rateMeans(peMeans<0), 'b', 'lineWidth', 2);
patch([peMeans(peMeans<0); flip(peMeans(peMeans<0))], [rateMeans(peMeans<0) + rateSems(peMeans<0); flip(rateMeans(peMeans<0) - rateSems(peMeans<0))], 'b', 'FaceAlpha', 0.4, 'EdgeColor', 'none');
title('lickRateNoRwd-pe')

subplot(4,6,16); hold on;
plot(peMeans(peMeans>0), rateMeans(peMeans>0), 'm', 'lineWidth', 2);
patch([peMeans(peMeans>0); flip(peMeans(peMeans>0))], [rateMeans(peMeans>0) + rateSems(peMeans>0); flip(rateMeans(peMeans>0) - rateSems(peMeans>0))], 'm', 'FaceAlpha', 0.4, 'EdgeColor', 'none');
title('lickRateRwd-pe')
%% bin choiceConf, plot lickLat and lickRate before

%% bin Qc, plot lickLat and lickRate before

%%


