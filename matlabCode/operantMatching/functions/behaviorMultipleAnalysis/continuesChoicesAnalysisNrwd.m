function continuesChoicesAnalysisNrwd(xlFile, sheet, col, varargin)
%% settings
p = inputParser;
% default parameters if none given
p.addParameter('modelName', '5params');
p.addParameter('len', 2:4);
p.addParameter('followLen', 15);
p.addParameter('trialBins',4);
p.addParameter('numBins',3);
p.parse(varargin{:});
modelName = '5params';
dayList = getDayList(xlFile, sheet, col);
%% without model
len = p.Results.len;
lickLats = cell(1,length(len));
lickRates = cell(1,length(len));
lickRatesRwd = cell(1,length(len));
lickRatesNoRwd = cell(1,length(len));
lcChoice = cell(1,length(len));
wsChoice = cell(1,length(len));
lickLatsAfterNoRwd = cell(1,length(len));
followChoices = cell(1,length(len));
followedKernel = cell(1,length(len));
followChoicesRef = cell(1,length(len));
followedKernelRef = cell(1,length(len));
distribution = cell(1,length(len));
swKernel = [];
kernelProp = [];
ws = [];
lc = [];
rwdDiff = [];
rwdSumDiff = [];
noRwdSumDiff = [];
choiceDiff = [];
wsDiff = [];
lcDiff = [];

% choice kernel 
myKernel = [1 2 4 2 1];
randChoice = rand(1,10000);
randChoice = double(randChoice>0.5);
randChoice = 2*randChoice-1;
randChoice = conv(randChoice, myKernel);
randChoice = abs(randChoice(0.5*length(myKernel)-0.5+1:end-(0.5*length(myKernel)-0.5)))/sum(myKernel);
kernelUniq = unique(randChoice(0.5*length(myKernel)+0.5:end-(0.5*length(myKernel)-0.5)));

for sess = 1:length(dayList)
    os = behAnalysisNoPlot_opMD(dayList{sess},'simpleFlag',1);
    switches = zeros(1, length(os.allChoices));
    switches(os.changeChoice_Inds) = 1;
    %% continued Choices
    lickInds = cell(1,length(len));
    refInds = cell(1, length(len));
    lickRwdInds = cell(1, length(len));
    lickNoRwdInds = cell(1, length(len));
%     myKernel(0.5*length(myKernel)+0.5) = 0;
    choiceKernel = conv(os.allChoices, myKernel);
    choiceKernel = abs(choiceKernel(0.5*length(myKernel)-0.5+1:end-(0.5*length(myKernel)-0.5)))/sum(myKernel);
%     figure;
%     histogram(choiceKernel);
    
    for i = 1:length(len)
        % generate choice and reward history
        choicesHis = conv(os.allChoices, ones(1,len(i)));
        choicesHis = choicesHis(1:end-(len(i)-1));
        rewardsHis = conv(os.allRewards, ones(1,len(i)-1));
        rewardsHis = rewardsHis(1:end-(len(i)-2));
        noRewardsHis = conv(os.allNoRewards, ones(1,len(i)-1));
        noRewardsHis = noRewardsHis(1:end-(len(i)-2));
        % detect len consecutive choices and len-1 consecutive no rewards
        conChoicesInds = find(abs(choicesHis)>=len(i));
        conPreRewardInds = find(abs(rewardsHis)>=(len(i)-1))+1;
        conPreNoRewardInds = find(abs(noRewardsHis)>=(len(i)-1))+1;
        lickInds{i} = intersect(conChoicesInds, conPreNoRewardInds); 
        lickNoRwdInds{i} = intersect(conChoicesInds, intersect(conPreNoRewardInds, os.nrwd_Inds)); 
        lickRwdInds{i} = intersect(conChoicesInds, intersect(conPreNoRewardInds, os.rwd_Inds)); 
         % take all stay trials as ref
        refInds{i} = find(abs(choicesHis)>=len(i));
    end
    for i = 1:(length(len)-1)
        lickInds{i} = setxor(lickInds{i},lickInds{i+1}); % take longer ones away from shorter groups
        lickRwdInds{i} = setxor(lickRwdInds{i},lickRwdInds{i+1});
        lickNoRwdInds{i} = setxor(lickNoRwdInds{i},lickNoRwdInds{i+1});
        refInds{i} = setxor(refInds{i+1}, refInds{i});
    end
    for i = 1:length(len)
        distribution{i} = [distribution{i} length(lickInds{i})/length(os.allChoices)];
        lickLats{i} = [lickLats{i} os.lickLatLogZ(lickInds{i})]; % take lickLatz of certain groups
        lickRates{i} = [lickRates{i} os.lickRateZ(lickInds{i})]; % take lick rate of certain groups
        lickRatesNoRwd{i} = [lickRatesNoRwd{i} os.lickRateRwdZ(lickNoRwdInds{i})]; % take cons lick Rate of certain groups
        lickRatesRwd{i} = [lickRatesRwd{i} os.lickRateRwdZ(lickRwdInds{i})]; % take cons lick Rate of certain groups
        % focus on nrwd on len(j) trial
        tempInds = lickNoRwdInds{i};
        tempInds = tempInds(tempInds<length(os.allChoices)); % take away the last trial
        lcChoice{i} = [lcChoice{i} switches(tempInds+1)]; % take next choice
        lickLatsAfterNoRwd{i} = [lickLatsAfterNoRwd{i} os.lickLatLogZ(tempInds+1)];
        % take all following trials
        tempFollowChoices = nan(length(tempInds), p.Results.followLen);
        tempFollowChoiceKernel = nan(length(tempInds), p.Results.followLen);
        for j = 1:length(tempInds)
            blEnd = min(length(os.allChoices), tempInds(j)+p.Results.followLen);
            tempFollowChoices(j,1:blEnd-tempInds(j)) = os.allChoices(tempInds(j)+1:blEnd);
            tempFollowChoices(j, tempFollowChoices(j,:)<0) = 0;
            tempFollowChoiceKernel(j,1:blEnd-tempInds(j)) = choiceKernel(tempInds(j)+1:blEnd);
            if os.allChoices(tempInds(j))<0
                tempFollowChoices(j,:) = 1 - tempFollowChoices(j,:);
            end
        end
        followChoices{i} = [followChoices{i}; tempFollowChoices];
        followedKernel{i} = [followedKernel{i}; tempFollowChoiceKernel];
         % focus on rwd trials
        tempInds = lickRwdInds{i};
        tempInds = tempInds(tempInds<length(os.allChoices)); % take away the last trial
        wsChoice{i} = [wsChoice{i} 1-switches(tempInds+1)]; % take next choice
        
        % take all stay trials as ref
        tempInds = refInds{i};
        tempInds = tempInds(tempInds<length(os.allChoices)); % take away the last trial
       % take all following trials
        tempFollowChoicesRef = nan(length(tempInds), p.Results.followLen);
        tempFollowChoiceKernelRef = nan(length(tempInds), p.Results.followLen);
        for j = 1:length(tempInds)
            blEnd = min(length(os.allChoices), tempInds(j)+p.Results.followLen);
            tempFollowChoicesRef(j,1:blEnd-tempInds(j)) = os.allChoices(tempInds(j)+1:blEnd);
            tempFollowChoicesRef(j, tempFollowChoicesRef(j,:)<0) = 0;
            tempFollowChoiceKernelRef(j,1:blEnd-tempInds(j)) = choiceKernel(tempInds(j)+1:blEnd);
            if os.allChoices(tempInds(j))<0
                tempFollowChoicesRef(j,:) = 1 - tempFollowChoicesRef(j,:);
            end
        end
        followChoicesRef{i} = [followChoicesRef{i}; tempFollowChoicesRef];
        followedKernelRef{i} = [followedKernelRef{i}; tempFollowChoiceKernelRef];        
    end

    %% 
    swKernelTemp = NaN(1,length(kernelUniq));
    kernelPropTemp = NaN(1,length(kernelUniq));
    for i = 1:length(kernelUniq)
        swKernelTemp(i) = mean(switches(choiceKernel==kernelUniq(i)));
        kernelPropTemp(i) = sum(choiceKernel==kernelUniq(i))/length(choiceKernel);
    end
    swKernel = [swKernel; swKernelTemp];
    kernelProp = [kernelProp; kernelPropTemp];
    %% long term learning

    wsTemp = NaN(1,p.Results.trialBins);
    lcTemp = NaN(1,p.Results.trialBins);
    wsLTemp = NaN(1,p.Results.trialBins);
    wsRTemp = NaN(1,p.Results.trialBins);
    lcLTemp = NaN(1,p.Results.trialBins);
    lcRTemp = NaN(1,p.Results.trialBins);
    rwdDiffTemp = NaN(1,p.Results.trialBins);
    rwdSumDiffTemp = NaN(1,p.Results.trialBins);
    nRwdSumDiffTemp = NaN(1,p.Results.trialBins);
    choiceDiffTemp = NaN(1,p.Results.trialBins);
    for i = 1:p.Results.trialBins
        if i == p.Results.trialBins
            currInds = floor(length(os.allChoices)/p.Results.trialBins)*(i-1)+1 : length(os.allChoices);
        else
            currInds = floor(length(os.allChoices)/p.Results.trialBins)*(i-1)+1 : floor(length(os.allChoices)/p.Results.trialBins)*i;
        end
        currChoices = os.allChoices(currInds);
        currRewards = os.allRewards(currInds);
        currsw = switches(currInds);

        wsTemp(i) = sum(abs(currRewards(1:end-1))==1 & currsw(2:end)==0)/sum(abs(currRewards(1:end-1))==1);
        lcTemp(i) = sum(abs(currRewards(1:end-1))==0 & currsw(2:end)>0)/sum(abs(currRewards(1:end-1))==0);
        wsLTemp(i) = sum(currRewards(1:end-1)==-1 & currsw(2:end)==0)/sum(currRewards(1:end-1)==-1);
        wsRTemp(i) = sum(currRewards(1:end-1)==1 & currsw(2:end)==0)/sum(currRewards(1:end-1)==1);
        lcLTemp(i) = sum(currRewards(1:end-1)==0 & currChoices(1:end-1)==-1 & currsw(2:end)>0)/sum(abs(currRewards(1:end-1))==0 & currChoices(1:end-1)==-1);
        lcRTemp(i) = sum(currRewards(1:end-1)==0 & currChoices(1:end-1)==1 & currsw(2:end)>0)/sum(abs(currRewards(1:end-1))==0 & currChoices(1:end-1)==1);
        rwdDiffTemp(i) = sum(currRewards==1)/sum(currChoices==1)-sum(currRewards==-1)/sum(currChoices==-1);
        rwdSumDiffTemp(i) = (sum(currRewards==1)-sum(currRewards==-1))/sum(currRewards~=0);
        nRwdSumDiffTemp(i) = 2*sum(currRewards==0 & currChoices==1)/sum(currRewards==0) - 1;
        choiceDiffTemp(i) = 2*sum(currChoices==1)/length(currChoices) - 1;
        
    end
    
        ws = [ws; wsTemp];
        lc = [lc; lcTemp];
        wsDiff = [wsDiff; wsRTemp-wsLTemp];
        lcDiff = [lcDiff; (lcRTemp-lcLTemp)];
        rwdDiff = [rwdDiff; rwdDiffTemp];
        rwdSumDiff = [rwdSumDiff; rwdSumDiffTemp];
        choiceDiff = [choiceDiff; choiceDiffTemp];
        noRwdSumDiff = [noRwdSumDiff; nRwdSumDiffTemp];
end

meanLats = cellfun(@(x) mean(x, 'omitnan'), lickLats);
meanRates = cellfun(@(x) mean(x, 'omitnan'), lickRates);
meanRatesNoRwd = cellfun(@(x) mean(x, 'omitnan'), lickRatesNoRwd);
meanRatesRwd = cellfun(@(x) mean(x, 'omitnan'), lickRatesRwd);
meanLc = cellfun(@(x) mean(x, 'omitnan'), lcChoice);
meanWs = cellfun(@(x) mean(x, 'omitnan'), wsChoice);
meanDist = cellfun(@(x) mean(x, 'omitnan'), distribution);
meanChoices = cellfun(@(x) mean(x, 'omitnan'), followChoices,'UniformOutput', false);
meanSw = cellfun(@(x) mean(x, 'omitnan'), followedKernel,'UniformOutput', false);
meanChoicesRef = cellfun(@(x) mean(x, 'omitnan'), followChoicesRef,'UniformOutput', false);
meanSwRef = cellfun(@(x) mean(x, 'omitnan'), followedKernelRef,'UniformOutput', false);
meanSwKernel = mean(swKernel, 'omitnan');
meanKernelProp = mean(kernelProp, 'omitnan');


semLats = cellfun(@(x) sem(x), lickLats);
semRates = cellfun(@(x) sem(x), lickRates);
semRatesNoRwd = cellfun(@(x) sem(x), lickRatesNoRwd);
semRatesRwd = cellfun(@(x) sem(x), lickRatesRwd);
semDiff = sqrt(semRatesNoRwd.^2 + semRatesRwd.^2);
semLc = cellfun(@(x) sem_bern(x), lcChoice);
semWs = cellfun(@(x) sem_bern(x), wsChoice);
semDist = cellfun(@(x) sem_bern(x), distribution);
semChoices = cellfun(@(x) sem_bern(x), followChoices, 'UniformOutput', false);
semSw = cellfun(@(x) sem(x), followedKernel, 'UniformOutput', false);
semChoicesRef = cellfun(@(x) sem_bern(x), followChoicesRef, 'UniformOutput', false);
semSwRef = cellfun(@(x) sem(x), followedKernelRef, 'UniformOutput', false);
semSwKernel = sem(swKernel);
semKernelProp = sem(kernelProp);

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

%% regression on bias and rwd history
% lc diff as bias
lcChange = -lcDiff(:,end) + lcDiff(:,1);

lm_bias_lc_choicesRwd = fitlm([mean(choiceDiff(:,2:end-1), 2),mean(rwdDiff(:,2:end-1),2)], -(lcDiff(:,end)-lcDiff(:,1)));
lm_bias_lc_choice = fitlm(mean(choiceDiff(:,2:end-1), 2), -(lcDiff(:,end)-lcDiff(:,1)));
lm_bias_lc_rwd = fitlm(mean(rwdDiff(:,2:end-1),2), -(lcDiff(:,end)-lcDiff(:,1)));
lm_bias_lc_rwdSum = fitlm(mean(rwdSumDiff(:,2:end-1),2), -(lcDiff(:,end)-lcDiff(:,1)));
lm_bias_lc_choicesRwdSum = fitlm([mean(choiceDiff(:,2:end-1), 2),mean(rwdSumDiff(:,2:end-1),2)], -(lcDiff(:,end)-lcDiff(:,1)));
lm_bias_lc_noRwdSum = fitlm(mean(noRwdSumDiff(:,2:end-1),2), lcChange);

biasMat = [];
biasVector = [];
for j = 1:length(dayList)
    temp = NaN(2,p.Results.trialBins-2);
    temp(1,:) = lcDiff(j, 2:end-1);
    temp(2,:) = lcDiff(j, 1:end-2);
    biasMat = [biasMat, temp];
    biasVector = [biasVector, lcDiff(j,end-1:end)];
end
lm_bias_lc_autoCorrN = fitlm(biasMat', biasVector);
% ws diff as bias
biasMat = [];
biasVector = [];
for j = 1:length(dayList)
    temp = NaN(2,p.Results.trialBins-2);
    temp(1,:) = wsDiff(j, 2:end-1);
    temp(2,:) = wsDiff(j, 1:end-2);
    biasMat = [biasMat, temp];
    biasVector = [biasVector, wsDiff(j,end-1:end)];
end
lm_bias_lc_autoCorrW = fitlm(biasMat', biasVector);
lm_bias_ws_choicesRwd = fitlm([mean(choiceDiff(:,2:end-1), 2),mean(rwdDiff(:,2:end-1), 2)], (wsDiff(:,end)-wsDiff(:,1)));


%%
figure;
screen = get(0,'Screensize');
screen(4) = screen(4) - 100;
set(gcf, 'Position', screen)
suptitle([sheet '  ' col])

subplot(4,6,1); hold on;
plot(len, meanLats, 'c', 'lineWidth', 2);
patch([len flip(len)], [meanLats + semLats, flip(meanLats - semLats)], 'c', 'FaceAlpha', 0.4, 'EdgeColor', 'none');
title('lickLat N st-rwd')
xlabel('No. st noRwds+1')

subplot(4,6,2); hold on;
plot(len, meanRates, 'color', [0.5 0.5 1], 'lineWidth', 2);
patch([len flip(len)], [meanRates + semRates, flip(meanRates - semRates)], [0.5 0.5 1], 'FaceAlpha', 0.4, 'EdgeColor', 'none');
title('lickRate')
title('lickRate N st-rwd')
xlabel('No. st nRwds+1 ')

subplot(4,6,3); hold on;
plot(len, meanRatesRwd, 'm', 'lineWidth', 2);
patch([len flip(len)], [meanRatesRwd + semRatesRwd, flip(meanRatesRwd - semRatesRwd)], 'm', 'FaceAlpha', 0.4, 'EdgeColor', 'none');
title('lickRateRwd')
title('lickRate rwd N st-rwd')
xlabel('No. st nRwds+1')

subplot(4,6,4); hold on;
plot(len, meanRatesNoRwd, 'b', 'lineWidth', 2);
patch([len flip(len)], [meanRatesNoRwd + semRatesNoRwd, flip(meanRatesNoRwd - semRatesNoRwd)], 'b', 'FaceAlpha', 0.4, 'EdgeColor', 'none');
title('lickRateNoRwd')
title('lickRate noRwd N st-rwd')
xlabel('No. st nRwds+1')

subplot(4,6,5); hold on;
plot(len, meanRatesRwd - meanRatesNoRwd, 'r', 'lineWidth', 2);
patch([len flip(len)], [meanRatesRwd - meanRatesNoRwd + semDiff, flip(meanRatesRwd - meanRatesNoRwd - semDiff)], 'r', 'FaceAlpha', 0.4, 'EdgeColor', 'none');
title('lickRateRwd-NoRwd')
xlabel('No. st nRwds+1')

subplot(4,6,6); hold on;
plot(len, meanLatsNoRwdStay, 'b', 'lineWidth', 2);
patch([len flip(len)], [meanLatsNoRwdStay + semLatsNoRwdsStay, flip(meanLatsNoRwdStay - semLatsNoRwdsStay)], 'b', 'FaceAlpha', 0.4, 'EdgeColor', 'none');
plot(len, meanLatsNoRwdSwitch, 'm', 'lineWidth', 2);
patch([len flip(len)], [meanLatsNoRwdSwitch + semLatsNoRwdsSwitch, flip(meanLatsNoRwdSwitch - semLatsNoRwdsSwitch)], 'm', 'FaceAlpha', 0.4, 'EdgeColor', 'none');

plot(len, meanLats, 'color', [0.5 0.5 0.5], 'lineWidth', 2);
patch([len flip(len)], [meanLats + semLats, flip(meanLats - semLats)], [0.5 0.5 0.5], 'FaceAlpha', 0.4, 'EdgeColor', 'none');
title('lickLatsAfterNoRwd-by choice')
legend({'stay', '', 'switch', '', 'pre'})
xlabel('No. st nRwds+1')

subplot(4,6,7); hold on;
plot(len, meanLc, 'b', 'lineWidth', 2);
patch([len flip(len)], [meanLc + semLc, flip(meanLc - semLc)], 'b', 'FaceAlpha', 0.4, 'EdgeColor', 'none');
title('P(sw|noRwd)')
xlabel('No. st nRwds+1')
subplot(4,6,8); hold on;
plot(len, meanWs, 'r', 'lineWidth', 2);
patch([len flip(len)], [meanWs + semWs, flip(meanWs - semWs)], 'r', 'FaceAlpha', 0.4, 'EdgeColor', 'none');
title('P(st|rwd)')
xlabel('No. st nRwds+1')
subplot(4,6,9); hold on;
plot(len, meanDist, 'b', 'lineWidth', 2);
patch([len flip(len)], [meanDist + semDist, flip(meanDist - semDist)], 'b', 'FaceAlpha', 0.4, 'EdgeColor', 'none');
title('prop cases')
xlabel('No. st nRwds')

subplot(4,6,10); hold on;
colors = cool(length(len));
% plot reference
for i = 1:length(len)
    plot(1:p.Results.followLen, meanChoicesRef{i}, 'Color' ,1-[i/length(len) i/length(len) i/length(len)], 'lineWidth', 2);
    patch([1:p.Results.followLen flip(1:p.Results.followLen)], [meanChoicesRef{i} + semChoicesRef{i}, flip(meanChoicesRef{i} - semChoicesRef{i})], 1-[i/length(len) i/length(len) i/length(len)], 'FaceAlpha', 0.4, 'EdgeColor', 'none');
end

for i = 1:length(len)
    plot(1:p.Results.followLen, meanChoices{i}, 'Color' ,colors(i,:), 'lineWidth', 2);
    patch([1:p.Results.followLen flip(1:p.Results.followLen)], [meanChoices{i} + semChoices{i}, flip(meanChoices{i} - semChoices{i})], colors(i,:), 'FaceAlpha', 0.4, 'EdgeColor', 'none');
end
title('following Choices')
xlabel('No. noRwd Trials')

subplot(4,6,11); hold on;
% plot reference
for i = 1:length(len)
    plot(1:p.Results.followLen, meanSwRef{i}, 'Color' ,1-[i/length(len) i/length(len) i/length(len)], 'lineWidth', 2);
    patch([1:p.Results.followLen flip(1:p.Results.followLen)], [meanSwRef{i} + semSwRef{i}, flip(meanSwRef{i} - semSwRef{i})], 1-[i/length(len) i/length(len) i/length(len)], 'FaceAlpha', 0.4, 'EdgeColor', 'none');
end

for i = 1:length(len)
    plot(1:p.Results.followLen, meanSw{i}, 'Color' ,colors(i,:), 'lineWidth', 2);
    patch([1:p.Results.followLen flip(1:p.Results.followLen)], [meanSw{i} + semSw{i}, flip(meanSw{i} - semSw{i})], colors(i,:), 'FaceAlpha', 0.4, 'EdgeColor', 'none');
end
title('following Kernel')
xlabel('No. noRwd Trials')

subplot(4,6,12); hold on;
yyaxis left
plot(kernelUniq, meanSwKernel, 'Color', [0 1 1], 'LineWidth', 2);
patch([kernelUniq flip(kernelUniq)], [meanSwKernel + semSwKernel, flip(meanSwKernel - semSwKernel)], 'c', 'FaceAlpha', 0.4, 'EdgeColor', 'none');
title('p(sw)/proportion-kernel')
ylabel('p(sw)')
yyaxis right
plot(kernelUniq, meanKernelProp, 'Color', [1 0.5 0], 'LineWidth', 2);
patch([kernelUniq flip(kernelUniq)], [meanKernelProp + semKernelProp, flip(meanKernelProp - semKernelProp)], [1 0.5 0], 'FaceAlpha', 0.4, 'EdgeColor', 'none');
ylabel('proportion')

subplot(4,6,13); hold on;
yyaxis left
plot(1:p.Results.trialBins, mean(ws), 'Color', [0 1 1], 'LineWidth', 2);
patch([1:p.Results.trialBins flip(1:p.Results.trialBins)], [mean(ws) + sem(ws), flip(mean(ws) - sem(ws))], 'c', 'FaceAlpha', 0.4, 'EdgeColor', 'none');
title('p(ws)/p(lc)')
ylabel('p(ws)')
yyaxis right
plot(1:p.Results.trialBins, mean(lc), 'Color', [1 0.5 0], 'LineWidth', 2);
patch([1:p.Results.trialBins flip(1:p.Results.trialBins)], [mean(lc) + sem(lc), flip(mean(lc) - sem(lc))], [1 0.5 0], 'FaceAlpha', 0.4, 'EdgeColor', 'none');
ylabel('p(lc)')

subplot(4,6,14); hold on;

yyaxis left
plot(1:p.Results.trialBins, mean(abs(wsDiff),'omitnan'), 'Color', [0 1 1], 'LineWidth', 2);
patch([1:p.Results.trialBins flip(1:p.Results.trialBins)], [mean(abs(wsDiff),'omitnan') + sem(abs(wsDiff)), flip(mean(abs(wsDiff),'omitnan') - sem(abs(wsDiff)))], 'c', 'FaceAlpha', 0.4, 'EdgeColor', 'none');
title('bias size')
ylabel('p(ws)')
yyaxis right
plot(1:p.Results.trialBins, mean(abs(lcDiff),'omitnan'), 'Color', [1 0.5 0], 'LineWidth', 2);
patch([1:p.Results.trialBins flip(1:p.Results.trialBins)], [mean(abs(lcDiff),'omitnan') + sem(abs(lcDiff)), flip(mean(abs(lcDiff),'omitnan') - sem(abs(lcDiff)))], [1 0.5 0], 'FaceAlpha', 0.4, 'EdgeColor', 'none');
ylabel('p(lc)')

% bin bias change by choiceDiff and rwdDiff
% bin lc by rwdDiff and choiceDiff
sumChoiceDiff = mean(choiceDiff(:,2:end-1),2);
binEdges = binEqualSize(sumChoiceDiff, p.Results.numBins);
meanlcChange = NaN(1, length(binEdges)-1);
semlcChange = NaN(1, length(binEdges)-1);
meanChoiceDiff = NaN(1, length(binEdges)-1);
for j = 1:length(binEdges)-1
    if j < length(binEdges)-1
        meanlcChange(j) = mean(lcChange(sumChoiceDiff>= binEdges(j) & sumChoiceDiff < binEdges(j+1)), 'omitnan');
        semlcChange(j) = sem(lcChange(sumChoiceDiff>= binEdges(j) & sumChoiceDiff< binEdges(j+1)));
        meanChoiceDiff(j) = mean(sumChoiceDiff(sumChoiceDiff>= binEdges(j) & sumChoiceDiff< binEdges(j+1)), 'omitnan');
    else
        meanlcChange(j) = mean(lcChange(sumChoiceDiff>= binEdges(j) & sumChoiceDiff <= binEdges(j+1)), 'omitnan');
        semlcChange(j) = sem(lcChange(sumChoiceDiff>= binEdges(j) & sumChoiceDiff <= binEdges(j+1)));
        meanChoiceDiff(j) = mean(sumChoiceDiff(sumChoiceDiff>= binEdges(j) & sumChoiceDiff<= binEdges(j+1)), 'omitnan');
    end
end

subplot(4,6,15); hold on
plot(meanChoiceDiff, meanlcChange, 'LineWidth',2, 'color', 'c');
patch([meanChoiceDiff flip(meanChoiceDiff)], [meanlcChange + semlcChange, flip(meanlcChange - semlcChange)], 'c', 'FaceAlpha', 0.4, 'EdgeColor', 'none');
xlabel('choiceDiff');
ylabel('bias')
% bin lc by rwdNum
sumRwdDiff = mean(rwdSumDiff(:,2:end-1),2);
binEdges = binEqualSize(sumRwdDiff, p.Results.numBins);
meanlcChange = NaN(1, length(binEdges)-1);
semlcChange = NaN(1, length(binEdges)-1);
meanRwdDiff = NaN(1, length(binEdges)-1);
for j = 1:length(binEdges)-1
    if j < length(binEdges)-1
        meanlcChange(j) = mean(lcChange(sumRwdDiff>= binEdges(j) & sumRwdDiff < binEdges(j+1)), 'omitnan');
        semlcChange(j) = sem(lcChange(sumRwdDiff>= binEdges(j) & sumRwdDiff< binEdges(j+1)));
        meanRwdDiff(j) = mean(sumRwdDiff(sumRwdDiff>= binEdges(j) & sumRwdDiff< binEdges(j+1)), 'omitnan');
    else
        meanlcChange(j) = mean(lcChange(sumRwdDiff>= binEdges(j) & sumRwdDiff < binEdges(j+1)), 'omitnan');
        semlcChange(j) = sem(lcChange(sumRwdDiff>= binEdges(j) & sumRwdDiff< binEdges(j+1)));
        meanRwdDiff(j) = mean(sumRwdDiff(sumRwdDiff>= binEdges(j) & sumRwdDiff< binEdges(j+1)), 'omitnan');
    end
end


subplot(4,6,16)
plot(meanRwdDiff, meanlcChange, 'LineWidth',2, 'color', 'm');
patch([meanRwdDiff flip(meanRwdDiff)], [meanlcChange + semlcChange, flip(meanlcChange - semlcChange)], 'm', 'FaceAlpha', 0.4, 'EdgeColor', 'none');
xlabel('rwdDiff')
ylabel('bias')

sumnoRwdSumDiff = mean(noRwdSumDiff(:,2:end-1),2);
binEdges = binEqualSize(sumnoRwdSumDiff, p.Results.numBins);
meanlcChange = NaN(1, length(binEdges)-1);
semlcChange = NaN(1, length(binEdges)-1);
meannoRwdDiff = NaN(1, length(binEdges)-1);
for j = 1:length(binEdges)-1
    if j < length(binEdges)-1
        meanlcChange(j) = mean(lcChange(sumnoRwdSumDiff>= binEdges(j) & sumnoRwdSumDiff < binEdges(j+1)), 'omitnan');
        semlcChange(j) = sem(lcChange(sumnoRwdSumDiff>= binEdges(j) & sumnoRwdSumDiff< binEdges(j+1)));
        meannoRwdDiff(j) = mean(sumnoRwdSumDiff(sumnoRwdSumDiff>= binEdges(j) & sumnoRwdSumDiff< binEdges(j+1)), 'omitnan');
    else
        meanlcChange(j) = mean(lcChange(sumnoRwdSumDiff>= binEdges(j) & sumnoRwdSumDiff < binEdges(j+1)), 'omitnan');
        semlcChange(j) = sem(lcChange(sumnoRwdSumDiff>= binEdges(j) & sumnoRwdSumDiff< binEdges(j+1)));
        meannoRwdDiff(j) = mean(sumnoRwdSumDiff(sumnoRwdSumDiff>= binEdges(j) & sumnoRwdSumDiff< binEdges(j+1)), 'omitnan');
    end
end


subplot(4,6,17)
plot(meannoRwdDiff, meanlcChange, 'LineWidth',2, 'color', 'b');
patch([meannoRwdDiff flip(meannoRwdDiff)], [meanlcChange + semlcChange, flip(meanlcChange - semlcChange)], 'b', 'FaceAlpha', 0.4, 'EdgeColor', 'none');
xlabel('noRwdDiff')
ylabel('bias')
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
    [params, modelName, ~, noSession] = getStanModelParams_sampsOnly(animalName, col, modelName, numSamps, 'sessionParamsFlag', 0, 'sessionName', dayList{sess});
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
%% long timescale
% frist half ws lc vs second half

% ws(R) - ws(L) second half ~ first half Left rwd - right rwd 

%% bin pe, plot lickRate after stay
numBins = 12;
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
subplot(4,6,19); hold on;
plot(peMeans(peMeans<0), rateMeans(peMeans<0), 'b', 'lineWidth', 2);
patch([peMeans(peMeans<0); flip(peMeans(peMeans<0))], [rateMeans(peMeans<0) + rateSems(peMeans<0); flip(rateMeans(peMeans<0) - rateSems(peMeans<0))], 'b', 'FaceAlpha', 0.4, 'EdgeColor', 'none');
title('lickRate-stay-noRwd')
xlabel('pe')

subplot(4,6,20); hold on;
plot(peMeans(peMeans>0), rateMeans(peMeans>0), 'm', 'lineWidth', 2);
patch([peMeans(peMeans>0); flip(peMeans(peMeans>0))], [rateMeans(peMeans>0) + rateSems(peMeans>0); flip(rateMeans(peMeans>0) - rateSems(peMeans>0))], 'm', 'FaceAlpha', 0.4, 'EdgeColor', 'none');
title('lickRate-stay-rwd')
xlabel('pe')


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
        peMeans(k) = mean(combinePeSwitch(combinePeSwitch >= edges(k) & combinePeSwitch < edges(k+1)), 'omitnan');
    else
        rateTemp = combineLickRatesRwdSwitch(combinePeSwitch >= edges(k) & combinePeSwitch <= edges(k+1));
        peMeans(k) = mean(combinePeSwitch(combinePeSwitch >= edges(k) & combinePeSwitch <= edges(k+1)), 'omitnan');
    end
    rateMeans(k) = mean(rateTemp, 'omitnan');
    rateSems(k) = sem(rateTemp);
end
%%
subplot(4,6,21); hold on;
plot(peMeans(peMeans<0), rateMeans(peMeans<0), 'b', 'lineWidth', 2);
patch([peMeans(peMeans<0); flip(peMeans(peMeans<0))], [rateMeans(peMeans<0) + rateSems(peMeans<0); flip(rateMeans(peMeans<0) - rateSems(peMeans<0))], 'b', 'FaceAlpha', 0.4, 'EdgeColor', 'none');
title('lickRate-sw-noRwd')
xlabel('pe')

subplot(4,6,22); hold on;
plot(peMeans(peMeans>0), rateMeans(peMeans>0), 'm', 'lineWidth', 2);
patch([peMeans(peMeans>0); flip(peMeans(peMeans>0))], [rateMeans(peMeans>0) + rateSems(peMeans>0); flip(rateMeans(peMeans>0) - rateSems(peMeans>0))], 'm', 'FaceAlpha', 0.4, 'EdgeColor', 'none');
title('lickRate-sw-rwd')
xlabel('pe')
%% bin choiceConf, plot lickLat and lickRate before

%% bin Qc, plot lickLat and lickRate before

%%


