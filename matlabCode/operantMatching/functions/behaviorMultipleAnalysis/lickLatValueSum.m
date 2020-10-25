function [combineQSum, combineLickLat] = lickLatValueSum(animal, category, varargin)
%task and model parameters
p = inputParser;
% default parameters if none given
p.addParameter('plotFlag', 1);
p.addParameter('maxTrials', 200);
p.addParameter('modelNameOld', 'fourParam_tF_oneside');
p.addParameter('modelName', '4params_tF_oneside');
p.addParameter('paramNames', {'aN', 'aP', 'tF', 'beta'});
p.addParameter('numBins', 10);
p.parse(varargin{:});

colors = cool(5);
paramNames = p.Results.paramNames;
maxTrial = p.Results.maxTrials;
%determine root for file location and load model file
[root, sep] = currComputer();
sampFile = [animal category '_', p.Results.modelNameOld];
path = [root animal sep animal 'sorted' sep 'stan' sep 'bernoulli' sep p.Results.modelNameOld sep];
load([path sampFile '.mat'], 'dayList');
samples = load([path sampFile '.mat'], sampFile);
samples = samples.(sampFile);


combineLickLat = [];
combineQSum = [];
combinQC = [];
combineQuitTrial = [];
combineQuitTrialFirst = [];
combineQuitSum = [];
combineQuitSumFirst = [];


if contains(p.Results.modelName,'tF')
    input = 'choice, outcome, ITI)';
else
    input = 'choice, outcome)';
end

rsq = zeros(size(dayList));
rmse = zeros(size(dayList));
for i = 1:length(dayList)
 %   disp(dayList{i})
    % load data
    behSessionData = loadBehavioralData([dayList{i} '.asc']);
    behavStruct = parseBehavioralData(behSessionData,maxTrial);
    choice = behavStruct.allChoices(1:min(maxTrial,length(behavStruct.allChoices)));
    choice(choice==-1)=0;
    outcome = behavStruct.allRewards(1:min(maxTrial,length(behavStruct.allChoices)));
    outcome = abs(outcome);
    ITI = behavStruct.timeBtwn(1:min(maxTrial,length(behavStruct.allChoices)));
    
    % compute session parameters and latent variables
    for j = 1:length(paramNames)
        tmp = samples.(paramNames{j})(:,i);
        [n,e] = histcounts(tmp, 50);
        [~, maxInd] = max(n);
        params(j) = median(tmp(tmp > e(maxInd) & tmp < e(maxInd+1)));
    end
    eval(['[~,~,Q,~,cQ] = qLearningModel_' p.Results.modelName '(params,' input ';'])
    
    % find quit trials
    A = cellfun(@(x) strcmp(x,'CSplus'),{behSessionData.trialType});
    B = cellfun(@isnan, {behSessionData.rewardTime});
    quitTrial = find(A>0 & B>0);
    lastLickTrial = ones(size(quitTrial));
    % if there's quitTrial within length, find the Qsums
    for j = 1:length(quitTrial) % find the last responding trial before quiting
       diffTrial = behavStruct.responseInds - quitTrial(j);
       if ~isnan(find(diffTrial<0))
           trialID = max(diffTrial(diffTrial<0));
           lastLickTrial(j) = behavStruct.responseInds(diffTrial==trialID);
       end
    end

    if length(behavStruct.responseInds)>p.Results.maxTrials % only keep those within session length
        quitTrial = quitTrial(lastLickTrial <= behavStruct.responseInds(200));
        lastLickTrial = lastLickTrial(lastLickTrial<= behavStruct.responseInds(200));
    end    
    
    if ~isempty(quitTrial) 
        quitSumQ = zeros(size(lastLickTrial));     
        % assign value to the quit trial
        for j = 1:length(quitTrial)
            if quitTrial(j)==1
                quitSumQ = 0;
            else
                nextTrial = find(behavStruct.responseInds == lastLickTrial(j))+1; % next respond trial in respondInds
                Qtemp = Q(nextTrial,:);
                % find the corresponding q with two different models
                if contains(p.Results.modelName,'tF')
                    if nextTrial~=maxTrial+1
                        timeToNext = (behSessionData(behavStruct.responseInds(nextTrial)).CSon - behSessionData(quitTrial(j)).CSon)/1000;
                        if choice(nextTrial-1)==1
                            Qtemp(1) = params(3)^(-timeToNext)*Qtemp(1);
                        else
                            Qtemp(2) = params(3)^(-timeToNext)*Qtemp(2);
                        end
                    else
                        timeFromMaxTrial = (behSessionData(quitTrial(j)).CSon - behSessionData(behavStruct.responseInds(maxTrial)).rewardTime)/1000;
                        if choice(nextTrial-1)==1
                            Qtemp(1) = params(3)^(timeFromMaxTrial)*Qtemp(1);
                        else
                            Qtemp(2) = params(3)^(timeFromMaxTrial)*Qtemp(2);
                        end
                    end
                end
                
                quitSumQ(j) = sum(Qtemp);
            end
        end
        % find the first quits in a sequence of quitting
        seqs = [1 find(diff(quitTrial)~=1) + 1];
        % project to respondInds
        for j = 1:length(quitTrial)
            if ~isempty(find(behavStruct.responseInds == lastLickTrial(j), 1))
             quitTrial(j) = find(behavStruct.responseInds == lastLickTrial(j));
            end
        end
        quitFirst = quitTrial(seqs);
        quitSumQFirst = quitSumQ(seqs);
            
        combineQuitTrial = [combineQuitTrial quitTrial];
        combineQuitSum = [combineQuitSum quitSumQ];
        combineQuitTrialFirst = [combineQuitTrialFirst quitFirst];
        combineQuitSumFirst = [combineQuitSumFirst quitSumQFirst];
    end
    lickLat = behavStruct.lickLatZ(1:min(length(behavStruct.lickLat), maxTrial));
    combineLickLat = [combineLickLat lickLat];
    combineQSum = [combineQSum sum(Q(1:end-1,:)')];
    combinQC = [combinQC cQ'];
    
    lmTemp = fitlm(sum(Q(1:end-1,:),2), lickLat');
    rsq(i) = lmTemp.Rsquared.Adjusted;
    rmse(i) = lmTemp.RMSE;
end
    

    

% group lickLats by rwd env
[group,~] = discretize(combineQSum,p.Results.numBins);

len = zeros(1,p.Results.numBins);

for i = 1:p.Results.numBins
    len(i) = length(find(group==i));
end

small = find(len<2);

if ~isempty(small)&&length(small)==1
    if small < 0.5*p.Results.numBins
        for j = small:1
            group(group==j) = j+1;
        end
    else 
        for j = small:p.Results.numBins
            group(group==j) = j-1;
        end
    end
end

lickLatM = zeros(1,max(group));
lickLatSEM = zeros(1,max(group));
binMean = zeros(1,max(group));

for i = 1:max(group)
    lickLatM(i) = mean(combineLickLat(group == i));
    lickLatSEM(i) = std(combineLickLat(group == i))/sqrt(len(i));
    binMean(i) = mean(combineQSum(group == i));
end

lm = fitlm(combineQSum,combineLickLat');
%% plot everything

figure2('Position', [1 1 900 600]); 
c = 4;
subplot(2,3,1); hold on;
% scatter(combineQSum, combineLickLat, 3, 'MarkerFaceColor',[.7 .7 .7], 'MarkerEdgeColor', 'none');
text(0.2, 0.7, sprintf('drug R^2 = %0.4f RMSE = %0.2f', lm.Rsquared.Adjusted, lm.RMSE));

errorbar(binMean, lickLatM, lickLatSEM, 'Color', colors(c,:), 'linewidth',1.5)
%ylim([200 1500]) 
ylabel('lickLat')
xlabel('Qsum')

subplot(2,3,2); hold on;
histogram(rsq, 10,'Normalization', 'Probability', 'FaceColor', colors(c,:))
title('adjustedRsquare')

subplot(2,3,3); hold on;
histogram(rmse, 10,'Normalization', 'Probability', 'FaceColor', colors(c,:))
title('RMSE')

subplot(2,3,4); hold on;
cdfplot(combineQuitTrialFirst)
xlabel('Trial')
ylabel('pStartQuit')

subplot(2,3,5); hold on;
cdfplot(combineQuitSumFirst)
xlabel('Qsum')
ylabel('pStartQuit')


subplot(2,3,6); hold on;
histogram(combineQuitSumFirst, 10, 'FaceColor', colors(c,:))
xlabel('Qsum')
title('numQuits')

suptitle([animal '  ' category])

























