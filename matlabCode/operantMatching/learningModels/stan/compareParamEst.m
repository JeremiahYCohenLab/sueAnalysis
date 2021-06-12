function compareParamEst(animalName, category, model, bin1, bin2, varargin)
%task and model parameters
p = inputParser;
p.addParameter('paramNames', {'aNmin', 'aP', 'aF', 'aPE', 'v', 'beta', 'bias'});
p.parse(varargin{:});
[root, sep] = currComputer();
if contains(model,'tF')
    input = 'choice, outcome, ITI)';
else
    input = 'choice, outcome)';
end
    
%% load model fitting results
sampFile = [animalName category '_', model];
path = [root animalName sep animalName 'sorted' sep 'stan' sep 'bernoulli' sep model sep category sep];
load([path sampFile '.mat'], 'dayList');
samples = load([path sampFile '.mat'], sampFile);
samples = samples.(sampFile);
paramNames = p.Results.paramNames;
paramEsts = zeros(length(dayList), length(paramNames));
paramEsts2 = zeros(length(dayList), length(paramNames));
aNSatRatio1 = zeros(length(dayList), 1);
aNSatRatio2 = zeros(length(dayList), 1);
for id = 1:length(dayList)
%     fprintf([dayList{id} '\n']);
    s = behAnalysisNoPlot_opMD(dayList{id});
    choice = s.allChoices';
    choice(choice<0) = 0;
    outcome = abs(s.allRewards);
    %generate best estimates of parameters
    allSamples = [];
    % first bins
    edges = cell(1,length(paramNames));
    for i = 1:length(paramNames)
        tmp = samples.(paramNames{i})(:,id);
%         if strcmp(paramNames{i}, 'aPE') || strcmp(paramNames{i}, 'v')
%             tmp = log(tmp);
%         end
        allSamples = [allSamples tmp];
        edges{i} = linspace(min(tmp), max(tmp),bin1+1);
    end
    n = histcnd(allSamples,edges); %bin samples by multiple dimensions
    [~, inds] = myMaxAll(n); %find the bin with max num in bin
    % second bins
    edges2 = cell(1,length(paramNames));
    for i = 1:length(paramNames) %use previous best bin as newbin
        edgeTmp = edges{i};
        edges2{i} = linspace(edgeTmp(inds(i)),edgeTmp(inds(i)+1),bin2+1);
    end
    n = histcnd(allSamples,edges2); %bin samples by multiple dimensions
    [~, inds] = myMaxAll(n); %find the bin with max num in bin   
    for i = 1:length(paramNames) %use median in bin as best estimate
        tmp = allSamples(:,i);
        edgeTmp = edges2{i};
        if inds(i) < bin2
            paramEsts(id,i) = median(tmp(tmp >= edgeTmp(inds(i)) & tmp < edgeTmp(inds(i)+1)));
        else
            paramEsts(id,i) = median(tmp(tmp >= edgeTmp(inds(i)) & tmp <= edgeTmp(inds(i)+1)));
        end
    end
    
%     paramEsts(id,contains(paramNames, 'aPE')) = exp(paramEsts(id,contains(paramNames, 'aPE')));
%     paramEsts(id,contains(paramNames, 'v')) = exp(paramEsts(id,contains(paramNames, 'v')));
    
    if contains(model, 'ord')
        modelSim = model(1:end-4);
    else
        modelSim = model;
    end
    eval(['[LL, probChoice, Q, pe, pePe, aN, peBar] = qLearningModel_' modelSim '(paramEsts(id,:),' input ';']);
    ll1(id) = LL/length(choice);
    aNSatRatio1(id) = sum(aN==1)/length(outcome);
    % version 2
    
    for i = 1:length(paramNames)
        tmp = samples.(paramNames{i})(:,id);
        if strcmp(paramNames{i}, 'aPE') || strcmp(paramNames{i}, 'v')
            tmp = log(tmp);
        end
        [n,e] = histcounts(tmp, bin1*bin2);
        [~, maxInd] = max(n);
        paramEsts2(id,i) = median(tmp(tmp > e(maxInd) & tmp < e(maxInd+1)));
    end
    paramEsts2(id,contains(paramNames, 'aPE')) = exp(paramEsts2(id,contains(paramNames, 'aPE')));
    paramEsts2(id,contains(paramNames, 'v')) = exp(paramEsts2(id,contains(paramNames, 'v')));
    eval(['[LL,probC, Q, pe,pePe, aN, peBar] = qLearningModel_' modelSim '(paramEsts2(id,:),' input ';']);
    ll2(id) = LL/length(outcome);
    aNSatRatio2(id) = sum(aN==1)/length(outcome);
end
figure2('position', [0 0 800 1200]); hold on;
colors = cool(length(ll1));
[~, ind] = sort(ll1);
row = ceil((length(paramNames)+2)/2);

subplot(row,2,1)
scatter(ll1(ind), ll2(ind), 15, colors);
xlabel('combine')
ylabel('single')
line([minmax(ll1')], [minmax(ll1')], 'color', [0.7, 0.7, 0.7]);
title(['trial ll ' num2str(bin1*bin2)])

subplot(row,2,2)
scatter(aNSatRatio1(ind), aNSatRatio2(ind), 15, colors);
xlabel('combine')
ylabel('single')
line([minmax(aNSatRatio1')], [minmax(aNSatRatio2')], 'color', [0.7, 0.7, 0.7]);
title(['aN saturation ' num2str(bin1*bin2)])

for i = 1:length(paramNames)
    subplot(row,2,i+2)
    if strcmp(paramNames{i}, 'aPE') || strcmp(paramNames{i}, 'v')
        scatter(log(paramEsts(ind,i)), log(paramEsts2(ind,i)), 15, colors);
        line(minmax(log(paramEsts(:,i)')), minmax(log(paramEsts2(:,i)')), 'color', [0.7, 0.7, 0.7]);
        title(['log(' paramNames{i} ') ' num2str(bin1*bin2)])
    else
        scatter(paramEsts(ind,i), paramEsts2(ind,i), 15, colors);
        line(minmax(paramEsts(:,i)'), minmax(paramEsts2(:,i)'), 'color', [0.7, 0.7, 0.7]);
        title([paramNames{i} ' ' num2str(bin1*bin2)])
    end
    xlabel('combine')
    ylabel('single')
end


 