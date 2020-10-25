function criterion = computeDIC(animal,category,modelNameOld,modelName,paramNames)
% animal = 'ZS040RwdDelay';
% category = 'good';
% modelNameOld = 'fourParam_tF_oneside';
% modelName = '4params_tF_oneside';
% paramNames = {'aN', 'aP', 'tF', 'beta'};
maxTrial = 200;
plotFlag = 1;

[root, sep] = currComputer();
sampFile = [animal category '_', modelNameOld];
path = [root animal sep animal 'sorted' sep 'stan' sep 'bernoulli' sep modelNameOld sep];
load([path sampFile '.mat'], 'dayList');
load([path sampFile '.mat'], 'paramEsts');
samples = load([path sampFile '.mat'], sampFile);
samples = samples.(sampFile);
%%
dic = zeros(length(dayList),1);
ell = zeros(length(dayList),1);
pd = zeros(length(dayList),1);
llhat = zeros(length(dayList),1);
llest = zeros(length(dayList),1);
llmed = zeros(length(dayList),1);
llmax = zeros(length(dayList),1);
maxll = max(samples.log_lik)';
llsesh = zeros(length(dayList),1);
for i = 1:length(dayList)
    ll = samples.log_lik(:,i);
    behSessionData = loadBehavioralData([dayList{i} '.asc']);
    behavStruct = parseBehavioralData(behSessionData, maxTrial);
    choice = behavStruct.allChoices(1:min(maxTrial,length(behavStruct.allChoices)));
    choice(choice==-1)=0;
    outcome = behavStruct.allRewards(1:min(maxTrial,length(behavStruct.allChoices)));
    outcome = abs(outcome);
    ITI = behavStruct.timeBtwn(1:min(maxTrial,length(behavStruct.allChoices)));
    if contains(modelName,'tF')
        input = 'choice, outcome, ITI)';
    else
        input = 'choice, outcome)';
    end
    % for dic
    params = zeros(length(paramNames),1);
    for j = 1:length(paramNames)
        temp = samples.(paramNames{j});
        params(j) = mean(temp(:,i));
    end
    eval(['llhat(i) = qLearningModel_' modelName '(params,' input ';']);
    dic(i) = 2*mean((-2)*ll) + 2*llhat(i);
    ell(i) = mean(ll);                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            
    pd(i) = mean((-2)*ll) + 2*llhat(i);
    % animal level
    eval(['llest(i) = qLearningModel_' modelName  '(paramEsts,' input ';'])
    % session level med
    params = zeros(length(paramNames),1);
    for j = 1:length(paramNames)
        temp = samples.(paramNames{j});
        params(j) = median(temp(:,i));
    end
    eval(['llmed(i) = qLearningModel_' modelName '(params,' input ';'])
    % session level max ll
    params = zeros(length(paramNames),1);
    [~, loc] = max(ll);
    for j = 1:length(paramNames)
        temp = samples.(paramNames{j});
        params(j) = temp(loc,i);
    end   
    eval(['llmax(i) = qLearningModel_' modelName '(params,' input ';'])
    % session level max params
    for j = 1:length(paramNames)
        tmp = samples.(paramNames{j})(:,i);
        [n,e] = histcounts(tmp, 50);
        [~, maxInd] = max(n);
        params(j) = median(tmp(tmp > e(maxInd) & tmp < e(maxInd+1)));
    end
    eval(['llsesh(i) = qLearningModel_' modelName '(params,' input ';'])
    
    criterion = struct;
    criterion.dic = dic;
    criterion.ell = ell;
    criterion.llhat = llhat;
    criterion.llani = llest;
    criterion.llmed = llmed;
    criterion.pd = pd;
    criterion.maxll = maxll;
    criterion.llmax = llmax;
    criterion.ll = samples.log_lik;
    criterion.llsesh = llsesh;
    
    save([path sampFile '.mat'],'criterion','-append')
end
%%