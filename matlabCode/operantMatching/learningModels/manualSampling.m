%% find infinite numbers
%% 
animalName = 'ZS061';
category = 'good';

[root, sep] = currComputer();
xlFile = [animalName '.xlsx'];
sheet = animalName;
[~, dayList, ~] = xlsread([root xlFile], sheet);
[~,col] = find(contains(dayList, category) == 1);

dayList = dayList(2:end,col);
endInd = find(cellfun(@isempty,dayList),1);
if ~isempty(endInd)
    dayList = dayList(1:endInd-1,:);
end
%% simulate values
sampleNum = 1000;
lambda = rand(sampleNum,1);
v0 = 2*rand(sampleNum,1);
omega = 5*rand(sampleNum,1);
beta = 5*rand(sampleNum,1);
%%
infM = cell(sampleNum,1);
parfor iter = 1:sampleNum  
    infMTemp = NaN(length(dayList),1);
    for ses = 1:length(dayList)
        behSessionData = loadBehavioralData([dayList{ses} '.asc']);
        behavStruct = parseBehavioralData(behSessionData, 600);
        choice = behavStruct.allChoices;
        choice(choice == -1) = 0;
        outcome = abs(behavStruct.allRewards);     
        
        [LH, probChosen, m, pe, v, w] = vkf_LL([lambda(iter) v0(iter) omega(iter) beta(iter)], choice, outcome);
        mdiff = m(:,2) - m(:,1);
        infMTemp(ses) = sum(abs(mdiff)>20);
    end
    infM{iter} = infMTemp;
end
%% simulation and recovery
% realParams
iteration = 200;
runs = 100;
maxTrial = 500;
paramNames = getParamNames_dF('vkf_fixV_kappa', 0);
% generate random paramters
params = zeros(iteration, 4);
params(:,1) = 0.9*betarnd(3, 5, iteration,1)+0.1; % lambda [0 1]
params(:,2) = 10*betarnd(5, 2, iteration,1); % v0 [0 10]
% params(:,2) = 2*3/7 * ones(iteration,1); % v0 [0 2]
params(:,3) = (9.9*betarnd(4, 2, iteration,1)+0.1); % omega [0 10]
params(:,4) = 0.7*betarnd(6, 3, iteration,1)+0.3; % beta [0 1]
% params(:,5) = 0.7*betarnd(6, 3, iteration,1)+0.3; % aF [0 1]
params(:,5) = 0.8*betarnd(6, 3, iteration,1)+0.2; % kappa [0 1]
% constraints for fitting
lambda_range = [0 1];
v0_range = [0 10];
omega_range = [0 10];
beta_range = [0 1];
aF_range = [0 1];
kappa_range = [0 1];
% start and fitted values for fitting
startValues = [];
startValues(:,1) = 0.9*rand(runs,1)+0.1;
startValues(:,2) = (10*rand(runs,1));
startValues(:,3) = (9.9*rand(runs,1)+0.1);
startValues(:,4) = 0.7*rand(runs,1)+0.3;
startValues(:,5) = (0.9*rand(runs,1)+0.1);

A=[eye(size(startValues, 2)); -eye(size(startValues, 2))];
b=[ lambda_range(2);  v0_range(2);  omega_range(2); beta_range(2); kappa_range(2);
   -lambda_range(1); -v0_range(1); -omega_range(1); -beta_range(1); -kappa_range(1)];

paramsEstimates = zeros(size(params));

% algorithm specs 
options = optimset('Algorithm', 'interior-point','ObjectiveLimit',...
        -1.000000000e+300,'TolFun',1e-15, 'Display','off', 'TolX', 1e-10);
    
figure;
for i = 1:length(paramNames)
subplot(2,3,i); hold on;
histogram(params(:,i),linspace(min(params(:,i)), max(params(:,i)), 20), 'Normalization', 'Probability','FaceColor', 'c');
histogram(startValues(:,i),linspace(min(startValues(:,i)), max(startValues(:,i)), 20), 'Normalization', 'Probability','FaceColor', 'm', 'FaceAlpha', 0.2, 'EdgeAlpha', 0.3);
title(paramNames{i})
if i == 1
    legend({'simulationParams', 'startValues'});
end
end
%%
msim = [];
mfit = [];
psim = [];
pfit = [];
choice = cell(iteration,1);
outcome = cell(iteration,1);
pRight = cell(iteration,1);
LHsim = ones(iteration,1);
LHest = ones(iteration,1);
for sim = 1:iteration
    %simulation
    [~, outcomeTemp, choiceTemp, pRightTemp] = vkfSim_kappa('params', params(sim,:),'randomSeed', sim,'maxTrials', maxTrial, 'plotFlag', 0);
    choiceTemp(choiceTemp<0) = 0;
    outcomeTemp = abs(outcomeTemp);
    choice{sim} = choiceTemp;
    outcome{sim} = outcomeTemp;
    pRight{sim} = pRightTemp;
    [LHsim(sim),~,mtemp,~,~,~,~,~,ptemp] = vkf_LL_kappa(params(sim,:), choiceTemp, abs(outcomeTemp), 0); 
    msim = [msim; mtemp];
    psim = [psim; ptemp];
    %%
    %model fitting

    % initialization
    allParams = zeros(size(startValues));
    LH = zeros(size(startValues, 1), 1);
    exitFl = zeros(size(startValues, 1), 1);
    hess = zeros(size(startValues, 1), size(startValues, 2), size(startValues, 2));

    % fmincon
    parfor r = 1:runs
        [allParams(r, :), LH(r, :), exitFl(r, :), ~, ~, ~, hess(r, :, :)] = ...
            fmincon(@vkf_LL_kappa, startValues(r,:), A, b, [], [], [], [], [], options, choiceTemp, outcomeTemp, 0);
    end
    [~,bestFit] = min(LH);
    [LHest(i),~,mtemp,~,~,~,~,~,ptemp] = vkf_LL_kappa(allParams(bestFit, :), choiceTemp, outcomeTemp, 0);
    mfit =[mfit; mtemp];
    pfit = [pfit; ptemp];
    paramsEstimates(sim,:) = allParams(bestFit, :);
end
%%
figure2;
colors = cool(iteration);
[~,int] = sort(params(:,5));
for i = 1:length(paramNames)
subplot(2,3,i); hold on;
scatter(params(int,i), paramsEstimates(int,i),15, colors, 'filled');
scatter(params(paramsEstimates(:,2)>3.95,i), paramsEstimates(paramsEstimates(:,2)>3.95,i), 15, 'r', 'filled');
% scatter(params(paramsEstimates(:,3)<0.05,i), paramsEstimates(paramsEstimates(:,3)<0.05,i), 15, 'm', 'filled');
line([minmax(paramsEstimates(:,i))], [minmax(paramsEstimates(:,i))], 'color', [0.5 0.5 0.5]);
[rho, pVal] = corr([params(:,i), paramsEstimates(:,i)]);
title(sprintf([paramNames{i} ' %.2f p=%.2f'], rho(1,2), pVal(1,2)))
end
suptitle(['maxTrial = ' num2str(maxTrial)])
figure2;
subplot(2,1,1)
colors = cool(size(msim,1));
[~,int] = sort(msim(:,2));
scatter(msim(int,2)-msim(int,1), mfit(int,2)-mfit(int,1), 15, colors, 'filled');
title('mdiff')
subplot(2,1,2)
scatter(psim(int), pfit(int), 15, colors, 'filled');
title('pRight')
suptitle(['maxTrial = ' num2str(maxTrial)])
 %%
int = find(paramsEstimates(:,2)>3.95);
int  = 1:iteration;
LH1 = zeros(length(int),1);
LH2 = zeros(length(int),1);
probChosen1 = cell(length(int),1);
probChosen2 = cell(length(int),1);
m1 = cell(length(int),1);
m2 = cell(length(int),1);
pe1 = cell(length(int),1);
pe2 = cell(length(int),1);
v1 = cell(length(int),1);
v2 = cell(length(int),1);
w1 = cell(length(int),1);
w2 = cell(length(int),1);
k1 = cell(length(int),1);
k2 = cell(length(int),1);
alpha1 = cell(length(int),1);
alpha2 = cell(length(int),1);
for i = 1:length(int)
    [LH1(i), probChosen1{i}, m1{i}, pe1{i}, v1{i}, w1{i}, k1{i}, alpha1] = vkf_LL_kappa(paramsEstimates(int(i),:), choice{int(i)}, abs(outcome{int(i)}), 0);
    [LH2(i), probChosen2{i}, m2{i}, pe2{i}, v2{i}, w2{i}, k2{i}, alpha2] = vkf_LL_kappa(params(int(i),:), choice{int(i)}, abs(outcome{int(i)}), 0);
end
%% plot pChosen sorted by params
[~,ints] = sort(paramsEstimates(:,5));
colors = cool(length(ints));
figure2; hold on;
for i = 1:length(ints)
    scatter(probChosen2{ints(i)}, probChosen1{ints(i)}, 8, colors(i,:), 'filled');
end
%% param gradient for simulation
i = 19;
vTest = linspace(0.1,1,20);
omegaTest = 10*linspace(0.1,1,20);
LHTest = zeros(size(vTest));
paramsTest = params(i,:);
probChosen = cell(length(vTest),1);
pRight = cell(length(vTest),1);
m = cell(length(vTest),1);
pe = cell(length(vTest),1);
v = cell(length(vTest),1);
w = cell(length(vTest),1);
alpha = cell(length(vTest),1);
k = cell(length(vTest),1);
[LHTestSim, probChosenSim,~,peSim,~,~,kSim,alphaSim,pRightSim] = vkf_LL_kappa(paramsTest, choice{i}, abs(outcome{i}), 0);
% LH, probChosen, m, pe, v, w, k, alpha, probChoice
for j = 1:length(vTest)
%     paramsTest(2) = vTest(j);
    paramsTest(3) = omegaTest(j);
    [LHTest(j), probChosen{j},m{j},pe{j},v{j},w{j},k{j},alpha{j},pRight{j}] = vkf_LL_kappa(paramsTest, choice{i}, abs(outcome{i}), 0);
end
%% plotting
figure2;
plot(omegaTest, LHTest)
figure2; hold on
colors = cool(length(vTest));
for j = 1:length(vTest)
    plot(pRight{j}, 'color', colors(j,:))
end
plot(pRightSim, 'color', 'r', 'LineWidth', 3)
xlim([1 100])
title('pRight')

figure2; hold on
colors = cool(length(vTest));
for j = 1:length(vTest)
    plot(probChosen{j}, 'color', colors(j,:))
end
plot(probChosenSim, 'color', 'r', 'LineWidth', 3)
xlim([1 100])
title('pChosen')

figure2; hold on
colors = cool(length(vTest));
for j = 1:length(vTest)
    plot(alpha{j}, 'color', colors(j,:))
end
plot(alphaSim, 'color', 'r', 'LineWidth', 3)
xlim([1 100])
title('alpha')

figure2; hold on
colors = cool(length(vTest));
for j = 1:length(vTest)
    plot(k{j}, 'color', colors(j,:))
end
plot(kSim, 'color', 'r', 'LineWidth', 3)
xlim([1 100])
title('k')

figure2; hold on
colors = cool(length(vTest));
for j = 1:length(vTest)
    plot(pe{j}, 'color', colors(j,:))
end
plot(peSim, 'color', 'r', 'LineWidth', 3)
xlim([1 100])
title('pe')
%% param gradient for real behavior
% load behavior and fit
model = 'vkf_fixV_kappa';
beh = 'good';
paramMatx = plotStanSessionParams({'ZS061'}, 'modelName', model, 'plotFlag', 0, 'saveFigFlag', 0, 'beh', beh);
i = 39;
xlFile = 'ZS61.xlsx';
filePath = [root 'ZS061' sep 'ZS061sorted' sep 'stan' sep 'bernoulli' sep model sep beh sep 'ZS061' beh '_' model '.mat'];
load(filePath, 'dayList');
session = dayList{i};
s = behAnalysisNoPlot_opMD(session);
paramsTest = paramMatx(i, :);
choice = s.allChoices;
choice(choice<0) = 0;
outcome = s.allRewards;
outcome = abs(outcome);
%%
% generate estimated values with gradient parameter
vTest = linspace(0.1,1,20);
omegaTest = 10*linspace(0.1,1,20);
LHTest = zeros(size(vTest));
probChosen = cell(length(vTest),1);
pRight = cell(length(vTest),1);
m = cell(length(vTest),1);
pe = cell(length(vTest),1);
v = cell(length(vTest),1);
w = cell(length(vTest),1);
alpha = cell(length(vTest),1);
k = cell(length(vTest),1);
[LHTestSim, probChosenSim,~,peSim,~,~,kSim,alphaSim,pRightSim] = vkf_LL(paramsTest, choice, outcome, 0);
% LH, probChosen, m, pe, v, w, k, alpha, probChoice
for j = 1:length(vTest)
%     paramsTest(2) = vTest(j);
    paramsTest(3) = omegaTest(j);
    [LHTest(j), probChosen{j},m{j},pe{j},v{j},w{j},k{j},alpha{j},pRight{j}] = vkf_LL(paramsTest, choice, outcome, 0);
end
%% plotting
figure2;
plot(omegaTest, LHTest)
figure2; hold on
colors = cool(length(vTest));
for j = 1:length(vTest)
    plot(pRight{j}, 'color', colors(j,:))
end
plot(pRightSim, 'color', 'r', 'LineWidth', 3)
xlim([1 100])
title('pRight')

figure2; hold on
colors = cool(length(vTest));
for j = 1:length(vTest)
    plot(probChosen{j}, 'color', colors(j,:))
end
plot(probChosenSim, 'color', 'r', 'LineWidth', 3)
xlim([1 100])
title('pChosen')

figure2; hold on
colors = cool(length(vTest));
for j = 1:length(vTest)
    plot(alpha{j}, 'color', colors(j,:))
end
plot(alphaSim, 'color', 'r', 'LineWidth', 3)
xlim([1 100])
title('alpha')

figure2; hold on
colors = cool(length(vTest));
for j = 1:length(vTest)
    plot(k{j}, 'color', colors(j,:))
end
plot(kSim, 'color', 'r', 'LineWidth', 3)
xlim([1 100])
title('k')

figure2; hold on
colors = cool(length(vTest));
for j = 1:length(vTest)
    plot(pe{j}, 'color', colors(j,:))
end
plot(peSim, 'color', 'r', 'LineWidth', 3)
xlim([1 100])
title('pe')
% plot
%%
m = 2.5;
sig = 0.4;
sigma = 2;
hpr = normrnd(0,2,10000,1);
hpr = m*ones(10000,1);
pr = normrnd(0,2,10000,1);
cauchy = makedist('tLocationScale','mu',0,'sigma',sigma);
halfcauchy = truncate(cauchy,0,inf);
hcs = random(halfcauchy, 10000, 1);
hcs = sig*ones(10000,1);
meow = 1./(1+exp(-(hpr + hcs .* pr)));
woof = normcdf(hpr + hcs .* pr);
figure2; subplot(1,2,1); histogram(meow, 0:0.02:1); xlim([0 1]);
subplot(1,2,2); histogram(woof,0:0.02:1); xlim([0 1]);
suptitle([num2str(m) '  ' num2str(sig) '  ' num2str(normpdf(m, 0, 3)*cauchypdf(sig, 0, 2))])
%%


