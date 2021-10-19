%% cat foraging task simulation
samp = 2000;
foodSize = rand(samp,1);
foodAmt = foodSize + 0.3*rand(samp,1);
foodQual = rand(samp,1);
foodMat = [ones(samp,1),foodSize,foodQual];
% speed = b0 + b1*foodSize + b2*foodQual + e; e ~ norm(0,sigma^2);
w = [1;1;1];
sigma = 1;
speed = foodMat*w + normrnd(0,sigma^2,samp,1);


%% lm
% wEstimate = (x'x)^(-1)*x'y; from mle
foodMat = [ones(samp,1),foodAmt,foodQual];
wEstimate = inv(foodMat'*foodMat)*foodMat'*speed;
speedEstimate = foodMat*wEstimate;
wSE = diag(sqrt(inv(foodMat'*foodMat)*norm(speedEstimate - speed)^2/(samp - size(foodMat,2))));
wTStats = wEstimate./wSE;

lm = fitlm(foodMat, speed, 'Intercept',false);

%% Correlated regressors
foodMat = [ones(samp,1),foodSize,foodQual];
lmCorr = fitlm(foodMat, speed, 'Intercept',false);
foodMat = [ones(samp,1),foodSize,foodAmt, foodQual];
lmCorrAll = fitlm(foodMat, speed, 'Intercept',false);

%% GLM
% choice logistic model (GLM)
% distribution and link function difference
load('~yourPath/s.mat')
%%
tMax = 10;
allChoices = s.allChoices;
allChoices(allChoices<0) = 0; % 1 or right choice, 0 for left
allRewards = s.allRewards; % 1 for right reward, -1 for left

rwdMatx = [];
for i = 1:tMax
    rwdMatx(i,:) = [NaN(1,i) allRewards(1:end-i)];
end

glm_choices = fitglm(rwdMatx', allChoices, 'distribution','binomial','link','logit');
% plotting
figure;
relevInds = 2:tMax+1;
coefVals = glm_choices.Coefficients.Estimate(relevInds);
CIbands = coefCI(glm_choices);
errorL = abs(coefVals - CIbands(relevInds,1));
errorU = abs(coefVals - CIbands(relevInds,2));
errorbar((1:tMax)+0.2,coefVals,errorL,errorU,'Color','c','linewidth',2)
line([0 tMax+1], [0 0],'color', [0.5 0.5 0.5])
xlabel('Outcome n trials back')
ylabel('\beta coeffs')

%% reinforcement learning models
% generative vs non-generative models 
iteration = 200;
runs = 200;
maxTrial = 300;
paramNames = {'aN','aP','aF','beta'};
% generate random paramters
params = zeros(iteration, 4);
params(:,1) = betarnd(3, 5, iteration,1); % aN
params(:,2) = betarnd(5, 2, iteration,1); % aP
params(:,3) = betarnd(4, 2, iteration,1); % aF
params(:,4) = normrnd(4, 2, iteration,1); % beta
params(params(:,4)<1,4) = 2; % beta


% constraints for fitting
alphaNPE_range = [0 1];
alphaPPE_range = [0 1];
alphaForget_range = [0 1];
beta_range = [0 10];

% start and fitted values for fitting
startValues = [];
startValues(:,1) = linspace(0,1,runs);
startValues(:,2) = linspace(0,1,runs);
startValues(:,3) = linspace(0,1,runs);
startValues(:,4) = linspace(0,1,runs)*range(beta_range);

A=[eye(size(startValues, 2)); -eye(size(startValues, 2))];
b=[ alphaNPE_range(2);  alphaPPE_range(2);  alphaForget_range(2); beta_range(2);
   -alphaNPE_range(1); -alphaPPE_range(1); -alphaForget_range(1); -beta_range(1)];

paramsEstimates = zeros(iteration, 4);

% algorithm specs 
options = optimset('Algorithm', 'interior-point','ObjectiveLimit',...
        -1.000000000e+300,'TolFun',1e-15, 'Display','off');
    
figure;
for i = 1:4
subplot(2,2,i); hold on;
histogram(params(:,i),linspace(min(params(:,i)), max(params(:,i)), 20), 'Normalization', 'Probability','FaceColor', 'c');
histogram(startValues(:,i),linspace(min(startValues(:,i)), max(startValues(:,i)), 20), 'Normalization', 'Probability','FaceColor', 'm', 'FaceAlpha', 0.2, 'EdgeAlpha', 0.3);
title(paramNames{i})
if i == 1
    legend({'simulationParams', 'startValues'});
end
end
%%
Qsim = [];
Qfit = [];
psim = [];
pfit = [];
for sim = 1:iteration
    %simulation
    [outcome, choice] = qLearningModel_simNoPlot('params', params(sim,:),'randomSeed', sim,'maxTrials', maxTrial);
    choice(choice<0) = 0;
    outcome = abs(outcome);
   [~,ptemp,Qtemp] = qLearningModel_4params(params(sim,:), choice, abs(outcome));
    Qsim = [Qsim; Qtemp];
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
            fmincon(@qLearningModel_4params, startValues(r,:), A, b, [], [], [], [], [], options, choice, outcome);
    end
    [~,bestFit] = min(LH);
    [~,ptemp,Qtemp] = qLearningModel_4params(allParams(bestFit, :), choice, outcome);
    Qfit =[Qfit; Qtemp];
    pfit = [pfit; ptemp];
    paramsEstimates(sim,:) = allParams(bestFit, :);
end
%%
figure2;
for i = 1:4
subplot(2,2,i); hold on;
scatter(params(:,i), paramsEstimates(:,i),15, [0.3 1 1], 'filled');
scatter(params(paramsEstimates(:,4)>9.5,i), paramsEstimates(paramsEstimates(:,4)>9.5,i), 15, 'r', 'filled');
scatter(params(paramsEstimates(:,3)<0.05,i), paramsEstimates(paramsEstimates(:,3)<0.05,i), 15, 'm', 'filled');
line([minmax(paramsEstimates(:,i))], [minmax(paramsEstimates(:,i))], 'color', [0.5 0.5 0.5]);
title(paramNames{i})
end
suptitle(['maxTrial = ' num2str(maxTrial)])
figure2;
subplot(2,1,1)
scatter(Qsim(:,2)-Qsim(:,1), Qfit(:,2)-Qfit(:,1), 15, 'c', 'filled');
title('Qdiff')
subplot(2,1,2)
scatter(psim, pfit, 15, 'm', 'filled');
title('pRight')
suptitle(['maxTrial = ' num2str(maxTrial)])
 %%