function qLearning_fitMLE(session, modelName)
[root, sep] = currComputerAcute();
pd = parseSessionString_allen(session, root, sep);
s = behAnalysisNoPlot_opMDNWB(session);
choice = 0.5 *(1 + s.allChoices);
outcome = abs(s.allRewards);
              
% Set up optimization problem
options = optimset('Algorithm', 'interior-point','ObjectiveLimit',...
    -1.000000000e+300,'TolFun',1e-15, 'Display','off');

runs = 500;

% parameter range
alphaLearn_range = [0 1];
alphaNPE_range = [0 1];
alphaPPE_range = [0 1];
alphaForget_range = [0 1];
beta_range = [0 10];
bias_range = [-10 10];
kappa_range = [-5 5];


% start and fitted values for fitting
startValues = [];
startValues(:,1) = rand(runs,1);
startValues(:,2) = rand(runs,1);
startValues(:,3) = rand(runs,1);
startValues(:,4) = rand(runs,1)*range(beta_range);
startValues(:,5) = rand(runs,1)*range(bias_range) + bias_range(1);

A=[eye(size(startValues, 2)); -eye(size(startValues, 2))];

% initialize output variables

allParams = zeros(size(startValues, 1), size(startValues, 2));
LH = zeros(size(startValues, 1), 1);
exitFl = zeros(size(startValues, 1), 1);
hess = zeros(size(startValues, 1), size(startValues, 2), size(startValues, 2));
numParam = size(startValues, 2);

if strcmp(modelName, '5params')
    b=[ alphaNPE_range(2);  alphaPPE_range(2);  alphaForget_range(2); beta_range(2); bias_range(2);
       -alphaNPE_range(1); -alphaPPE_range(1); -alphaForget_range(1); -beta_range(1); -bias_range(1)];
    parfor r = 1:runs
        [allParams(r, :), LH(r, :), exitFl(r, :), ~, ~, ~, hess(r, :, :)] = ...
            fmincon(@qLearningModel_5params, startValues(r,:), A, b, [], [], [], [], [], options, choice, outcome);
    end
    [minLH,bestFit] = min(LH);
    paramEstimates = allParams(bestFit, :);
end

paramNames = getParamNames_dF(modelName, 1);

BIC = log(length(outcome))*numParam - 2*LH;

bestHess = squeeze(hess(bestFit, :, :));
CIvals = sqrt(diag(inv(bestHess)))'*1.96;
exitFl = exitFl(bestFit, :);
modelVars = getModelVariables_dF(modelName, paramEstimates, choice, outcome);
figure2;
scatterAll(allParams, {'aN', 'aP', 'aF', 'beta', 'bias'}, 10, 'k')
sgtitle(session, 'interpreter', 'none')    
save([pd.sortedFolder session '_' modelName  'MLE.mat'], "paramNames", "paramEstimates", 'bestFit', 'modelVars', 'CIvals')
