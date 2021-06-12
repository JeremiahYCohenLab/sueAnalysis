function model = idealObserverSuccessFit(modelNames, varargin)

p = inputParser;
% default parameters if none given
p.addParameter('testFlag', 0)   %limits fitting to 1 run
p.addParameter('numVals', 40)   %number of starting values and runs
p.addParameter('rSeed', 241657)
p.parse(varargin{:});
    
%initialize model start value ranges
startValueRanges = [];
for currM = 1:length(modelNames)
    switch modelNames{currM}
        case 'vkf'                    %lambda   v0   omega  beta
            startValueRanges{currM} = [0.1 1; 0.1 1; 0.1 1; 1 20];
    end
end
              
%set up optimization problem
% options = optimset('Algorithm', 'interior-point', 'ObjectiveLimit', 0,...
%     'TolFun',1e-15, 'Display','iter');
ms = MultiStart('FunctionTolerance',2e-4,'UseParallel',true, 'StartPointsToRun', 'bounds-ineqs');

%set ranges for parameter estimates

lambda_range = [0 1];  % volatility update rate
v0_range = [0 1];     % initial volatility
omega_range = [0 1];   % observation noise
beta_range = [0 20];    % inverse-temp-like parameter for softmax decision function

%fit models one at a time
for currMod = 1:length(modelNames)
    %generate start values from defined ranges
    startValues = createStartValues(startValueRanges{currMod}, p.Results.numVals);
    if p.Results.testFlag == 1
        startValues = startValues(1, :);    %only one run if testing
    end
    
    %initialize output variables
    runs = size(startValues, 1);
    allParams = zeros(size(startValues, 1), size(startValues, 2));
    succMin = zeros(size(startValues, 1), 1);
    exitFl = zeros(size(startValues, 1), 1);
    hess = zeros(size(startValues, 1), size(startValues, 2), size(startValues, 2));
    numParam = size(startValues, 2);

    A=[eye(size(startValues, 2)); -eye(size(startValues, 2))];
    if strcmp(modelNames{currMod}, 'vkf')
        lb = [lambda_range(1); v0_range(1); omega_range(1); beta_range(1)];
        ub = [lambda_range(2); v0_range(2); omega_range(2); beta_range(2)];
        for r = 1:runs
            problem = createOptimProblem('fmincon', 'x0', startValues(r, :), 'objective', @vkf, 'lb', lb, 'ub', ub);
            [allParams(r, :), succMin(r, :)] = run(ms, problem, 20);
%             [allParams(r, :), succMin(r, :), exitFl(r, :), ~, ~, ~, hess(r, :, :)] = ...
%                 fmincon(@vkf, startValues(r, :), A, b, [], [], [], [], [], options, p.Results.rSeed);
        end
        [~, bestFit] = min(succMin);
        model.(modelNames{currMod}).bestParams = allParams(bestFit, :);
        [~, model.(modelNames{currMod}).succ] = ...
            vkf(model.(modelNames{currMod}).bestParams);
    else 
        error('model name not found')
    end

    model.(modelNames{currMod}).succMin = succMin(bestFit, :);
    
    bestHess = squeeze(hess(bestFit, :, :));
    model.(modelNames{currMod}).CIvals = sqrt(diag(inv(bestHess)))'*1.96;
    model.(modelNames{currMod}).exitFl = exitFl(bestFit, :);
end


end
