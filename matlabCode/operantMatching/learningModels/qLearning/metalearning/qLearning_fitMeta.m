function [model] = qLearning_fitMeta(filename, varargin)

p = inputParser;
% default parameters if none given
p.addParameter('figFlag', 0)
p.addParameter('testFlag', 0)
p.addParameter('revForFlag',0)
p.addParameter('mouse',[])
p.addParameter('category', [])
p.parse(varargin{:});


if ~isempty(p.Results.mouse)
    behavStruct = parseBehavioralData_multiple(filename, p.Results.mouse, p.Results.category, p.Results.revForFlag);
else
    [behSessionData, unCorrectedBlockSwitch, out] = loadBehavioralData(filename, p.Results.revForFlag);
    behavStruct = parseBehavioralData(behSessionData, unCorrectedBlockSwitch);
end


outcome = abs([behavStruct.allReward_R; behavStruct.allReward_L])';
choice = abs([behavStruct.allChoice_R; behavStruct.allChoice_L])';
ITI = [behavStruct.timeBtwn]';

% Initialize models
modelNames = {'fourParams_twoLearnRates_alphaForget', 'fiveParams_rBeta_scale'};
startValueCSVs = {'qLearningModel_4params_2learnRates_alphaForget_startValues.csv',...
                    'qLearningModel_5params_rBeta_scale_startValues.csv',};



% Set up optimization problem
options = optimset('Algorithm', 'interior-point','ObjectiveLimit',...
    -1.000000000e+300,'TolFun',1e-15, 'Display','off');

alpha_range = [0 1];
alphaNPE_range = [0 1];
alphaPPE_range = [0 1];
alphaForget_range = [0 1];
beta_range = [0 200];
betaRate_range = [0 1];
betaMin_range = [0 5];
betaMax_range = [6 30];
v_range = [0 1];
rBarStart_range = [0 1]; 
tForget_range = [0 1];
bias_range = [-5 5];

if p.Results.figFlag
    figure; hold on;
end
for currMod = 1:length(modelNames)
    startValues = csvread(startValueCSVs{currMod});
    if p.Results.testFlag == 1
        startValues = startValues(1, :);
    end
   
    % initialize output variables
    runs = size(startValues, 1);
    allParams = zeros(size(startValues, 1), size(startValues, 2));
    LH = zeros(size(startValues, 1), 1);
    exitFl = zeros(size(startValues, 1), 1);
    hess = zeros(size(startValues, 1), size(startValues, 2), size(startValues, 2));
    numParam = size(startValues, 2);
    A=[eye(size(startValues, 2)); -eye(size(startValues, 2))];
    

    if strcmp(modelNames{currMod}, 'fourParams_twoLearnRates_alphaForget')
        b=[ alphaNPE_range(2);  alphaPPE_range(2);  alphaForget_range(2); beta_range(2);
           -alphaNPE_range(1); -alphaPPE_range(1); -alphaForget_range(1); -beta_range(1)];
        parfor r = 1:runs
            [allParams(r, :), LH(r, :), exitFl(r, :), ~, ~, ~, hess(r, :, :)] = ...
                fmincon(@qLearningModel_4params_2learnRates_alphaForget, startValues(r,:), A, b, [], [], [], [], [], options, choice, outcome);
        end
        [~,bestFit] = min(LH);
        model.(modelNames{currMod}).bestParams = allParams(bestFit, :);
        [~, model.(modelNames{currMod}).probChoice, model.(modelNames{currMod}).Q, model.(modelNames{currMod}).pe] = ...
            qLearningModel_4params_2learnRates_alphaForget(model.(modelNames{currMod}).bestParams, choice, outcome);
    end
    if strcmp(modelNames{currMod}, 'fiveParams_opponency')
        b=[ alphaNPE_range(2);  alphaPPE_range(2);  alphaForget_range(2); beta_range(2); v_range(2);
           -alphaNPE_range(1); -alphaPPE_range(1); -alphaForget_range(1); -beta_range(1); -v_range(1)];
        parfor r = 1:runs
            [allParams(r, :), LH(r, :), exitFl(r, :), ~, ~, ~, hess(r, :, :)] = ...
                fmincon(@qLearningModel_5params_opponency, startValues(r,:), A, b, [], [], [], [], [], options, choice, outcome);
        end
        [~,bestFit] = min(LH);
        model.(modelNames{currMod}).bestParams = allParams(bestFit, :);
        [~, model.(modelNames{currMod}).probChoice, model.(modelNames{currMod}).Q, model.(modelNames{currMod}).pe, model.(modelNames{currMod}).rBar] = ...
            qLearningModel_5params_opponency(model.(modelNames{currMod}).bestParams, choice, outcome);
    end
    if strcmp(modelNames{currMod}, 'fiveParams_peBeta')
        b=[ alpha_range(2);  alphaForget_range(2); betaRate_range(2); betaMin_range(2); betaMax_range(2);
           -alpha_range(1); -alphaForget_range(1); -betaRate_range(1); -betaMin_range(1); -betaMax_range(1)];
        parfor r = 1:runs
            [allParams(r, :), LH(r, :), exitFl(r, :), ~, ~, ~, hess(r, :, :)] = ...
                fmincon(@qLearningModel_5params_peBeta, startValues(r,:), A, b, [], [], [], [], [], options, choice, outcome);
        end
        [~,bestFit] = min(LH);
        model.(modelNames{currMod}).bestParams = allParams(bestFit, :);
        [~, model.(modelNames{currMod}).probChoice, model.(modelNames{currMod}).Q, model.(modelNames{currMod}).pe,...
            model.(modelNames{currMod}).beta, model.(modelNames{currMod}).R ] = ...
            qLearningModel_5params_peBeta(model.(modelNames{currMod}).bestParams, choice, outcome);
    end
    if strcmp(modelNames{currMod}, 'sixParams_peBeta')
        b=[ alphaNPE_range(2);  alphaPPE_range(2);  alphaForget_range(2); betaRate_range(2); betaMin_range(2); betaMax_range(2);
           -alphaNPE_range(1); -alphaPPE_range(1); -alphaForget_range(1); -betaRate_range(1); -betaMin_range(1); -betaMax_range(2)];
        parfor r = 1:runs
            [allParams(r, :), LH(r, :), exitFl(r, :), ~, ~, ~, hess(r, :, :)] = ...
                fmincon(@qLearningModel_6params_peBeta, startValues(r,:), A, b, [], [], [], [], [], options, choice, outcome);
        end
        [~,bestFit] = min(LH);
        model.(modelNames{currMod}).bestParams = allParams(bestFit, :);
        [~, model.(modelNames{currMod}).probChoice, model.(modelNames{currMod}).Q, model.(modelNames{currMod}).pe,...
            model.(modelNames{currMod}).beta, model.(modelNames{currMod}).R ] = ...
            qLearningModel_6params_peBeta(model.(modelNames{currMod}).bestParams, choice, outcome);
    end
    if strcmp(modelNames{currMod}, 'fiveParams_rBeta_scale')
        b=[ alphaNPE_range(2);  alphaPPE_range(2);  alphaForget_range(2); v_range(2); betaRate_range(2);
           -alphaNPE_range(1); -alphaPPE_range(1); -alphaForget_range(1); -v_range(1); -betaRate_range(1)];
        parfor r = 1:runs
            [allParams(r, :), LH(r, :), exitFl(r, :), ~, ~, ~, hess(r, :, :)] = ...
                fmincon(@qLearningModel_5params_rBeta_scale, startValues(r,:), A, b, [], [], [], [], [], options, choice, outcome);
        end
        [~,bestFit] = min(LH);
        model.(modelNames{currMod}).bestParams = allParams(bestFit, :);
        [~, model.(modelNames{currMod}).probChoice, model.(modelNames{currMod}).Q, model.(modelNames{currMod}).pe,...
            model.(modelNames{currMod}).beta, model.(modelNames{currMod}).R ] = ...
            qLearningModel_5params_rBeta_scale(model.(modelNames{currMod}).bestParams, choice, outcome);
    end

    model.(modelNames{currMod}).LH = -1 * LH(bestFit, :);
    model.(modelNames{currMod}).BIC = log(length(outcome))*numParam - 2*model.(modelNames{currMod}).LH;
    
    bestHess = squeeze(hess(bestFit, :, :));
    model.(modelNames{currMod}).CIvals = sqrt(diag(inv(bestHess)))'*1.96;
    model.(modelNames{currMod}).exitFl = exitFl(bestFit, :);
    
    if p.Results.figFlag
        plot(model.(modelNames{currMod}).probChoice(:,1)/max(model.(modelNames{currMod}).probChoice(:,1)), 'linewidth', 2);
    end
end

model.choice = choice;
model.outcome = outcome;

if p.Results.figFlag
    for i = 1:length(choice)
        if choice(i,1) == 1
            if outcome(i,1) == 1
                plot([i i], [0.8 1], '-k')
            else
                plot([i i], [0.8 0.9], '-k')
            end
        else
            if outcome(i,2) == 1
                plot([i i], [0 0.2], '-k')
            else
                plot([i i], [0.1 0.2], '-k')
            end
        end
    end
    str = [modelNames];
    str = regexprep(str, '\_', ' ');
    legend(str)
end