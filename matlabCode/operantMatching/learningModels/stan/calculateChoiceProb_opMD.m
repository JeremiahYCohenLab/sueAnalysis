function [model] = calculateChoiceProb_opMD(sessionName,varargin)

%task and model parameters
p = inputParser;
% default parameters if none given
p.addParameter('revForFlag', false);
p.addParameter('modelNames', {'fiveParamO', 'twoParams'})
p.addParameter('params', [0.0596149,0.305917,0.642195,3.31916,0.1]);
p.addParameter('plotFlag', 0);
p.parse(varargin{:});

modelNames = p.Results.modelNames;

filename = [sessionName '.asc'];
if p.Results.revForFlag == 1
    [behSessionData, unCorrectedBlockSwitch, out] = loadBehavioralData_revFor(filename);
    behavStruct = parseBehavioralData_revFor(behSessionData, unCorrectedBlockSwitch);
else
    [behSessionData, unCorrectedBlockSwitch, out] = loadBehavioralData(filename);
    behavStruct = parseBehavioralData(behSessionData, unCorrectedBlockSwitch);
end

outcome = abs([behavStruct.allReward_R; behavStruct.allReward_L])';
choice = abs([behavStruct.allChoice_R; behavStruct.allChoice_L])';
ITI = [behavStruct.timeBtwn]';


if p.Results.plotFlag
    normKern = normpdf(-15:15,0,4);
    normKern = normKern / sum(normKern);
    normKern = fliplr(normKern(1:16));
    allChoices = choice(:,1) - choice(:,2);
    figure; hold on;
    yyaxis left; plot(conv(allChoices',normKern)/max(conv(allChoices',normKern)),'k','linewidth',1);
    ylim([-1 1]); yticks([-1 0 1]); yticklabels([0 0.5 1])
    set(gca, 'Ycolor', 'k');
    
    numParams = length(modelNames);
    blue = [0 1 1];
    purp = [0.7 0 1];
    colors = [linspace(blue(1),purp(1),numParams)', linspace(blue(2),purp(2),numParams)', linspace(blue(3),purp(3),numParams)'];
end

for currMod = 1:length(modelNames)
    if strcmp(modelNames{currMod}, 'fiveParamO_rBarStart')
        [model.(modelNames{currMod}).LH, model.(modelNames{currMod}).probChoice, model.(modelNames{currMod}).Q, model.(modelNames{currMod}).pe, model.(modelNames{currMod}).rBar] = ...
            qLearningModel_5params_opponency(p.Results.params{currMod}, choice, outcome);
    end
    
    if strcmp(modelNames{currMod}, 'twoParams')
        [model.(modelNames{currMod}).LH, model.(modelNames{currMod}).probChoice, model.(modelNames{currMod}).Q, model.(modelNames{currMod}).pe] = ...
            qLearningModel_2params(p.Results.params{currMod}, choice, outcome);
    end
    
    if strcmp(modelNames{currMod}, 'fourParam')
        [model.(modelNames{currMod}).LH, model.(modelNames{currMod}).probChoice, model.(modelNames{currMod}).Q, model.(modelNames{currMod}).pe] = ...
            qLearningModel_4params_2learnRates_alphaForget(p.Results.params{currMod}, choice, outcome);
    end
    
    if strcmp(modelNames{currMod}, 'fiveParam_peBeta')
        [model.(modelNames{currMod}).LH, model.(modelNames{currMod}).probChoice, model.(modelNames{currMod}).Q, model.(modelNames{currMod}).pe,...
            model.(modelNames{currMod}).beta, model.(modelNames{currMod}).R] = ...
            qLearningModel_5params_peBeta(p.Results.params{currMod}, choice, outcome);
    end
    if strcmp(modelNames{currMod}, 'fourParam_rBeta_confQ')
        [model.(modelNames{currMod}).LH, model.(modelNames{currMod}).probChoice, model.(modelNames{currMod}).Q, model.(modelNames{currMod}).pe,...
            model.(modelNames{currMod}).beta, model.(modelNames{currMod}).R] = ...
            qLearningModel_4params_rBeta_confQ(p.Results.params{currMod}, choice, outcome);
    end
    if strcmp(modelNames{currMod}, 'fourParam_rBeta_scale')
        [model.(modelNames{currMod}).LH, model.(modelNames{currMod}).probChoice, model.(modelNames{currMod}).Q, model.(modelNames{currMod}).pe,...
            model.(modelNames{currMod}).beta, model.(modelNames{currMod}).R] = ...
            qLearningModel_4params_rBeta_scale(p.Results.params{currMod}, choice, outcome);
    end
    if strcmp(modelNames{currMod}, 'fiveParam_rBeta_scale')
        [model.(modelNames{currMod}).LH, model.(modelNames{currMod}).probChoice, model.(modelNames{currMod}).Q, model.(modelNames{currMod}).pe,...
            model.(modelNames{currMod}).beta, model.(modelNames{currMod}).R] = ...
            qLearningModel_5params_rBeta_scale(p.Results.params{currMod}, choice, outcome);
    end
    if strcmp(modelNames{currMod}, 'fiveParam_rBeta_scale_min')
        [model.(modelNames{currMod}).LH, model.(modelNames{currMod}).probChoice, model.(modelNames{currMod}).Q, model.(modelNames{currMod}).pe,...
            model.(modelNames{currMod}).beta, model.(modelNames{currMod}).R] = ...
            qLearningModel_6params_rBeta_scale_min(p.Results.params{currMod}, choice, outcome);
    end
    if strcmp(modelNames{currMod}, 'sixParam_rBeta_scale_initR')
        [model.(modelNames{currMod}).LH, model.(modelNames{currMod}).probChoice, model.(modelNames{currMod}).Q, model.(modelNames{currMod}).pe,...
            model.(modelNames{currMod}).beta, model.(modelNames{currMod}).R] = ...
            qLearningModel_6params_rBeta_scale_initR(p.Results.params{currMod}, choice, outcome);
    end
    if strcmp(modelNames{currMod}, 'sixParam_rBeta_scale_kappa')
        [model.(modelNames{currMod}).LH, model.(modelNames{currMod}).probChoice, model.(modelNames{currMod}).Q, model.(modelNames{currMod}).pe,...
            model.(modelNames{currMod}).beta, model.(modelNames{currMod}).R] = ...
            qLearningModel_6params_rBeta_scale_kappa(p.Results.params{currMod}, choice, outcome);
    end
    if strcmp(modelNames{currMod}, 'sixParam_peBeta')
        [model.(modelNames{currMod}).LH, model.(modelNames{currMod}).probChoice, model.(modelNames{currMod}).Q, model.(modelNames{currMod}).pe,...
            model.(modelNames{currMod}).beta, model.(modelNames{currMod}).R] = ...
            qLearningModel_6params_peBeta(p.Results.params{currMod}, choice, outcome);
    end
    
    if strcmp(modelNames{currMod}, 'sixParam_rBeta_rRPE')
        [model.(modelNames{currMod}).LH, model.(modelNames{currMod}).probChoice, model.(modelNames{currMod}).Q, model.(modelNames{currMod}).pe,...
            model.(modelNames{currMod}).beta, model.(modelNames{currMod}).R] = ...
            qLearningModel_6params_rBeta_rRPE(p.Results.params{currMod}, choice, outcome);
    end
    
    if strcmp(modelNames{currMod}, 'sixParam_rBeta_rV')
        [model.(modelNames{currMod}).LH, model.(modelNames{currMod}).probChoice, model.(modelNames{currMod}).Q, model.(modelNames{currMod}).pe,...
            model.(modelNames{currMod}).beta, model.(modelNames{currMod}).R] = ...
            qLearningModel_6params_rBeta_rV(p.Results.params{currMod}, choice, outcome);
    end
    
   if strcmp(modelNames{currMod}, 'sixParam_absPePeAN_bi')
        [model.(modelNames{currMod}).LH, model.(modelNames{currMod}).probChoice, model.(modelNames{currMod}).Q, model.(modelNames{currMod}).pe,...
            model.(modelNames{currMod}).pePe, model.(modelNames{currMod}).aN, model.(modelNames{currMod}).peBar] = ...
            qLearningModel_6params_absPePeAN_bi(p.Results.params{currMod}, choice, outcome);
    end
    
    
    if p.Results.plotFlag
        yyaxis right; 
        eval(['plot(model.' modelNames{currMod} '.probChoice(:,1), ''-'', ''Color'', [' num2str(colors(currMod,:)) '],''linewidth'', 1)']);
        ylim([0 1]); yticks([0 0.5 1]);
        set(gca, 'Ycolor', 'k');
    end

end

if p.Results.plotFlag
    ax = gca;
    ax.TickDir = 'out';
    set(gcf, 'Position', [-1910 159 1905 750]);
    set(gcf, 'Renderer', 'Painters')
    legend(['choices', modelNames], 'Interpreter', 'none')
    suptitle(sessionName)
end



