function [t] = generateStanModelTerms_opMD(modelType, modelPath, sessionName, sessionFlag, revForFlag)

if nargin < 4
    sessionFlag = 1;
end

if nargin < 5
    revForFlag = 0;
end


%get model params
load(modelPath);
if sessionFlag
    sessionInd = find(~cellfun(@isempty,strfind(dayList,sessionName)));
end
tmp = whos;
samples = eval(tmp(1).name);
mdlFields = fields(samples);

if any(~cellfun(@isempty,strfind(mdlFields,'bias')))
    paramInds = find(~cellfun(@isempty,strfind(mdlFields,'mu_')));
    paramInds = paramInds(2:end);
    paramInds = paramInds - length(paramInds) - 1;
else
    paramInds = find(~cellfun(@isempty,strfind(mdlFields,'mu_')));
    paramInds = paramInds(2:end);
    paramInds = paramInds - length(paramInds);
end

for j = 1:length(paramInds)
    if sessionFlag
        tmp = eval(['samples.' mdlFields{paramInds(j)}]);
        startValues(j) = median(tmp(:,sessionInd));
    else
        startValues(j) = median(eval(['samples.mu_' mdlFields{paramInds(j)}]));
    end
end 

if sessionFlag && isfield(samples, 'bias')
    startValues = [startValues median(samples.bias(:,sessionInd))];
end
    
t = struct;
t.params = startValues;

if ~isempty(sessionName)
    %get session behavior
    [behSessionData,blockSwitch,~] = loadBehavioralData([sessionName '.asc'], revForFlag);
    o = parseBehavioralData(behSessionData, blockSwitch);
    outcome = abs([o.allReward_R; o.allReward_L])';
    choice = abs([o.allChoice_R; o.allChoice_L])';


    switch modelType
        case 'twoParam'
            [t.LH, t.probChoice, t.Q, t.pe] = qLearningModel_2params(startValues, choice, outcome);
        case 'threeParam'
            [t.LH, t.probChoice, t.Q, t.pe] = qLearningModel_3params_alphaForget(startValues, choice, outcome);
        case 'fourParam'
            [t.LH, t.probChoice, t.Q, t.pe, t.cQ] = qLearningModel_4params_2learnRates_alphaForget(startValues, choice, outcome);
        case 'fiveParamO'
            [t.LH, t.probChoice, t.Q, t.pe, t.rBar] = qLearningModel_5params_opponency(startValues, choice, outcome);
        case 'fiveParamO_peUpdate'
            [t.LH, t.probChoice, t.Q, t.pe, t.rBar] = qLearningModel_5params_opponency_peUpdate(startValues, choice, outcome);
        case 'fiveParam_bilal'
            [t.LH, t.probChoice, t.Q, t.pe] = qLearningModel_5params_bilal(startValues, choice, outcome);
         case 'fiveParam_peBeta_avg'
            [t.LH, t.probChoice, t.Q, t.pe, t.beta, t.R] = qLearningModel_5params_peBeta_avg(startValues, choice, outcome);
        case 'sixParam_peBeta'
            [t.LH, t.probChoice, t.Q, t.pe, t.beta, t.R] = qLearningModel_6params_peBeta(startValues, choice, outcome);
        case 'sixParam_t.peBeta_diff'
            [t.LH, t.probChoice, t.Q, t.pe, t.beta, t.R] = qLearningModel_6params_peBeta_diff(startValues, choice, outcome);
        case 'sevenParam_peBeta_k'
            [t.LH, t.probChoice, t.Q, t.pe, t.beta, t.R] = qLearningModel_7params_peBeta_k(startValues, choice, outcome);
        case 'fourParam_rBeta_confQ'
            [t.LH, t.probChoice, t.Q, t.pe, t.beta, t.R] = qLearningModel_4params_rBeta_confQ(startValues, choice, outcome);
        case 'fiveParam_rBeta_scale'
            [t.LH, t.probChoice, t.Q, t.pe, t.beta, t.R] = qLearningModel_5params_rBeta_scale(startValues, choice, outcome);
        case 'sixParam_rBeta_scale_min'
            [t.LH, t.probChoice, t.Q, t.pe, t.beta, t.R] = qLearningModel_6params_rBeta_scale_min(startValues, choice, outcome);
        case 'fiveParam_rBeta_confQ'
            [t.LH, t.probChoice, t.Q, t.pe, t.beta, t.R] = qLearningModel_5params_rBeta_scale(startValues, choice, outcome);
        case 'sixParam_rBeta_rRPE'
            [t.LH, t.probChoice, t.Q, t.pe, t.beta, t.R] = qLearningModel_6params_rBeta_rRPE(startValues, choice, outcome);
        case 'sixParam_rBeta_rV'
            [t.LH, t.probChoice, t.Q, t.pe, t.beta, t.R] = qLearningModel_6params_rBeta_rV(startValues, choice, outcome);
        case 'sixParam_peAN'
            [t.LH, t.probChoice, t.Q, t.pe, t.peBar, t.aN] = qLearningModel_6params_peAN(startValues, choice, outcome);
        case 'sixParam_pePeAN'
            [t.LH, t.probChoice, t.Q, t.pe, t.pePe, t.aN, t.peBar] = qLearningModel_6params_pePeAN(startValues, choice, outcome);
        case 'sixParam_pePeAN_lag'
            [t.LH, t.probChoice, t.Q, t.pe, t.pePe, t.aN, t.peBar] = qLearningModel_6params_pePeAN_lag(startValues, choice, outcome);
        case 'fiveParam_absPePeAN'
            [t.LH, t.probChoice, t.Q, t.pe, t.pePe, t.aN, t.peBar] = qLearningModel_5params_absPePeAN(startValues, choice, outcome);
        case 'sixParam_absPePeAN'
            [t.LH, t.probChoice, t.Q, t.pe, t.pePe, t.aN, t.peBar] = qLearningModel_6params_absPePeAN(startValues, choice, outcome);
        case 'sixParam_absPePeAN_bi'
            [t.LH, t.probChoice, t.Q, t.pe, t.pePe, t.aN, t.peBar] = qLearningModel_6params_absPePeAN_bi(startValues, choice, outcome);
        case 'sixParam_absPePeAN_biSep'
            [t.LH, t.probChoice, t.Q, t.pe, t.pePe, t.aN, t.peBar, t.peBar_L, t.peBar_R] = qLearningModel_6params_absPePeAN_biSep(startValues, choice, outcome);
        case 'sevenParam_absPePeAN_biSep_f'
            [t.LH, t.probChoice, t.Q, t.pe, t.pePe, t.aN, t.peBar, t.peBar_L, t.peBar_R] = qLearningModel_7params_absPePeAN_biSep_f(startValues, choice, outcome);
        case 'sevenParam_absPePeAN_bi_k'
            [t.LH, t.probChoice, t.Q, t.pe, t.pePe, t.aN, t.peBar] = qLearningModel_7params_absPePeAN_bi_k(startValues, choice, outcome);
        case 'sevenParam_absPePeLR'
            [t.LH, t.probChoice, t.Q, t.pe, t.pePe, t.peBar] = qLearningModel_7params_absPePeLR(startValues, choice, outcome);
        case 'sevenParam_peLR'
            [t.LH, t.probChoice, t.Q, t.pe, t.R, t.aN, t.aP] = qLearningModel_7params_peLR(startValues, choice, outcome);
        case 'sixParam_rAN'
            [t.LH, t.probChoice, t.Q, t.pe, t.R, t.aN] = qLearningModel_6params_rAN(startValues, choice, outcome);
        case 'eightParam_rBeta_peAN'
            [t.LH, t.probChoice, t.Q, t.pe, t.R, t.aN, t.peBar] = qLearningModel_8params_rBeta_peAN(startValues, choice, outcome);
        case 'eightParam_rBeta_pePeAN'
            [t.LH, t.probChoice, t.Q, t.pe, t.R, t.aN, t.peBar] = qLearningModel_8params_rBeta_pePeAN(startValues, choice, outcome);
    end 
    
end

end