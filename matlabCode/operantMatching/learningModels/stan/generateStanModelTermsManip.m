function [t] = generateStanModelTermsManip(modelType, modelPath, sessionName, preFlag, sessionFlag, revForFlag)

if nargin < 5
    sessionFlag = 1;
end

if nargin < 6
    revForFlag = 0;
end

if strfind(modelType, 'four')
    numParams = 8;
elseif strfind(modelType, 'five')
    numParams = 5;
elseif strfind(modelType, 'six')
    numParams = 6;
end

%get model params
load(modelPath);
tmp = whos;
samples = eval(tmp(1).name);
mdlFields = fields(samples);
tmpInd = find(~cellfun(@isempty,strfind(mdlFields,'log_lik')),1);
if sessionFlag
    sessionInd = find(~cellfun(@isempty,strfind(dayList,sessionName)));
    paramInds(1,:) = [tmpInd-numParams*4:tmpInd-numParams*3-1];
    paramInds(2,:) = paramInds(1,:) + numParams;
else
    paramInds(1,:) = [tmpInd-numParams*2:tmpInd-numParams-1];
    paramInds = [paramInds; paramInds(1,:) + numParams];
end

if preFlag
    paramInds = paramInds(1,:);
else
    paramInds = paramInds(2,:);
    sessionInd = sessionInd - size(dayList,1);
end

for j = 1:length(paramInds)
    if sessionFlag
        tmp = eval(['samples.' mdlFields{paramInds(1,j)}]);
        startValues(j) = median(tmp(:,sessionInd(1)));
    else
        startValues(j) = median(eval(['samples.mu_' mdlFields{paramInds(j)}]));
    end
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
    
    t.o = abs(o.allRewards)';
    
end

end