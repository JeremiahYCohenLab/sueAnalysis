function [t] = generateStanModelTermsManip_opMD(modelType, modelPath, sessionName, sessionFlag, revForFlag)

if nargin < 4
    sessionFlag = 1;
end

if nargin < 5
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
    [sessionInd(1),~] = find(~cellfun(@isempty,strfind(dayList,sessionName{1})));
    [sessionInd(2),~] = find(~cellfun(@isempty,strfind(dayList,sessionName{2})));
    paramInds(1,:) = [tmpInd-numParams*4:tmpInd-numParams*3-1];
    paramInds(2,:) = paramInds(1,:) + numParams;
else
    paramInds(1,:) = [tmpInd-numParams*2:tmpInd-numParams-1];
    paramInds = [paramInds; paramInds(1,:) + numParams];
end

for j = 1:length(paramInds)
    if sessionFlag
        tmp = eval(['samples.' mdlFields{paramInds(1,j)}]);
        startValues{1}(j) = median(tmp(:,sessionInd(1)));
        tmp = eval(['samples.' mdlFields{paramInds(2,j)}]);
        startValues{2}(j) = median(tmp(:,sessionInd(2)));
    else
        startValues(j) = median(eval(['samples.mu_' mdlFields{paramInds(j)}]));
    end
end 

t = struct;
t.params = startValues;

for i = 1:2
    %get session behavior
    [behSessionData,blockSwitch,~] = loadBehavioralData([sessionName{1} '.asc'], revForFlag);
    o = parseBehavioralData(behSessionData, blockSwitch);
    outcome = abs([o.allReward_R; o.allReward_L])';
    choice = abs([o.allChoice_R; o.allChoice_L])';


    switch modelType
        case 'fourParam'
            [t.LH{i}, t.probChoice{i}, t.Q{i}, t.pe{i}] = qLearningModel_4params_2learnRates_alphaForget(startValues{i}, choice, outcome);
        case 'fiveParamO'
            [t.LH{i}, t.probChoice{i}, t.Q{i}, t.pe{i}, t.rBar{i}] = qLearningModel_5params_opponency(startValues{i}, choice, outcome);
        case 'fiveParamO_peUpdate'
            [t.LH{i}, t.probChoice{i}, t.Q{i}, t.pe{i}, t.rBar{i}] = qLearningModel_5params_opponency_peUpdate(startValues{i}, choice, outcome);
         case 'fiveParam_peBeta_avg'
            [t.LH{i}, t.probChoice{i}, t.Q{i}, t.pe{i}, t.beta{i}, t.R{i}] = qLearningModel_5params_peBeta_avg(startValues{i}, choice, outcome);
        case 'sixParam_peBeta'
            [t.LH{i}, t.probChoice{i}, t.Q{i}, t.pe{i}, t.beta{i}, t.R{i}] = qLearningModel_6params_peBeta(startValues{i}, choice, outcome);
        case 'sixParam_t.pe{i}Beta_diff'
            [t.LH{i}, t.probChoice{i}, t.Q{i}, t.pe{i}, t.beta{i}, t.R{i}] = qLearningModel_6params_peBeta_diff(startValues{i}, choice, outcome);
        case 'fourParam_rBeta_confQ'
            [t.LH{i}, t.probChoice{i}, t.Q{i}, t.pe{i}, t.beta{i}, t.R{i}] = qLearningModel_4params_rBeta_confQ(startValues{i}, choice, outcome);
        case 'fiveParam_rBeta_scale'
            [t.LH{i}, t.probChoice{i}, t.Q{i}, t.pe{i}, t.beta{i}, t.R{i}] = qLearningModel_5params_rBeta_scale(startValues{i}, choice, outcome);
        case 'sixParam_rBeta_scale_min'
            [t.LH{i}, t.probChoice{i}, t.Q{i}, t.pe{i}, t.beta{i}, t.R{i}] = qLearningModel_6params_rBeta_scale_min(startValues{i}, choice, outcome);
        case 'fiveParam_rBeta_confQ'
            [t.LH{i}, t.probChoice{i}, t.Q{i}, t.pe{i}, t.beta{i}, t.R{i}] = qLearningModel_5params_rBeta_scale(startValues{i}, choice, outcome);
        case 'sixParam_rBeta_rRPE'
            [t.LH{i}, t.probChoice{i}, t.Q{i}, t.pe{i}, t.beta{i}, t.R{i}] = qLearningModel_6params_rBeta_rRPE(startValues{i}, choice, outcome);
        case 'sixParam_rBeta_rV'
            [t.LH{i}, t.probChoice{i}, t.Q{i}, t.pe{i}, t.beta{i}, t.R{i}] = qLearningModel_6params_rBeta_rV(startValues{i}, choice, outcome);
        case 'sixParam_peAN_avg'
            [t.LH{i}, t.probChoice{i}, t.Q{i}, t.pe{i}, t.R{i}, t.aN{i}] = qLearningModel_6params_peAN_avg(startValues{i}, choice, outcome);
        case 'sixParam_absPePeAN_bi'
            [t.LH{i}, t.probChoice{i}, t.Q{i}, t.pe{i}, t.pePe{i}, t.aN{i}, t.peBar{i}] = qLearningModel_6params_absPePeAN_bi(startValues{i}, choice, outcome);
        case 'sevenParam_peLR_avg'
            [t.LH{i}, t.probChoice{i}, t.Q{i}, t.pe{i}, t.R{i}, t.aN{i}, t.aP{i}] = qLearningModel_7params_peLR_avg(startValues{i}, choice, outcome);
        case 'sixParam_rAN'
            [t.LH{i}, t.probChoice{i}, t.Q{i}, t.pe{i}, t.R{i}, t.aN{i}] = qLearningModel_6params_rAN(startValues{i}, choice, outcome);
    end 
    
end

end