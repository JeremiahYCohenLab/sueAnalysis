function t = getModelVariablesLaser_dF(modelName, params, choice, outcome, laser)

t =  struct;
switch modelName
    case '5params'
        [t.LH, t.probChoice, t.Q, t.pe] = qLearningModel_5params(params, choice, outcome);
    case '5params_k_bias'
        [t.LH, t.probChoice, t.Q, t.pe] = qLearningModel_5params_k(params, choice, outcome);
    case '5paramsLaserNegRPE'
        [t.LH, t.probChoice, t.Q, t.pe, t.peChange] = qLearningModel_5paramsLaserNegRPE(params, choice, outcome, laser);
    case '5paramsLaserNegRPERotation'
        [t.LH, t.probChoice, t.Q, t.pe] = qLearningModel_5paramsLaserNegRPERotation(params, choice, outcome, laser);
    case '5params_k_bias_LaserNegRPE'
        [t.LH, t.probChoice, t.Q, t.pe, t.peOri] = qLearningModel_5params_k_bias_LaserNegRPE(params, choice, outcome, laser);
    case '5params_k_bias_LaserNegOnlyRPE'
        [t.LH, t.probChoice, t.Q, t.pe] = qLearningModel_5params_k_bias_LaserNegOnlyRPE(params, choice, outcome, laser);
    case '5params_k_bias_LaserDisengageScale'
        [t.LH, t.probChoice, t.Q, t.pe] = qLearningModel_5params_k_bias_LaserDisengageScale(params, choice, outcome, laser);
end 