function t = getModelVariables_dF(modelName, params, choice, outcome)

t =  struct;
switch modelName
    case 'twoParam'
        [t.LH, t.probChoice, t.Q, t.pe] = qLearningModel_2params(params, choice, outcome);
    case 'threeParam'
        [t.LH, t.probChoice, t.Q, t.pe] = qLearningModel_3params_alphaForget(params, choice, outcome);
    case 'threeParam_bias'
        [t.LH, t.probChoice, t.Q, t.pe] = qLearningModel_3params_bias(params, choice, outcome);
    case 'fourParam'
        [t.LH, t.probChoice, t.Q, t.pe, t.cQ] = qLearningModel_4params_2learnRates_alphaForget(params, choice, outcome);
    case 'fiveParamO'
        [t.LH, t.probChoice, t.Q, t.pe, t.rBar] = qLearningModel_5params_opponency(params, choice, outcome);
    case 'fiveParamO_peUpdate'
        [t.LH, t.probChoice, t.Q, t.pe, t.rBar] = qLearningModel_5params_opponency_peUpdate(params, choice, outcome);
    case 'fiveParam_bilal'
        [t.LH, t.probChoice, t.Q, t.pe] = qLearningModel_5params_bilal(params, choice, outcome);
     case 'fiveParam_peBeta_avg'
        [t.LH, t.probChoice, t.Q, t.pe, t.beta, t.R] = qLearningModel_5params_peBeta_avg(params, choice, outcome);
    case 'sixParam_peBeta'
        [t.LH, t.probChoice, t.Q, t.pe, t.beta, t.R] = qLearningModel_6params_peBeta(params, choice, outcome);
    case 'sixParam_t.peBeta_diff'
        [t.LH, t.probChoice, t.Q, t.pe, t.beta, t.R] = qLearningModel_6params_peBeta_diff(params, choice, outcome);
    case 'sevenParam_peBeta_k'
        [t.LH, t.probChoice, t.Q, t.pe, t.beta, t.R] = qLearningModel_7params_peBeta_k(params, choice, outcome);
    case 'fourParam_rBeta_confQ'
        [t.LH, t.probChoice, t.Q, t.pe, t.beta, t.R] = qLearningModel_4params_rBeta_confQ(params, choice, outcome);
    case 'fiveParam_rBeta_scale'
        [t.LH, t.probChoice, t.Q, t.pe, t.beta, t.R] = qLearningModel_5params_rBeta_scale(params, choice, outcome);
    case 'sixParam_rBeta_scale_min'
        [t.LH, t.probChoice, t.Q, t.pe, t.beta, t.R] = qLearningModel_6params_rBeta_scale_min(params, choice, outcome);
    case 'fiveParam_rBeta_confQ'
        [t.LH, t.probChoice, t.Q, t.pe, t.beta, t.R] = qLearningModel_5params_rBeta_scale(params, choice, outcome);
    case 'sixParam_rBeta_rRPE'
        [t.LH, t.probChoice, t.Q, t.pe, t.beta, t.R] = qLearningModel_6params_rBeta_rRPE(params, choice, outcome);
    case 'sixParam_rBeta_rV'
        [t.LH, t.probChoice, t.Q, t.pe, t.beta, t.R] = qLearningModel_6params_rBeta_rV(params, choice, outcome);
    case 'sixParam_peAN'
        [t.LH, t.probChoice, t.Q, t.pe, t.peBar, t.aN] = qLearningModel_6params_peAN(params, choice, outcome);
    case 'sixParam_pePeAN'
        [t.LH, t.probChoice, t.Q, t.pe, t.pePe, t.aN, t.peBar] = qLearningModel_6params_pePeAN(params, choice, outcome);
    case 'sixParam_pePeAN_bi'
        [t.LH, t.probChoice, t.Q, t.pe, t.pePe, t.aN, t.peBar] = qLearningModel_6params_pePeAN_bi(params, choice, outcome);
    case 'sixParam_pePeAN_lag'
        [t.LH, t.probChoice, t.Q, t.pe, t.pePe, t.aN, t.peBar] = qLearningModel_6params_pePeAN_lag(params, choice, outcome);
    case 'fiveParam_absPePeAN'
        [t.LH, t.probChoice, t.Q, t.pe, t.pePe, t.aN, t.peBar] = qLearningModel_5params_absPePeAN(params, choice, outcome);
    case 'sixParam_absPePeAN'
        [t.LH, t.probChoice, t.Q, t.pe, t.pePe, t.aN, t.peBar] = qLearningModel_6params_absPePeAN(params, choice, outcome);
    case 'sixParam_absPePeAN_bi'
        [t.LH, t.probChoice, t.Q, t.pe, t.pePe, t.aN, t.peBar] = qLearningModel_6params_absPePeAN_bi(params, choice, outcome);
    case 'sixParam_absPePeAN_bi_tW'
        [t.LH, t.probChoice, t.Q, t.pe, t.pePe, t.aN, t.peBar] = qLearningModel_6params_absPePeAN_bi(params, choice, outcome);
    case 'sevenParam_absPePeAN_bi_bias'
        [t.LH, t.probChoice, t.Q, t.pe, t.pePe, t.aN, t.peBar] = qLearningModel_7params_absPePeAN_bi_bias(params, choice, outcome);
    case 'eightParam_absPePeAN_bi_bias_k'
        [t.LH, t.probChoice, t.Q, t.pe, t.pePe, t.aN, t.peBar] = qLearningModel_8params_absPePeAN_bi_bias_k(params, choice, outcome);
    case 'eightParam_absPePeAN_bi_int_bias'
        [t.LH, t.probChoice, t.Q, t.pe, t.pePe, t.aN, t.peBar] = qLearningModel_8params_absPePeAN_bi_int_bias(params, choice, outcome);
    case 'sevenParam_absPePeAN_bi_bias_bothF'
        [t.LH, t.probChoice, t.Q, t.pe, t.pePe, t.aN, t.peBar] = qLearningModel_7params_absPePeAN_bi_bias_bothF(params, choice, outcome);
    case 'sevenParam_absPePeAN_bi_bias_tW'
        [t.LH, t.probChoice, t.Q, t.pe, t.pePe, t.aN, t.peBar] = qLearningModel_7params_absPePeAN_bi_bias(params, choice, outcome);
    case 'sixParam_absPePeAN_biSep'
        [t.LH, t.probChoice, t.Q, t.pe, t.pePe, t.aN, t.peBar, t.peBar_L, t.peBar_R] = qLearningModel_6params_absPePeAN_biSep(params, choice, outcome);
    case 'sevenParam_absPePeAN_biSep_f'
        [t.LH, t.probChoice, t.Q, t.pe, t.pePe, t.aN, t.peBar, t.peBar_L, t.peBar_R] = qLearningModel_7params_absPePeAN_biSep_f(params, choice, outcome);
    case 'sevenParam_absPePeAN_bi_k'
        [t.LH, t.probChoice, t.Q, t.pe, t.pePe, t.aN, t.peBar] = qLearningModel_7params_absPePeAN_bi_k(params, choice, outcome);
    case 'sixParam_absPePeAN_exp'
        [t.LH, t.probChoice, t.Q, t.pe, t.pePe, t.aN, t.peBar] = qLearningModel_6params_absPePeAN_exp(params, choice, outcome);
    case 'sixParam_absPePeAN_exp_tW'
        [t.LH, t.probChoice, t.Q, t.pe, t.pePe, t.aN, t.peBar] = qLearningModel_6params_absPePeAN_exp(params, choice, outcome);
    case 'sixParam_absPePeAN_exp_bias'
        [t.LH, t.probChoice, t.Q, t.pe, t.pePe, t.aN, t.peBar] = qLearningModel_6params_absPePeAN_exp_bias(params, choice, outcome);
    case 'sevenParam_absPePeAN_exp_bias'
        [t.LH, t.probChoice, t.Q, t.pe, t.pePe, t.aN, t.peBar] = qLearningModel_7params_absPePeAN_exp_bias(params, choice, outcome);
    case 'eightParam_absPePeAN_exp_bias_k'
        [t.LH, t.probChoice, t.Q, t.pe, t.pePe, t.aN, t.peBar] = qLearningModel_8params_absPePeAN_exp_bias_k(params, choice, outcome);
    case 'sevenParam_absPePeAN_exp_bias_tW'
        [t.LH, t.probChoice, t.Q, t.pe, t.pePe, t.aN, t.peBar] = qLearningModel_7params_absPePeAN_exp_bias(params, choice, outcome);
    case 'sixParam_absPePeAN_log'
        [t.LH, t.probChoice, t.Q, t.pe, t.pePe, t.aN, t.peBar] = qLearningModel_6params_absPePeAN_log(params, choice, outcome);
    case 'sevenParam_absPePeAN_log'
        [t.LH, t.probChoice, t.Q, t.pe, t.pePe, t.aN, t.peBar] = qLearningModel_7params_absPePeAN_log(params, choice, outcome); 
    case 'eightParam_absPePeAN_thresh_bias'
        [t.LH, t.probChoice, t.Q, t.pe, t.pePe, t.peBar] = qLearningModel_8params_absPePeAN_thresh_bias(params, choice, outcome);
    case 'sixParam_absPePe_scale_bias'
        [t.LH, t.probChoice, t.Q, t.pe, t.pePe, t.peBar] = qLearningModel_6params_absPePe_scale_bias(params, choice, outcome);
    case 'sixParam_absPePe_scaleBoth_bias'
        [t.LH, t.probChoice, t.Q, t.pe, t.pePe, t.peBar] = qLearningModel_6params_absPePe_scaleBoth_bias(params, choice, outcome); 
    case 'sixParam_absPePe_scaleSep_bias'
        [t.LH, t.probChoice, t.Q, t.pe, t.pePe, t.peBar] = qLearningModel_6params_absPePe_scaleSep_bias(params, choice, outcome);  
    case 'sevenParam_absPePe_scaleBoth_bias'
        [t.LH, t.probChoice, t.Q, t.pe, t.pePe, t.peBar] = qLearningModel_7params_absPePe_scaleBoth_bias(params, choice, outcome); 
    case 'sixParam_absPePeAN_scale_bias_tight'
        [t.LH, t.probChoice, t.Q, t.pe, t.pePe, t.peBar] = qLearningModel_6params_absPePe_scale_bias(params, choice, outcome);
    case 'sevenParam_absPePeAN_scale_bias'
        [t.LH, t.probChoice, t.Q, t.pe, t.pePe, t.aN, t.peBar] = qLearningModel_7params_absPePeAN_scale_bias(params, choice, outcome);
    case 'eightParam_absPePeAN_scale_int_bias'
        [t.LH, t.probChoice, t.Q, t.pe, t.pePe, t.aN, t.peBar] = qLearningModel_8params_absPePeAN_scale_int_bias(params, choice, outcome);
    case 'eightParam_absPePeAN_int_bias'
        [t.LH, t.probChoice, t.Q, t.pe, t.pePe, t.aN, t.peBar] = qLearningModel_8params_absPePeAN_int_bias(params, choice, outcome);
    case 'sevenParam_absPePeLR'
        [t.LH, t.probChoice, t.Q, t.pe, t.pePe, t.peBar] = qLearningModel_7params_absPePeLR(params, choice, outcome);
    case 'sevenParam_peLR'
        [t.LH, t.probChoice, t.Q, t.pe, t.R, t.aN, t.aP] = qLearningModel_7params_peLR(params, choice, outcome);
    case 'sixParam_rAN'
        [t.LH, t.probChoice, t.Q, t.pe, t.R, t.aN] = qLearningModel_6params_rAN(params, choice, outcome);
    case 'eightParam_rBeta_peAN'
        [t.LH, t.probChoice, t.Q, t.pe, t.R, t.aN, t.peBar] = qLearningModel_8params_rBeta_peAN(params, choice, outcome);
    case 'eightParam_rBeta_pePeAN'
        [t.LH, t.probChoice, t.Q, t.pe, t.R, t.aN, t.peBar] = qLearningModel_8params_rBeta_pePeAN(params, choice, outcome);

%% most recent
    case '5params'
        [t.LH, t.probChoice, t.Q, t.pe] = qLearningModel_5params(params, choice, outcome);
    case '5paramsLaserNegRPE'
        [t.LH, t.probChoice, t.Q, t.pe, t.peChange] = qLearningModel_5paramsLaserNegRPE(params, choice, outcome, laser);
    case '5paramsLaserNegRPERotation'
        [t.LH, t.probChoice, t.Q, t.pe, t.peChange] = qLearningModel_5paramsLaserNegRPERotation_simNoPlot(params, choice, outcome, laser);
    case '5params_k_bias'
        [t.LH, t.probChoice, t.Q, t.pe] = qLearningModel_5params_k(params, choice, outcome);
    case 'fiveParam_ph_bias'
        [t.LH, t.probChoice, t.Q, t.pe, t.alpha] = qLearningModel_5params_ph_bias(params, choice, outcome);
    case 'sixParam_ph_bias'
        [t.LH, t.probChoice, t.Q, t.pe, t.alpha] = qLearningModel_6params_ph_bias(params, choice, outcome);
    case 'sixParam_absPePeAN_bi_bias'
        [t.LH, t.probChoice, t.Q, t.pe, t.pePe, t.aN, t.peBar] = qLearningModel_6params_absPePeAN_bi_bias(params, choice, outcome);
    case 'sixParam_absPePeAN_bi_scale_bias'
        [t.LH, t.probChoice, t.Q, t.pe, t.pePe, t.aN, t.peBar] = qLearningModel_6params_absPePeAN_bi_scale_bias(params, choice, outcome);
    case 'sevenParam_absPePeAN_int_bias'
        [t.LH, t.probChoice, t.Q, t.pe, t.pePe, t.aN, t.peBar] = qLearningModel_7params_absPePeAN_int_bias(params, choice, outcome);  
    case 'sevenParam_absPePeAN_int_bias_ord'
        [t.LH, t.probChoice, t.Q, t.pe, t.pePe, t.aN, t.peBar] = qLearningModel_7params_absPePeAN_int_bias(params, choice, outcome);  
    case 'sevenParam_absPePeAN_scale_int_bias'
        [t.LH, t.probChoice, t.Q, t.pe, t.pePe, t.aN, t.peBar] = qLearningModel_7params_absPePeAN_scale_int_bias(params, choice, outcome);
    case '7params_absPePeAN_scale_int_bias_ord'
        [t.LH, t.probChoice, t.Q, t.pe, t.pePe, t.aN, t.peBar] = qLearningModel_7params_absPePeAN_scale_int_bias(params, choice, outcome);

    case 'sixParam_rTrace'
        [t.LH, t.probChoice, t.Q, t.pe, t.rBar] = qLearningModel_6params_rTrace(params, choice, outcome);
    case 'sevenParam_rTrace_k'
        [t.LH, t.probChoice, t.Q, t.pe, t.rBar] = qLearningModel_7params_rTrace_k(params, choice, outcome);
end 
    