function paramNames = getParamNames_dF(modelName, biasFlag)


switch modelName
    case 'fourParam'
        paramNames = [{'aN', 'aP', 'aF', 'beta'}];
    case '4params_bias_LaserNegRPE'
        paramNames = [{'a', 'aF', 'beta', 'diff'}];
    case '5params'
        paramNames = [{'aN', 'aP', 'aF', 'beta'}];     
    case '5params_biForget'
        paramNames = [{'aN', 'aP', 'aF', 'beta'}]; 
    case '5paramsHmm'
        paramNames = [{'aN', 'aP', 'aF', 'betaE', 'betaT'}];
    case '5paramsLaserNegRPE'
        paramNames = [{'aN', 'aP', 'aF', 'beta', 'diff'}]; 
    case '5paramsLaserNegRPERotation'
        paramNames = [{'aN', 'aP', 'aF', 'beta', 'diff'}]; 
    case '5paramsLaserUpdate'
        paramNames = [{'aN', 'aP', 'aF', 'beta', 'diff'}]; 
    case '5paramsExpForgetHmm'
        paramNames = [{'aN', 'aP', 'aF', 'betaE', 'betaT'}];
    case '5paramsNoPriorForgetHmm'
        paramNames = [{'aN', 'aP', 'aF', 'betaE', 'betaT'}];
    case '5paramsNoPriorForget'
        paramNames = [{'aN', 'aP', 'aF', 'beta'}];
    case '5params_tF'
        paramNames = [{'aN', 'aP', 'tF', 'beta'}];
    case '5params_biaF'
        paramNames = [{'aN', 'aP', 'aF', 'beta'}];
    case '5params_expF'
        paramNames = [{'aN', 'aP', 'aF', 'beta'}];
    case '5params_k_bias'
        paramNames = [{'a', 'aF', 'beta', 'k'}];
    case '5params_k_bias_biForget'
        paramNames = [{'a', 'aF', 'beta', 'k'}];
    case '5params_k_bias_LaserNegRPE'
        paramNames = [{'a', 'aF', 'beta', 'k', 'diff'}];
    case '5params_k_bias_LaserNegRPE_expPrior'
        paramNames = [{'a', 'aF', 'beta', 'k', 'diff'}];
    case '5params_k_bias_LaserDisengage'
        paramNames = [{'a', 'aF', 'beta', 'k', 'diff'}];
    case '5params_k_bias_LaserDisengageScale'
        paramNames = [{'a', 'aF', 'beta', 'k', 'scale'}];
    case '5params_k_bias_LaserDisengageScale_expPrior'
        paramNames = [{'a', 'aF', 'beta', 'k', 'scale'}];        
    case '5params_k_bias_LaserNegOnlyRPE'
        paramNames = [{'a', 'aF', 'beta', 'k', 'diff'}];
    case '5params_k_bias_LaserNegRPERotation'
        paramNames = [{'a', 'aF', 'beta', 'k', 'diff'}];
    case '5params_2LR_k_bias'
        paramNames = [{'aN', 'aP', 'beta', 'k'}];
    case '5params_k_biForget_bias'
        paramNames = [{'a', 'aF', 'beta', 'k'}];
    case '5params_kExp_bias'
        paramNames = [{'a', 'aF', 'aChoice', 'beta', 'k'}];
    case '5params_kExp_bias_LaserNegRPERotation'
        paramNames = [{'a', 'aF', 'aChoice', 'beta', 'k', 'diff'}];
    case '5params_inv'
        paramNames = [{'aN', 'aP', 'aF', 'beta'}];
    case 'fiveParam_kappa'
        paramNames = [{'aN', 'aP', 'aF', 'beta', 'k'}];
    case 'fiveParam_ph_bias'
        paramNames = [{'eta', 'kappa', 'aF', 'beta'}];
    case 'sixParam_ph_bias'
        paramNames = [{'eta', 'kappaN', 'kappaP', 'aF', 'beta'}];
    case 'sixParam_absPePe_scale_bias'
        paramNames = [{'aN', 'aP', 'aF', 'aPE', 'beta'}];
    case 'sixParam_absPePe_scaleBoth_bias'
        paramNames = [{'aN', 'aP', 'aF', 'aPE', 'beta'}];
    case 'sixParam_absPePeAN_bi_bias'
        paramNames = [{'aNmin', 'aP', 'aF', 'aPE', 'beta'}];
    case 'sixParam_absPePeAN_bi_scale_bias'
        paramNames = [{'aNmin', 'aP', 'aF', 'aPE', 'beta'}];
    case 'sixParam_absPePeAN_bi_bias_noF'
        paramNames = [{'aNmin', 'aNscale', 'aP', 'aPE', 'beta'}];
    case 'sevenParam_absPePeAN_bi_bias'
        paramNames = [{'aNmin', 'aNscale', 'aP', 'aF', 'aPE', 'beta'}];
    case 'sevenParam_absPePeAN_int_bias'
        paramNames = [{'aNmin', 'aP', 'aF', 'aPE', 'v', 'beta'}];
    case 'sevenParam_absPePeAN_scale_int_bias'
        paramNames = [{'aNmin', 'aP', 'aF', 'aPE', 'v', 'beta'}];
    case '7params_absPePeAN_scale_int_bias_ord'
        paramNames = [{'aNmin', 'aP', 'aF', 'aPE', 'v', 'beta'}];
    case 'sevenParam_absPePeAN_int_bias_ord'
        paramNames = [{'aNmin', 'aP', 'aF', 'aPE', 'v', 'beta'}];
    case 'eightParam_absPePeAN_scale_int_bias'
        paramNames = [{'aNmin', 'aNscale', 'aP', 'aF', 'aPE', 'v', 'beta'}];
    case 'sixParam_rTrace'
        paramNames = [{'a', 'aF', 'beta', 'v', 'w'}];
    case 'sevenParam_rTrace_k'
        paramNames = [{'a', 'aF', 'beta', 'v', 'w', 'k'}];
    case 'sixParam_ph_bias'
        paramNames = [{'eta', 'kappaN', 'kappaP', 'aF', 'beta'}];
        
    % delta models
    case 'delta_sevenParam_absPePeAN_scale_int_bias_ord'
        paramNames = [{'aNmin', 'aP', 'aF', 'aPE', 'v', 'beta', 'dPeBar'}];

    % bayesian models
    case 'fbm_tsPrior'
        paramNames = [{'a', 'b'}];
    case 'fbm_softmax'
        paramNames = [{'a', 'b', 'beta'}];
    case 'dbm_tsPrior'
        paramNames = [{'gamma', 'a', 'b'}];
    case 'dbm_softmax'
        paramNames = [{'gamma', 'a', 'b', 'beta'}];
    case 'dbm_softmax_bias'
        paramNames = [{'gamma', 'a', 'b', 'beta', 'bias'}];
    case 'dbm_softmax_bias_probs'
        paramNames = [{'gamma', 'beta', 'bias'}];
    case 'dbm_softmax_ab'
        paramNames = [{'gamma', 'ab', 'beta'}];
    % vkf models
    case 'vkf'
        paramNames = [{'lambda', 'vInit', 'omega', 'beta'}];
    case 'vkf_fixV'
        paramNames = [{'lambda', 'vInit', 'omega', 'beta'}];        
    case 'vkf_fixV_aF'
        paramNames = [{'lambda', 'vInit', 'omega', 'beta', 'aF'}];
    case 'vkf_fixV_kappa'
        paramNames = [{'lambda', 'vInit', 'omega', 'beta', 'kappa'}];
end

if biasFlag
    paramNames = [paramNames {'bias'}];
end