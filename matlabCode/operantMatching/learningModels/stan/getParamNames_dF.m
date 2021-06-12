function paramNames = getParamNames_dF(modelName, biasFlag)


switch modelName
    case 'fourParam'
        paramNames = [{'aN', 'aP', 'aF', 'beta'}];
    case '5params'
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
end

if biasFlag
    paramNames = [paramNames {'bias'}];
end