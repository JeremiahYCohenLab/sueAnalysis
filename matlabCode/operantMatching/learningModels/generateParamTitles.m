function [titles] = generateParamTitles(paramNames)


for i = 1:length(paramNames)
    switch paramNames{i}
        case 'aN'
            titles{i} = '\alpha_N_P_E';
        case 'aNscale'
            titles{i} = '\alpha_N_P_E scale';
        case 'aNmin'
            titles{i} = '\alpha_N_P_E minimum';
        case 'aP'
            titles{i} = '\alpha_P_P_E';
        case 'aF'
            titles{i} = '\xi';
        case 'aPE'
            titles{i} = '\alpha_\sigma';
        case 'beta'
            titles{i} = '\beta';
        case 'betaScale'
            titles{i} = '\beta scale';
        case 'betaMin'
            titles{i} = '\beta minimum';
        case 'k'
            titles{i} = '\kappa';
        case 'peBar'
            titles{i} = '$\bar{\upsilon}$';
        otherwise
            titles{i} = paramNames{i};
    end
end


end