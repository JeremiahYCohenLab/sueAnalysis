 tMax =12;
figure;
    subplot(1,2,1); hold on;
    relevInds = 2:tMax+1;
    coefVals = CG14_preMdl.Coefficients.Estimate(relevInds);
    CIbands = coefCI(CG14_preMdl);
    errorL = abs(coefVals - CIbands(relevInds,1));
    errorU = abs(coefVals - CIbands(relevInds,2));
    errorbar((1:tMax)+0.2,coefVals,errorL,errorU,'Color',[0.7 0 1],'linewidth',2)
    coefVals = CG14_postMdl.Coefficients.Estimate(relevInds);
    CIbands = coefCI(CG14_postMdl);
    errorL = abs(coefVals - CIbands(relevInds,1));
    errorU = abs(coefVals - CIbands(relevInds,2));
    errorbar((1:tMax)+0.2,coefVals,errorL,errorU,'Color',[0.7 0.5 1],'linewidth',2)
    legend('pre', 'no rBar')
    xlabel('Reward n Trials Back')
    ylabel('\beta Coefficient')
    xlim([0.5 tMax+0.5])

    subplot(1,2,2); hold on;
    relevInds = tMax+2:length(CG14_preMdl.Coefficients.Estimate);
    coefVals = CG14_preMdl.Coefficients.Estimate(relevInds);
    CIbands = coefCI(CG14_preMdl);
    errorL = abs(coefVals - CIbands(relevInds,1));
    errorU = abs(coefVals - CIbands(relevInds,2));
    errorbar((1:tMax)+0.2,coefVals,errorL,errorU,'b','linewidth',2)
    coefVals = CG14_postMdl.Coefficients.Estimate(relevInds);
    CIbands = coefCI(CG14_postMdl);
    errorL = abs(coefVals - CIbands(relevInds,1));
    errorU = abs(coefVals - CIbands(relevInds,2));
    errorbar((1:tMax)+0.2,coefVals,errorL,errorU,'Color',[0.5 0.5 1],'linewidth',2)
    legend('pre', 'no rBar')
    xlabel('No Reward n Trials Back')
    ylabel('\beta Coefficient')
    xlim([0.5 tMax+0.5])