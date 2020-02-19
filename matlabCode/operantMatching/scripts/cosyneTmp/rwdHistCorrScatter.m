sigInds_maxFR = find(pVals_u.maxFR_rwdHist < 0.05);
nonsigInds_maxFR = find(pVals_u.maxFR_rwdHist > 0.05);
sigInds_preCS = find(pVals_u.preCS_rwdHist < 0.05);
nonsigInds_preCS = find(pVals_u.preCS_rwdHist > 0.05);
bothSigInds = intersect(sigInds_preCS, sigInds_maxFR);
sigInds_maxFR(ismember(sigInds_maxFR, bothSigInds)) = [];
sigInds_preCS(ismember(sigInds_preCS, bothSigInds)) = [];
bothNonSigInds = intersect(nonsigInds_preCS, nonsigInds_maxFR);

figure;
subplot(1,2,1);  hold on;
scatter(rhos_u.maxFR_rwdHist(bothSigInds), rhos_u.preCS_rwdHist(bothSigInds), 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'r')
scatter(rhos_u.maxFR_rwdHist(sigInds_maxFR), rhos_u.preCS_rwdHist(sigInds_maxFR), 'MarkerFaceColor', [0 1 1], 'MarkerEdgeColor', [0 1 1])
scatter(rhos_u.maxFR_rwdHist(sigInds_preCS), rhos_u.preCS_rwdHist(sigInds_preCS), 'MarkerFaceColor', 'g', 'MarkerEdgeColor', 'g')
scatter(rhos_u.maxFR_rwdHist(bothNonSigInds), rhos_u.preCS_rwdHist(bothNonSigInds), 'MarkerFaceColor', 'b', 'MarkerEdgeColor', 'b')
legend('both significant', 'max FR significant', 'pre CS significant', 'neither significant')
title('unidentified')
xlabel('post-CS max FR')
ylabel('pre-CS spike counts')
set(gca, 'TickDir', 'out')


sigInds_maxFR = find(pVals.maxFR_rwdHist < 0.05);
nonsigInds_maxFR = find(pVals.maxFR_rwdHist > 0.05);
sigInds_preCS = find(pVals.preCS_rwdHist < 0.05);
nonsigInds_preCS = find(pVals.preCS_rwdHist > 0.05);
bothSigInds = intersect(sigInds_preCS, sigInds_maxFR);
sigInds_maxFR(ismember(sigInds_maxFR, bothSigInds)) = [];
sigInds_preCS(ismember(sigInds_preCS, bothSigInds)) = [];
bothNonSigInds = intersect(nonsigInds_preCS, nonsigInds_maxFR);

subplot(1,2,2); hold on;
scatter(rhos.maxFR_rwdHist(bothSigInds), rhos.preCS_rwdHist(bothSigInds), 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'r')
scatter(rhos.maxFR_rwdHist(sigInds_maxFR), rhos.preCS_rwdHist(sigInds_maxFR), 'MarkerFaceColor', [0 1 1], 'MarkerEdgeColor', [0 1 1])
scatter(rhos.maxFR_rwdHist(sigInds_preCS), rhos.preCS_rwdHist(sigInds_preCS), 'MarkerFaceColor', 'g', 'MarkerEdgeColor', 'g')
scatter(rhos.maxFR_rwdHist(bothNonSigInds), rhos.preCS_rwdHist(bothNonSigInds), 'MarkerFaceColor', 'b', 'MarkerEdgeColor', 'b')
xlabel('post-CS max FR')
ylabel('pre-CS spike counts')
title('5ht')
suptitle('reward history spearman correlation coefficients')
set(gca, 'TickDir', 'out')

set(gcf, 'Renderer', 'Painters')

