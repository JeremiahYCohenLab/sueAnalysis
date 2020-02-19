figure;
subplot(2,3,1); hold on;
histogram(CG14preS_fiveParamO_fourStart.mu_aN, 50, 'Normalization', 'Probability', 'FaceColor', 'r')
histogram(CG14effect_fiveParamOFixed_aN.mu_aN, 50, 'Normalization', 'Probability', 'FaceColor', [1 0.7 0.7])
ax = gca;
ax.TickDir = 'out';
legend('pre-lesion', 'post-lesion')
title('\alpha NPE')
set(gca, 'FontSize',18)

subplot(2,3,2); hold on;
histogram(CG14preS_fiveParamO_fourStart.mu_aP, 50,  'Normalization', 'Probability', 'FaceColor', 'b')
histogram(CG14effect_fiveParamOFixed_aP.mu_aP, 50,  'Normalization', 'Probability', 'FaceColor', [0.7 0.7 1])
ax = gca;
ax.TickDir = 'out';
%legend('pre-lesion', 'post-lesion')
title('\alpha PPE')
set(gca, 'FontSize',18)

subplot(2,3,3); hold on;
histogram(CG14preS_fiveParamO_fourStart.mu_aF, 50, 'Normalization', 'Probability', 'FaceColor', 'c')
histogram(CG14effect_fiveParamOFixed_aF.mu_aF, 50, 'Normalization', 'Probability', 'FaceColor', [0.7 1 1])
ax = gca;
ax.TickDir = 'out';
%legend('pre-lesion', 'post-lesion')
title('\alpha forget')
set(gca, 'FontSize',18)

subplot(2,3,4); hold on;
histogram(CG14preS_fiveParamO_fourStart.mu_beta, 50, 'Normalization', 'Probability', 'FaceColor', 'm')
histogram(CG14effect_fiveParamOFixed_beta.mu_beta, 50, 'Normalization', 'Probability', 'FaceColor', [1 0.7 1])
ax = gca;
ax.TickDir = 'out';
%legend('pre-lesion', 'post-lesion')
title('\beta')
set(gca, 'FontSize',18)

subplot(2,3,5); hold on;
histogram(CG14preS_fiveParamO_fourStart.mu_v, 50, 'Normalization', 'Probability', 'FaceColor', [0.6 0 1])
histogram(CG14effect_fiveParamOFixed_v.mu_v, 200, 'Normalization', 'Probability', 'FaceColor', [0.9 0.7 1])
ax = gca;
ax.TickDir = 'out';
%legend('pre-lesion', 'post-lesion')
title('v')
set(gca, 'FontSize',18)


ax = gca;
ax.TickDir = 'out';
%set(gcf, 'Position', [       -1911         490        1904         419]);
set(gcf, 'Renderer', 'Painters')