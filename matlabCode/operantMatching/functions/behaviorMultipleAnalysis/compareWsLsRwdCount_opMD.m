function compareWsLsRwdCount_opMD(xlFile, sheet, pre, post, tMax)

if nargin < 5
    tMax = 10;
end

[wS_rwdHx_pre, lS_rwdHx_pre, ~, ~] = wslsRwdCount_opMD(xlFile, sheet, pre, 'tMax', tMax);
[wS_rwdHx_post, lS_rwdHx_post, ~, ~] = wslsRwdCount_opMD(xlFile, sheet, post, 'tMax', tMax);


figure;
subplot(1,2,1); hold on;
scatter([0:tMax], wS_rwdHx_pre, [], 'b', 'filled')
scatter([0:tMax], wS_rwdHx_post, [], 'r', 'filled')
xlim([-0.5 tMax+0.5])
xlabel('number of rewards in last 10 trials')
ylabel('probability')
title('win-stay')
legend('pre', 'post')
set(gca, 'tickdir', 'out')

subplot(1,2,2); hold on;
scatter([0:tMax], lS_rwdHx_pre, [], 'b', 'filled')
scatter([0:tMax], lS_rwdHx_post, [], 'r', 'filled')
xlim([-0.5 tMax+0.5])
xlabel('number of rewards in last 10 trials')
ylabel('probability')
title('lose-shift')
legend('pre', 'post')
set(gca, 'tickdir', 'out')

set(gcf, 'renderer', 'painters', 'position', [-1919 41 1920 963])