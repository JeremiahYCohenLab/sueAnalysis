function [preGLM, postGLM] = compareLinRegLick_opMD(file, animal, pre, post, revForFlag, trialFlag)

if nargin < 6
    trialFlag = 0;
end
if nargin < 5
    revForFlag = 0;
end

%run function to generate lrm 
if trialFlag
    [preGLM, preStayLickLat, preSwitchLickLat, tMax, preITIlicks]= combineLinRegLickLat_opMD(file, animal, pre, 'revForFlag', revForFlag);
    [postGLM, postStayLickLat, postSwitchLickLat, ~, postITIlicks]= combineLinRegLickLat_opMD(file, animal, post, 'revForFlag', revForFlag);
else
    [preGLM, preStayLickLat, preSwitchLickLat, binSize, timeMax, preITIlicks]= combineLinRegLickLatTime_opMD(file, animal, pre, 'revForFlag', revForFlag);
    [postGLM, postStayLickLat, postSwitchLickLat, ~, ~, postITIlicks]= combineLinRegLickLatTime_opMD(file, animal, post, 'revForFlag', revForFlag);
    timeBinEdges = [1000:binSize:timeMax];
    tMax = length(timeBinEdges) - 1;
end


%% plot beta coeffs for multiple covariate type model
figure;
subplot(1,2,1); hold on;
relevInds = 2:tMax+1;
coefVals = preGLM.Coefficients.Estimate(relevInds);
CIbands = coefCI(preGLM);
errorL = abs(coefVals - CIbands(relevInds,1));
errorU = abs(coefVals - CIbands(relevInds,2));
if trialFlag
    errorbar([1:tMax],coefVals,errorL,errorU,'b','linewidth',2)
else
    errorbar(((1:tMax)*binSize/1000),coefVals,errorL,errorU,'b','linewidth',2)
end

coefVals = postGLM.Coefficients.Estimate(relevInds);
CIbands = coefCI(postGLM);
errorL = abs(coefVals - CIbands(relevInds,1));
errorU = abs(coefVals - CIbands(relevInds,2));
if trialFlag
    errorbar([1:tMax],coefVals,errorL,errorU,'Color', [0.7 0 1],'linewidth',2)
    xlabel('Reward n trials back')
else
    errorbar(((1:tMax)*binSize/1000),coefVals,errorL,errorU,'Color', [0.7 0 1],'linewidth',2)
    xlim([0 (tMax*binSize/1000 + binSize/1000)])
    xlabel('Reward n seconds back')
end
y = ylim;
title('LRM - Rewards on Z-scored Licks')
legend('pre', 'post')
ylabel('\beta Coefficient')
set(gca, 'tickdir', 'out')


subplot(1,2,2); hold on;
relevInds = tMax+2:2*tMax+1;
coefVals = preGLM.Coefficients.Estimate(relevInds);
CIbands = coefCI(preGLM);
errorL = abs(coefVals - CIbands(relevInds,1));
errorU = abs(coefVals - CIbands(relevInds,2));
if trialFlag
    errorbar([1:tMax],coefVals,errorL,errorU,'b','linewidth',2)
else
    errorbar(((1:tMax)*binSize/1000),coefVals,errorL,errorU,'b','linewidth',2)
end

coefVals = postGLM.Coefficients.Estimate(relevInds);
CIbands = coefCI(postGLM);
errorL = abs(coefVals - CIbands(relevInds,1));
errorU = abs(coefVals - CIbands(relevInds,2));
if trialFlag
    errorbar([1:tMax],coefVals,errorL,errorU,'Color', [0.7 0 1],'linewidth',2)
    xlabel('No reward n trials back')
else
    errorbar(((1:tMax)*binSize/1000),coefVals,errorL,errorU,'Color', [0.7 0 1],'linewidth',2)
    xlim([0 (tMax*binSize/1000 + binSize/1000)])
    xlabel('No reward n seconds back')
end
ylim([y])
title('LRM - No rewards on Z-scored Licks')
ylabel('\beta Coefficient')
set(gca, 'tickdir', 'out')

suptitle(animal)

%% lick latency by stay and switch figure
figure;
mag = [1 0 1];
cyan = [0 1 1];
set(gcf,'defaultAxesColorOrder',[mag; cyan]);

subplot(2,2,1)
yyaxis left; histogram(preStayLickLat, 30, 'Normalization', 'probability')
ylabel('Probability')
yyaxis right; histogram(postStayLickLat, 30, 'Normalization', 'probability')
legend('pre', 'post')
title('stay lick latency')
set(gca, 'tickdir', 'out')


subplot(2,2,2)
yyaxis left; histogram(preSwitchLickLat, 30, 'Normalization', 'probability')
yyaxis right; histogram(postSwitchLickLat, 30, 'Normalization', 'probability')
legend('pre', 'post')
title('switch lick latency')
set(gca, 'tickdir', 'out')


subplot(2,2,3)
yyaxis left; histogram(preStayLickLat, 30, 'Normalization', 'probability')
ylabel('Probability')
yyaxis right; histogram(preSwitchLickLat, 30, 'Normalization', 'probability')
legend('stay', 'switch')
title('pre lick latency')
set(gca, 'tickdir', 'out')


subplot(2,2,4)
yyaxis left; histogram(postStayLickLat, 30, 'Normalization', 'probability')
yyaxis right; histogram(postSwitchLickLat, 30, 'Normalization', 'probability')
legend('stay', 'switch')
title('post lick latency')
set(gca, 'tickdir', 'out')

suptitle(animal)

figure; hold on;
histogram(preITIlicks, 100, 'Normalization', 'probability', 'FaceColor', 'c')
ylabel('Probability')
histogram(postITIlicks, 100, 'Normalization', 'probability' , 'FaceColor', 'm')
legend('pre', 'post')
title(['ITI lick rate - ' animal])
set(gca, 'tickdir', 'out')

