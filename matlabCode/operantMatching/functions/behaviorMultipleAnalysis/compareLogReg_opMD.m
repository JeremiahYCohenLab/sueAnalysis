function [preAll, postAll] = compareLogReg_opMD(file, animal, pre, post, varargin)

%task and model parameters
p = inputParser;
% default parameters if none given
p.addParameter('trialFlag', 0)
p.addParameter('revForFlag', 0)
p.addParameter('binSize', 6000)
p.addParameter('numBins', 10)
p.addParameter('plotFlag', 1);
p.parse(varargin{:});


%run function to generate lrm
if p.Results.trialFlag
    [preAll, tMax]= combineLogReg_opMD(file, animal, pre, 'revForFlag', p.Results.revForFlag,...
       'numBins', p.Results.numBins);
    [postAll,~]= combineLogReg_opMD(file, animal, post, 'revForFlag', p.Results.revForFlag,...
       'numBins', p.Results.numBins);
else
    [preAll, s]= combineLogRegTime_opMD(file, animal, pre, 'revForFlag', p.Results.revForFlag,...
       'binSize', p.Results.binSize, 'numBins', p.Results.numBins);
    [postAll,~]= combineLogRegTime_opMD(file, animal, post, 'revForFlag',  p.Results.revForFlag,...
       'binSize', p.Results.binSize, 'numBins', p.Results.numBins);
    tMax = s.tMax;
end


%plot beta coeffs for multiple covariate type model
figure;
subplot(1,3,1); hold on;
relevInds = 2:tMax+1;
coefVals = preAll.Coefficients.Estimate(relevInds);
CIbands = coefCI(preAll);
errorL = abs(coefVals - CIbands(relevInds,1));
errorU = abs(coefVals - CIbands(relevInds,2));
if p.Results.trialFlag
    errorbar([1:tMax],coefVals,errorL,errorU,'b','linewidth',2)
else
    errorbar(((1:s.tMax)*s.binSize/1000),coefVals,errorL,errorU,'b','linewidth',2)
end

coefVals = postAll.Coefficients.Estimate(relevInds);
CIbands = coefCI(postAll);
errorL = abs(coefVals - CIbands(relevInds,1));
errorU = abs(coefVals - CIbands(relevInds,2));
if p.Results.trialFlag
    errorbar([1:tMax],coefVals,errorL,errorU,'Color', [0.7 0 1],'linewidth',2)
    xlabel('Reward n trials back')
else
    errorbar(((1:s.tMax)*s.binSize/1000),coefVals,errorL,errorU,'Color', [0.7 0 1],'linewidth',2)
    xlim([0 (s.tMax*s.binSize/1000 + s.binSize/1000)])
    xlabel('Reward n seconds back')
end
set(gca, 'TickDir', 'Out')
title('Combined Model - Reward')
legend(['pre | intercept: ' num2str(preAll.Coefficients.Estimate(1))], ['post | intercept: ' num2str(postAll.Coefficients.Estimate(1))])
ylabel('\beta Coefficient')


subplot(1,3,2); hold on;
relevInds = tMax+2:tMax*2+1;
coefVals = preAll.Coefficients.Estimate(relevInds);
CIbands = coefCI(preAll);
errorL = abs(coefVals - CIbands(relevInds,1));
errorU = abs(coefVals - CIbands(relevInds,2));
if p.Results.trialFlag
    errorbar([1:tMax],coefVals,errorL,errorU,'b','linewidth',2)
else
    errorbar(((1:s.tMax)*s.binSize/1000),coefVals,errorL,errorU,'b','linewidth',2)
end

coefVals = postAll.Coefficients.Estimate(relevInds);
CIbands = coefCI(postAll);
errorL = abs(coefVals - CIbands(relevInds,1));
errorU = abs(coefVals - CIbands(relevInds,2));
if p.Results.trialFlag
    errorbar([1:tMax],coefVals,errorL,errorU,'Color', [0.7 0 1],'linewidth',2)
    xlabel('No reward n trials back')
else
    errorbar(((1:s.tMax)*s.binSize/1000),coefVals,errorL,errorU,'Color', [0.7 0 1],'linewidth',2)
    xlim([0 (s.tMax*s.binSize/1000 + s.binSize/1000)])
    xlabel('No reward n seconds back')
end
set(gca, 'TickDir', 'Out')
title('Combined Model - No Reward')
legend('pre', 'post')
ylabel('\beta Coefficient')

subplot(1,3,3); hold on;
relevInds = tMax*2+1:length(preAll.Coefficients.Estimate);
coefVals = preAll.Coefficients.Estimate(relevInds);
CIbands = coefCI(preAll);
errorL = abs(coefVals - CIbands(relevInds,1));
errorU = abs(coefVals - CIbands(relevInds,2));
errorbar([1:length(relevInds)],coefVals,errorL,errorU,'b','linewidth',2)

relevInds = tMax*2+1:length(postAll.Coefficients.Estimate);
coefVals = postAll.Coefficients.Estimate(relevInds);
CIbands = coefCI(postAll);
errorL = abs(coefVals - CIbands(relevInds,1));
errorU = abs(coefVals - CIbands(relevInds,2));
errorbar([1:length(relevInds)],coefVals,errorL,errorU,'Color', [0.7 0 1],'linewidth',2)
xlabel('Session number')

set(gca, 'TickDir', 'Out')
title('Combined Model - Session bias')
legend('pre', 'post')
ylabel('\beta Coefficient')
suptitle(animal)