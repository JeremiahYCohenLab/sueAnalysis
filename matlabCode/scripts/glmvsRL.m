% parameters range
tarInd = 1;
stepNum = 5;
maxTrials = 800;
sessNum = 150;
binNum = 10;
paramNames = getParamNames_dF('5params', 0);
paramsFix = [0.6 0.7 0.4 4.5];
paramsRange = [[0.1 0.9]; [0.2, 0.95]; [0.4 0.8]; [3.0 7.0]];
params = zeros(stepNum, length(paramsFix));

for i = 1:length(paramsFix)
    if i ~= tarInd
        params(:,i) = paramsFix(i)*ones(stepNum,1);
    else
        params(:,i) = linspace(paramsRange(i,1), paramsRange(i,2),stepNum);
    end
end

%% simulation
clear allGLM
clear rwdGLM
for i = 1:stepNum
    combinedChoices = [];
    combinedRwdMat = [];
    combinedNrwdMat = [];
    combinedChoiceMat = [];
    % simulation
    for j = 1:sessNum
        [~, allRewards, allChoices] = qLearningModel_simNoPlot('maxTrials', maxTrials, 'params', params(i,:), 'randomSeed', i*10+j);
        rwdMatxTmp = [];
        choiceMatxTmp = [];
        noRwdMatxTmp = [];
        allChoices(allChoices==-1) = 0;
        allNoRewards = allChoices;
        allNoRewards(allRewards~=0) = 0;
        for j = 1:binNum
            rwdMatxTmp(j,:) = [NaN(1,j) allRewards(1:end-j)];
            choiceMatxTmp(j,:) = [NaN(1,j) allChoices(1:end-j)];
            noRwdMatxTmp(j,:) = [NaN(1,j) allNoRewards(1:end-j)];
        end
        combinedChoiceMat = [combinedChoiceMat, choiceMatxTmp];
        combinedRwdMat = [combinedRwdMat, rwdMatxTmp];
        combinedNrwdMat = [combinedNrwdMat, noRwdMatxTmp];
        combinedChoices = [combinedChoices, allChoices];
    end
    % model fitting
     glmTemp = fitglm([combinedRwdMat' combinedNrwdMat'], combinedChoices,'distribution','binomial','Link','logit', 'Intercept', true); 
     coeffTemp = glmTemp.Coefficients.Estimate;
     CIbands = coefCI(glmTemp);
     allGLM(i).coeff = coeffTemp;
     allGLM(i).CI = CIbands; 
     
     glmTemp = fitglm([combinedRwdMat'], combinedChoices,'distribution','binomial','Link','logit', 'Intercept', true); 
     coeffTemp = glmTemp.Coefficients.Estimate;
     CIbands = coefCI(glmTemp);
     rwdGLM(i).coeff = coeffTemp;
     rwdGLM(i).CI = CIbands; 
end
%% plotting
colorsRwd = [ones(stepNum,1), linspace(0.7, 0, stepNum)', linspace(0.7, 0, stepNum)'];
colorsNrwd = [linspace(0.7, 0, stepNum)', linspace(0.7, 0, stepNum)', ones(stepNum,1)];

% compare simulated glm plots
figure2;  
subplot(1,2,1); hold on; title([paramNames(tarInd) ' Rwd']); 
for i = 1:stepNum
    relevInds = 2:binNum+1;
    coefVals = allGLM(i).coeff;
    CIbands = allGLM(i).CI;
    errorL = abs(coefVals(relevInds) - CIbands(relevInds,1));
    errorU = abs(coefVals(relevInds) - CIbands(relevInds,2));
    errorbar((1:binNum),coefVals(relevInds),errorL,errorU,'Color', colorsRwd(i,:),'linewidth',2)   
end
line([0 binNum], [0 0], 'color', [0.7 0.7 0.7], 'LineStyle', '--')

subplot(1,2,2); hold on; title([paramNames(tarInd) ' NoRwd']); 
for i = 1:stepNum
    relevInds = binNum+2:2*binNum+1;
    coefVals = allGLM(i).coeff;
    CIbands = allGLM(i).CI;
    errorL = abs(coefVals(relevInds) - CIbands(relevInds,1));
    errorU = abs(coefVals(relevInds) - CIbands(relevInds,2));
    errorbar((1:binNum),coefVals(relevInds),errorL,errorU,'Color',colorsNrwd(i,:),'linewidth',2)   
end
line([0 binNum], [0 0], 'color', [0.7 0.7 0.7], 'LineStyle', '--')
%% plot rwd
% compare simulated glm plots
figure2;  hold on; title([paramNames(tarInd) ' Rwd']); 
for i = 1:stepNum
    relevInds = 2:binNum+1;
    coefVals = rwdGLM(i).coeff;
    CIbands = rwdGLM(i).CI;
    errorL = abs(coefVals(relevInds) - CIbands(relevInds,1));
    errorU = abs(coefVals(relevInds) - CIbands(relevInds,2));
    errorbar((1:binNum),coefVals(relevInds),errorL,errorU,'Color', colorsRwd(i,:),'linewidth',2)   
end
line([0 binNum], [0 0], 'color', [0.7 0.7 0.7], 'LineStyle', '--')
%%

