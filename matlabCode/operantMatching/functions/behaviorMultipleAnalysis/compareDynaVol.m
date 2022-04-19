function compareDynaVol(xlFile, sheet, col, varargin)
    p = inputParser;
    % default parameters if none given
    p.addParameter('revForFlag', 0)
    p.addParameter('numBins', 10)
    p.addParameter('plotFlag', 1);
    p.addParameter('maxTrials', 600);
    p.addParameter('halves', 0);
    p.parse(varargin{:});
    
    switch p.Results.halves 
        case 0
            position = 'all';
        case 1
            position = 'early';
        case 2
            position = 'late';
    end


    
    dayList = getDayList(xlFile, sheet, col);
    
    combinedChoicesMatxH = []; 
    combinedRewardsMatxH = [];
    combinedNoRewardsMatxH = [];
    combinedAllChoice_RH = [];
    
    combinedChoicesMatxL = []; 
    combinedRewardsMatxL = [];
    combinedNoRewardsMatxL = [];
    combinedAllChoice_RL = [];
    
    tMax = p.Results.numBins;
 

for i = 1: length(dayList)
    sessionName = dayList{i};
    s = behAnalysisNoPlot_opMD(sessionName, 'simpleFlag', 1);
    allRewards = s.allRewards;
    allChoices = s.allChoices;
    allNoRewards = allChoices;
    allNoRewards(allRewards~=0) = 0;
    allChoice_R = allChoices;
    allChoice_R(allChoice_R<0) = 0;
 
    
    rwdMatxTmp = [];
    choiceMatxTmp = [];
    noRwdMatxTmp = [];
    for j = 1:tMax
        rwdMatxTmp(j,:) = [NaN(1,j) allRewards(1:end-j)];
        choiceMatxTmp(j,:) = [NaN(1,j) allChoices(1:end-j)];
        noRwdMatxTmp(j,:) = [NaN(1,j) allNoRewards(1:end-j)];
    end
    
%    allRewards(allRewards == -1) = 1;
    if p.Results.halves==0
        combinedRewardsMatxH = [combinedRewardsMatxH NaN(tMax,100) rwdMatxTmp(:,s.vol==1)];
        combinedNoRewardsMatxH = [combinedNoRewardsMatxH NaN(tMax,100) noRwdMatxTmp(:,s.vol==1)];
        combinedChoicesMatxH = [combinedChoicesMatxH NaN(tMax,100) choiceMatxTmp(:,s.vol==1)];
        combinedAllChoice_RH = [combinedAllChoice_RH NaN(1,100) allChoice_R(1,s.vol==1)];


        combinedRewardsMatxL = [combinedRewardsMatxL NaN(tMax,100) rwdMatxTmp(:,s.vol==0)]; 
        combinedNoRewardsMatxL = [combinedNoRewardsMatxL NaN(tMax,100) noRwdMatxTmp(:,s.vol==0)];
        combinedChoicesMatxL  = [combinedChoicesMatxL NaN(tMax,100) choiceMatxTmp(:,s.vol==0)];
        combinedAllChoice_RL = [combinedAllChoice_RL NaN(1,100) allChoice_R(1,s.vol==0)];
    else
        if p.Results.halves==1
           if s.vol(1) == 1
                combinedRewardsMatxH = [combinedRewardsMatxH NaN(tMax,100) rwdMatxTmp(:,s.vol==1)];
                combinedNoRewardsMatxH = [combinedNoRewardsMatxH NaN(tMax,100) noRwdMatxTmp(:,s.vol==1)];
                combinedChoicesMatxH = [combinedChoicesMatxH NaN(tMax,100) choiceMatxTmp(:,s.vol==1)];
                combinedAllChoice_RH = [combinedAllChoice_RH NaN(1,100) allChoice_R(1,s.vol==1)];               
           else
                combinedRewardsMatxL = [combinedRewardsMatxL NaN(tMax,100) rwdMatxTmp(:,s.vol==0)];
                combinedNoRewardsMatxL = [combinedNoRewardsMatxL NaN(tMax,100) noRwdMatxTmp(:,s.vol==0)];
                combinedChoicesMatxL = [combinedChoicesMatxL NaN(tMax,100) choiceMatxTmp(:,s.vol==0)];
                combinedAllChoice_RL = [combinedAllChoice_RL NaN(1,100) allChoice_R(1,s.vol==0)];                              
           end
        else
            if s.vol(1) == 1
                combinedRewardsMatxL = [combinedRewardsMatxL NaN(tMax,100) rwdMatxTmp(:,s.vol==0)];
                combinedNoRewardsMatxL = [combinedNoRewardsMatxL NaN(tMax,100) noRwdMatxTmp(:,s.vol==0)];
                combinedChoicesMatxL = [combinedChoicesMatxL NaN(tMax,100) choiceMatxTmp(:,s.vol==0)];
                combinedAllChoice_RL = [combinedAllChoice_RL NaN(1,100) allChoice_R(1,s.vol==0)];   
            else
                combinedRewardsMatxH = [combinedRewardsMatxH NaN(tMax,100) rwdMatxTmp(:,s.vol==1)];
                combinedNoRewardsMatxH = [combinedNoRewardsMatxH NaN(tMax,100) noRwdMatxTmp(:,s.vol==1)];
                combinedChoicesMatxH = [combinedChoicesMatxH NaN(tMax,100) choiceMatxTmp(:,s.vol==1)];
                combinedAllChoice_RH = [combinedAllChoice_RH NaN(1,100) allChoice_R(1,s.vol==1)]; 
            end
        end
    end

end
    glm_rwdNoRwdH = fitglm([combinedRewardsMatxH' combinedNoRewardsMatxH'], combinedAllChoice_RH,'distribution','binomial','link','logit'); 
    glm_rwdNoRwdL = fitglm([combinedRewardsMatxL' combinedNoRewardsMatxL'], combinedAllChoice_RL,'distribution','binomial','link','logit');
    
    glm_rwdChoiceH = fitglm([combinedRewardsMatxH' combinedChoicesMatxH'], combinedAllChoice_RH,'distribution','binomial','link','logit'); 
    glm_rwdChoiceL = fitglm([combinedRewardsMatxL' combinedChoicesMatxL'], combinedAllChoice_RL,'distribution','binomial','link','logit');    
%% rwd and norwd    
    figure;
    subplot(1,2,1); hold on;
    relevInds = 2:tMax+1;
    coefVals = glm_rwdNoRwdH.Coefficients.Estimate(relevInds);
    CIbands = coefCI(glm_rwdNoRwdH);
    errorL = abs(coefVals - CIbands(relevInds,1));
    errorU = abs(coefVals - CIbands(relevInds,2));
    errorbar((1:tMax),coefVals,errorL,errorU,'Color', [0.7 0 1],'linewidth',2)
    
    coefVals = glm_rwdNoRwdL.Coefficients.Estimate(relevInds);
    CIbands = coefCI(glm_rwdNoRwdL);
    errorL = abs(coefVals - CIbands(relevInds,1));
    errorU = abs(coefVals - CIbands(relevInds,2));
    errorbar((1:tMax),coefVals,errorL,errorU,'Color', 'b','linewidth',2)    
    line([0.5 tMax+0.5], [0 0], 'Color',[0.5 0.5 0.5],'LineStyle','--')
    
    legend({'high', 'low'});
    xlabel('Reward n Trials Back')
    ylabel('\beta Coefficient')
    xlim([0.5 tMax+0.5])
    
    title('rwd')

    subplot(1,2,2); hold on;
    
    relevInds = tMax+2:length(glm_rwdNoRwdH.Coefficients.Estimate);
    coefVals = glm_rwdNoRwdH.Coefficients.Estimate(relevInds);
    CIbands = coefCI(glm_rwdNoRwdH);
    errorL = abs(coefVals - CIbands(relevInds,1));
    errorU = abs(coefVals - CIbands(relevInds,2));
    errorbar((1:tMax),coefVals,errorL,errorU,'Color', [0.7 0 1],'linewidth',2)

    relevInds = tMax+2:length(glm_rwdNoRwdL.Coefficients.Estimate);
    coefVals = glm_rwdNoRwdL.Coefficients.Estimate(relevInds);
    CIbands = coefCI(glm_rwdNoRwdL);
    errorL = abs(coefVals - CIbands(relevInds,1));
    errorU = abs(coefVals - CIbands(relevInds,2));
    errorbar((1:tMax),coefVals,errorL,errorU,'Color', 'b','linewidth',2)
    
    line([0.5 tMax+0.5], [0 0], 'Color',[0.5 0.5 0.5],'LineStyle','--')

    xlabel('Reward n Trials Back')
    xlim([0.5 tMax+0.5])
    legend('high', 'low')
    title('no rwd')

    suptitle([sheet ' ' col ' ' position])
    
    
    %% choice and rwd
    figure;
    subplot(1,2,1); hold on;
    relevInds = 2:tMax+1;
    coefVals = glm_rwdChoiceH.Coefficients.Estimate(relevInds);
    CIbands = coefCI(glm_rwdChoiceH);
    errorL = abs(coefVals - CIbands(relevInds,1));
    errorU = abs(coefVals - CIbands(relevInds,2));
    errorbar((1:tMax),coefVals,errorL,errorU,'Color', [0.7 0 1],'linewidth',2)
    
    coefVals = glm_rwdChoiceL.Coefficients.Estimate(relevInds);
    CIbands = coefCI(glm_rwdChoiceL);
    errorL = abs(coefVals - CIbands(relevInds,1));
    errorU = abs(coefVals - CIbands(relevInds,2));
    errorbar((1:tMax),coefVals,errorL,errorU,'Color', 'b','linewidth',2)    
    line([0.5 tMax+0.5], [0 0], 'Color',[0.5 0.5 0.5],'LineStyle','--')
    
    legend({'high', 'low'});
    xlabel('Reward n Trials Back')
    ylabel('\beta Coefficient')
    xlim([0.5 tMax+0.5])
    
    title('rwd')

    subplot(1,2,2); hold on;
    
    relevInds = tMax+2:length(glm_rwdChoiceH.Coefficients.Estimate);
    coefVals = glm_rwdChoiceH.Coefficients.Estimate(relevInds);
    CIbands = coefCI(glm_rwdChoiceH);
    errorL = abs(coefVals - CIbands(relevInds,1));
    errorU = abs(coefVals - CIbands(relevInds,2));
    errorbar((1:tMax),coefVals,errorL,errorU,'Color', [0.7 0 1],'linewidth',2)

    relevInds = tMax+2:length(glm_rwdChoiceL.Coefficients.Estimate);
    coefVals = glm_rwdChoiceL.Coefficients.Estimate(relevInds);
    CIbands = coefCI(glm_rwdChoiceL);
    errorL = abs(coefVals - CIbands(relevInds,1));
    errorU = abs(coefVals - CIbands(relevInds,2));
    errorbar((1:tMax),coefVals,errorL,errorU,'Color', 'b','linewidth',2)
    
    line([0.5 tMax+0.5], [0 0], 'Color',[0.5 0.5 0.5],'LineStyle','--')

    xlabel('Choices n Trials Back')
    xlim([0.5 tMax+0.5])
    legend('high', 'low')
    title('choice')

    suptitle([sheet ' ' col ' ' position])

    

    
    
