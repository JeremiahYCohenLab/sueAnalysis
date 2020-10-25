%% summary
% lickLat analysis of drug inhibition
%% load data from drug and saline 
    xlFile = 'combineDrug';
    animal = 'combineDrugSaline';
    cmbSaline = combineSesh(xlFile, 'combineDrugSaline', 'saline', 'maxTrials',200);
    cmbDrug = combineSesh(xlFile, 'combineDrugDrug', 'drug', 'maxTrials',200);
    %% lickLats
    responseInds = find(~isnan([cmbSaline.respondTime]))';
    allReward_R = [cmbSaline(responseInds).rewardR]; 
    allReward_L = [cmbSaline(responseInds).rewardL]; 
    rewardInd = ((allReward_L == 1) + (allReward_R == 1)) > 0;
    lickLatSaline = [cmbSaline(responseInds).respondTime] - [cmbSaline(responseInds).CSon];
    lickLatSalineR = lickLatSaline(logical([0,rewardInd(1:end-1)]));
    lickLatSalineN = lickLatSaline(logical([0,~rewardInd(1:end-1)]));

    responseInds = find(~isnan([cmbDrug.respondTime]))';
    allReward_R = [cmbDrug(responseInds).rewardR]; 
    allReward_L = [cmbDrug(responseInds).rewardL]; 
    rewardInd = ((allReward_L == 1) + (allReward_R == 1)) > 0;
    lickLatDrug = [cmbDrug(responseInds).respondTime] - [cmbDrug(responseInds).CSon];
    lickLatDrugR = lickLatDrug(logical([0,rewardInd(1:end-1)]));
    lickLatDrugN = lickLatDrug(logical([0,~rewardInd(1:end-1)]));
    %% plot lickLat
    figure; 
    subplot(1,3,1);hold on;
    histogram(lickLatSaline,0:50:1500,'Normalization','probability', 'FaceColor', 'm'); 
    histogram(lickLatDrug,0:50:1500,'Normalization','probability', 'FaceColor', 'c');
    title('allLicks')
    legend('saline', 'drug')
    subplot(1,3,2);hold on;
    histogram(lickLatSalineR,0:50:1500,'Normalization','probability', 'FaceColor', 'm'); 
    histogram(lickLatDrugR,0:50:1500,'Normalization','probability', 'FaceColor', 'c');    
    title('afterRwd')
    legend('saline', 'drug')
    subplot(1,3,3);hold on;
    histogram(lickLatSalineN,0:50:1500,'Normalization','probability', 'FaceColor', 'm'); 
    histogram(lickLatDrugN,0:50:1500,'Normalization','probability', 'FaceColor', 'c');   
    title('afterNrwd')
    legend('saline', 'drug')
    %% anovan
    lickLat = [lickLatDrugN, lickLatDrugR, lickLatSalineN, lickLatSalineR];
    g1 = [ones(size([lickLatDrugN, lickLatDrugR])), zeros(size([lickLatSalineN, lickLatSalineR]))]';
    g2 = [ones(size(lickLatDrugN)), zeros(size(lickLatDrugR)), ones(size(lickLatSalineN)), zeros(size(lickLatSalineR))]';
    p = anovan(lickLat,{g1 g2},'model','interaction','varnames',{'drug','reward'});
    %%
    
    
    
    %% nrwd history  25s for drug:rwd, 30s for drug:nrwd; 15s for no pre; 25 second for pre and remove pre:rwd1
    [~, rwdMatxD, noRwdMatxD, combinedPreLickD, combinedLickLatD] = combineLinRegLickLatTime_opMD('combineDrug','combineDrugDrug','drug','binSize', 30000, 'plotFlag', 0,'maxTrials', 200);
    [~, rwdMatxS, noRwdMatxS, combinedPreLickS, combinedLickLatS] = combineLinRegLickLatTime_opMD('combineDrug','combineDrugSaline','saline','binSize',30000, 'plotFlag', 0,'maxTrials', 200);
    %[~, ~, ~, combinedLickLat] = combineLinRegLickLatTime_opMD('combineDrug','combineDrugDrug','drug','binSize',5000,'maxTrials', 350);
    lickLat = [combinedLickLatD, combinedLickLatS];
    preLick = [combinedPreLickD, combinedPreLickS];
    rwdMatx = [rwdMatxD, rwdMatxS];
    noRwdMatx = [noRwdMatxD, noRwdMatxS];
    drug = [ones(size(combinedLickLatD)), (-1)*ones(size(combinedLickLatS))];
 %  tbl = table(drug', rwdMatx(1,:)', rwdMatx(2,:)', rwdMatx(3,:)', noRwdMatx(1,:)',noRwdMatx(2,:)',noRwdMatx(3,:)',lickLat', 'VariableNames', {'drug', 'rwd1','rwd2','rwd3', 'nRwd1', 'nRwd2','nRwd3','lickLat'});
 %  tbl = table(drug', noRwdMatx(1,:)',noRwdMatx(2,:)',noRwdMatx(3,:)',lickLat', 'VariableNames', {'drug','nRwd1', 'nRwd2','nRwd3','lickLat'});
 %  tbl = table(drug', rwdMatx(1,:)', rwdMatx(2,:)', rwdMatx(3,:)',lickLat', 'VariableNames', {'drug', 'rwd1','rwd2','rwd3','lickLat'});
  %  tbl = table(drug', rwdMatx(1,:)', noRwdMatx(1,:)', lickLat', 'VariableNames', {'drug', 'rwd1', 'nRwd1', 'lickLat'});
 %   tbl = table(preLick', drug', rwdMatx(1,:)', noRwdMatx(1,:)', lickLat', 'VariableNames', {'pre', 'drug', 'rwd1', 'nRwd1', 'lickLat'});
      tbl = table(preLick', drug', rwdMatx(1,:)', lickLat', 'VariableNames', {'pre', 'drug', 'rwd1', 'lickLat'});
  %  tbl = table(preLick', drug', rwdMatx(1,:)', [noRwdMatx(1,:)'+rwdMatx(1,:)'], lickLat', 'VariableNames', {'pre', 'drug', 'rwd', 'num', 'lickLat'});
 %  tbl = table(drug', noRwdMatx(1,:)',lickLat', 'VariableNames', {'drug', 'nRwd1','lickLat'});
    mdl = stepwiselm(tbl,'interactions');
    %% kernel rwdHis interaction
    [combinedLickHxD, combinedLickHxTrialD,combinedLickLatD] = lickLatRwdHis('combineDrug','combineDrugdrug','drug','binSize', 10000, 'plotFlag', 0,'maxTrials', 200);
    [combinedLickHxS, combinedLickHxTrialS,combinedLickLatS] = lickLatRwdHis('combineDrug','combineDrugSaline','saline','binSize', 10000, 'plotFlag', 0,'maxTrials', 200);
    lickLat = [combinedLickLatD, combinedLickLatS];
    QMatx = [combinedLickHxD, combinedLickHxS];
    hisMatxT = [combinedLickHxTrialD, combinedLickHxTrialS];
    drug = [ones(size(combinedLickLatD)), (-1)*ones(size(combinedLickLatS))];    
    tbl = table(drug', QMatx', lickLat', 'VariableNames', {'drug', 'histbyTime', 'lickLat'});
%    tbl = table(preLick', drug', noRwdMatx(1,:)', lickLat', 'VariableNames', {'pre', 'drug', 'nrwd1', 'lickLat'});
    mdl = stepwiselm(tbl,'interactions');
    
    %% value sum interaction
    [combineQSumS, combinedLickLatS] = lickLatValueSum('combineDrugSaline', 'saline');
    [combineQSumD, combinedLickLatD] = lickLatValueSum('combineDrugdrug', 'drug');
    lickLat = [combinedLickLatD, combinedLickLatS];
    QMatx = [combineQSumD, combineQSumS];
    drug = [ones(size(combinedLickLatD)), zeros(size(combinedLickLatS))]; 
    tbl = table(zscore(drug'), zscore(QMatx'), zscore(lickLat'), 'VariableNames', {'drug', 'Qsum', 'lickLat'});
%    tbl = table(preLick', drug', noRwdMatx(1,:)', lickLat', 'VariableNames', {'pre', 'drug', 'nrwd1', 'lickLat'});
    mdl = stepwiselm(tbl,'interactions');
    %% plot results
    figure; hold on;
    
    coefVals = mdl.Coefficients.Estimate(2:end);
    CIbands = coefCI(mdl);
    errorL = abs(coefVals - CIbands(2:end,1));
    errorU = abs(coefVals - CIbands(2:end,2));
    in = 1/length(coefVals);
    height = max(abs(CIbands)');
    xlim([0 length(coefVals)+1])
    ylim([min(0, min(1.5*CIbands(2:end,1))) max(0, max(1.5*CIbands(2:end,2)))])
    for i = 1:length(coefVals)
        bar(i,coefVals(i),'FaceColor',[0.5+0.49*in*i 0.5 1-0.49*in*i],'EdgeColor',[0.3+0.69*in*i .2 .8-0.79*in*i],'LineWidth',1.5); hold on;
        errorbar(i,coefVals(i),errorL(i),errorU(i),'.','Color',[0.3+0.69*in*i .2 .8-0.79*in*i],'LineWidth',1.5);
        text(i-0.2,1.2*(sign(coefVals(i))*height(i+1)), mdl.CoefficientNames{i+1})
    end
    title('lrm: on lickLat')
    ylabel('\beta Coefficient')
    text(length(coefVals)-0.5,0.8*max(CIbands(2:end,2)),sprintf('R^2 = %d',mdl.Rsquared.Adjusted))
    %% lick analysis
    % combineLogRegLickTime_opMD generate matrics
    tbl = table(combinedRewardsMatx(1,:)', combinedRewardsMatx(2,:)',combinedRewardsMatx(3,:)',...
                 combinedNoRewardsMatx(1,:)', combinedNoRewardsMatx(2,:)',combinedNoRewardsMatx(3,:)',...
                 combineAntiLicksMatx(1,:)', combineAntiLicksMatx(2,:)', combineAntiLicksMatx(3,:)',...
                 combinedAllChoice_R', 'VariableNames', {'r1', 'r2', 'r3', 'nr1', 'nr2','nr3', 'l1','l2','l3','choice'});
     mdl = stepwiselm(tbl,'interactions');
             %%
    


