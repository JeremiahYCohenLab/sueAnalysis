function combineLickCatAnalysis(animalName, model, col, varargin)
p = inputParser;
p.addParameter('paramNames', {'aN', 'aP', 'aF', 'beta', 'bias'});
p.parse(varargin{:})

[root, sep] = currComputer();

%% load model fitting results and calculate DVs
sampFile = [animalName col '_', model];
path = [root animalName sep animalName 'sorted' sep 'stan' sep 'bernoulli' sep model sep col sep];
load([path sampFile '.mat'], 'dayList');
samples = load([path sampFile '.mat'], sampFile);
samples = samples.(sampFile);
paramNames = p.Results.paramNames;
%decide if time forget
if contains(model,'tF')
    input = 'choice, outcome, ITI)';
else
    input = 'choice, outcome)';
end
% intializing for combining stats
combineLatCat = [];
combineSvS = [];
combineConf = [];
combinePreITI = [];
combinePrePe = [];
combineQsum = [];
combineLick = [];
combineLickRaw = [];
combineRwdHis =[];
combineRwdTimeMat = [];
combineNoRwdTimeMat = [];
combineChoice = [];
combineChoiceMat = [];
combineRwdMat = [];
combineBias = [];
combineSessionID = [];
stayedRwdedProp = zeros(length(dayList),2); % stay when rwded in all long licks or short licks
stayedNoRwdedProp = zeros(length(dayList),2); % stay when rwded in all long licks or short licks
switchedRwdedProp = zeros(length(dayList),2); % switch when norwded in all long licks or short licks
switchedNoRwdProp = zeros(length(dayList),2); % switch when norwded in all long licks or short licks
RwdedStayProp = zeros(length(dayList),2); % stay when rwded in all long licks or short licks
NoRwdedStayProp = zeros(length(dayList),2); % stay when rwded in all long licks or short licks
RwdedSwitchProp = zeros(length(dayList),2); % switch when norwded in all long licks or short licks
NoRwdedSwitchProp = zeros(length(dayList),2); % switch when norwded in all long licks or short licks
paramEsts = zeros(length(dayList), length(paramNames));
logLH = zeros(length(dayList),1);
longCatProp = zeros(length(dayList), 1);


% model 
for ses = 1:length(dayList)
    s = behAnalysisNoPlot_opMD(dayList{ses});
    %generate best estimates of parameters
    allSamples = [];
    edges = cell(1,length(paramNames));
    for i = 1:length(paramNames)
        tmp = samples.(paramNames{i})(:,ses);
        allSamples = [allSamples tmp];
        edges{i} = linspace(min(tmp), max(tmp),50);
    end
    n = histcnd(allSamples,edges); %bin samples by multiple dimensions
    [~, inds] = myMaxAll(n); %find the bin with max num in bin
    for i = 1:length(paramNames) %use median in bin as best estimate
        tmp = allSamples(:,i);
        edgeTmp = edges{i};
        if inds(i) < 50
            paramEsts(ses, i) = median(tmp(tmp >= edgeTmp(inds(i)) & tmp < edgeTmp(inds(i)+1)));
        else
            paramEsts(ses, i) = edgeTmp(inds(i));
        end
    end
    
    choice = s.allChoices';
    choice(choice<0) = 0;
    outcome = abs(s.allRewards);
    % calculate model varaibles
    eval(['[LL,probC,Q,pe] = qLearningModel_' model '(paramEsts(ses,:),' input ';'])
    logLH(ses) = LL/length(s.responseInds);
    % diff value
    Qdiff = abs(Q(:,2)-Q(:,1));
    % total value
    Qsum = sum(Q,2);
    % prepe
    prePe = [NaN pe(1:end-1)];
    % choice confsesence
    choiceConf = 2.*probC - 1;


    % kmeans c = 2;
    ind = s.lickInds;
    svs = [NaN; - ones(length(ind)-1,1)];
    hmm = zeros(length(ind),1);
    svs(s.changeChoice_Inds) = 1;
    hmm(s.hmmStates == 1) = 1;
    preRwd = abs(s.allRewards);
    preRwd = [NaN; preRwd(1:end-1)'];
    preQsum = Qsum;
    preQsum = [NaN; preQsum(1:end-1)];
    bias = zeros(size(ind));
    biasInds = contains(p.Results.paramNames, 'bias');
    if paramEsts(ses,biasInds)>0
        bias(s.lickR_Inds)=1;
    else
        bias(s.lickL_Inds)=1;
    end    

    sessionID = zeros(length(ind), length(dayList));
    sessionID(:,ses) = 1;
    combinePrePe = [combinePrePe; [NaN, zscore(abs(prePe(2:end)))]'];
    combineQsum = [combineQsum; zscore(Qsum)];
    combineConf = [combineConf; choiceConf];
    combinePreITI = [combinePreITI; s.timeBtwn'];
    combineSvS = [combineSvS; svs];
    combineLatCat = [combineLatCat; ind];
    combineLickRaw = [combineLickRaw; s.lickLat'];
    combineRwdTimeMat = [combineRwdTimeMat; s.rwdMatxForLick'];
    combineNoRwdTimeMat = [combineNoRwdTimeMat; s.noRwdMatxForLick'];
    combineChoice = [combineChoice; choice];
    combineChoiceMat = [combineChoiceMat; s.choiceMatx'];
    combineRwdMat = [combineRwdMat; s.rwdMatx'];
    %combineSessionID = [combineSessionID; sessionID]
    combineBias = [combineBias; bias];
    % zscore within each group
    lickLatCatZ = zeros(length(s.lickLat),1);
    lickLatCatZ(ind == 0) = zscore(s.lickLat(ind == 0));
    lickLatCatZ(ind == 1) = zscore(s.lickLat(ind == 1));
    combineLick = [combineLick; lickLatCatZ];
    
    % calculate rwd History
    combineRwdHis = [combineRwdHis; zscore(s.rwdHx')];
    % prop of long cat
    longCatProp(ses) = sum(ind == 1)/length(ind);
    
    % winStay/loseShift
    rwdInds = find(outcome > 0);
    noRwdInds = find(outcome < 1);
    rwdstayInds = intersect(find(svs < 0), find(outcome > 0)+1);
    noRwdstayInds = intersect(find(svs < 0), find(outcome < 1)+1);
    rwdswitchInds = intersect(find(svs > 0), find(outcome > 0)+1);
    noRwdSwitchInds = intersect(find(svs > 0), find(outcome < 1)+1);
    
    stayedRwdedProp(ses,1) = length(intersect(rwdstayInds, find(ind==0)))/length(intersect(find(ind==0), rwdInds+1)); % stay when rwded in all short licks
    stayedRwdedProp(ses,2) = length(intersect(rwdstayInds, find(ind==1)))/length(intersect(find(ind==1), rwdInds+1)); % stay when rwded in all long licks

    stayedNoRwdedProp(ses,1) = length(intersect(noRwdstayInds, find(ind==0)))/length(intersect(find(ind==0), noRwdInds+1)); % stay when nonrwded in all short licks
    stayedNoRwdedProp(ses,2) = length(intersect(noRwdstayInds, find(ind==1)))/length(intersect(find(ind==1), noRwdInds+1)); % stay when nonrwded in all long licks 
    
    switchedRwdedProp(ses,1) = length(intersect(rwdswitchInds, find(ind==0)))/length(intersect(find(ind==0), rwdInds+1));  % switch when rwded in all short licks
    switchedRwdedProp(ses,2) = length(intersect(rwdswitchInds, find(ind==1)))/length(intersect(find(ind==1), rwdInds+1));  % switch when rwded in all long licks
    
    switchedNoRwdProp(ses,1) = length(intersect(noRwdSwitchInds, find(ind==0)))/length(intersect(find(ind==0), noRwdInds+1)); % switch when norwded in all short licks
    switchedNoRwdProp(ses,2) = length(intersect(noRwdSwitchInds, find(ind==1)))/length(intersect(find(ind==1), noRwdInds+1)); % switch when norwded in all long licks
    
    RwdedStayProp(ses,1) = length(intersect(rwdstayInds-1, find(ind==0)))/length(intersect(find(ind==0), rwdInds)); % stay in next trial when rwded in all short licks
    RwdedStayProp(ses,2) = length(intersect(rwdstayInds-1, find(ind==1)))/length(intersect(find(ind==1), rwdInds)); % stay in next trial when rwded in all long licks

    NoRwdedStayProp(ses,1) = length(intersect(noRwdstayInds-1, find(ind==0)))/length(intersect(find(ind==0), noRwdInds)); % stay in next trial when nonrwded in all short licks
    NoRwdedStayProp(ses,2) = length(intersect(noRwdstayInds-1, find(ind==1)))/length(intersect(find(ind==1), noRwdInds)); % stay in next trial when nonrwded in all long licks 
    
    RwdedSwitchProp(ses,1) = length(intersect(rwdswitchInds-1, find(ind==0)))/length(intersect(find(ind==0), rwdInds));  % switch in next trial when rwded in all short licks
    RwdedSwitchProp(ses,2) = length(intersect(rwdswitchInds-1, find(ind==1)))/length(intersect(find(ind==1), rwdInds));  % switch in next trial when rwded in all long licks
    
    NoRwdedSwitchProp(ses,1) = length(intersect(noRwdSwitchInds-1, find(ind==0)))/length(intersect(find(ind==0), noRwdInds)); % switch in next trial when norwded in all short licks
    NoRwdedSwitchProp(ses,2) = length(intersect(noRwdSwitchInds-1, find(ind==1)))/length(intersect(find(ind==1), noRwdInds)); % switch in next trial when norwded in all long licks
    
end

% glm with interactions
tblCat = table(combineRwdHis, combinePrePe, combineQsum, combineConf, combinePreITI, combineSvS, combineBias, combineLatCat, 'VariableNames', {'rwdHis', 'absprePe', 'Qsum', 'confidence', 'preITI', 'svs', 'bias', 'lickCat'});
tblLatShort = table(combineRwdHis(combineLatCat == 0), combinePrePe(combineLatCat == 0), combineQsum(combineLatCat == 0), combineConf(combineLatCat == 0), combinePreITI(combineLatCat == 0), combineSvS(combineLatCat == 0), combineBias(combineLatCat == 0), combineLick(combineLatCat == 0), 'VariableNames', {'rwdHis', 'absprePe', 'Qsum', 'confidence', 'preITI', 'svs', 'bias', 'lickLat'});
tblLatLong = table(combineRwdHis(combineLatCat == 1), combinePrePe(combineLatCat == 1), combineQsum(combineLatCat == 1), combineConf(combineLatCat == 1), combinePreITI(combineLatCat == 1), combineSvS(combineLatCat == 1), combineBias(combineLatCat == 1),combineLick(combineLatCat == 1), 'VariableNames', {'rwdHis', 'absprePe', 'Qsum', 'confidence', 'preITI', 'svs', 'bias', 'lickLat'});

mdlCat = stepwiseglm(tblCat,'interactions','Distribution','binomial', 'Link','logit');
mdlShort = fitlm(tblLatShort,'lickLat ~ rwdHis + Qsum + confidence + preITI + svs + svs*bias');
mdlLong = fitlm(tblLatLong,'lickLat ~ rwdHis + Qsum + confidence + preITI + svs + svs*bias');

% glm for rwd history & lick Lat without interactions
shortlm = fitlm([combineRwdTimeMat(combineLatCat == 0,:) combineNoRwdTimeMat(combineLatCat == 0,:)], combineLick(combineLatCat == 0));
longlm = fitlm([combineRwdTimeMat(combineLatCat == 1,:) combineNoRwdTimeMat(combineLatCat == 1,:)], combineLick(combineLatCat == 1));

% choice autocorrelation 
glm_choiceShort = fitglm([combineChoiceMat(combineLatCat == 0,:)], combineChoice(combineLatCat == 0), 'distribution','binomial','link','logit');
glm_choiceLong = fitglm([combineChoiceMat(combineLatCat == 1,:)], combineChoice(combineLatCat == 1), 'distribution','binomial','link','logit');

% % rwd glm
% glm_rwdShort = fitglm([combineRwdMat(combineLatCat == 0,1:0.5*size(combineRwdMat,2))], combineChoice(combineLatCat == 0), 'distribution','binomial','link','logit');
% glm_rwdLong = fitglm([combineRwdMat(combineLatCat == 1,1:0.5*size(combineRwdMat,2))], combineChoice(combineLatCat == 1), 'distribution','binomial','link','logit');
% 
% %no rwd glm
% glm_noRwdShort = fitglm([combineRwdMat(combineLatCat == 0,(0.5*size(combineRwdMat,2)+1):end)], combineChoice(combineLatCat == 0), 'distribution','binomial','link','logit');
% glm_noRwdLong = fitglm([combineRwdMat(combineLatCat == 1,(0.5*size(combineRwdMat,2)+1):end)], combineChoice(combineLatCat == 1), 'distribution','binomial','link','logit');

% rwd and no rwd glm
glm_rwdShort = fitglm(combineRwdMat(combineLatCat == 0,:), combineChoice(combineLatCat == 0), 'distribution','binomial','link','logit');
glm_rwdLong = fitglm(combineRwdMat(combineLatCat == 1,:), combineChoice(combineLatCat == 1), 'distribution','binomial','link','logit');


%% plot everything
figure;
screenSize = get(0,'Screensize');
screenSize(4) = screenSize(4) - 100;
set(gcf, 'Position', screenSize)
suptitle(animalName)
% check clustering
subplot(4,6,7); hold on;
histogram(combineLickRaw(combineLatCat == 0),0:20:800,'FaceColor', 'c');
histogram(combineLickRaw(combineLatCat == 1),0:20:800,'FaceColor', 'm');
legend({'cluster1', 'cluster2'})
xlabel('lickLat /ms')

subplot(4,6,1);hold on;
colors = cool(length(dayList));
scatter(longCatProp,logLH,15,colors, 'filled');
xlabel('longProp')
ylabel('logLH/trial')

subplot(4,6,2);hold on;
scatter(longCatProp,paramEsts(:,1),15, colors, 'filled');
xlabel('longProp')
ylabel('aN')


subplot(4,6,3);hold on;
scatter(longCatProp,paramEsts(:,3),15, colors, 'filled');
xlabel('longProp')
ylabel('aF')

subplot(4,6,4);hold on;
scatter(longCatProp,paramEsts(:,4),15, colors, 'filled');
xlabel('longProp')
ylabel('beta')

subplot(4,6,5);hold on;
scatter(longCatProp,paramEsts(:,2),15, colors, 'filled');
xlabel('longProp')
ylabel('aP')

subplot(4,6,6);hold on;
scatter(longCatProp,abs(paramEsts(:,5)),15, colors, 'filled');
xlabel('longProp')
ylabel('abs(bias)')


subplot(4,6,8:10);hold on;

coefVals = mdlCat.Coefficients.Estimate(2:end);
CIbands = coefCI(mdlCat);
errorL = abs(coefVals - CIbands(2:end,1));
errorU = abs(coefVals - CIbands(2:end,2));
in = 1/length(coefVals);
height = max(abs(CIbands)');
xlim([0 length(coefVals)+1])
ylim([min(0, min(1.5*CIbands(2:end,1))) max(0, max(1.5*CIbands(2:end,2)))])
for i = 1:length(coefVals)
    bar(i,coefVals(i),'FaceColor',[0.5+0.49*in*i 0.5 1-0.49*in*i],'EdgeColor',[0.3+0.69*in*i .2 .8-0.79*in*i],'LineWidth',1.5); hold on;
    errorbar(i,coefVals(i),errorL(i),errorU(i),'.','Color',[0.3+0.69*in*i .2 .8-0.79*in*i],'LineWidth',1.5);
    text(i-0.4,1.2*(sign(coefVals(i))*height(i+1)), mdlCat.CoefficientNames{i+1})
end
title('glm: on lickGroup')
text(0.5,0.8*max(CIbands(2:end,2)),sprintf('R^2 = %d',mdlCat.Rsquared.Adjusted))
ylabel('\beta Coefficient')

subplot(4,6,11);hold on;

coefVals = mdlShort.Coefficients.Estimate(2:end);
CIbands = coefCI(mdlShort);
errorL = abs(coefVals - CIbands(2:end,1));
errorU = abs(coefVals - CIbands(2:end,2));
in = 1/length(coefVals);
height = max(abs(CIbands)');
xlim([0 length(coefVals)+1])
ylim([min(0, min(1.5*CIbands(2:end,1))) max(0, max(1.5*CIbands(2:end,2)))])
for i = 1:length(coefVals)
    bar(i,coefVals(i),'FaceColor',[0.5+0.49*in*i 0.5 1-0.49*in*i],'EdgeColor',[0.3+0.69*in*i .2 .8-0.79*in*i],'LineWidth',1.5); hold on;
    errorbar(i,coefVals(i),errorL(i),errorU(i),'.','Color',[0.3+0.69*in*i .2 .8-0.79*in*i],'LineWidth',1.5);
    text(i-0.4,1.2*(sign(coefVals(i))*height(i+1)), mdlShort.CoefficientNames{i+1})
end
title('lm: on lickLat in short group')
ylabel('\beta Coefficient')

subplot(4,6,12);hold on;

coefVals = mdlLong.Coefficients.Estimate(2:end);
CIbands = coefCI(mdlLong);
errorL = abs(coefVals - CIbands(2:end,1));
errorU = abs(coefVals - CIbands(2:end,2));
in = 1/length(coefVals);
height = max(abs(CIbands)');
xlim([0 length(coefVals)+1])
ylim([min(0, min(1.5*CIbands(2:end,1))) max(0, max(1.5*CIbands(2:end,2)))])
for i = 1:length(coefVals)
    bar(i,coefVals(i),'FaceColor',[0.5+0.49*in*i 0.5 1-0.49*in*i],'EdgeColor',[0.3+0.69*in*i .2 .8-0.79*in*i],'LineWidth',1.5); hold on;
    errorbar(i,coefVals(i),errorL(i),errorU(i),'.','Color',[0.3+0.69*in*i .2 .8-0.79*in*i],'LineWidth',1.5);
    text(i-0.4,1.2*(sign(coefVals(i))*height(i+1)), mdlLong.CoefficientNames{i+1})
end
title('lm: on lickLat in long group')
ylabel('\beta Coefficient')

subplot(4,6,[17,18]);hold on;
tMax = size(s.rwdMatxForLick,1);
binSize = s.timeBinSize;
relevInds = 2:tMax+1;
coefVals = shortlm.Coefficients.Estimate(relevInds);
CIbands = coefCI(shortlm);
errorL = abs(coefVals - CIbands(relevInds,1));
errorU = abs(coefVals - CIbands(relevInds,2));
errorbar(((1:tMax)*binSize/1000),coefVals,errorL,errorU,'Color', colors(end,:),'linewidth',2)

relevInds = tMax+2:2*tMax+1;
coefVals = shortlm.Coefficients.Estimate(relevInds);
CIbands = coefCI(shortlm);
errorL = abs(coefVals - CIbands(relevInds,1));
errorU = abs(coefVals - CIbands(relevInds,2));
errorbar(((1:tMax)*binSize/1000),coefVals,errorL,errorU,'Color', colors(1,:),'linewidth',2)
line([0 tMax*binSize/1000], [0 0], 'Color',[0.5 0.5 0.5],'LineStyle','--')

legend('reward', 'no reward')
xlabel('reward n seconds back')
ylabel('\beta Coefficient')
xlim([0 (tMax*binSize/1000 + 5)])
set(gca, 'tickdir', 'out')
title('RewardHis on lickLat in short group')

subplot(4,6,[23,24]);hold on;
relevInds = 2:tMax+1;
coefVals = longlm.Coefficients.Estimate(relevInds);
CIbands = coefCI(longlm);
errorL = abs(coefVals - CIbands(relevInds,1));
errorU = abs(coefVals - CIbands(relevInds,2));
errorbar(((1:tMax)*binSize/1000),coefVals,errorL,errorU,'Color', colors(end,:),'linewidth',2)

relevInds = tMax+2:2*tMax+1;
coefVals = longlm.Coefficients.Estimate(relevInds);
CIbands = coefCI(longlm);
errorL = abs(coefVals - CIbands(relevInds,1));
errorU = abs(coefVals - CIbands(relevInds,2));
errorbar(((1:tMax)*binSize/1000),coefVals,errorL,errorU,'Color', colors(1,:),'linewidth',2)
line([0 tMax*binSize/1000], [0 0], 'Color',[0.5 0.5 0.5],'LineStyle','--')

legend('reward', 'no reward')
xlabel('reward n seconds back')
ylabel('\beta Coefficient')
xlim([0 (tMax*binSize/1000 + 5)])
set(gca, 'tickdir', 'out')
title('RewardHis on lickLat in long group')

subplot(4,6,13); hold on 
scatter(switchedRwdedProp(:,1), switchedNoRwdProp(:,1), 20, 'c', 'filled'); 
scatter(switchedRwdedProp(:,2), switchedNoRwdProp(:,2), 20, 'm', 'filled')
legend({'short', 'long'})
xlabel('switch from rwd')
ylabel('switch from no rwd')
title('stay or switch on current trial')
rotate3d on

subplot(4,6,14); hold on 
scatter(RwdedSwitchProp(:,1), NoRwdedSwitchProp(:,1), 20, 'c', 'filled'); 
scatter(RwdedSwitchProp(:,2), NoRwdedSwitchProp(:,2), 20, 'm', 'filled')
legend({'short', 'long'})
xlabel('switch from rwd')
ylabel('switch from no rwd')
title('stay or switch on next trial')
rotate3d on

subplot(4,6,[15,16]); hold on
tMax = 0.5*size(combineRwdMat,2);
relevInds = 2:tMax+1;
coefVals = glm_rwdShort.Coefficients.Estimate(relevInds);
CIbands = coefCI(glm_rwdShort);
errorL = abs(coefVals - CIbands(relevInds,1));
errorU = abs(coefVals - CIbands(relevInds,2));
errorbar((1:tMax)+0.2,coefVals,errorL,errorU,'Color', colors(1,:),'linewidth',2)

coefVals = glm_rwdLong.Coefficients.Estimate(relevInds);
CIbands = coefCI(glm_rwdLong);
errorL = abs(coefVals - CIbands(relevInds,1));
errorU = abs(coefVals - CIbands(relevInds,2));
errorbar((1:tMax)+0.2,coefVals,errorL,errorU,'Color', colors(end,:),'linewidth',2)

line([0 tMax+0.2], [0 0], 'Color',[0.5 0.5 0.5],'LineStyle','--')

legend({'short', 'long'})
title('glm: rwd on choice')

subplot(4,6,[21,22]); hold on
relevInds = 2+tMax:2*tMax+1;
coefVals = glm_rwdShort.Coefficients.Estimate(relevInds);
CIbands = coefCI(glm_rwdShort);
errorL = abs(coefVals - CIbands(relevInds,1));
errorU = abs(coefVals - CIbands(relevInds,2));
errorbar((1:tMax)+0.2,coefVals,errorL,errorU,'Color', colors(1,:),'linewidth',2)

coefVals = glm_rwdLong.Coefficients.Estimate(relevInds);
CIbands = coefCI(glm_rwdLong);
errorL = abs(coefVals - CIbands(relevInds,1));
errorU = abs(coefVals - CIbands(relevInds,2));
errorbar((1:tMax)+0.2,coefVals,errorL,errorU,'Color', colors(end,:),'linewidth',2)

line([0 tMax+0.2], [0 0], 'Color',[0.5 0.5 0.5],'LineStyle','--')

legend({'short', 'long'})
title('glm: norwd on choice')

% subplot(4,6,[19,20]); hold on
% relevInds = 2:tMax+1;
% coefVals = glm_choiceShort.Coefficients.Estimate(relevInds);
% CIbands = coefCI(glm_choiceShort);
% errorL = abs(coefVals - CIbands(relevInds,1));
% errorU = abs(coefVals - CIbands(relevInds,2));
% errorbar((1:tMax)+0.2,coefVals,errorL,errorU,'Color', colors(1,:),'linewidth',2)
% 
% coefVals = glm_choiceLong.Coefficients.Estimate(relevInds);
% CIbands = coefCI(glm_choiceLong);
% errorL = abs(coefVals - CIbands(relevInds,1));
% errorU = abs(coefVals - CIbands(relevInds,2));
% errorbar((1:tMax)+0.2,coefVals,errorL,errorU,'Color', colors(end,:),'linewidth',2)
% 
% line([0 tMax+0.2], [0 0], 'Color',[0.5 0.5 0.5],'LineStyle','--')
% 
% legend({'short', 'long'})
% title('glm: choice autoCorr')

end
