% random simulation for laser effect
clear all;
model = '5paramsLaserNegRPERotation';
paramNames = getParamNames_dF(model, 1);
sessionNum = 80;
maxTrial = 700;
allParams = [0.5+0.2*betarnd(1, 3, sessionNum, 1), ...
    0.5+0.3*betarnd(7, 8, sessionNum, 1), ...
    betarnd(2, 4, sessionNum, 1), ...
    3*betarnd(5, 5, sessionNum, 1)+4, ...
    0.3*betarnd(5, 5, sessionNum, 1)-0.3,...
    zeros(sessionNum, 1)];
figure2;
for i = 1:size(allParams,2)
    edges = linspace(min(allParams(:,i))-0.01, max(allParams(:,i))+0.01, 10);
    subplot(1, size(allParams,2), i);
    histogram(allParams(:,i), edges);
    title(paramNames{i});
end
%% simulation with behavior
clear all;
ani = 'allGt';
sheet = 'inhibitionGt';
col = 'cueOnShamLate';
modelName = '5params_k_bias_LaserDisengageScale';
maxTrial = 400;
dayList = getDayList(sheet, ani, col);
[root, sep] =  currComputer();
aniList = cellfun(@(x) x(2:6), dayList, 'UniformOutput', 0);
aniList = unique(aniList);
sessionNum = 10;
params = [];
for i = 1:length(dayList)
    pd = parseSessionString_df(dayList{i}, root, sep);
    paramsTemp = getStanModelParams_sampsOnly(pd.aniName, col, modelName, sessionNum, 'sessionName', [], 'biasFlag',0, 'sessionParamsFlag', 0);
    params = [params; paramsTemp];
end
allParams = [params, zeros(size(params,1),1)];

%% check model consistency
clear all;
i = 2;
[root, sep] = currComputer();
col = 'cueOnGood';
ani = 'ZS085';
modelName = '5params_k_bias_LaserNegRPE';
filename = [ani col '_' modelName];
load([root ani sep ani 'sorted' sep 'stan\bernoulli' sep modelName sep col sep filename '.mat']);
eval(['samps =' filename ';'])
eval(['clear ' filename]);
sampNum = 2;
inds = randperm(length(samps.mu_a), sampNum);
paramNames = getParamNames_dF(modelName, 1);
sessionInd = 1;

for currP = 1:length(paramNames)
    tmp = eval(['samps.' paramNames{currP}]);
    tmp = tmp(:,sessionInd);

    params(:,currP) = tmp(inds);
end

ll = samps.log_lik(inds,sessionInd);


for m = 1:size(params,1)
    s = behAnalysisNoPlot_opMD(dayList{sessionInd}, 'simpleFlag', 1);
    choice = s.allChoices;
    choice(choice==-1)= 0;
    LH(m) = qLearningModel_5params_k_bias_LaserNegRPE(params(m,:), choice, abs(s.allRewards), s.laser);
end
%%
for i = 1:size(allParams,1)
    [t(i), allRewards{i}, allChoices{i}, ~, ~, laser{i}] = qLearningModel_5paramsLaserNegRPERotation_simNoPlot(allParams(i,:), maxTrial, i+150*i);
end
%% LaserNegRPE
for i = 1:size(allParams,1)
    [t(i), allRewards{i}, allChoices{i}, ~, ~, laser{i}] = qLearningModel_5params_k_bias_LaserNegRPE_simNoPlot(allParams(i,:), maxTrial,  16+150*i);
end
%% laserDisengage
for i = 1:size(allParams,1)
    [t(i), allRewards{i}, allChoices{i}, ~, ~, laser{i}] = qLearningModel_5params_k_bias_LaserDisengage_simNoPlot(allParams(i,:), maxTrial,  16+150*i);
end
%%
for i = 1:size(allParams,1)
    [t(i), allRewards{i}, allChoices{i}, ~, ~, laser{i}] = qLearningModel_5params_k_bias_LaserNegOnlyRPE_simNoPlot(allParams(i,:), maxTrial,  2+150*i);
end
%%
for i = 1:size(allParams,1)
    [t(i), allRewards{i}, allChoices{i}, ~, ~, laser{i}] = qLearningModel_5params_k_bias_LaserNeg(allParams(i,:), maxTrial,  2+150*i);
end
%% check simulation
for i = 1:size(allParams,1)
    Q = t(i).Q;
    choiceInd = allChoices{i}*0.5 + 1.5;
    Qdiff = allParams(i,4)*(Q(:,2) - Q(:,1));
    Qdiff = Qdiff(2:end);
    Qdiff(allChoices{i}==-1) = -Qdiff(allChoices{i}==-1);
    s = behAnalysisNoPlotSim_opMD(allChoices{i}, allRewards{i}, laser{i});
    pStay = t(i).pRight;
    pStay = [pStay(1:end-1);NaN];
    pStay(allChoices{i}'== -1) = 1-pStay(allChoices{i}'== -1);
    figure2; 
    subplot(3,1,1);
    plot([1:length(s.allChoices); 1:length(s.allChoices)], 0.5*[zeros(1, length(s.allChoices)); s.allChoices + s.allRewards], 'k');
%     subplot(1,2,1); hold on;
%     edges = linspace(min(Qdiff)-0.001, max(Qdiff)+0.001, 15);
%     histogram(Qdiff(intersect(find(s.laser==1), s.rwd_Inds)), edges, 'FaceColor', [0 0 1], 'Normalization', 'probability');
%     histogram(Qdiff(intersect(find(s.laser==0), s.rwd_Inds)), edges, 'FaceColor', 'k', 'Normalization', 'probability');
%     subplot(1,2,2); hold on;
%     edges = linspace(min(pStay)-0.001, max(pStay)+0.001, 15);
%     histogram(pStay(intersect(find(s.laser==1), s.rwd_Inds)), edges, 'FaceColor', [0 0 1], 'Normalization', 'probability');
%     histogram(pStay(intersect(find(s.laser==0), s.rwd_Inds)), edges, 'FaceColor', 'k', 'Normalization', 'probability');    
end
%%
%% scatter plots pSwitch
lcInhi = zeros(size(allParams,1), 1);
lcCtrl = zeros(size(allParams,1), 1);
wsInhi = zeros(size(allParams,1), 1);
wsCtrl = zeros(size(allParams,1), 1);
for i = 1:length(t)
    s = behAnalysisNoPlotSim_opMD(allChoices{i}, allRewards{i}, laser{i});
    lcInhi(i) = length(mintersect(s.nrwd_Inds, find(s.laser ==1), s.changeChoice_Inds-1))/length(intersect(s.nrwd_Inds, find(s.laser ==1)));
    lcCtrl(i) = length(mintersect(s.nrwd_Inds, find(s.laser ==0), s.changeChoice_Inds-1))/length(intersect(s.nrwd_Inds, find(s.laser ==0)));
    wsInhi(i) = length(mintersect(s.rwd_Inds, find(s.laser ==1), s.stayChoice_Inds-1))/length(mintersect(s.rwd_Inds', find(s.laser' ==1), 1:length(allChoices{i})-1 ));
    wsCtrl(i) = length(mintersect(s.rwd_Inds, find(s.laser ==0), s.stayChoice_Inds-1))/length(mintersect(s.rwd_Inds', find(s.laser' ==0), 1:length(allChoices{i})-1 ));
end
%%
figure2; 
subplot(1,2,1); hold on;
scatter(lcCtrl, lcInhi, 20, 'm', 'Filled');
plot([0 0.6], [0 0.6], 'LineStyle', '--', 'Color', [0.6 0.6 0.6]);
set(gca, 'TickDir', 'out');
p = signrank(lcInhi, lcCtrl);
title(['loss switch' num2str(p)])
xlabel('no laser')
ylabel('laser')

subplot(1,2,2); hold on;
scatter(wsCtrl, wsInhi, 20, 'm', 'Filled');
plot([0.7 1], [0.7 1], 'LineStyle', '--', 'Color', [0.6 0.6 0.6]);
xlim([0.8 1.02])
ylim([0.8 1.02])
set(gca, 'TickDir', 'out');
p = signrank(wsInhi, wsCtrl);
title(['win stay' num2str(p)])
xlabel('no laser')

%% scatter plots estimated swt
aniNames = cellfun(@(x) x(2:6), dayList, 'UniformOutput', false);
anis = unique(aniNames);
lcLaserMean = NaN(size(anis));
lcLaserSem = NaN(size(anis));
lcControlMean = NaN(size(anis));
lcControlSem = NaN(size(anis));
wsLaserMean = NaN(size(anis));
wsLaserSem = NaN(size(anis));
wsControlMean = NaN(size(anis));
wsControlSem = NaN(size(anis));

for i = 1:length(anis)
    sessionInd = cellfun(@(x) strcmp(x, anis{i}), aniNames);
    sessionInd = ((find(sessionInd==1, 1)-1)*sessionNum+1) : (find(sessionInd==1, 1, 'last'))*sessionNum;
    lcLaserMean(i) = mean(lcInhi(sessionInd));
    lcControlMean(i) = mean(lcCtrl(sessionInd));
    lcLaserSem(i) = sem(lcInhi(sessionInd));
    lcControlSem(i) = sem(lcCtrl(sessionInd));

    wsLaserMean(i) = mean(wsInhi(sessionInd));
    wsControlMean(i) = mean(wsCtrl(sessionInd));
    wsLaserSem(i) = sem(wsInhi(sessionInd));
    wsControlSem(i) = sem(wsCtrl(sessionInd));
end
%% 
figure2;
subplot(1,2,1); hold on;
errorbar(lcControlMean, lcLaserMean, lcLaserSem, lcLaserSem, lcControlSem, lcControlSem, "LineStyle","none", 'Color', 'k', 'Marker','diamond', 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'none')
plot([0 0.6], [0 0.6], 'LineStyle', '--', 'Color', [0.6 0.6 0.6]);
subplot(1,2,2); hold on;
errorbar(wsControlMean, wsLaserMean, wsLaserSem, wsLaserSem, wsControlSem, wsControlSem, "LineStyle","none", 'Color', 'k', 'Marker','diamond', 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k')
plot([0.7 1], [0.7 1], 'LineStyle', '--', 'Color', [0.6 0.6 0.6]);

%%
figure2; 
subplot(1,2,1); hold on;
scatter(lcCtrl, lcInhi, 20, 'm', 'Filled');
plot([0 0.6], [0 0.6], 'LineStyle', '--', 'Color', [0.6 0.6 0.6]);
set(gca, 'TickDir', 'out');
p = signrank(lcInhi, lcCtrl);
title(['loss switch' num2str(p)])
xlabel('no laser')
ylabel('laser')

subplot(1,2,2); hold on;
scatter(wsCtrl, wsInhi, 20, 'm', 'Filled');
plot([0.7 1], [0.7 1], 'LineStyle', '--', 'Color', [0.6 0.6 0.6]);
xlim([0.8 1.02])
ylim([0.8 1.02])
set(gca, 'TickDir', 'out');
p = signrank(wsInhi, wsCtrl);
title(['win stay' num2str(p)])
xlabel('no laser')
%% peri-stimulus analysis
% pre-locate all params
ani = 'allGt';
sheet = 'inhibitionGt';
col = 'cueOnGood';
dayList = getDayList(sheet, ani, col);
[root, sep] = currComputer();
trialB = 3;
trialF = 3;
periStiSwitchesNextInhi = [];
periStiRwdsInhi = [];
periStiSwitchesNextCtrl = [];
periStiRwdsCtrl = [];
periStiChoiceCtrl = [];
periStiChoiceInhi = [];
% find all lasers
for i = 1:length(t)
    % getChoices, Switches and Hmm
    s = behAnalysisNoPlotSim_opMD(allChoices{i}', allRewards{i}', laser{i}');
%     session = dayList{i};
%     s = behAnalysisNoPlot_opMD(session, 'simpleFlag', 1);
%     hmmInd = cellfun(@(x) strcmp(session, x), allSession);
%     states = allStates{hmmInd};
%     transEInd = find(states(2:end-1)==1 & states(1:end-2)~=states(3:end) & states(1:end-2)~= 1 & states(3:end)~= 1)+1;
%     states(transEInd) = NaN;
    choices = s.allChoices';
    choices(choices==-1) = 0;
    svs = NaN(size(choices));
    svs(s.changeChoice_Inds) = 1;
    svs(s.stayChoice_Inds) = 0;
    svsNext = NaN(size(choices));
    svsNext(s.changeChoice_Inds-1) = 1;
    svsNext(s.stayChoice_Inds-1) = 0;
    % find trial Inds with laser
    laserInd = find(s.laser==1);
    noLaserInd = find(s.laser==0);
    
    
    choiceTemp = [NaN(1, trialB), choices, NaN(1, trialF)];
    laserTemp = [NaN(1, trialB), s.laser', NaN(1, trialF)];
    svsTemp = [NaN(1, trialB), svs, NaN(1, trialF)];
    svsNextTemp = [NaN(1, trialB), svsNext, NaN(1, trialF)];
    rwdTemp = [NaN(1, trialB), abs(s.allRewards'), NaN(1, trialF)];

    % focus on inhibition trials
    focus = laserInd;
    currInd = mintersect(focus, s.nrwd_Inds, s.changeChoice_Inds-1);
    currInd = currInd(currInd<=length(choices));
    currInd = currInd + trialB;
    for j = 1:length(currInd)
        % rwd
        currSeq = rwdTemp(currInd(j)-trialB:currInd(j)+trialF);
        periStiRwdsInhi = [periStiRwdsInhi; currSeq];

        % Switch
        currSeq = svsNextTemp(currInd(j)-trialB:currInd(j)+trialF);
        periStiSwitchesNextInhi = [periStiSwitchesNextInhi; currSeq];
        
        % choice
        currSeq = choiceTemp(currInd(j)-trialB:currInd(j)+trialF);
        if choiceTemp(currInd(j))==0
            currSeq = 1-currSeq;
        end
        periStiChoiceInhi = [periStiChoiceInhi; currSeq];
    end

    % focus on control trials
    focus = noLaserInd;
    currInd = mintersect(focus, s.nrwd_Inds, s.changeChoice_Inds-1);  
    currInd = currInd(currInd<=length(choices));
    currInd = currInd + trialB;
    for j = 1:length(currInd)
        % rwd
        currSeq = rwdTemp(currInd(j)-trialB:currInd(j)+trialF);
        periStiRwdsCtrl = [periStiRwdsCtrl; currSeq];

        % Switch
        currSeq = svsNextTemp(currInd(j)-trialB:currInd(j)+trialF);
        periStiSwitchesNextCtrl = [periStiSwitchesNextCtrl; currSeq];
        
        % choice
        currSeq = choiceTemp(currInd(j)-trialB:currInd(j)+trialF);
        if choiceTemp(currInd(j))==0
            currSeq = 1-currSeq;
        end
        periStiChoiceCtrl = [periStiChoiceCtrl; currSeq];
    end
end
%% P(switch|nrwd)
nrwdInd = periStiRwdsCtrl ==0;
nrwdSwNextMat = periStiSwitchesNextCtrl;
nrwdSwNextMat(~nrwdInd) = 0;
nrwdNum = sum(nrwdInd, 1, 'omitnan');
nrwdSwNum = sum(nrwdSwNextMat, 1, 'omitnan');

meanLCCtrl = nrwdSwNum./nrwdNum;
semLCCtrl = sem_bernoulli(nrwdSwNum, nrwdNum);

nrwdInd = periStiRwdsInhi ==0;
nrwdSwNextMat = periStiSwitchesNextInhi;
nrwdSwNextMat(~nrwdInd) = 0;
nrwdNum = sum(nrwdInd, 1, 'omitnan');
nrwdSwNum = sum(nrwdSwNextMat, 1, 'omitnan');

meanLCInhi = nrwdSwNum./nrwdNum;
semLCInhi = sem_bernoulli(nrwdSwNum, nrwdNum);

colorCtrl = [0.4 0.4 0.4];
colorInhi = [0.2 0.2 1];
figure2; hold on;
errorbar(-trialB:1:trialF,  meanLCCtrl, semLCCtrl, 'Color', colorCtrl, 'LineWidth', 2)
errorbar(-trialB:1:trialF,  meanLCInhi, semLCInhi, 'Color', colorInhi, 'LineWidth', 2)
ylabel('P(switch|nrwd)')
%% Choice
meanChoiceCtrl = mean(periStiChoiceCtrl, 'omitnan');
meanChoiceInhi = mean(periStiChoiceInhi, 'omitnan');
semChoiceCtrl = sem(periStiChoiceCtrl);
semChoiceInhi = sem(periStiChoiceInhi);

colorCtrl = [0.4 0.4 0.4];
colorInhi = [0.2 0.2 1];
figure2; hold on;
errorbar(-trialB:1:trialF,  meanChoiceCtrl, semChoiceCtrl, 'Color', colorCtrl, 'LineWidth', 2)
errorbar(-trialB:1:trialF,  meanChoiceInhi, semChoiceInhi, 'Color', colorInhi, 'LineWidth', 2)
ylabel('P(choice)')
%%


% title(currTitle)
%%

%% GLM analysis
combineRwdMat = [];
combineNrwdMat = [];
combineChoices = [];
combineLaser = [];
combineLaserMat = [];
maxBack = 10;

for i = 1:length(allChoices)
    s = behAnalysisNoPlotSim_opMD(allChoices{i}', allRewards{i}', laser{i}');
    rwdMatTmp = NaN(maxBack, length(allChoices{i}));
    nRwdMatTmp = NaN(maxBack, length(allChoices{i}));
    laserMatTmp = NaN(maxBack, length(allChoices{i}));
    choices = s.allChoices';
    choices(choices<0) = 0;
   
    for j = 1:maxBack
        tmp = [NaN(1,j) s.allRewards(1:end-j)'];
        rwdMatTmp(j,:) = tmp;
        tmp = [NaN(1,j) s.allNoRewards(1:end-j)'];
        nRwdMatTmp(j,:) = tmp;
        tmp = [NaN(1,j) s.laser(1:end-j)'];
        laserMatTmp(j,:) = tmp;
    end

    combineRwdMat = [combineRwdMat NaN(maxBack, 100) rwdMatTmp];
    combineNrwdMat = [combineNrwdMat NaN(maxBack, 100) nRwdMatTmp];
    combineChoices = [combineChoices NaN(1, 100) choices];
    combineLaser = [combineLaser NaN(1,100) s.laser'];
    combineLaserMat = [combineLaserMat NaN(maxBack, 100) laserMatTmp];
end
%% GLMs
% raw choices
glm = fitglm([combineRwdMat' combineNrwdMat'], combineChoices,'distribution','binomial','Link','logit');

coeffs = glm.Coefficients.Estimate(2:end);
CIbands = coefCI(glm);
CIbands = CIbands(2:end,:);
errorL = abs(coeffs - CIbands(:,2));
errorU = abs(coeffs - CIbands(:,1));

figure2; hold on;
currInd = 1:maxBack;
currNrwdInd = (maxBack+1):2*maxBack;
errorbar(currInd,coeffs(currInd),errorL(currInd),errorU(currInd),'Color', 'r','linewidth',2);
errorbar(currInd,coeffs(currNrwdInd),errorL(currNrwdInd),errorU(currNrwdInd),'Color', 'k','linewidth',2);
line([0 maxBack], [0 0], 'Color',[0.5 0.5 0.5],'LineStyle','--')
title('no separation')
%%
% raw choices + laser on nrwd

combineLaserNrwd = combineNrwdMat;
combineLaserNrwd(combineLaserMat == 0) = 0;

glm = fitglm([combineRwdMat' combineNrwdMat' combineLaserNrwd'], combineChoices,'distribution','binomial','Link','logit');

coeffs = glm.Coefficients.Estimate(2:end);
CIbands = coefCI(glm);
CIbands = CIbands(2:end,:);
errorL = abs(coeffs - CIbands(:,2));
errorU = abs(coeffs - CIbands(:,1));

figure2; hold on;
currInd = 1:maxBack;
currNrwdInd = (maxBack+1):2*maxBack;
currNrwdLaserInd = 2*maxBack+1:3*maxBack;
errorbar(currInd,coeffs(currInd),errorL(currInd),errorU(currInd),'Color', 'r','linewidth',2);
errorbar(currInd,coeffs(currNrwdInd),errorL(currNrwdInd),errorU(currNrwdInd),'Color', 'k','linewidth',2);
errorbar(currInd,coeffs(currNrwdLaserInd),errorL(currNrwdLaserInd),errorU(currNrwdLaserInd),'Color', 'b','linewidth',2);
line([0 maxBack], [0 0], 'Color',[0.5 0.5 0.5],'LineStyle','--')
title('allNrwd+laserNrwd')
%% rwd choices, nrwd laser/no laser
% raw choices + laser on nrwd

combineLaserNrwd = combineNrwdMat;
combineLaserNrwd(combineLaserMat == 0) = 0;

combineNoLaserNrwd = combineNrwdMat;
combineNoLaserNrwd(combineLaserMat == 1) = 0;

glm = fitglm([combineRwdMat' combineNoLaserNrwd' combineLaserNrwd'], combineChoices,'distribution','binomial','Link','logit');

coeffs = glm.Coefficients.Estimate(2:end);
CIbands = coefCI(glm);
CIbands = CIbands(2:end,:);
errorL = abs(coeffs - CIbands(:,2));
errorU = abs(coeffs - CIbands(:,1));

figure2; hold on;
currInd = 1:maxBack;
currNrwdNoLaserInd = (maxBack+1):2*maxBack;
currNrwdLaserInd = 2*maxBack+1:3*maxBack;
errorbar(currInd,coeffs(currInd),errorL(currInd),errorU(currInd),'Color', 'r','linewidth',2);
errorbar(currInd,coeffs(currNrwdNoLaserInd),errorL(currNrwdNoLaserInd),errorU(currNrwdNoLaserInd),'Color', 'k','linewidth',2);
errorbar(currInd,coeffs(currNrwdLaserInd),errorL(currNrwdLaserInd),errorU(currNrwdLaserInd),'Color', 'b','linewidth',2);
title('noLaserNrwd+laserNrwd')
%%
combineLaserRwd = combineRwdMat;
combineLaserRwd(combineLaserMat == 0) = 0;

combineNoLaserRwd = combineRwdMat;
combineNoLaserRwd(combineLaserMat == 1) = 0;

glm = fitglm([combineLaserRwd' combineNoLaserRwd' combineNrwdMat'], combineChoices,'distribution','binomial','Link','logit');

coeffs = glm.Coefficients.Estimate(2:end);
CIbands = coefCI(glm);
CIbands = CIbands(2:end,:);
errorL = abs(coeffs - CIbands(:,2));
errorU = abs(coeffs - CIbands(:,1));

figure2; hold on;
currNrwdInd = 2*maxBack+1:3*maxBack; 
currRwdNoLaserInd = (maxBack+1):2*maxBack;
currRwdLaserInd = 1:maxBack;
errorbar(currInd,coeffs(currInd),errorL(currInd),errorU(currInd),'Color', 'r','linewidth',2);
errorbar(currInd,coeffs(currRwdNoLaserInd),errorL(currRwdNoLaserInd),errorU(currRwdNoLaserInd),'Color', 'k','linewidth',2);
errorbar(currInd,coeffs(currRwdLaserInd),errorL(currRwdLaserInd),errorU(currRwdLaserInd),'Color', 'b','linewidth',2);
title('noLaserRwd+laserRwd')
%%









%% stan fitting and comparison
choiceTmp = {};
outcomeTmp = {};
laserTmp = {};
Tsesh = [];
for i = 1:length(allChoices)
    behavStruct = behAnalysisNoPlotSim_opMD(allChoices{i}, allRewards{i}, laserComb{i});
    choiceTmp{i} = behavStruct.allChoices;
    choiceTmp{i}(choiceTmp{i} == -1) = 0;
    outcomeTmp{i} = abs(behavStruct.allRewards); 
    laserTmp{i} = [behavStruct.laser];
    Tsesh(i,1) = length(outcomeTmp{i});
end

T = max(Tsesh);
N = length(allChoices);
choice = zeros(N, T);
outcome = zeros(N, T);
laserComb = zeros(N, T);

for i = 1:N
    choice(i, 1:Tsesh(i)) = choiceTmp{i};
    outcome(i, 1:Tsesh(i)) = outcomeTmp{i};
    laserComb(i, 1:Tsesh(i)) = laserTmp{i};
end

session_dat = struct('N',N,'T',T, 'Tsesh', Tsesh, 'choice', choice, 'outcome', outcome, 'laser', laserComb);
%%
filePath = 'C:\Users\zhixi\Documents\gitRepositories\sueAnalysis\matlabCode\operantMatching\learningModels\stan\bernoulli\';
fullName = 'stan_qLearning_5paramsLaserNegRPE.stan';
savePath = 'F:\tmpData';
fit = stan('file',[filePath fullName],'data',session_dat,'verbose',true,...
            'iter', 10000, 'warmup', 5000, 'working_dir', savePath, 'chains', 8, 'refresh', 200);
%%
params = allParams;

%%
i = 20;
choices = allChoices{i};
choices(choices==-1) = 0;
outcomes = abs(allRewards{i});
probChosenSim = t(i).pRight;
probChosenSim(choices == 0) = 1 - probChosenSim(choices == 0);
[LH, probChosen, Q, pe] = qLearningModel_5paramsLaserNegRPE(params(i,:), choices, outcomes, laser{i});
figure; subplot(3,1,1); hold on
plot([1:length(choices); 1:length(choices)], 0.5*[zeros(1, length(choices)); allChoices{i}' + allRewards{i}'], 'k')
subplot(3,1,1); hold on;
plot(1:length(choices), probChosen, 'LineWidth', 2, 'Color', 'r');
plot(1:length(choices), probChosenSim, 'LineWidth', 2, 'Color', 'b', 'LineStyle', '--');

%% stan fitting
clear all;
i = 2;
col = 'cueOnGood';
ani = 'ZS085';
modelName = '5paramsLaserNegRPE';
filename = 'ZS085cueOnGood_5paramsLaserNegRPE';
load(['F:\ZS085\ZS085sorted\stan\bernoulli\5paramsLaserNegRPE\cueOnGood\' filename '.mat']);
eval(['samps =' filename ';'])
eval(['clear ' filename]);
sampNum = 20;
inds = randperm(length(samps.mu_aN), sampNum);
paramNames = getParamNames_dF(modelName, 1);
sessionInd = 1;

for currP = 1:length(paramNames)
    tmp = eval(['samps.' paramNames{currP}]);
    tmp = tmp(:,sessionInd);

    params(:,currP) = tmp(inds);
end

ll = samps.log_lik(inds,sessionInd);
%%
s = behAnalysisNoPlot_opMD(dayList{i}, 'simpleFlag', 1);
choices = s.allChoices;
choices(choices==-1) = 0;
outcomes = abs(s.allRewards);
for j = 1:size(params,1)
    [LH(j), probChosen{j}, Q{j}, pe{j}] = qLearningModel_5paramsLaserNegRPE(params(j,:), choices, outcomes, s.laser);
end
figure; subplot(3,1,1); hold on
plot([1:length(choices); 1:length(choices)], 0.5*[zeros(1, length(choices)); s.allChoices + s.allRewards], 'k');
scatter(find(s.laser==1), 1.3*s.allChoices(s.laser==1), 10, 'b', 'filled');
colors = cool(sampNum);
subplot(3,1,1)
[~, sortInd] = sort(params(:,5));
for j = 1:size(params,1)
    plot(1:length(choices), probChosen{sortInd(j)}, 'LineWidth', 1, 'Color', colors(j,:));
end
%% change of prediction by change of diff
params = allParams;
for i = 1:length(dayList)
    clear LH;
    clear probChosen;
    clear Q;
    clear pe;

    s = behAnalysisNoPlot_opMD(dayList{i}, 'simpleFlag', 1);
    choices = s.allChoices;
    choices(choices==-1) = 0;
    outcomes = abs(s.allRewards);
    randDiff = sort(0.3*(rand(20,1)-0.5));
    colors = cool(length(randDiff));
    for j = 1:3
        paramsTemp = params(j,:);
        figure; subplot(3,1,1); hold on
        plot([1:length(choices); 1:length(choices)], 0.5*[zeros(1, length(choices)); s.allChoices + s.allRewards], 'k');
        scatter(find(s.laser==1), 1.3*s.allChoices(s.laser==1), 10, 'b', 'filled');
        for m = 1:length(randDiff)
            paramsTemp(5) = randDiff(m);
            [LH, probChosen, Q, pe] = qLearningModel_5params_k_bias_LaserNegRPE(paramsTemp, choices, outcomes, s.laser);
            subplot(3,1,1); hold on;
            plot(1:length(choices), probChosen, 'LineWidth', 0.5, 'Color', colors(m,:));
        end
        sgtitle(dayList{i})
        [LH, probChosen, Q, pe] = qLearningModel_5params_k_bias_LaserNegRPE(params(j,:), choices, outcomes, s.laser);
        plot(1:length(choices), probChosen, 'LineWidth', 0.5, 'Color', 'k');
        screen = get(0,'Screensize');
        screen(4) = screen(4) - 100;
        set(gcf, 'Position', screen);
    end

end
%%



%% model comparison
modelNames = {'5params_k_bias', '5params_k_bias_LaserNegRPE', '5params_k_bias_LaserNegOnlyRPE'};
colors = hsv(length(modelNames));
aniNames = {'ZS082', 'ZS083', 'ZS085', 'ZS086'};
[root, sep] = currComputer();
col = 'cueOnGood';
allLL = cell(length(modelNames), 1);
nbins = 50;
for i = 1:length(aniNames)
    ani = aniNames{i};
    figure2; hold on;
    
    for m = 1:length(modelNames)
        modelCurr = modelNames{m};
        filename = [ani col '_' modelCurr];
        load([root ani sep ani 'sorted' sep 'stan\bernoulli' sep modelCurr sep col sep filename '.mat']);
        eval(['samps =' filename ';'])
        eval(['clear ' filename]);
        ll = sum([samps.log_lik],2);
        noDInd = samps.divergent__ ~= 1;
        ll = ll(noDInd);
        histogram(ll, 100, 'FaceColor', colors(m,:), 'EdgeColor', 'none', 'Normalization', 'probability');
        for j = 1:length(dayList)
            tmpLL = [samps.log_lik(:,j)];
            tmpLL = tmpLL(noDInd);
            [N,edges] = histcounts(tmpLL,nbins);
            [~, winInd] = max(N);
            if winInd == nbins
                sessLL = mean(tmpLL(tmpLL>=edges(winInd) & tmpLL<=edges(winInd+1)));
            else
                sessLL = mean(tmpLL(tmpLL>=edges(winInd) & tmpLL<edges(winInd+1)));
            end
            allLL{m} = [allLL{m}; sessLL];
        end
        clear samps; 
    end
    legend(modelNames, 'Interpreter', 'none')
    title(ani);
end
%% plot each model's prediction
modelNames = {'5params_k_bias', '5params_k_bias_LaserNegRPE', '5params_k_bias_LaserNegOnlyRPE'};
colors = hsv(length(modelNames));
xlFile = 'inhibitionGt';
sheet = 'allGt';
col = 'cueOnGood';
[root, sep] = currComputer();
dayList = getDayList(xlFile, sheet, col);
sampNum = 200;

for i = 1:length(dayList)
    session = dayList{i};
    pd = parseSessionString_df(session, root, sep);
    s = behAnalysisNoPlot_opMD(session, 'simpleFlag', 1);
    figure;
    clear pPred;
    for m = 1:length(modelNames)
        modelCurr = modelNames{m};
        params = getStanModelParams_sampsOnly(pd.aniName, col, modelCurr, sampNum, 'sessionName', dayList{i}, 'biasFlag',1);
        t = inferModelVar(dayList{i}, params, modelCurr, 'perturb', 1);
        pPred(m,:) = t.probChoice';
        subplot(3,1,1); hold on;
        if m == 2
            plot(1:length(s.allChoices), pPred(m, :), 'LineWidth', 2, 'Color', colors(m,:), 'LineStyle', '--');
        else
            plot(1:length(s.allChoices), pPred(m, :), 'LineWidth', 2, 'Color', colors(m,:));
        end
    end
    subplot(3,1,1); hold on;
    legend(modelNames, 'Interpreter', 'none')
    title(session);
    subplot(3,1,2); hold on;
    plot([1:length(s.allChoices); 1:length(s.allChoices)], 0.5*[zeros(1, length(s.allChoices)); s.allChoices + s.allRewards], 'k');
    scatter(find(s.laser==1), 1.1*s.allChoices(s.laser==1), 12, 'b', 'filled');
    subplot(3,1,3); hold on;
    for m = 1:length(modelNames)
        plot(1:length(s.allChoices), pPred(m,:)-pPred(1,:), 'Color', colors(m,:));
    end
    legend(modelNames, 'Interpreter', 'none')
    screen = get(0,'Screensize');
    screen(4) = screen(4) - 100;
    set(gcf, 'Position', screen)
end
%% find pPred around time of laser
% plot each model's prediction
modelNames = {'5params_k_bias', '5params_k_bias_LaserNegRPE', '5params_k_bias_LaserNegOnlyRPE'};
colors = hsv(length(modelNames));
xlFile = 'inhibitionGt';
sheet = 'allGt';
col = 'cueOnGood';
[root, sep] = currComputer();
dayList = getDayList(xlFile, sheet, col);
sampNum = 200;
tb = 3;
tf = 3;
laserMatCombine = cell(1,length(modelNames));
noLaserMatCombine = cell(1,length(modelNames));
for i = 1:length(dayList)
    session = dayList{i};
    pd = parseSessionString_df(session, root, sep);
    s = behAnalysisNoPlot_opMD(session, 'simpleFlag', 1);
    clear pPred;
    clear pPredAligned;
    clear noLaserMat;
    clear laserMat;
    for m = 1:length(modelNames)
        modelCurr = modelNames{m};
        params = getStanModelParams_sampsOnly(pd.aniName, col, modelCurr, sampNum, 'sessionName', dayList{i}, 'biasFlag',1);
        t = inferModelVar(dayList{i}, params, modelCurr, 'perturb', 1);
        pPred = t.probChoice';
        pPred = [NaN(1,tb) pPred NaN(1, tf)];
        pPredAlignedTmp = NaN(length(s.allChoices), tb+tf+1);
        for j = 1:(tb+tf+1)
            pPredAlignedTmp(:,j) = pPred(j:j+length(s.allChoices)-1)';
        end
        focusCase = mintersect(1:length(s.allChoices), s.nrwd_Inds, s.changeChoice_Inds-1);
        
        laserMat{m} = pPredAlignedTmp(mintersect(focusCase, find(s.laser==1)),:);
        noLaserMat{m} = pPredAlignedTmp(mintersect(focusCase, find(s.laser==0)),:);
        pPredAligned{m} = pPredAlignedTmp;
        laserMatCombine{m} = [laserMatCombine{m}; laserMat{m}];
        noLaserMatCombine{m} = [noLaserMatCombine{m}; noLaserMat{m}];
  
    end
    figure2Wide;
    subplot(1,3,1); hold on;
    for m = 1:length(modelNames)
      plotFilled(-tb:tf, laserMat{m}, colors(m,:));
    end
    emp = cell(size(modelNames));
    emp{i} = '';
    title('laser')
%     legends = [modelNames; emp];
%     legends = legends(1:2*length(modelNames));
%     legend(legends, 'Interpreter', 'none')
%     
    subplot(1,3,2); hold on;
    for m = 1:length(modelNames)
      plotFilled(-tb:tf, noLaserMat{m}, colors(m,:));
    end
    title('no laser')
    
    subplot(1,3,3); hold on;
    for m = 1:length(modelNames)
      plotFilled(-tb:tf, pPredAligned{m}, colors(m,:));
    end
    title('all')
    sgtitle(session);
    
%     legends = [modelNames; cell(size(modelNames))];
%     legends = legends(1:2*length(modelNames));
%     legend(legends, 'Interpreter', 'none')
%     screen = get(0,'Screensize');
%     screen(4) = screen(4) - 100;
%     set(gcf, 'Position', screen)
end
%% plot posteriors animal by animal or session by session
% animal by animal
col = 'cueOnShamLate';
color = 'k';
xlFile = 'inhibitionGt';
modelName = '5params_k_bias_LaserDisengageScale';
sheet = 'allGt';
dayList = getDayList(xlFile, sheet, col);
aniList = cellfun(@(x) x(2:6), dayList, 'UniformOutput', false);
aniList =  unique(aniList);
[root, sep] = currComputer();
%%
figure2; hold on;

for i = 1:length(aniList)
    
    ani = aniList{i};
    filename = [ani col '_' modelName];
    load([root ani sep ani 'sorted' sep 'stan\bernoulli' sep modelName sep col sep filename '.mat']);
    eval(['samps =' filename ';'])
    eval(['clear ' filename]);
%     histogram(samps.mu_diff, 50, 'FaceColor', color, 'EdgeColor', 'none', 'FaceAlpha', 0.25);
    noDInd = samps.divergent__ ~= 1;
    [N, edges] = histcounts(samps.mu_scale(noDInd), 50, 'Normalization', 'probability');
    hPos = 0.5*edges(1:end-1) + 0.5*edges(2:end);
    vPos = 0.25*(i-1);
    patch([hPos, flip(hPos)], [vPos+N, flip(vPos*ones(size(N)))], color, 'edgeColor', 'none', 'FaceAlpha',0.25);
    p5__ = summary{{'mu_scale'}, {'p5_'}};
    p95__ = summary{{'mu_scale'}, {'p95_'}};
    line([p5__ p5__], [vPos vPos+0.05], 'LineStyle', '--', 'Color', 'k','LineWidth',2)
    line([p95__ p95__], [vPos vPos+0.05], 'LineStyle', '--', 'Color', 'k','LineWidth',2)
end

plot([1 1], [0 0.3*(i-1) + 0.15], 'LineStyle', '--', 'Color', 'r','LineWidth',2)
%%


















%%
