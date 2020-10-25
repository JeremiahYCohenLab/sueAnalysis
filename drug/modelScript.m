
%% compare 
% goal:
% fitting of models in normal behavior (find a better one)
% fitting of models in drug vs saline behavior (should be comparable)
stan_qLearningFit('combineRwdDelay', 'ZS040RwdDelay', 'good','modelName','fourParam','paramNames', {'aN', 'aP', 'aF', 'beta'});
stan_qLearningFit('combineRwdDelay', 'ZS040RwdDelay', 'good','modelName','fourParam_tF','paramNames', {'aN', 'aP', 'tF', 'beta'});
stan_qLearningFit('combineRwdDelay', 'ZS040RwdDelay', 'good','modelName','fourParam_tF_oneside','paramNames', {'aN', 'aP', 'tF', 'beta'});

stan_qLearningFit('combineRwdDelay', 'ZS041RwdDelay', 'good','modelName','fourParam','paramNames', {'aN', 'aP', 'aF', 'beta'});
stan_qLearningFit('combineRwdDelay', 'ZS041RwdDelay', 'good','modelName','fourParam_tF','paramNames', {'aN', 'aP', 'tF', 'beta'});
stan_qLearningFit('combineRwdDelay', 'ZS041RwdDelay', 'good','modelName','fourParam_tF_oneside','paramNames', {'aN', 'aP', 'tF', 'beta'});

stan_qLearningFit('combineDrug', 'combineDrugSaline', 'saline','modelName','fourParam','paramNames', {'aN', 'aP', 'aF', 'beta'});
stan_qLearningFit('combineDrug', 'combineDrugSaline', 'saline','modelName','fourParam_tF','paramNames', {'aN', 'aP', 'tF', 'beta'});
stan_qLearningFit('combineDrug', 'combineDrugSaline', 'saline','modelName','fourParam_tF_oneside','paramNames', {'aN', 'aP', 'tF', 'beta'});

stan_qLearningFit('combineDrug', 'combineDrugSaline', 'salineLate','modelName','fourParam','paramNames', {'aN', 'aP', 'aF', 'beta'});
stan_qLearningFit('combineDrug', 'combineDrugSaline', 'salineLate','modelName','fourParam_tF','paramNames', {'aN', 'aP', 'tF', 'beta'});
stan_qLearningFit('combineDrug', 'combineDrugSaline', 'salineLate','modelName','fourParam_tF_oneside','paramNames', {'aN', 'aP', 'tF', 'beta'});

stan_qLearningFit('combineDrug', 'combineDrugdrug', 'drug','modelName','fourParam','paramNames', {'aN', 'aP', 'aF', 'beta'});
stan_qLearningFit('combineDrug', 'combineDrugdrug', 'drug','modelName','fourParam_tF','paramNames', {'aN', 'aP', 'tF', 'beta'});
stan_qLearningFit('combineDrug', 'combineDrugdrug', 'drug','modelName','fourParam_tF_oneside','paramNames', {'aN', 'aP', 'tF', 'beta'});

stan_qLearningFit('combineDrug', 'combineDrug', 'all','modelName','fourParam','paramNames', {'aN', 'aP', 'aF', 'beta'});
stan_qLearningFit('combineDrug', 'combineDrug', 'all','modelName','fourParam_tF_oneside','paramNames', {'aN', 'aP', 'tF', 'beta'});

stan_qLearningFit('combineDrug', 'combineDrugHabituation', 'habituation','modelName','4Params','paramNames', {'aN', 'aP', 'aF', 'beta'});
stan_qLearningFit('combineDrug', 'combineDrugHabituation', 'habituation','modelName','4Params_tF_oneside','paramNames', {'aN', 'aP', 'tF', 'beta'});

%% fixed params
load('C:\Users\zhixi\Documents\data\combineDrug\combineDrugsorted\stan\bernoulli\fourParam\combineDrugall_4params.mat')
load('C:\Users\zhixi\Documents\data\combineDrug\combineDrugsorted\stan\bernoulli\fourParam_tF_oneside\combineDrugall_fourParam_tF_oneside.mat')

samples = combineDrugall_fourParam_tF_oneside;
params = zeros(0.5*length(dayList),length(paramNames));
for i = 1:0.5*length(dayList)
    for j = 1:length(paramNames)                                                          
            tmp = samples.(paramNames{j})(:,i+6);
            [n,e] = histcounts(tmp, 50);
            [~, maxInd] = max(n);
            params(i,j) = median(tmp(tmp > e(maxInd) & tmp < e(maxInd+1)));
    end
end
%%
stan_qLearningFit('combineDrug', 'combineDrugSaline', 'salineLate','modelName','4params_tF_oneside','paramNames', {'aN', 'aP', 'tF', 'beta'}, 'nonfixedParams', 'aN', 'fixedParams', paramEsts);
stan_qLearningFit('combineDrug', 'combineDrugdrug', 'drug','modelName','4params_tF_oneside','paramNames', {'aN', 'aP', 'tF', 'beta'}, 'nonfixedParams', 'aN', 'fixedParams', paramEsts);


%% compare model
% goal:
% fitting of models in normal behavior (find a better one)
% Four compare ways:
criterion_ZS40 = computeDIC('ZS040RwdDelay','good','fourParam','4params',{'aN', 'aP', 'aF', 'beta'});
criterion_ZS40_tF = computeDIC('ZS040RwdDelay','good','fourParam_tF','4params_tF',{'aN', 'aP', 'tF', 'beta'});
criterion_ZS40_tF_one = computeDIC('ZS040RwdDelay','good','fourParam_tF_oneside','4params_tF_oneside',{'aN', 'aP', 'tF', 'beta'});

criterion_ZS041 = computeDIC('ZS041RwdDelay','good','fourParam','4params',{'aN', 'aP', 'aF', 'beta'});
criterion_ZS041_tF = computeDIC('ZS041RwdDelay','good','fourParam_tF','4params_tF',{'aN', 'aP', 'tF', 'beta'});
criterion_ZS041_tF_one = computeDIC('ZS041RwdDelay','good','fourParam_tF_oneside','4params_tF_oneside',{'aN', 'aP', 'tF', 'beta'});
%% compare in drug and saline
% for criterion, still use 'fourparams'
criterion_drug = computeDIC('combineDrugdrug', 'drug','fourParam','4params',{'aN', 'aP', 'aF', 'beta'});
criterion_drug_tF = computeDIC('combineDrugdrug', 'drug','fourParam_tF','4params_tF',{'aN', 'aP', 'tF', 'beta'});
criterion_drug_tF_one = computeDIC('combineDrugdrug', 'drug','fourParam_tF_oneside','4params_tF_oneside',{'aN', 'aP', 'tF', 'beta'});

criterion_saline = computeDIC('combineDrugSaline', 'saline','fourParam','4params',{'aN', 'aP', 'aF', 'beta'});
criterion_saline_tF = computeDIC('combineDrugSaline', 'saline','fourParam_tF','4params_tF',{'aN', 'aP', 'tF', 'beta'});
criterion_saline_tF_one = computeDIC('combineDrugSaline', 'saline','fourParam_tF_oneside','4params_tF_oneside',{'aN', 'aP', 'tF', 'beta'});

criterion_saline = computeDIC('combineDrugSaline', 'salineLate','fourParam','4params',{'aN', 'aP', 'aF', 'beta'});
criterion_saline_tF = computeDIC('combineDrugSaline', 'salineLate','fourParam_tF','4params_tF',{'aN', 'aP', 'tF', 'beta'});
criterion_saline_tF_one = computeDIC('combineDrugSaline', 'salineLate','fourParam_tF_oneside','4params_tF_oneside',{'aN', 'aP', 'tF', 'beta'});


% get each session's ll for three models.
% ll: all 80000 samples, histogram them for all each animal and all
% animals.
%% color
blue = [0 1 1];
purp = [0.7 0 1];
num = 3;
colors = [linspace(blue(1),purp(1),num)', linspace(blue(2),purp(2),num)', linspace(blue(3),purp(3),num)'];
r = [0 0 0];
b = [0 0 1];
map = [linspace(b(1),r(1), 100)', linspace(b(2),r(2), 100)',linspace(b(3),r(3), 100)'];
%% ZS040 all ll
figure2; hold on;    
histogram(reshape(criterion_ZS40.ll,[],1),100,'Normalization', 'Probability', 'FaceColor', colors(1,:), 'EdgeColor', 'none');
histogram(reshape(criterion_ZS40_tF.ll,[],1),100,'Normalization', 'Probability', 'FaceColor', colors(2,:), 'EdgeColor', 'none');
histogram(reshape(criterion_ZS40_tF_one.ll,[],1),100,'Normalization', 'Probability', 'FaceColor', colors(3,:), 'EdgeColor', 'none');
legend('aF','tF','tF-one');
title('ZS040-rwdDelay')
%% ZS040 sum ll
figure2; hold on;    
histogram(sum(criterion_ZS40.ll,2),100,'Normalization', 'Probability', 'FaceColor', colors(1,:), 'EdgeColor', 'none');
histogram(sum(criterion_ZS40_tF.ll,2),100,'Normalization', 'Probability', 'FaceColor', colors(2,:), 'EdgeColor', 'none');
histogram(sum(criterion_ZS40_tF_one.ll,2),100,'Normalization', 'Probability', 'FaceColor', colors(3,:), 'EdgeColor', 'none');
legend('aF','tF','tF-one');
title('ZS040-rwdDelay')

%% ZS040 correlation
figure2('position',[0 0 900 300]);

subplot(1,3,1)
[rho, p] = corr(criterion_ZS40.ll);
rho = rho - diag(rho);
imagesc(rho);
title('aF')
subplot(1,3,2)
[rho, p] = corr(criterion_ZS40_tF.ll);
rho = rho - diag(rho);
imagesc(rho);
title('tF')
subplot(1,3,3)
[rho, p] = corr(criterion_ZS40_tF_one.ll);
rho = rho - diag(rho);
imagesc(rho);
title('tF-one')
suptitle('ZS040-rwdDelay')

%% ZS040 dic
figure2('position', [0 0 900 600]);

blue = [0 1 1];
purp = [1 0 1];
num = length(criterion_ZS40.dic);
colors = [linspace(blue(1),purp(1),num)', linspace(blue(2),purp(2),num)', linspace(blue(3),purp(3),num)'];

subplot(2,3,1); hold on;
scatter(criterion_ZS40.dic, criterion_ZS40_tF.dic, 10, colors,'filled');
line([0 300],[0 300], 'Color', [1 0 0], 'LineStyle','--');
xlabel('aF');
ylabel('tF');

subplot(2,3,4);hold on;
histogram(criterion_ZS40.dic-criterion_ZS40_tF.dic,20);
line([0 0],[0 8], 'Color', [1 0 0], 'LineStyle','--');
title('aF-tF')

subplot(2,3,2); hold on;
scatter(criterion_ZS40.dic, criterion_ZS40_tF_one.dic, 10, colors,'filled');
line([0 300],[0 300], 'Color', [1 0 0], 'LineStyle','--');
xlabel('aF');
ylabel('tF-one');

subplot(2,3,5); hold on;
histogram(criterion_ZS40.dic-criterion_ZS40_tF_one.dic,20);
line([0 0],[0 8], 'Color', [1 0 0], 'LineStyle','--');
title('aF-tF-one')

subplot(2,3,3); hold on;
scatter(criterion_ZS40_tF_one.dic, criterion_ZS40_tF.dic, 10, colors,'filled');
line([0 300],[0 300], 'Color', [1 0 0], 'LineStyle','--');
xlabel('tF-one');
ylabel('tF');

subplot(2,3,6);hold on;
histogram(criterion_ZS40_tF_one.dic-criterion_ZS40_tF.dic,20);
line([0 0],[0 8], 'Color', [1 0 0], 'LineStyle','--');
title('tF-one-tF')
suptitle('ZS040-DIC')
%% ZS040 max
figure2('position', [0 0 900 600]);
subplot(2,3,1); hold on;
scatter(criterion_ZS40.llmax, criterion_ZS40_tF.llmax, 10, colors,'filled');
line([-100 0],[-100 0], 'Color', [1 0 0], 'LineStyle','--');
xlabel('aF');
ylabel('tF');

subplot(2,3,4);hold on;
histogram(criterion_ZS40.llmax-criterion_ZS40_tF.llmax,20);
line([0 0],[0 8], 'Color', [1 0 0], 'LineStyle','--');
title('aF-tF')

subplot(2,3,2); hold on;
scatter(criterion_ZS40.llmax, criterion_ZS40_tF_one.llmax, 10, colors,'filled');
line([-100 0],[-100 0], 'Color', [1 0 0], 'LineStyle','--');
xlabel('aF');
ylabel('tF-one');

subplot(2,3,5); hold on;
histogram(criterion_ZS40.llmax-criterion_ZS40_tF_one.llmax,20);
line([0 0],[0 8], 'Color', [1 0 0], 'LineStyle','--');
title('aF-tF-one')

subplot(2,3,3); hold on;
scatter(criterion_ZS40_tF_one.llmax, criterion_ZS40_tF.llmax, 10, colors,'filled');
line([-100 0],[-100 0], 'Color', [1 0 0], 'LineStyle','--');
xlabel('tF-one');
ylabel('tF');

subplot(2,3,6);hold on;
histogram(criterion_ZS40_tF_one.llmax-criterion_ZS40_tF.llmax,20);
line([0 0],[0 8], 'Color', [1 0 0], 'LineStyle','--');
title('tF-one-tF')
suptitle('ZS040-max')
%% ZS040 ani
figure2('position', [0 0 900 600]);
subplot(2,3,1); hold on;
scatter(criterion_ZS40.llani, criterion_ZS40_tF.llani, 10, colors,'filled');
line([-100 0],[-100 0], 'Color', [1 0 0], 'LineStyle','--');
xlabel('aF');
ylabel('tF');

subplot(2,3,4);hold on;
histogram(criterion_ZS40.llani-criterion_ZS40_tF.llani,20);
line([0 0],[0 8], 'Color', [1 0 0], 'LineStyle','--');
title('aF-tF')

subplot(2,3,2); hold on;
scatter(criterion_ZS40.llani, criterion_ZS40_tF_one.llani, 10, colors,'filled');
line([-100 0],[-100 0], 'Color', [1 0 0], 'LineStyle','--');
xlabel('aF');
ylabel('tF-one');

subplot(2,3,5); hold on;
histogram(criterion_ZS40.llani-criterion_ZS40_tF_one.llani,20);
line([0 0],[0 8], 'Color', [1 0 0], 'LineStyle','--');
title('aF-tF-one')

subplot(2,3,3); hold on;
scatter(criterion_ZS40_tF_one.llani, criterion_ZS40_tF.llani, 10, colors,'filled');
line([-100 0],[-100 0], 'Color', [1 0 0], 'LineStyle','--');
xlabel('tF-one');
ylabel('tF');

subplot(2,3,6);hold on;
histogram(criterion_ZS40_tF_one.llani-criterion_ZS40_tF.llani,20);
line([0 0],[0 8], 'Color', [1 0 0], 'LineStyle','--');
title('tF-one-tF')
suptitle('ZS040-ani')

%% ZS040 sesh
figure2('position', [0 0 900 300]);
subplot(2,3,1); hold on;
scatter(criterion_ZS40.llsesh, criterion_ZS40_tF.llsesh, 10, colors,'filled');
line([-100 0],[-100 0], 'Color', [1 0 0], 'LineStyle','--');
xlabel('aF');
ylabel('tF');

subplot(2,3,4);hold on;
histogram(criterion_ZS40.llsesh-criterion_ZS40_tF.llsesh,20);
line([0 0],[0 8], 'Color', [1 0 0], 'LineStyle','--');
title('aF-tF')

subplot(2,3,2); hold on;
scatter(criterion_ZS40.llsesh, criterion_ZS40_tF_one.llsesh, 10, colors,'filled');
line([-100 0],[-100 0], 'Color', [1 0 0], 'LineStyle','--');
xlabel('aF');
ylabel('tF-one');

subplot(2,3,5); hold on;
histogram(criterion_ZS40.llsesh-criterion_ZS40_tF_one.llsesh,20);
line([0 0],[0 8], 'Color', [1 0 0], 'LineStyle','--');
title('aF-tF-one')

subplot(2,3,3); hold on;
scatter(criterion_ZS40_tF_one.llsesh, criterion_ZS40_tF.llsesh, 10, colors,'filled');
line([-100 0],[-100 0], 'Color', [1 0 0], 'LineStyle','--');
xlabel('tF-one');
ylabel('tF');

subplot(2,3,6);hold on;
histogram(criterion_ZS40_tF_one.llsesh-criterion_ZS40_tF.llsesh,20);
line([0 0],[0 8], 'Color', [1 0 0], 'LineStyle','--');
title('tF-one-tF')
suptitle('ZS040-sesh')


% find best fit and compare:
% 1. best fit by ll itself;
% 2. best fit by best parameters for each session;
% 3. best fit by best parameters for each animal/group;
% 4. dic




%% fitting of models in drug vs saline behavior (should be comparable)


% Four compare ways:

% get each session's ll for three models.
% ll: all 80000 samples, histogram them for all each animal and all
% animals.

% find best fit and compare:
% 1. best fit by ll itself;
% 2. best fit by best parameters for each session;
% 3. best fit by best parameters for each animal/group;

%% color
blue = [0 1 1];
purp = [0.7 0 1];
num = 2;
colors = [linspace(blue(1),purp(1),num)', linspace(blue(2),purp(2),num)', linspace(blue(3),purp(3),num)'];
r = [0 0 0];
b = [0 0 1];
map = [linspace(b(1),r(1), 100)', linspace(b(2),r(2), 100)',linspace(b(3),r(3), 100)'];
%% all lls
figure2; 
subplot(1,2,1);
hold on;    
histogram(reshape(criterion_drug.ll,[],1),100,'Normalization', 'Probability', 'FaceColor', colors(1,:), 'EdgeColor', 'none');
histogram(reshape(criterion_saline.ll,[],1),100,'Normalization', 'Probability', 'FaceColor', colors(2,:), 'EdgeColor', 'none');
legend('drug','saline');
title('aF')

subplot(1,2,2);
hold on;    
histogram(reshape(criterion_drug_tF_one.ll,[],1),100,'Normalization', 'Probability', 'FaceColor', colors(1,:), 'EdgeColor', 'none');
histogram(reshape(criterion_saline_tF_one.ll,[],1),100,'Normalization', 'Probability', 'FaceColor', colors(2,:), 'EdgeColor', 'none');
legend('drug','saline');
title('tF-one')
%% all mean lls
figure2; 
subplot(1,2,1); hold on;    
histogram(mean(criterion_drug.ll,2),100,'Normalization', 'Probability', 'FaceColor', colors(1,:), 'EdgeColor', 'none');
histogram(mean(criterion_saline.ll,2),100,'Normalization', 'Probability', 'FaceColor', colors(2,:), 'EdgeColor', 'none');
legend('drug','saline');
title('aF');

subplot(1,2,2); hold on;
histogram(mean(criterion_drug_tF_one.ll,2),100,'Normalization', 'Probability', 'FaceColor', colors(1,:), 'EdgeColor', 'none');
histogram(mean(criterion_saline_tF_one.ll,2),100,'Normalization', 'Probability', 'FaceColor', colors(2,:), 'EdgeColor', 'none');
legend('drug','saline');
title('tF-one');

%% best lls
figure2; 
subplot(1,2,1); hold on;    
histogram(criterion_drug.llmax,10,'FaceColor', colors(1,:), 'EdgeColor', 'none');
histogram(criterion_saline.llmax,10,'FaceColor', colors(2,:), 'EdgeColor', 'none');
legend('drug','saline');
title('aF-llmax');

subplot(1,2,2); hold on;
histogram(criterion_drug_tF_one.llmax,10,'FaceColor', colors(1,:), 'EdgeColor', 'none');
histogram(criterion_saline_tF_one.llmax,10,'FaceColor', colors(2,:), 'EdgeColor', 'none');
legend('drug','saline');
title('tF-one-llmax');
%% dic
figure2; 
subplot(1,2,1); hold on;    
histogram(criterion_drug.dic,10,'FaceColor', colors(1,:), 'EdgeColor', 'none');
histogram(criterion_saline.dic,10,'FaceColor', colors(2,:), 'EdgeColor', 'none');
legend('drug','saline');
title('aF-dic');

subplot(1,2,2); hold on;
histogram(criterion_drug_tF_one.dic,10,'FaceColor', colors(1,:), 'EdgeColor', 'none');
histogram(criterion_saline_tF_one.dic,10,'FaceColor', colors(2,:), 'EdgeColor', 'none');
legend('drug','saline');
title('tF-one-dic');
%% ani lls
figure2; 
subplot(1,2,1); hold on;    
histogram(criterion_drug.llani,10,'FaceColor', colors(1,:), 'EdgeColor', 'none');
histogram(criterion_saline.llani,10,'FaceColor', colors(2,:), 'EdgeColor', 'none');
legend('drug','saline');
title('aF-llani');

subplot(1,2,2); hold on;
histogram(criterion_drug_tF_one.llani,10,'FaceColor', colors(1,:), 'EdgeColor', 'none');
histogram(criterion_saline_tF_one.llani,10,'FaceColor', colors(2,:), 'EdgeColor', 'none');
legend('drug','saline');
title('tF-one-llani');

%% sesh lls
figure2; 
subplot(1,2,1); hold on;    
histogram(criterion_drug.llsesh,10, 'FaceColor', colors(1,:), 'EdgeColor', 'none');
histogram(criterion_saline.llsesh,10, 'FaceColor', colors(2,:), 'EdgeColor', 'none');
legend('drug','saline');
title('aF-llsesh');

subplot(1,2,2); hold on;
histogram(criterion_drug_tF_one.llsesh,10, 'FaceColor', colors(1,:), 'EdgeColor', 'none');
histogram(criterion_saline_tF_one.llsesh,10, 'FaceColor', colors(2,:), 'EdgeColor', 'none');
legend('drug','saline');
title('tF-one-llsesh');

%% drug analysis
% repeat what has been done with kernel with prefered model
% need to point to the model output construct
% log(lickLat)
%% compare parameters
samples = combineDrugdrugdrug_fourParam_tF_oneside;
paramNames = {'aN','aP','tF','beta'};
numParams = length(paramNames);
blue = [0.3 1 0.7];
purp = [1 0 0.7];
colors = [linspace(blue(1),purp(1),numParams)', linspace(blue(2),purp(2),numParams)', linspace(blue(3),purp(3),numParams)'];
figure2;
for i = 1:numParams
    subplot(1,numParams,i); hold on;
    histogram(eval(['samples.mu_' paramNames{i}]) , 100,...
        'Normalization', 'Probability', 'FaceColor', colors(i,:), 'EdgeColor', 'none')
    set(gca,'tickdir', 'out') 
    title(paramNames{i})
end
%%



