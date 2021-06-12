function lickCatAnalysis(session, model, col, varargin)
p = inputParser;
p.addParameter('paramNames', {'aN', 'aP', 'aF', 'beta', 'bias'});
p.parse(varargin{:})

[root, sep] = currComputer();
[animalName, date] = strtok(session, 'd'); 
animalName = animalName(2:end);
date = date(1:9);
s = behAnalysisNoPlot_opMD(session);
%% load model fitting results and calculate DVs
sampFile = [animalName col '_', model];
path = [root animalName sep animalName 'sorted' sep 'stan' sep 'bernoulli' sep model sep col sep];
load([path sampFile '.mat'], 'dayList');
samples = load([path sampFile '.mat'], sampFile);
samples = samples.(sampFile);
paramNames = p.Results.paramNames;
% model 
id = find(strcmp(dayList,session),1);
%generate best estimates of parameters
paramEsts = [];
allSamples = [];
edges = cell(1,length(paramNames));
for i = 1:length(paramNames)
    tmp = samples.(paramNames{i})(:,id);
    allSamples = [allSamples tmp];
    edges{i} = linspace(min(tmp), max(tmp),50);
end
n = histcnd(allSamples,edges); %bin samples by multiple dimensions
[~, inds] = myMaxAll(n); %find the bin with max num in bin
for i = 1:length(paramNames) %use median in bin as best estimate
    tmp = allSamples(:,i);
    edgeTmp = edges{i};
    if inds(i) < 50
        paramEsts(i) = median(tmp(tmp >= edgeTmp(inds(i)) & tmp < edgeTmp(inds(i)+1)));
    else
        paramEsts(i) = edgeTmp(inds(i));
    end
end
% decide if input includes time forget
if contains(model,'tF')
    input = 'choice, outcome, ITI)';
else
    input = 'choice, outcome)';
end
choice = s.allChoices';
choice(choice<0) = 0;
outcome = abs(s.allRewards);
% calculate model varaibles
eval(['[LL,probC,Q,pe] = qLearningModel_' model '(paramEsts,' input ';'])
% diff value
Qdiff = abs(Q(:,2)-Q(:,1));
% total value
Qsum = sum(Q,2);
% prepe
prePe = [NaN pe(1:end-1)];
% choice confidence
choiceConf = 2.*probC - 1;


% kmeans c = 2;
ind = kmeans(log(s.lickLat'),2,'Start', [4.5; 5.5], 'OnlinePhase', 'on');
ind(ind == 3) = 2;
if mean(s.lickLat(ind==1))>mean(s.lickLat(ind==2))
    ind = 3-ind;
end
svs = [NaN; - ones(length(ind)-1,1)];
hmm = zeros(length(ind),1);
svs(s.changeChoice_Inds) = 1;
hmm(s.hmmStates == 1) = 1;
preRwd = abs(s.allRewards);
preRwd = [NaN; preRwd(1:end-1)'];
preQsum = Qsum;
preQsum = [NaN; preQsum(1:end-1)];
% glm with interactions

tbl = table(abs(prePe)', zscore(Qsum), choiceConf, s.timeBtwn', svs, ind-1, 'VariableNames', {'absprePe', 'Qsum', 'confidence', 'preITI', 'svs', 'lickLat'});
%tbl = table(choiceConf, svs, ind-1, 'VariableNames', {'confidence','switch', 'lickLat'});


%tbl = table(zscore(Qsum), preRwd, 0.5*(-sign(choiceConf)+1), s.timeBtwn', svs, zscore(s.lickLat'), 'VariableNames', {'Qsum', 'preRwd', 'confidence', 'preITI', 'switch', 'lickLat'});
%tbl = table(preRwd, zscore(choiceConf), s.timeBtwn', svs, zscore(s.lickLat'), 'VariableNames', {'preRwd', 'confidence', 'preITI', 'switch', 'lickLat'});
%tbl = table(zscore(Qsum), zscore(choiceConf), s.timeBtwn', svs, zscore(s.lickLat'), 'VariableNames', {'Qsum', 'confidence', 'preITI', 'switch', 'lickLat'});
%tbl = table(zscore(choiceConf), svs, zscore(s.lickLat'), 'VariableNames', {'confidence','switch', 'lickLat'});

% mdl = stepwiselm(tbl,'interactions');

mdl = stepwiseglm(tbl,'interactions','Distribution','binomial', 'Link','logit');
%% plot everything
% check clustering
figure;
screenSize = get(0,'Screensize');
screenSize(4) = screenSize(4) - 100;
set(gcf, 'Position', screenSize)
suptitle(session)
subplot(4,3,1:3); hold on;
% trial plot
    rMag = 1;
    nrMag = rMag/2;
for i = 1:length(s.lickLat)
    if ind(i) == 2
       fill([i-0.5, i+0.5, i+0.5, i-0.5],[-1*rMag, -1*rMag, rMag, rMag],[1 0 0],'LineStyle','none');
    else
       fill([i-0.5, i+0.5, i+0.5, i-0.5],[-1*rMag, -1*rMag, rMag, rMag],[0.7 1 1],'LineStyle','none');
    end
end

for i = 1:length(ind)
    if s.allChoices(i) == 1
        if s.allRewards(i) == 1 % R side rewarded
            plot([i i],[0 rMag],'k')
        else
            plot([i i],[0 nrMag],'k') % R side not rewarded
        end
    else
        if s.allRewards(i) == -1 % L side rewarded
            plot([i i],[-1*rMag 0],'k')
        else
            plot([i i],[-1*nrMag 0],'k')
        end
    end
end

text(0,1.5,'L/R');
xlim([0.5 length(ind)+0.5]);
ylabel('<-- L       R  -->')

subplot(4,3,4); hold on;
histogram(s.lickLat(ind == 1),0:20:800,'FaceColor', 'c');
histogram(s.lickLat(ind == 2),0:20:800,'FaceColor', 'r');
legend({'cluster 1', 'cluster2'})
xlabel('lickLat /ms')

subplot(4,3,5);hold on;
histogram(Qsum(ind == 1),linspace(min(Qsum), max(Qsum), 15),'Normalization','probability', 'FaceColor', 'c');
histogram(Qsum(ind == 2),linspace(min(Qsum), max(Qsum), 15),'Normalization','probability', 'FaceColor', 'r');
xlabel('Qsum')

subplot(4,3,6);hold on;
histogram(Qdiff(ind == 1),linspace(min(Qdiff), max(Qdiff), 15),'Normalization','probability', 'FaceColor', 'c');
histogram(Qdiff(ind == 2),linspace(min(Qdiff), max(Qdiff), 15),'Normalization','probability', 'FaceColor', 'r');
xlabel('Qdiff')


subplot(4,3,7);hold on;
histogram(choiceConf(ind == 1),linspace(min(choiceConf), max(choiceConf), 15),'Normalization','probability', 'FaceColor', 'c');
histogram(choiceConf(ind == 2),linspace(min(choiceConf), max(choiceConf), 15),'Normalization','probability', 'FaceColor', 'r');
xlabel('confidence')

subplot(4,3,8);hold on;
histogram(s.lickLat(s.stayChoice_Inds),0:20:800,'Normalization','probability', 'FaceColor', 'c');
histogram(s.lickLat(s.changeChoice_Inds),0:20:800,'Normalization','probability', 'FaceColor', 'r');
legend({'stay','switch' })
xlabel('switch/stay')

subplot(4,3,9);hold on;
histogram(s.lickLat(s.hmmStates~=1),0:20:800,'Normalization','probability', 'FaceColor', 'c');
histogram(s.lickLat(s.hmmStates==1),0:20:800,'Normalization','probability', 'FaceColor', 'r');
legend({'exploit','explore'})
xlabel('hmm')

subplot(4,3,10);hold on;
histogram(prePe(ind == 1),linspace(min(prePe), max(prePe), 15),'Normalization','probability', 'FaceColor', 'c');
histogram(prePe(ind == 2),linspace(min(prePe), max(prePe), 15),'Normalization','probability', 'FaceColor', 'r');
xlabel('prePE')

subplot(4,3,11);hold on;
currInds = intersect(s.stayChoice_Inds, find(choiceConf<0.1));
currIndsInv = intersect(s.changeChoice_Inds, find(choiceConf>0.3));
histogram(s.lickLat,0:20:800,'FaceColor', [0.7,0.7,0.7])
histogram(s.lickLat(currInds),0:20:800,'Normalization','probability', 'FaceColor', 'c');
histogram(s.lickLat(currIndsInv),0:20:800,'Normalization','probability', 'FaceColor', 'm')
legend({'','stay in unpreferred','change to preferred'})

subplot(4,3,12); hold on

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
    text(i-0.4,1.2*(sign(coefVals(i))*height(i+1)), mdl.CoefficientNames{i+1})
end
title('glm: on lickGroup')
ylabel('\beta Coefficient')


end
