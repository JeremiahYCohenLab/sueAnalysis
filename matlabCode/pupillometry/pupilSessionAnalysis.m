function pupilSessionAnalysis(session, category, varargin)
%task and model parameters
p = inputParser;
% default parameters if none given
p.addParameter('plotFlag', 1);
p.addParameter('maxTrial', 1000);
p.addParameter('modelNameOld', '5params');
p.addParameter('modelName', '5params');
p.addParameter('paramNames', {'aN', 'aP', 'aF', 'beta', 'bias'});
p.addParameter('binSizeFrame', 6)% need to be even number
p.addParameter('stepSize', 250)
p.addParameter('startTime', -2000)
p.addParameter('endTime',6000)
p.addParameter('saveFigFlag', 1)
p.parse(varargin{:});

paramNames = p.Results.paramNames;
maxTrial = p.Results.maxTrial;
time = p.Results.startTime:p.Results.stepSize:p.Results.endTime;
numBins = length(time);


if contains(p.Results.modelName,'tF')
    input = 'choice, outcome, ITI)';
else
    input = 'choice, outcome)';
end

% load model fitting results 
[animalName, date] = strtok(session, 'd'); 
animalName = animalName(2:end);
date = date(1:9);
[root, sep] = currComputer();
sampFile = [animalName category '_', p.Results.modelNameOld];
path = [root animalName sep animalName 'sorted' sep 'stan' sep 'bernoulli' sep p.Results.modelNameOld sep category sep ];
load([path sampFile '.mat'], 'dayList');
samples = load([path sampFile '.mat'], sampFile);
id = contains(dayList, session);
samples = samples.(sampFile);

%% load files, first check if there's pupil data
sessionFolder = ['m' animalName date];  
% paths
videopath = [root animalName sep sessionFolder sep 'pupil'];  
if isstrprop(session(end), 'alpha')
    behSessionDataPath = [root animalName sep sessionFolder sep 'sorted' sep 'session ' session(end) sep session '_sessionData_behav.mat'];
    pupilAlignPath = [root animalName sep sessionFolder sep 'sorted' sep 'session ' session(end) sep session '_pupil.mat'];
    pupilPath = [root animalName sep sessionFolder sep 'pupil'];
else
    behSessionDataPath = [root animalName sep sessionFolder sep 'sorted' sep 'session' sep session '_sessionData_behav.mat'];
    pupilAlignPath = [root animalName sep sessionFolder sep 'sorted' sep 'session' sep session '_pupil.mat'];
    pupilPath = [root animalName sep sessionFolder sep 'pupil'];
    savePath = [root animalName sep sessionFolder sep 'figures' sep];
end
% behavior
if exist(behSessionDataPath,'file')
    load(behSessionDataPath)
else
    sessionData = generateSessionData_operantMatchingDecoupledRwdDelay(session);
end
% pupil
if ~exist(pupilPath, 'dir')
   fprintf([session ' no pupil data'])
    return
end
if exist(pupilAlignPath, 'file')
    load(pupilAlignPath);
else
    fprintf([session ' pupil not aligned yet'])
    return
end
% jump the current loop if not well aligned
if errorProp > 0.2 || isnan(errorProp)
    fprintf([session ' pupil not well aligned'])
    return
end
%% preparation 
% parse behavior
os = behAnalysisNoPlot_opMD(session);
choice = os.allChoices';
choice(choice<0) = 0;
outcome = abs(os.allRewards)';
choice = choice(1:min(length(choice), maxTrial));
outcome = outcome(1:length(choice));
responseInds = os.responseInds(1:min(length(choice), maxTrial)); 
ITI = os.timeBtwn(1:length(choice));

% load diameter
list = dir(videopath);
expression = ['^' session 'DLC' '\w*' 'shuffle1_100000.csv' '$'];
expressionSkeleton = ['^' session 'DLC' '\w*' 'shuffle1_100000_skeleton.csv' '$'];
position = list(~cellfun(@isempty, cellfun(@(x) regexp(x, expression), {list.name}, 'UniformOutput', false))).name;
skeleton = list(~cellfun(@isempty, cellfun(@(x) regexp(x, expressionSkeleton), {list.name}, 'UniformOutput', false))).name;
positionRaw = csvread([videopath sep position], 3, 0);
diaRaw = csvread([videopath sep skeleton], 2, 0);
qualF = positionRaw(:,4) > 0.99 &  positionRaw(:,7) > 0.99;   % quality control for good frames with good diameter estimation. 
dia = interp1(find(qualF>0),diaRaw(qualF,2),1:length(diaRaw(:,2)));

% zscore based on reliability
% m = sum(dia'.* positionRaw(:,4).*positionRaw(:,7))/sum(positionRaw(:,4).*positionRaw(:,7));
% m = mean(dia);
% sd = sum((dia' - m).^2.*positionRaw(:,4).*positionRaw(:,7))/(sum(positionRaw(:,4).*positionRaw(:,7))-1);
% sd = sqrt(sd);
% diaZ = (dia-m)/sd;
diaZ = dia;
diaZ(~isnan(dia)) = zscore(dia(~isnan(dia)));
% autocorrelation
%         diaAutoMat = NaN(20, length(diaZ));
%         for j = 1:20
%             diaAutoMat(j,j+1:end) = diaZ(1:end-j);            
%         end
%         autolm = fitlm(diaAutoMat',diaZ');
%         diaZ = autolm.Residuals.Raw;

%% create diameter matrix 
cueMat = zeros(numBins,length(responseInds));
steps = round(0.001*(p.Results.startTime:p.Results.stepSize:p.Results.endTime)*FR);
% only for responded csplus
cueFT = cueFT(os.responseInds);
resMat = zeros(numBins,length(responseInds));
resFT = cueFT + round(0.001*FR*os.lickLat);    
% calcualte cueMat
for j = 1:length(responseInds)
    if cueFT(j)==0
        cueMat(:,j)= NaN;
        continue
    end
    for s = 1:numBins
    firstFr = cueFT(j)+steps(s)- 0.5*p.Results.binSizeFrame;
    if j < length(cueFT) && cueFT(j+1)>0 % find last frame
        [lastFr,I] = min([cueFT(j)+steps(s)+ 0.5*p.Results.binSizeFrame,length(diaZ),cueFT(j+1)-1]);
    else
        [lastFr,I] = min([cueFT(j)+steps(s)+ 0.5*p.Results.binSizeFrame,length(diaZ)]);
    end
    % calculate mean if calculatable
    if firstFr<=lastFr
        if firstFr < 1
            fprintf(num2str(j));
        end
        cueMat(s,j) = mean(diaZ(firstFr:lastFr));
    else
        cueMat(s,j) = NaN;
    end     
    % break the loop for current trial if meet the next trial's start
    if (I==3 || I==2) && j<length(responseInds)
       cueMat(s+1:end,j) = NaN;
       break
    end
    end
end
% calcualte resMat
for j = 1:length(responseInds)
    if cueFT(j)==0
        resMat(:,j)= NaN;
        continue
    end
    for s = 1:numBins
    firstFr = min(resFT(j)+steps(s)- 0.5*p.Results.binSizeFrame,length(diaZ));
    if j < length(resFT) && cueFT(j+1)>0
        [lastFr,I] = min([resFT(j)+steps(s)+ 0.5*p.Results.binSizeFrame,length(diaZ),cueFT(j+1)-1]);
    else
        [lastFr,I] = min([resFT(j)+steps(s)+ 0.5*p.Results.binSizeFrame,length(diaZ)]);
    end
    % calculate mean is calculatable
    if firstFr<=lastFr
        resMat(s,j) = mean(diaZ(firstFr:lastFr));
    else
        resMat(s,j) = NaN;
    end     
    % break the loop for current trial if meet the next trial's start
    if (I==3 || I==2) && j<length(responseInds)
       resMat(s+1:end,j) = NaN;
       break
    end
    end
end
%% behavior
% lickLat
% rwd   
% switch
svsTemp = find(choice(2:end) ~= choice(1:end-1)) + 1;
svs = zeros(1,length(responseInds));
svs(svsTemp) = 1;
% model 
paramEsts = [];
allSamples = [];
edges = cell(1,length(paramNames));
for i = 1:length(paramNames)
    tmp = samples.(p.Results.paramNames{i})(:,id);
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
eval(['[LL,probC,Q,pe] = qLearningModel_' p.Results.modelName '(paramEsts,' input ';'])
% diff value
Qdiff = abs(Q(:,2)-Q(:,1));
% total value
Qsum = sum(Q,2);
% pe
pe = pe';
% prepe
prePe = [NaN; pe(1:end-1)];
% conf
choiceConf = 2.*probC - 1;
%time in session
timeTemp = [behSessionData(responseInds).CSon] - behSessionData(responseInds(1)).CSon;

% side bias
biasSide = zeros(size(choice));
biasInds = contains(paramNames,'bias');
if paramEsts(biasInds)>0
    biasSide(os.lickR_Inds)=1;
else
    biasSide(os.lickL_Inds)=1;
end

% autocorrelation
% for j = 1:numBins
%     diaAutoMat = NaN(20, size(combineCueMat,1));
%     for a = 1:20
%         diaAutoMat(a,a+1:end) = combineCueMat(1:end-a,j);
%     end
%     autolm = fitlm(diaAutoMat',combineCueMat(:,j));
%     combineCueMatAuto(:,j) = autolm.Residuals.Raw;
% end
% 
% for j = 1:numBins
%     diaAutoMat = NaN(20, size(combineResMat(:,j),1));
%     for a = 1:20
%         diaAutoMat(a,a+1:end) = combineResMat(1:end-a,j);
%     end
%     autolm = fitlm(diaAutoMat',combineResMat(:,j));
%     combineResMat(:,j) = autolm.Residuals.Raw;
% end
%%
combineMat = [biasSide pe biasSide.*pe];
regressors = {'bias', 'pe', 'bias*pe'};
% combineMat = [os.lickLatZ(qualInd),...
%                 combineRwd(qualInd)-0.5,...
%                combinePreRwd(combineQualInd>0)-0.5,...
%                (combineRwd(combineQualInd>0)-0.5).*(combinePreRwd(combineQualInd>0)-0.5),...
%                 sign(combinePe(combineQualInd>0)),...
%                 abs(combinePrepe(combineQualInd>0)),...
%                 abs(combinePe(combineQualInd>0))-abs(combinePrepe(combineQualInd>0)),...
%                 (combineQSum(combineQualInd>0)),...
%                abs(combineQDiff(combineQualInd>0)),...                
%                 combineSvS(combineQualInd>0),...
%                 combineTimeInSession(combineQualInd>0)];
%             
% 
% combineMat = [combineLickLat(combineQualInd>0),...
%                combineTimeInSession(combineQualInd>0),...
%                 combineRwd(combineQualInd>0)-0.5,...
%                combinePreRwd(combineQualInd>0)-0.5,...
%                 combineSvS(combineQualInd>0)];
% 
%             
% combineMat = [combineLickLat(combineQualInd>0),...
%                 combineTimeInSession(combineQualInd>0),...
%                 combineQSum(combineQualInd>0),...
%                 abs(combinePrepe(combineQualInd>0)),...
%                 abs(combinePe(combineQualInd>0)),...
%                 combineSvS(combineQualInd>0)];
%              
%             
% combineMat = [sign(combinePe(combineQualInd>0)),...
%                 combineQSum(combineQualInd>0),...
%                 abs(combinePe(combineQualInd>0)),...
%                 abs(combinePrepe(combineQualInd>0)),...
%                 combineSvS(combineQualInd>0)];
%             
% combineMat = [combineLickLat(combineQualInd>0),...
%                 combineQSum(combineQualInd>0),...
%                 combineTimeInSession(combineQualInd>0),...
%                 sign(combinePe(combineQualInd>0)),...
%                 abs(combinePe(combineQualInd>0)),...
%                 combineSvS(combineQualInd>0)];
%             
%             
%             
% combineMat = [
%                 (combineQSum(combineQualInd>0)),...
%                abs(combineQDiff(combineQualInd>0)),...
%                 combineSvS(combineQualInd>0)];

%% zscore and fit model
for i = 1:size(combineMat,2)
    combineMat(~isnan(combineMat(:,i)),i) = zscore(combineMat(~isnan(combineMat(:,i)),i),0,1);
end

coeff = zeros(numBins,size(combineMat,2)+1,3);
sigs = zeros(numBins, size(combineMat,2)+1);
rsq = zeros(numBins,1);
resMat = resMat';
for i = 1:numBins
    lm = fitlm(combineMat(qualInd(responseInds), :),resMat(qualInd(responseInds),i));
    for j = 1:size(combineMat,2)+1
    coeff(i,j,1) = lm.Coefficients.Estimate(j);
    ci = coefCI(lm);
    coeff(i,j,2:3) = ci(j,:);
    sigs(i,j) = double(lm.Coefficients.pValue(j)<0.05);
    end
    rsq(i)=lm.Rsquared.Adjusted;
end


%% bin trials by pe
numBins = 10; % bins in pe
width = 4; % in num of bins
% decide start of the time windows
solenoidTimeEnd = max([behSessionData(responseInds).CSon] + 3100 - [behSessionData(responseInds).respondTime]); % 3600 ms is 500ms after sol coming back after 3100ms
[~,soleEndInd] = min(abs(time - solenoidTimeEnd)); %bin ind with solenoid end
solenoidTimeStart = min([behSessionData(responseInds).CSon] + 3100 - [behSessionData(responseInds).respondTime]);
[~,soleStartInd] = min(abs(time - solenoidTimeStart)); % bin ind with solenoid start
[~,rwdInd] = min(abs(time - os.rwdDelay)); % bin ind with solenoid start
peInd = find(contains(regressors,'pe'));
if ~isempty(peInd) % if pe is regressor
    sigIndEarly = find(sum(sigs(:,peInd+1),2)>0 & time'<solenoidTimeStart & time'>os.rwdDelay); 
    if ~isempty(sigIndEarly) % if any significance
        peEffect = sum(squeeze(abs(coeff(sigIndEarly,peInd(1:end-1)+1,1))),2);
        [~,maxInd] = max(peEffect);
        startTimeInd = max(sigIndEarly(maxInd) - 0.5*width, rwdInd);
    else 
        startTimeInd = rwdInd;
    end
    endTimeInd  = min(startTimeInd+width, soleStartInd);

    sigIndLate = find(sum(sigs(:,peInd+1),2)>0 & time' > solenoidTimeEnd & time'<p.Results.endTime);
    if ~isempty(sigIndLate)
        peLateEffect = sum(squeeze(abs(coeff(sigIndLate,peInd(1:end-1)+1,1))),2);
        [~,maxInd] = max(peLateEffect);
        startTimeLateInd = max(sigIndLate(maxInd) - 0.5*width, soleEndInd);
    else 
        startTimeLateInd = soleEndInd;
    end
    endTimeLateInd  = min(startTimeLateInd+width, length(time));
else
    startTimeInd = rwdInd;
    endTimeInd  = min(startTimeInd+width, soleStartInd);
    startTimeLateInd = soleEndInd;
    endTimeLateInd  = min(startTimeLateInd+width, length(time));
end
diaEarly = nanmean(resMat(:,startTimeInd:endTimeInd),2);
diaLate = nanmean(resMat(:,startTimeLateInd:endTimeLateInd),2);
%% caculate on both sides
edges = binEqualSize(pe, numBins);
diaMeans = zeros(numBins,1);
diaSems = zeros(numBins,1);
peMeans = zeros(numBins,1);
diaMeansLate = zeros(numBins,1);
diaSemsLate = zeros(numBins,1);
for k = 1:numBins
    if k < numBins
        diaTemp = diaEarly(pe >= edges(k) & pe < edges(k+1));
        diaTempLate = diaLate(pe >= edges(k) & pe < edges(k+1));
        peMeans(k) = mean(pe(pe >= edges(k) & pe < edges(k+1)));
    else
        diaTemp = diaEarly(pe >= edges(k) & pe <= edges(k+1));
        diaTempLate = diaLate(pe >= edges(k) & pe <= edges(k+1));
        peMeans(k) = mean(pe(pe >= edges(k) & pe <= edges(k+1)));
    end
    diaMeans(k) = nanmean(diaTemp);
    diaMeansLate(k) = nanmean(diaTempLate);
    diaSems(k) = sem(diaTemp(~isnan(diaTemp)));
    diaSemsLate(k) = sem(diaTempLate(~isnan(diaTempLate)));
end

%% separate left
peL = pe(os.lickL_Inds);
diaEarlyL = diaEarly(os.lickL_Inds);
diaLateL = diaLate(os.lickL_Inds);
edgesL = binEqualSize(peL, numBins);
diaMeansL = zeros(numBins,1);
diaSemsL = zeros(numBins,1);
diaMeansLLate = zeros(numBins,1);
diaSemsLLate = zeros(numBins,1);
peMeansL = zeros(numBins,1);
for k = 1:numBins
    if k < numBins
        diaTemp = diaEarlyL(peL >= edgesL(k) & peL < edgesL(k+1));
        diaTempLate = diaLateL(peL >= edgesL(k) & peL < edgesL(k+1));
        peMeansL(k) = mean(peL(peL >= edgesL(k) & peL < edgesL(k+1)));
    else
        diaTemp = diaEarlyL(peL >= edgesL(k) & peL <= edgesL(k+1));
        diaTempLate = diaLateL(peL >= edgesL(k) & peL <= edgesL(k+1));
        peMeansL(k) = mean(peL(peL >= edgesL(k) & peL <= edgesL(k+1)));
    end
    diaMeansL(k) = nanmean(diaTemp);
    diaMeansLLate(k) = nanmean(diaTempLate);
    diaSemsL(k) = sem(diaTemp(~isnan(diaTemp)));
    diaSemsLLate(k) = sem(diaTempLate(~isnan(diaTempLate)));
end

%% separate right
peR = pe(os.lickR_Inds);
diaEarlyR = diaEarly(os.lickR_Inds);
diaLateR = diaLate(os.lickR_Inds);
edgesR = binEqualSize(peR, numBins);
diaMeansR = zeros(numBins,1);
diaSemsR = zeros(numBins,1);
peMeansR = zeros(numBins,1);
diaMeansRLate = zeros(numBins,1);
diaSemsRLate = zeros(numBins,1);
for k = 1:numBins
    if k < numBins
        diaTemp = diaEarlyR(peR >= edgesR(k) & peR < edgesR(k+1));
        diaTempLate = diaLateR(peR >= edgesR(k) & peR < edgesR(k+1));
        peMeansR(k) = mean(peR(peR >= edgesR(k) & peR < edgesR(k+1)));
    else
        diaTemp = diaEarlyR(peR >= edgesR(k) & peR <= edgesR(k+1));
        diaTempLate = diaLateR(peR >= edgesR(k) & peR <= edgesR(k+1));
        peMeansR(k) = mean(peR(peR >= edgesR(k) & peR <= edgesR(k+1)));
    end
    diaMeansR(k) = nanmean(diaTemp);
    diaMeansRLate(k) = nanmean(diaTempLate);
    diaSemsR(k) = sem(diaTemp(~isnan(diaTemp)));
    diaSemsRLate(k) = sem(diaTempLate(~isnan(diaTempLate)));       
end

%% plots 
GLM=figure2('position',[0 0 800 1600]);
colors = cool(size(combineMat,2));

if ismember('bias', p.Results.paramNames)
    biasInd = contains(p.Results.paramNames, 'bias');
    if paramEsts(biasInd) > 0
        colorR = [1 0 1];
        colorL = [0 1 1];
    else
        colorR = [0 1 1];
        colorL = [1 0 1];    
    end
else
    colorR = [0 1 1];
    colorR = [0 1 1];
end
    
subplot(4,2,[1 2]); hold on;
x = 0.001*time;
yyaxis right
plot(x, coeff(:,1,1), 'Color', [0.5 0.5 0.5], 'LineStyle','-', 'Marker','none', 'linewidth', 2);
ylim([1.2*min(coeff(:,:,2),[],'all') 1.2*max(coeff(:,:,3),[],'all')])
ylabel('intercept')

yyaxis left 
for i = 1:size(combineMat,2)
    plot(x, coeff(:,i+1,1), 'Color', colors(i,:), 'LineStyle','-', 'Marker','none', 'linewidth', 2); 
end

for i = 1:size(combineMat,2)
   fill([x fliplr(x)], [coeff(:,i+1,2)' fliplr(coeff(:,i+1,3)')], colors(i,:), 'facealpha', 0.25, 'edgecolor', 'none')
end

yyaxis right
fill([x fliplr(x)], [coeff(:,1,2)' fliplr(coeff(:,1,3)')], [0.5 0.5 0.5], 'facealpha', 0.25, 'edgecolor', 'none')
plt = gca;
plt.YAxis(2).Color = [0.2 0.2 0.2];

yyaxis left
line(minmax(x), [0 0],'Color', [0.2 0.2 0.2], 'LineStyle','--')
line([0.3 0.3], [-0.2 0.4],'Color', [1 0 0], 'LineStyle','--')
line([0 0], [-0.2 0.4],'Color', [0.2 0.2 0.2], 'LineStyle','--')
line([x(startTimeInd) x(endTimeInd)], [max(coeff(:,2:end,3),[],'all'),max(coeff(:,2:end,3),[],'all')],'Color', [1 0 0], 'LineWidth',4)
line([x(startTimeLateInd) x(endTimeLateInd)],[max(coeff(:,2:end,3),[],'all'),max(coeff(:,2:end,3),[],'all')],'Color', [0 0 1], 'LineWidth',4)
ylim([1.2*min(coeff(:,2:end,2),[],'all') 1.2*max(coeff(:,2:end,3),[],'all')]);
xlim(minmax(x))
ylabel('\beta coefficients')
xlabel('respond time /s')
legend(regressors, 'Location', 'best')
titileText = [session ': pupil diameter aligned to choice'];
title(titileText,'interpreter','none');
subplot(4,2,3); hold on;
errorbar(peMeans, diaMeans, diaSems,'linewidth', 2)

title('pupil dia change with pe')

subplot(4,2,4); hold on;
errorbar(peMeans, diaMeansLate, diaSemsLate,'linewidth', 2)
title('pupil dia after solenoid change with pe')

subplot(4,2,5); hold on;
errorbar(peMeansL, diaMeansL, diaSemsL,'linewidth', 2, 'color', colorL)
title(['pupil dia change with left pe' sprintf(('bias = %.2f left = %d'), paramEsts(biasInd), length(os.lickL_Inds))])

subplot(4,2,6); hold on;
errorbar(peMeansL, diaMeansLLate, diaSemsLLate,'linewidth', 2, 'color', colorL);
title(['pupil dia change with right pe' sprintf(('bias = %.2f right = %d'), paramEsts(biasInd), length(os.lickR_Inds))])

subplot(4,2,7); hold on;
errorbar(peMeansR, diaMeansR, diaSemsR,'linewidth', 2, 'color', colorR)
title('pupil dia change with right pe')

subplot(4,2,8); hold on;
errorbar(peMeansR, diaMeansRLate, diaSemsRLate,'linewidth', 2, 'color', colorR)
title('pupil dia after solenoid change with right pe')

if p.Results.saveFigFlag == 1 
    saveFigurePDF(GLM,[savePath session '_ pupilGLM'])
end



