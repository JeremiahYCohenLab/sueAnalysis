function pupilAnalysisRaw(animal, category, varargin)
%task and model parameters
p = inputParser;
% default parameters if none given
p.addParameter('plotFlag', 1);
p.addParameter('maxTrial', 1000);
p.addParameter('modelNameOld', '4params');
p.addParameter('modelName', '4params');
p.addParameter('paramNames', {'aN', 'aP', 'aF', 'beta'});
p.addParameter('binSizeFrame', 6)% need to be even number
p.addParameter('stepSize', 250)
p.addParameter('startTime', -2000)
p.addParameter('endTime',3500)
p.parse(varargin{:});

paramNames = p.Results.paramNames;
maxTrial = p.Results.maxTrial;
numBins = (p.Results.endTime - p.Results.startTime)/p.Results.stepSize +1;


combineCueMat = [];
combineResMat = [];
combineRwd = [];
combinePe = [];
combineQSum = [];
combineQDiff = [];
combineLickLat = [];
combineQualInd = [];
combineSvS = [];
combinePrepe = [];
combineTimeInSession = [];

if contains(p.Results.modelName,'tF')
    input = 'choice, outcome, ITI)';
else
    input = 'choice, outcome)';
end


for z = 1:length(animal)
% load model fitting results 
[root, sep] = currComputer();
sampFile = [animal{z} category '_', p.Results.modelNameOld];
path = [root animal{z} sep animal{z} 'sorted' sep 'stan' sep 'bernoulli' sep p.Results.modelNameOld sep];
load([path sampFile '.mat'], 'dayList');
samples = load([path sampFile '.mat'], sampFile);
samples = samples.(sampFile);

    for i = 1:length(dayList)
        %% load files, first check if there's pupil data
        fprintf([dayList{i} '\n'])
        session = dayList{i};
        [animalName, date] = strtok(session, 'd'); 
        animalName = animalName(2:end);
        date = date(1:9);
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
        end
        % behavior
        if exist(behSessionDataPath,'file')
            load(behSessionDataPath)
        else
            behSessionData = generateSessionData_operantMatchingDecoupledRwdDelay(session);
        end
        % pupil
        if exist(pupilAlignPath, 'file')
            load(pupilAlignPath)
        else
            [errorProp, cueFT, qualInd, FR] = timeAlign(session, 1);
        end
        % jump the current loop if not well aligned
        if errorProp > 0.2 || isnan(errorProp)
            continue
        end
        %% preparation 
        % parse behavior
        os = parseBehavioralData(behSessionData, maxTrial);
        choice = os.allChoices;
        choice(choice<0) = 0;
        outcome = abs(os.allRewards);
        choice = choice(1:min(length(choice), maxTrial));
        outcome = outcome(1:length(choice));
        responseInds = os.responseInds(1:min(length(choice), maxTrial)); 
        ITI = os.timeBtwn(1:length(choice));

        % load diameter
        list = dir(videopath);
        expression = ['^' session 'DLC' '\w*' 'shuffle2_480000.csv' '$'];
        expressionSkeleton = ['^' session 'DLC' '\w*' 'shuffle2_480000_skeleton.csv' '$'];
        position = list(~cellfun(@isempty, cellfun(@(x) regexp(x, expression), {list.name}, 'UniformOutput', false))).name;
        skeleton = list(~cellfun(@isempty, cellfun(@(x) regexp(x, expressionSkeleton), {list.name}, 'UniformOutput', false))).name;
        positionRaw = csvread([videopath sep position], 3, 0);
        diaRaw = csvread([videopath sep skeleton], 2, 0);
        qualF = positionRaw(:,4) > 0.99 &  positionRaw(:,7) > 0.99;   
        dia = interp1(find(qualF>0),diaRaw(qualF,2),1:length(diaRaw(:,2)));

        % zscore based on reliability
%         m = sum(dia'.* positionRaw(:,4).*positionRaw(:,7))/sum(positionRaw(:,4).*positionRaw(:,7));
%         sd = sum((dia' - m).^2.*positionRaw(:,4).*positionRaw(:,7))/(sum(positionRaw(:,4).*positionRaw(:,7))-1);
%         sd = sqrt(sd);
%         diaZ = (dia-m)/sd;
        diaZ = dia;
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
            if j < length(cueFT) && cueFT(j+1)>0
                [lastFr,I] = min([cueFT(j)+steps(s)+ 0.5*p.Results.binSizeFrame,length(diaZ),cueFT(j+1)-1]);
            else
                [lastFr,I] = min([cueFT(j)+steps(s)+ 0.5*p.Results.binSizeFrame,length(diaZ)]);
            end
            % calculate mean is calculatable
            if firstFr<=lastFr && firstFr>0
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
            % calculate mean
            if firstFr<=lastFr && firstFr>0
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
        for j = 1:length(paramNames)
            tmp = samples.(paramNames{j})(:,i);
            [n,e] = histcounts(tmp, 50);
            [~, maxInd] = max(n);
            params(j) = median(tmp(tmp > e(maxInd) & tmp < e(maxInd+1)));
        end

        eval(['[~,~,Qtemp,petemp,cQtemp] = qLearningModel_' p.Results.modelName '(params,' input ';'])
        % diff value
        Qdiff = abs(Qtemp(1:end-1,2)-Qtemp(1:end-1,1));
        % total value
        Qsum = sum(Qtemp(1:end-1,:),2);
        % chosen value
        cQ = cQtemp;
        % pe
        pe = petemp;
        % prepe
        prePe = [NaN pe(1:end-1)'];
        
        %time in session
        timeTemp = [behSessionData(responseInds).CSon] - behSessionData(responseInds(1)).CSon;
        %% consolidaterep
        combineCueMat = [combineCueMat; NaN(5,size(cueMat,1)); cueMat'];
        combineResMat = [combineResMat; NaN(5,size(resMat,1)); resMat'];
        combineRwd = [combineRwd; NaN(5,1); outcome'];
        combinePe = [combinePe; zeros(5,1); pe];
        combineQDiff = [combineQDiff; NaN(5,1); Qdiff];
        combineQSum = [combineQSum; NaN(5,1); Qsum];
        combineLickLat = [combineLickLat; NaN(5,1); os.lickLatLogZ'];
        combineQualInd = [combineQualInd; zeros(5,1); qualInd(responseInds)'];
        combineSvS = [combineSvS; NaN(5,1); svs'];
        combinePrepe = [combinePrepe; NaN(5,1); prePe'];
        combineTimeInSession = [combineTimeInSession; NaN(5,1); timeTemp'];

    end
end
combinePreRwd = [0;combineRwd(1:end-1)];
combineCueMatAuto = NaN(size(combineCueMat));
% autocorrelation
for j = 1:numBins
    diaAutoMat = NaN(20, size(combineCueMat,1));
    for a = 1:20
        diaAutoMat(a,a+1:end) = combineCueMat(1:end-a,j);
    end
    autolm = fitlm(diaAutoMat',combineCueMat(:,j));
    combineCueMatAuto(:,j) = autolm.Residuals.Raw;
end

for j = 1:numBins
    diaAutoMat = NaN(20, size(combineResMat(:,j),1));
    for a = 1:20
        diaAutoMat(a,a+1:end) = combineResMat(1:end-a,j);
    end
    autolm = fitlm(diaAutoMat',combineResMat(:,j));
    combineResMat(:,j) = autolm.Residuals.Raw;
end
%%
combineMat = [combineLickLat(combineQualInd>0),...
                combineRwd(combineQualInd>0)-0.5,...
               combinePreRwd(combineQualInd>0)-0.5,...
               (combineRwd(combineQualInd>0)-0.5).*(combinePreRwd(combineQualInd>0)-0.5),...
                sign(combinePe(combineQualInd>0)),...
                abs(combinePrepe(combineQualInd>0)),...
                abs(combinePe(combineQualInd>0))-abs(combinePrepe(combineQualInd>0)),...
                (combineQSum(combineQualInd>0)),...
               abs(combineQDiff(combineQualInd>0)),...                
                combineSvS(combineQualInd>0),...
                combineTimeInSession(combineQualInd>0)];
            

combineMat = [combineLickLat(combineQualInd>0),...
               combineTimeInSession(combineQualInd>0),...
                combineRwd(combineQualInd>0)-0.5,...
               combinePreRwd(combineQualInd>0)-0.5,...
                combineSvS(combineQualInd>0)];

            
combineMat = [combineLickLat(combineQualInd>0),...
                combineTimeInSession(combineQualInd>0),...
                combineQSum(combineQualInd>0),...
                abs(combinePrepe(combineQualInd>0)),...
                abs(combinePe(combineQualInd>0)),...
                combineSvS(combineQualInd>0)];
             
            
combineMat = [sign(combinePe(combineQualInd>0)),...
                combineQSum(combineQualInd>0),...
                abs(combinePe(combineQualInd>0)),...
                abs(combinePrepe(combineQualInd>0)),...
                combineSvS(combineQualInd>0)];
            
combineMat = [combineLickLat(combineQualInd>0),...
                combineQSum(combineQualInd>0),...
                combineTimeInSession(combineQualInd>0),...
                sign(combinePe(combineQualInd>0)),...
                abs(combinePe(combineQualInd>0)),...
                combineSvS(combineQualInd>0)];
            
            
            
combineMat = [
                (combineQSum(combineQualInd>0)),...
               abs(combineQDiff(combineQualInd>0)),...
                combineSvS(combineQualInd>0)];

%% zscore and fit model
for i = 1:size(combineMat,2)
    combineMat(~isnan(combineMat(:,i)),i) = zscore(combineMat(~isnan(combineMat(:,i)),i),0,1);
end

coeff = zeros(numBins,size(combineMat,2)+1,3);
rsq = zeros(numBins,1);
for i = 1:numBins
    lm = fitlm(combineMat,combineResMat(combineQualInd>0,i));
    for j = 1:size(combineMat,2)+1
    coeff(i,j,1) = lm.Coefficients.Estimate(j);
    ci = coefCI(lm);
    coeff(i,j,2:3) = ci(j,:);
    end
    rsq(i)=lm.Rsquared.Adjusted;
end
%% plots 
figure2('position',[0 0 800 1200]);
colors = jet(size(combineMat,2));
subplot(3,1,[1 2]); hold on;
x = 0.001*(p.Results.startTime:p.Results.stepSize:p.Results.endTime);
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
ylim([1.2*min(coeff(:,2:end,2),[],'all') 1.2*max(coeff(:,2:end,3),[],'all')]);
xlim(minmax(x))
ylabel('\beta coefficients')
xlabel('respond time')

subplot(3,1,3); hold on;
plot(x, rsq, 'linewidth', 2)
ylim([0 1.2*max(rsq)])
title('adjusted Rsquare')
%%
legend('lickLat','Rwd','PreRwd','interaction','SvS')

legend('lickLat','Qsum','pe', 'inter', 'SvS')

legend('lickLat','Pe','sign','interaction','SvS')

legend('sum','diff','SvS')
%% analysis

% load behavior



