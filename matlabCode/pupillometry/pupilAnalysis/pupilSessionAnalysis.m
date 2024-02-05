function pupilSessionAnalysis(session, varargin)
%task and model parameters
p = inputParser;
% default parameters if none given
p.addParameter('maxTrial', 1000);
p.addParameter('startTime', -2000)
p.addParameter('endTime',10000) % no longer than 10
p.addParameter('saveFigFlag', 1)
p.parse(varargin{:});

% load model fitting results 
[animalName, date] = strtok(session, 'd'); 
animalName = animalName(2:end);
date = date(1:9);
[root, sep] = currComputer();

%% load files, first check if there's pupil data
sessionFolder = ['m' animalName date];  
% paths
if isstrprop(session(end), 'alpha')
    behSessionDataPath = [root animalName sep sessionFolder sep 'sorted' sep 'session ' session(end) sep session '_sessionData_behav.mat'];
    pupilAlignPath = [root animalName sep sessionFolder sep 'sorted' sep 'session ' session(end) sep session '_pupil.mat'];
%     pupilPath = [root animalName sep sessionFolder sep 'pupil'];
    savePath = [root animalName sep sessionFolder sep 'figures' sep];
else
    behSessionDataPath = [root animalName sep sessionFolder sep 'sorted' sep 'session' sep session '_sessionData_behav.mat'];
    pupilAlignPath = [root animalName sep sessionFolder sep 'sorted' sep 'session' sep session '_pupil.mat'];
    savePath = [root animalName sep sessionFolder sep 'figures' sep];
end
% behavior
if exist(behSessionDataPath,'file')
    load(behSessionDataPath)
end
if ~exist('behSessionData', 'var')
    behSessionData = sessionData;
end


% pupil   
if exist(pupilAlignPath, 'file')
    load(pupilAlignPath);
else
    fprintf([session ' pupil not aligned yet \n'])
    return
end
% jump the current loop if not well aligned
if exist('eyeBlurryLate')
    if eyeBlurryLate && errorProp > 0.8
        fprintf([session ' pupil not well aligned \n'])
        return
    else
        if ~eyeBlurryLate && errorProp > 0.7
            fprintf([session ' pupil not well aligned \n'])
            return
        end
    end
else
    if errorProp > 0.3
        fprintf([session ' pupil not well aligned \n'])
        return       
    end
end
%% preparation 
% parse behavior
os = behAnalysisNoPlot_opMD(session, 'simpleFlag', 0);
origBlockSwitch = blockSwitch(blockSwitch<=length(behSessionData));




%% plot everything with behavior
figure;
screenSize = get(0,'Screensize');
screenSize(4) = screenSize(4) - 100;
set(gcf, 'Position', screenSize)
sgtitle(session)
%% baseline vs dilation
baseline = mean(sessionPupilCue(:,1:round(2*FR)),2,'omitnan');
dilation = max(sessionPupilCue(:,round(2*FR)+1:end),[],2,'omitnan')-baseline;
subplot(6,8,17)
scatter(baseline(qualInd>0),dilation(qualInd>0),6,'c','filled');
[R,P] = corrcoef(baseline(qualInd),dilation(qualInd));
title(sprintf('%0.2f p:%0.2f', R(1,2), P(1,2)));
xlabel('baseline')
ylabel('dilation')
subplot(6,8,18)
scatter(baseline(qualInd),dilation(qualInd)./baseline(qualInd),6,'m','filled');
[R,P] = corrcoef(baseline(qualInd),dilation(qualInd)./baseline(qualInd));
title(sprintf('%0.2f p:%0.2f', R(1,2), P(1,2)));
xlabel('baseline')
ylabel('dilation/baseline')
%%
% time plot
subplot(6,8,[1:8]); hold on
xlabel('Time (min)')
nrMag = 0.5;
rMag = 1;

for i = 1:length(behSessionData)
    currTime = (behSessionData(i).CSon - behSessionData(1).CSon)/1000/60; %convert to min
    if strcmp(behSessionData(i).trialType,'CSplus')
        if ~isnan(behSessionData(i).rewardR)
            if behSessionData(i).rewardR == 1 % R side rewarded
                plot([currTime currTime],[0 rMag],'k')
            else
                plot([currTime currTime],[0 nrMag],'k') % R side not rewarded
            end
        elseif ~isnan(behSessionData(i).rewardL)
            if behSessionData(i).rewardL == 1 % L side rewarded
                plot([currTime currTime],[-1*rMag 0],'k')
            else
                plot([currTime currTime],[-1*nrMag 0],'k')
            end
            
        else % CSplus trial but no rewardL or rewardR
            plot([currTime currTime],[-rMag rMag],'r')
        end
    else % CS minus trial
        plot([currTime currTime],0,'ko','markersize',4,'linewidth',2)
        if ~isempty(behSessionData(i).licksL) & behSessionData(i).licksL(1) - behSessionData(i).CSon < 2000
            plot([currTime currTime],[-1*nrMag 0],'k')
        elseif ~isempty(behSessionData(i).licksR) & behSessionData(i).licksR(1) - behSessionData(i).CSon < 2000
            plot([currTime currTime],[nrMag 0],'k')
        end 
    end
    if any(i == origBlockSwitch)
        plot([currTime currTime],[-1*rMag rMag],'--','linewidth',1,'Color',[30 144 255]./255)
    end
%     if i > responseInds(1) & lickLat(i) > 250 & ~isnan(lickLat(i))
%         plot([currTime currTime], [0 normRespLat(j)], '-', 'Color', [0.85 0.325 0.098])
%         j = j + 1;
%     end
end
xlim([0 currTime]);
%% overlay baseline pupil diameter
for i = 1:length(cueFT)
    if qualInd(i)
        if i==1
            currTime = [behSessionData(i).CSon - 2000, behSessionData(i).CSon] - behSessionData(1).CSon;
        else
            currTime = [behSessionData(i-1).CSon, behSessionData(i).CSon] - behSessionData(1).CSon;

        end
        currTime = currTime/60000;
        currPupil = (baseline(i)-min(baseline))*[1,1]/(max(baseline)-min(baseline));
        fill([currTime flip(currTime)], [currPupil -flip(currPupil)], 'c', 'edgeColor', 'none','FaceAlpha', '0.3'); hold on;
    end
end
%% pupil dilation
subplot(6,8,[9:16]); hold on;
for i = 1:length(behSessionData)
    currTime = (behSessionData(i).CSon - behSessionData(1).CSon)/1000/60; %convert to min
    if strcmp(behSessionData(i).trialType,'CSplus')
        if ~isnan(behSessionData(i).rewardR)
            if behSessionData(i).rewardR == 1 % R side rewarded
                plot([currTime currTime],[0 rMag],'k')
            else
                plot([currTime currTime],[0 nrMag],'k') % R side not rewarded
            end
        elseif ~isnan(behSessionData(i).rewardL)
            if behSessionData(i).rewardL == 1 % L side rewarded
                plot([currTime currTime],[-1*rMag 0],'k')
            else
                plot([currTime currTime],[-1*nrMag 0],'k')
            end
            
        else % CSplus trial but no rewardL or rewardR
            plot([currTime currTime],[-rMag rMag],'r')
        end
    else % CS minus trial
        plot([currTime currTime],0,'ko','markersize',4,'linewidth',2)
        if ~isempty(behSessionData(i).licksL) & behSessionData(i).licksL(1) - behSessionData(i).CSon < 2000
            plot([currTime currTime],[-1*nrMag 0],'k')
        elseif ~isempty(behSessionData(i).licksR) & behSessionData(i).licksR(1) - behSessionData(i).CSon < 2000
            plot([currTime currTime],[nrMag 0],'k')
        end 
    end
    if any(i == origBlockSwitch)
        plot([currTime currTime],[-1*rMag rMag],'--','linewidth',1,'Color',[30 144 255]./255)
    end
%     if i > responseInds(1) & lickLat(i) > 250 & ~isnan(lickLat(i))
%         plot([currTime currTime], [0 normRespLat(j)], '-', 'Color', [0.85 0.325 0.098])
%         j = j + 1;
%     end
end
xlim([0 currTime]);
%% overlay baseline pupil dilation
for i = 1:length(cueFT)
    if qualInd(i)
        if i==length(cueFT)
            currTime = [behSessionData(i).CSon, behSessionData(i).CSon+5000] - behSessionData(1).CSon;
        else
            currTime = [behSessionData(i).CSon, behSessionData(i+1).CSon] - behSessionData(1).CSon;
        end
        currTime = currTime/60000;
        currPupil = (dilation(i)-min(dilation))*[1,1]/(max(dilation)-min(dilation));
        fill([currTime flip(currTime)], [currPupil -flip(currPupil)], 'r', 'edgeColor', 'none','FaceAlpha', '0.3'); hold on;
    end
end
%% autoCorrlation
subplot(6,8,[19:21])
autocorr(baseline(qualInd),'NumLags',min([sum(qualInd)-1, 50]),'NumSTD',3);
title('baseline autoCorr')
xlabel('lags in trial')
subplot(6,8,[22:24]);
autocorr(dilation(qualInd),'NumLags',min([sum(qualInd)-1, 50]),'NumSTD',3);
title('dilation autoCorr')
xlabel('lags in trial')
%% preparation
time = linspace(-2,10,size(sessionPupilCue,2));
minP = min(sessionPupilCue,[],'all');
maxP = max(sessionPupilCue,[],'all');
%% compare rwd vs non-rwd
subplot(6,8,[25:26]); hold on;
my_SDF_1 = sessionPupilCue(intersect(os.responseInds(os.rwd_Inds), find(qualInd>0)),:);
my_SDF_2 = sessionPupilCue(intersect(os.responseInds(os.nrwd_Inds), find(qualInd>0)),:);
plotFilled(time, my_SDF_1,'r');
plotFilled(time, my_SDF_2,'b');
legend({'rwd', '', 'nrwd', ''})

%% compare swtich vs stay
subplot(6,8,[27:28]); hold on;
my_SDF_1 = sessionPupilCue(intersect(os.responseInds(os.changeChoice_Inds), find(qualInd>0)),:);
my_SDF_2 = sessionPupilCue(intersect(os.responseInds(os.stayChoice_Inds), find(qualInd>0)),:);
plotFilled(time, my_SDF_1,'r');
plotFilled(time, my_SDF_2,'b');
legend({'switch', '', 'stay', ''})
%% compare swtich/stay rwd/non-rwd
subplot(6,8,[29:30]); hold on;
my_SDF_1 = sessionPupilCue(intersect(os.responseInds(intersect(os.changeChoice_Inds, os.rwd_Inds)), find(qualInd>0)),:);
my_SDF_2 = sessionPupilCue(intersect(os.responseInds(intersect(os.changeChoice_Inds, os.nrwd_Inds)), find(qualInd>0)),:);
my_SDF_3 = sessionPupilCue(intersect(os.responseInds(intersect(os.stayChoice_Inds, os.rwd_Inds)), find(qualInd>0)),:);
my_SDF_4 = sessionPupilCue(intersect(os.responseInds(intersect(os.stayChoice_Inds, os.nrwd_Inds)), find(qualInd>0)),:);

plotFilled(time, my_SDF_1,'r');
plotFilled(time, my_SDF_2,'b');
plotFilled(time, my_SDF_3,'m');
plotFilled(time, my_SDF_4,'c');
legend({'c-r', '', 'c-n', '','s-r', '', 's-n', ''})
%% compare rwd/non-rwd then swtich/stay 
subplot(6,8,[31:32]); hold on;
my_SDF_1 = sessionPupilCue(intersect(os.responseInds(intersect(os.changeChoice_Inds, os.rwd_Inds+1)), find(qualInd>0)),:);
my_SDF_2 = sessionPupilCue(intersect(os.responseInds(intersect(os.changeChoice_Inds, os.nrwd_Inds+1)), find(qualInd>0)),:);
my_SDF_3 = sessionPupilCue(intersect(os.responseInds(intersect(os.stayChoice_Inds, os.rwd_Inds+1)), find(qualInd>0)),:);
my_SDF_4 = sessionPupilCue(intersect(os.responseInds(intersect(os.stayChoice_Inds, os.nrwd_Inds+1)), find(qualInd>0)),:);

plotFilled(time, my_SDF_1,'r');
plotFilled(time, my_SDF_2,'b');
plotFilled(time, my_SDF_3,'m');
plotFilled(time, my_SDF_4,'c');
legend({'r-c', '', 'n-c', '','r-s', '', 'n-s', ''})
%% compare CS+ CS-
subplot(6,8,[33:34]); hold on;
my_SDF_1 = sessionPupilCue(os.CSplus&qualInd,:);
my_SDF_2 = sessionPupilCue(os.CSminus&qualInd,:);
plotFilled(time, my_SDF_1,'r');
plotFilled(time, my_SDF_2,'b');
legend({'CS+', '', 'CS-',''})
%% split by rwd history
% Generate smoothed reward-history values over trials

[~,rwdHx_Inds] = sort(os.rwdHx);

%for tercile analysis
tercile = floor(length(rwdHx_Inds)/3);
rwdHxI_Inds = rwdHx_Inds(1:tercile);
rwdHxII_Inds = rwdHx_Inds(tercile+1:tercile*2);
rwdHxIII_Inds = rwdHx_Inds(tercile*2+1:end);

subplot(6,8,[35:36]); hold on;
my_SDF_1 = sessionPupilCue(intersect(os.responseInds(intersect(rwdHxI_Inds, os.rwd_Inds)), find(qualInd>0)),:);
my_SDF_2 = sessionPupilCue(intersect(os.responseInds(intersect(rwdHxII_Inds, os.rwd_Inds)), find(qualInd>0)),:);
my_SDF_3 = sessionPupilCue(intersect(os.responseInds(intersect(rwdHxIII_Inds, os.rwd_Inds)), find(qualInd>0)),:);
plotFilled(time, my_SDF_1, [1 0.6 0.6]);
plotFilled(time, my_SDF_2, [1 0.3 0.3]);
plotFilled(time, my_SDF_3, [1 0 0]);
legend({'low-r', '', 'medium-r','', 'high-r', ''})

subplot(6,8,[43:44]); hold on;
my_SDF_1 = sessionPupilCue(intersect(os.responseInds(intersect(rwdHxI_Inds, os.nrwd_Inds)), find(qualInd>0)),:);
my_SDF_2 = sessionPupilCue(intersect(os.responseInds(intersect(rwdHxII_Inds, os.nrwd_Inds)), find(qualInd>0)),:);
my_SDF_3 = sessionPupilCue(intersect(os.responseInds(intersect(rwdHxIII_Inds, os.nrwd_Inds)), find(qualInd>0)),:);
plotFilled(time, my_SDF_1, [0.6 0.6 1]);
plotFilled(time, my_SDF_2, [0.3 0.3 1]);
plotFilled(time, my_SDF_3, [0 0 1]);
legend({'low-nr', '', 'medium-nr','', 'high-nr', ''})
%% split by lick latency
[~,lickLat_Inds] = sort(os.lickLat);

%for tercile analysis
tercile = floor(length(lickLat_Inds)/3);
lickLatI_Inds = lickLat_Inds(1:tercile);
lickLatII_Inds = lickLat_Inds(tercile+1:tercile*2);
lickLatIII_Inds = lickLat_Inds(tercile*2+1:end);

subplot(6,8,[41:42]); hold on;
my_SDF_1 = sessionPupilCue(intersect(os.responseInds(lickLatI_Inds), find(qualInd>0)),:);
my_SDF_2 = sessionPupilCue(intersect(os.responseInds(lickLatII_Inds), find(qualInd>0)),:);
my_SDF_3 = sessionPupilCue(intersect(os.responseInds(lickLatIII_Inds), find(qualInd>0)),:);
plotFilled(time, my_SDF_1, [1 0 0]);
plotFilled(time, my_SDF_2, [1 0.3 0.3]);
plotFilled(time, my_SDF_3, [1 0.6 0.6]);
legend({'fast', '', 'medium','',  'slow', ''})
%% rwd and switch glm
% rwdMat
rwdOrnMat = nan(5, size(os.rwdMatx,2));
rwdOrn = zeros(1,length(os.responseInds));
rwdOrn(os.rwd_Inds) = 1;
for i = 1:size(rwdOrnMat,1)
    rwdOrnMat(i, i:end) = rwdOrn(1:end+1-i);
end
% switchMat
switchMat = nan(size(rwdOrnMat));
switches = zeros(1,length(os.responseInds));
switches(os.changeChoice_Inds) = 1;
for i = 1:size(switchMat,1)
    switchMat(i, i:end) = switches(1:end+1-i);
end

glm = fitlm([rwdOrnMat(:,ismember(os.responseInds, find(qualInd>0)))', switchMat(:,ismember(os.responseInds, find(qualInd>0)))'], dilation(intersect(os.responseInds, find(qualInd>0))));

%% rwdHist, switch, current outcome, interaction
stepSize = 3; % in frames
binSize = 6;  % in frames

midpoints = round(0.5*binSize+1):stepSize:size(sessionPupilCue,2)-0.5*binSize+1;
pupilSlide = zeros(size(sessionPupilCue,1), length(midpoints));

for i = 1:length(midpoints)
    pupilSlide(:,i) = mean(sessionPupilCue(:,midpoints(i)-0.5*binSize:midpoints(i)+0.5*binSize-1),2,'omitnan');
end

% regressors
if ~isempty(find(os.laser)==1)
    combineMat = [rwdOrn', switches', zscore(os.rwdHx)', [os.laser]'];
else
    combineMat = [rwdOrn', switches', zscore(os.rwdHx)'];
end

coeff = zeros(length(midpoints),size(combineMat,2)+1,3);
sigs = zeros(length(midpoints), size(combineMat,2)+1);
rsq = zeros(length(midpoints),1);
for i = 1:length(midpoints)
    lm = fitlm(combineMat(ismember(os.responseInds, find(qualInd>0)),:),pupilSlide(os.responseInds(ismember(os.responseInds, find(qualInd>0))),i));
    for j = 1:length(lm.CoefficientNames)
    coeff(i,j,1) = lm.Coefficients.Estimate(j);
    ci = coefCI(lm);
    coeff(i,j,2:3) = ci(j,:);
    sigs(i,j) = double(lm.Coefficients.pValue(j)<0.05);
    end
    rsq(i)=lm.Rsquared.Adjusted;
end
%plot
subplot(6,8,[45, 46, 37, 38]); hold on
colors = cool(size(coeff,2));
slideTime = linspace(-2,10, length(midpoints));
for i = 1:size(coeff,2)-1
    plot(slideTime, coeff(:,i+1,1), 'Color', colors(i,:), 'LineStyle','-', 'Marker','none', 'linewidth', 2); 
end

for i = 1:size(coeff,2)-1
   fill([slideTime fliplr(slideTime)], [coeff(:,i+1,2)' fliplr(coeff(:,i+1,3)')], colors(i,:), 'facealpha', 0.25, 'edgecolor', 'none')
end
line([slideTime(1) slideTime(end)], [0 0], 'color', [0.4 0.4 0.4], 'LineStyle','--');
if ~isempty(find(os.laser)==1)
    legend({'rwd', 'switch', 'hist', 'laser'});
else
    legend({'rwd', 'switch', 'hist'});
end

% plot laser difference

laser = [behSessionData.laser];
LRInd = intersect(os.responseInds(os.rwd_Inds), find(laser>0));
LNInd = intersect(os.responseInds(os.nrwd_Inds), find(laser>0));
NRInd = intersect(os.responseInds(os.rwd_Inds), find(laser==0));
NNInd = intersect(os.responseInds(os.nrwd_Inds), find(laser==0));

if ~isempty(LRInd) && ~isempty(NRInd)
    subplot(6,4,20); hold on;
    my_SDF_1 = sessionPupilCue(intersect(LRInd, find(qualInd>0)),:);
    my_SDF_2 = sessionPupilCue(intersect(NRInd, find(qualInd>0)),:);
    plotFilled(time, my_SDF_1, 'r');
    plotFilled(time, my_SDF_2, 'b');
    legend({'Laser', '', 'no-Laser', ''})
    title('Rwd')
end

if ~isempty(LNInd) && ~isempty(NNInd)
    subplot(6,4,24); hold on;
    my_SDF_1 = sessionPupilCue(intersect(LNInd, find(qualInd>0)),:);
    my_SDF_2 = sessionPupilCue(intersect(NNInd, find(qualInd>0)),:);
    plotFilled(time, my_SDF_1, 'r');
    plotFilled(time, my_SDF_2, 'b');
    legend({'Laser', '', 'no-Laser', ''})
    title('noRwd')
end
%% save
if p.Results.saveFigFlag == 1
    if isempty(dir(savePath))
        mkdir(savePath)
    end
    saveFigurePDF(gcf,[savePath session '_pupil.pdf'])
end










