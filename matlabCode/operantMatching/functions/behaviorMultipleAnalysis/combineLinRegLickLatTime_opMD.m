function [glm_rwdLick, rwdMatx, noRwdMatx, combinedPreLick, combinedLickLatZ] = combineLinRegLickLatTime_opMD(xlFile, animal, category, varargin)

p = inputParser;
% default parameters if none given
p.addParameter('revForFlag',0)
p.addParameter('plotFlag', 0)
p.addParameter('binSize', 10000);
p.addParameter('numBins', 20);
p.addParameter('maxTrials', 1000);
p.parse(varargin{:});

[root, sep] = currComputer();

[~, dayList, ~] = xlsread([root xlFile], animal);
[~,col] = find(~cellfun(@isempty,strfind(dayList, category)) == 1);
dayList = dayList(2:end,col);
endInd = find(cellfun(@isempty,dayList),1);
if ~isempty(endInd)
    dayList = dayList(1:endInd-1,:);
end

binSize = p.Results.binSize;
timeMax = binSize * p.Results.numBins + 1000;
timeBinEdges = [1000:binSize:timeMax];  %no trials shorter than 1s between outcome and CS on
tMax = length(timeBinEdges) - 1;
rwdMatx = [];
noRwdMatx = [];
combinedLickLat = [];
combinedLickLatZ = [];
stayLickLat = [];
switchLickLat = [];
exploreLickLat = [];
exploitLickLat = [];
combinedITIlicks = [];
combinedTimeInSesh = [];
combinedChangeChoice = [];
combinedPreLick = [];
preITI = [];
latAndTime = [];


for i = 1: length(dayList)
    sessionName = dayList{i};
    clear behSessionData
    [animalName, date] = strtok(sessionName, 'd'); 
    animalName = animalName(2:end);
    date = date(1:9);
    sessionFolder = ['m' animalName date];

    if isstrprop(sessionName(end), 'alpha')
        sessionDataPath = [root animalName sep sessionFolder sep 'sorted' sep 'session ' sessionName(end) sep sessionName '_sessionData_behav.mat'];
    else
        sessionDataPath = [root animalName sep sessionFolder sep 'sorted' sep 'session' sep sessionName '_sessionData_behav.mat'];
    end

    if exist(sessionDataPath,'file')
%         fprintf([sessionName '\n']);
        load(sessionDataPath)
    else
        fprintf([sessionName '\n']);
        [behSessionData, ~] = generateSessionData_operantMatchingDecoupledRwdDelay(sessionName);        
    end
    
    if ~exist('behSessionData', 'var')
        behSessionData = sessionData;
    end
    responseInds = find(~isnan([behSessionData.rewardTime])); 
    responseInds = responseInds(1:min(length(responseInds),p.Results.maxTrials));   
    behSessionData = behSessionData(1:responseInds(end));
    %%generate reward matrix for tMax trials
    allReward_R = [behSessionData(responseInds).rewardR]; 
    allReward_L = [behSessionData(responseInds).rewardL]; 
    allChoices = NaN(1,length(behSessionData(responseInds)));
    allChoices(~isnan(allReward_R)) = 1;
    allChoices(~isnan(allReward_L)) = -1;
    allReward_R(isnan(allReward_R)) = 0;
    allReward_L(isnan(allReward_L)) = 0;
    allRewards = zeros(1,length(allChoices));
    allRewards(logical(allReward_R)) = 1;
    allRewards(logical(allReward_L)) = 1;
    timeInSesh = ([behSessionData(responseInds).CSon] - behSessionData(1).CSon) / (1000 * behSessionData(responseInds(end)).CSon - behSessionData(responseInds(1)).CSon);
    changeChoice = [0 abs(diff(allChoices)) > 0];
%     fprintf([sessionName '\n'])
    s = behAnalysisNoPlot_opMD(sessionName);
    hmmStates = s.hmmStates(1:length(responseInds));
    timeBtwn = s.timeBtwn(1:length(responseInds));
    %% determine lick latency distributions for each spout
    lickLat = [behSessionData(responseInds).respondTime] - [behSessionData(responseInds).CSon];
    lickLatLog = log(lickLat);
    indsR = find(allChoices == 1);
    indsL = find(allChoices == -1);
    lickLat_R = zscore(lickLat(indsR));
    lickLat_L = zscore(lickLat(indsL));
    lickLatZ = NaN(1, length(allChoices));
    lickLatZ(indsR) = lickLat_R;
    lickLatZ(indsL) = lickLat_L;
    
    %combinedLickLat = [combinedLickLat NaN(1,101) lickLatZ(2:end)];
    combinedLickLat = [combinedLickLat NaN(1,100) lickLat];
    combinedLickLatZ = [combinedLickLatZ NaN(1,100) lickLatZ];
    preITI = [preITI NaN(1,100) timeBtwn];
    %% create binned outcome matrices
    rwdTmpMatx = zeros(tMax, length(responseInds));     %initialize matrices for number of response trials x number of time bins
    noRwdTmpMatx = zeros(tMax, length(responseInds));
    preLickTmp = NaN(1,length(responseInds));
    preC = NaN(1,length(responseInds));
    for j = 2:length(responseInds)          
        k = 1;
        %find time between "current" choice and previous rewards, up to timeMax in the past 
        timeTmp = [];
        timeTmpNoRwd = [];
        while j-k > 0 & behSessionData(responseInds(j)).respondTime - behSessionData(responseInds(j-k)).rewardTime < timeMax
            if behSessionData(responseInds(j-k)).rewardL == 1 || behSessionData(responseInds(j-k)).rewardR == 1
                timeTmp = [timeTmp (behSessionData(responseInds(j)).respondTime - behSessionData(responseInds(j-k)).rewardTime)];
            end
            if behSessionData(responseInds(j-k)).rewardL == 0 || behSessionData(responseInds(j-k)).rewardR == 0
                timeTmpNoRwd = [timeTmpNoRwd (behSessionData(responseInds(j)).respondTime - behSessionData(responseInds(j-k)).rewardTime)];
            end
            k = k + 1;
        end
        %bin outcome times and use to fill matrices
        if ~isempty(timeTmp)
            binnedRwds = discretize(timeTmp,timeBinEdges);
            for k = 1:tMax
                if ~isempty(find(binnedRwds == k, 1))
                    rwdTmpMatx(k,j) = sum(binnedRwds == k);
                end
            end
        end
        if ~isempty(timeTmpNoRwd)
            binnedNoRwds = discretize(timeTmpNoRwd,timeBinEdges);
            for k = 1:tMax
                if ~isempty(find(binnedNoRwds == k, 1))
                    noRwdTmpMatx(k,j) = sum(binnedNoRwds == k);
                end
            end
        end
        % The last lick binSize ago. 
        m = 1;
        while j-m > 1 && behSessionData(responseInds(j)).respondTime - behSessionData(responseInds(j-m)).rewardTime < p.Results.binSize
            m = m + 1;
        end
        if allChoices(j-m) == 1
            preC(j) = 1;
        else
            preC(j) = -1;
        end
        preLickTmp(j) = lickLatZ(j-m);
    end
    
    %fill in NaNs at beginning of session
    j = 2;
    while behSessionData(responseInds(j)).respondTime - behSessionData(responseInds(1)).respondTime < timeMax
        tmpDiff = behSessionData(responseInds(j)).respondTime - behSessionData(responseInds(1)).respondTime;
        binnedDiff = discretize(tmpDiff, timeBinEdges);
        rwdTmpMatx(binnedDiff+1:tMax,j) = NaN;
        noRwdTmpMatx(binnedDiff+1:tMax,j) = NaN;
        j = j+1;
    end
    %concatenate temp matrix with combined matrix
    rwdTmpMatx(:,1) = NaN;
    rwdMatx = [rwdMatx NaN(length(timeBinEdges)-1, 100) rwdTmpMatx];
    noRwdTmpMatx(:,1) = NaN;
    noRwdMatx = [noRwdMatx NaN(length(timeBinEdges)-1, 100) noRwdTmpMatx];
    combinedChangeChoice = [combinedChangeChoice NaN(1, 100) changeChoice];
    combinedTimeInSesh = [combinedTimeInSesh NaN(1, 100) timeInSesh];
    combinedPreLick = [combinedPreLick NaN(1, 100) preLickTmp];
    %% determine lick latency for stay v switch trials; explore vs exploit trails
    changeChoice = [false abs(diff(allChoices)) > 0];
    stayLickLat = [stayLickLat lickLat(~changeChoice)]; 
    switchLickLat = [switchLickLat lickLat(changeChoice)];
    exploitLickLat = [exploitLickLat lickLat(hmmStates~=1)];
    exploreLickLat = [exploreLickLat lickLat(hmmStates==1)];
    
    %% combine lick latency and time in session
    tmp = [lickLatZ; timeInSesh];
    latAndTime = [latAndTime tmp];
    
    %% determine iti lick rates on each trial
    itiLicks = zeros(1,length(behSessionData));
    for currT = 1:length(behSessionData)
        itiLicksTmp = [behSessionData(currT).licksL(behSessionData(currT).licksL > behSessionData(currT).CSon + 2000)...
                        behSessionData(currT).licksR(behSessionData(currT).licksR > behSessionData(currT).CSon + 2000)];
        if ~isempty(itiLicksTmp)
            itiLicks(currT) = length(itiLicksTmp) / (behSessionData(currT).trialEnd - behSessionData(currT).CSon - 2000);
        end
    end
    combinedITIlicks = [combinedITIlicks itiLicks];
end

%linear regression model
glm_rwdLick = fitlm([rwdMatx'], combinedLickLatZ);
glm_rwdLickAll = fitlm([rwdMatx' noRwdMatx' combinedTimeInSesh' combinedChangeChoice' preITI'], combinedLickLatZ);
tbl = table(combinedPreLick', combinedChangeChoice', rwdMatx(1,:)', noRwdMatx(1,:)', preITI',combinedLickLatZ', 'VariableNames', {'pre', 'sw', 'rwd1', 'nRwd1', 'preITI', 'lickLat'});
mdl = stepwiselm(tbl,'interactions');

if p.Results.plotFlag
    figure2('Position', [1 1 800 800]); suptitle([animal ' ' category])
    colors = cool(5);
    subplot(2,2,1); hold on
    relevInds = 2:tMax+1;
    coefVals = glm_rwdLickAll.Coefficients.Estimate(relevInds);
    CIbands = coefCI(glm_rwdLickAll);
    errorL = abs(coefVals - CIbands(relevInds,1));
    errorU = abs(coefVals - CIbands(relevInds,2));
    errorbar(((1:tMax)*binSize/1000),coefVals,errorL,errorU,'Color', colors(1,:),'linewidth',2)
    
    relevInds = tMax+2:2*tMax+1;
    coefVals = glm_rwdLickAll.Coefficients.Estimate(relevInds);
    CIbands = coefCI(glm_rwdLickAll);
    errorL = abs(coefVals - CIbands(relevInds,1));
    errorU = abs(coefVals - CIbands(relevInds,2));
    errorbar(((1:tMax)*binSize/1000),coefVals,errorL,errorU,'Color', colors(4,:),'linewidth',2)
    line([0 tMax*binSize/1000], [0 0], 'Color',[0.5 0.5 0.5],'LineStyle','--')
    
    legend('reward', 'no reward')
    xlabel('reward n seconds back')
    ylabel('\beta Coefficient')
    xlim([0 (tMax*binSize/1000 + 5)])
    set(gca, 'tickdir', 'out')
        
    subplot(4,2,2); hold on
    histogram(stayLickLat,min(lickLat)-25:25:max(lickLat)+25,'FaceColor', colors(2,:), 'Normalization', 'Probability')
    histogram(switchLickLat,min(lickLat)-25:25:max(lickLat)+25,'FaceColor', colors(5,:), 'Normalization', 'Probability')
    ylabel('probability')
    xlabel('lick latency')
    legend('stay', 'switch')

    subplot(4,2,4); hold on
    histogram(exploitLickLat,min(lickLat)-25:25:max(lickLat)+25,'FaceColor', colors(2,:), 'Normalization', 'Probability')
    histogram(exploreLickLat,min(lickLat)-25:25:max(lickLat)+25,'FaceColor', colors(5,:), 'Normalization', 'Probability')
    ylabel('probability')
    xlabel('lick latency')
    legend('exploit', 'explore')
    
    subplot(2,2,3); hold on
    numBins = 6;
    sortInds = discretize(latAndTime(2,:), numBins);
    for currBin = 1:numBins
        tmp = latAndTime(1, sortInds==currBin);
        meanLat(currBin) = mean(tmp);
        semLat(currBin) = std(tmp) / sqrt(length(tmp));
    end
    errorbar([1:numBins], meanLat, semLat, 'Color', colors(3,:),'linewidth',2);
    xlim([0 numBins+1]);
    xticks([])
    xlabel('normalized time in session')
    ylabel('z-scored lick latency')
    
    subplot(2,2,4); hold on;
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
    title('lrm: on lickLat')
    ylabel('\beta Coefficient')
    text(length(coefVals)-0.5,0.8*max(CIbands(2:end,2)),sprintf('R^2 = %d',mdl.Rsquared.Adjusted))
    hold off
    
end