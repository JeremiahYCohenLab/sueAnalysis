function [glm_rwdLick, stayLickLat, switchLickLat, binSize, timeMax, combinedITIlicks] = combineLogRegNrwdSwitchTime_opMD(animal, category, varargin)

p = inputParser;
% default parameters if none given
p.addParameter('revForFlag',0)
p.addParameter('plotFlag', 0)
p.parse(varargin{:});

[root, sep] = currComputer();
workbookFile = [root 'combineAnimals'];
[~, dayList, ~] = xlsread(workbookFile, animal);
col = contains(dayList(1,:),category);
dayList = dayList(2:end,col);
endInd = find(cellfun(@isempty,dayList),1);
if ~isempty(endInd)
    dayList = dayList(1:endInd-1,:);
end

timeMax = 121000;
binSize = 10000;
timeBinEdges = [1000:binSize:timeMax];  %no trials shorter than 1s between outcome and CS on
tMax = length(timeBinEdges) - 1;
rwdMatx = [];
svs = [];
eve = [];
preITIs = [];

for i = 1: length(dayList)
    sessionName = dayList{i};
    [animalName, date] = strtok(sessionName, 'd'); 
    animalName = animalName(2:end);
    date = date(1:9);
    sessionFolder = ['m' animalName date];

    if isstrprop(sessionName(end), 'alpha')
        sessionDataPath = [root animalName sep sessionFolder sep 'sorted' sep 'session' sep sessionName(end) sep sessionName '_sessionData_behav.mat'];
    else
        sessionDataPath = [root animalName sep sessionFolder sep 'sorted' sep 'session' sep sessionName '_sessionData_behav.mat'];
    end

    if exist(sessionDataPath,'file')
        load(sessionDataPath)
        if p.Results.revForFlag
            behSessionData = sessionData;
        end
    else
        [behSessionData, ~] = generateSessionData_operantMatchingDecoupled(sessionName);
    end
    
    %%generate reward matrix for tMax trials
    responseInds = find(~isnan([behSessionData.rewardTime])); % find CS+ trials with a response in the lick window
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

    %create binned outcome matrices
    rwdTmpMatx = zeros(tMax, length(responseInds));     %initialize matrices for number of response trials x number of time bins
    for j = 2:length(responseInds)          
        k = 1;
        %find time between "current" choice and previous rewards, up to timeMax in the past 
        timeTmp = [];
        while j-k > 0 & behSessionData(responseInds(j)).rewardTime - behSessionData(responseInds(j-k)).rewardTime < timeMax
            if behSessionData(responseInds(j-k)).rewardL ~= 1 && behSessionData(responseInds(j-k)).rewardR ~= 1
                timeTmp = [timeTmp (behSessionData(responseInds(j)).rewardTime - behSessionData(responseInds(j-k)).rewardTime)];
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
    end
    
    %fill in NaNs at beginning of session
    j = 2;
    while behSessionData(responseInds(j)).rewardTime - behSessionData(responseInds(1)).rewardTime < timeMax
        tmpDiff = behSessionData(responseInds(j)).rewardTime - behSessionData(responseInds(1)).rewardTime;
        binnedDiff = discretize(tmpDiff, timeBinEdges);
        rwdTmpMatx(binnedDiff:tMax,j) = NaN;
        j = j+1;
    end
    %concatenate temp matrix with combined matrix
    rwdTmpMatx(:,1) = NaN;
    rwdMatx = [rwdMatx NaN(length(timeBinEdges)-1, 100) rwdTmpMatx];
    
    
%% stay v switch trials
    svstemp = [NaN 0.5*abs(diff(allChoices))];
    svs = [svs NaN(1, 100) svstemp];
  % ore vs oit
  %  [behSessionData, states,] = fitHmmOpt(dayList{i});
    hmm = [behSessionData(responseInds).hmm];
    evetemp = [behSessionData(responseInds).hmm];
    evetemp(evetemp ~= 1) = 0;
    eve = [eve NaN(1, 100) evetemp];
    
%% ITI
    allITItemp = [behSessionData(responseInds(1:end-1) + 1).CSon] - [behSessionData(responseInds(1:end-1)).CSon];
    allITItemp = [allITItemp  20000];
    preITItemp = [NaN, allITItemp(1:end-1)]; 
    preITIs = [preITIs NaN(1, 100) preITItemp];
end
%%
preITIs = 0.0005 * preITIs;

%regression model
glm_all = fitglm([rwdMatx',preITIs'], svs, 'distribution','binomial','link','logit');
%glm_all = fitglm(rwdMatx', svs, 'distribution','binomial','link','logit');
if p.Results.plotFlag
    figure; hold on;
    relevInds = 2:tMax+1;
    coefVals = glm_all.Coefficients.Estimate(relevInds);
    CIbands = coefCI(glm_all);
    errorL = abs(coefVals - CIbands(relevInds,1));
    errorU = abs(coefVals - CIbands(relevInds,2));
    errorbar(((1:tMax)*binSize/1000),coefVals,errorL,errorU,'Color',[0.2 0.2 1],'linewidth',2)
    bar((tMax*binSize/1000 + 20),glm_all.Coefficients.Estimate(end), 5, 'FaceColor',[0.2 0.2 1]);
    errorbar((tMax*binSize/1000 + 20),glm_all.Coefficients.Estimate(end), abs(glm_all.Coefficients.Estimate(end) - CIbands(end,1)), abs(glm_all.Coefficients.Estimate(end) - CIbands(end,2)),'Color', [0.2 0.2 0.2],'linewidth',1)

    xlabel('nreward n seconds back')
    ylabel('\beta Coefficient')
    xlim([0 (tMax*binSize/1000 + 25)])
    line([0 (tMax*binSize/1000 + 15)], [0 0], 'Color','k','LineStyle','--')
    suptitle([animal ' ' category])
    set(gca, 'tickdir', 'out')
    
end