function [glm_seshInd, t] = combineLogRegTimeLickRate_opMD(xlFile, animal, category, varargin)

%task and model parameters
p = inputParser;
% default parameters if none given
p.addParameter('revForFlag', 0)
p.addParameter('binSize', 5000)
p.addParameter('numBins', 10)
p.addParameter('plotFlag', 0);
p.addParameter('maxTrials', 350);
p.parse(varargin{:});

%determine root for file location
[root, sep] = currComputer();

%import behavior session titles for desired category
[~, dayList, ~] = xlsread([root xlFile], animal);
[~,col] = find(~cellfun(@isempty,strfind(dayList, category)) == 1);
dayList = dayList(2:end,col);
endInd = find(cellfun(@isempty,dayList),1);
if ~isempty(endInd)
    dayList = dayList(1:endInd-1,:);
end

numBins = p.Results.numBins;
binSize = p.Results.binSize;  %6300 for trial length
timeMax = binSize * numBins + 1000;
timeBinEdges = [1000:binSize:timeMax];  %no trials shorter than 1s between outcome and CS on
tMax = length(timeBinEdges) - 1;
rwdMatx =[];                            %initialize matrices for combining session data
noRwdMatx = [];
combinedAllChoice_R = [];
seshLength = [];

%loop for each session in the list
for i = 1: length(dayList)              
    sessionName = dayList{i};                       %extract relevant info from session title
    [animalName, date] = strtok(sessionName, 'd'); 
    animalName = animalName(2:end);
    date = date(1:9);
    sessionFolder = ['m' animalName date];

    if isstrprop(sessionName(end), 'alpha')         %define appropriate data path
        sessionDataPath = [root animalName sep sessionFolder sep 'sorted' sep 'session ' sessionName(end) sep sessionName '_sessionData_behav.mat'];
    else
        sessionDataPath = [root animalName sep sessionFolder sep 'sorted' sep 'session' sep sessionName '_sessionData_behav.mat'];
    end

    if exist(sessionDataPath,'file')        %load preprocessed struct if there is one
        load(sessionDataPath)
        if p.Results.revForFlag
            behSessionData = sessionData;
        end
    elseif p.Results.revForFlag                                    %otherwise generate the struct
        [behSessionData, ~] = generateSessionData_operantMatching(sessionName);
    else
        [behSessionData, ~, ~, ~] = generateSessionData_operantMatchingDecoupled(sessionName);
    end
    behSessionData = behSessionData(1:min(length(behSessionData),p.Results.maxTrials));
    %create arrays for choices and rewards
    responseInds = find(~isnan([behSessionData.rewardTime])); % find CS+ trials with a response in the lick window
    allReward_R = [behSessionData(responseInds).rewardR]; 
    allReward_L = [behSessionData(responseInds).rewardL]; 
    allChoices = NaN(1,length(behSessionData(responseInds)));
    allChoices(~isnan(allReward_R)) = 1;
    allChoices(~isnan(allReward_L)) = -1;
    allChoice_R = double(allChoices == 1);
    allChoice_L = double(allChoices == -1);
    
    allReward_R(isnan(allReward_R)) = 0;
    allReward_L(isnan(allReward_L)) = 0;
    allRewards = zeros(1,length(allChoices));
    allRewards(logical(allReward_R)) = 1;
    allRewards(logical(allReward_L)) = 1;
    
    

    %create binned outcome matrices
    rwdTmpMatx = NaN(tMax, length(responseInds));     %initialize matrices for number of response trials x number of time bins
    noRwdTmpMatx = NaN(tMax, length(responseInds));
    
    for j = 2:length(responseInds)          
        k = 1;
        %find time between "current" choice and previous rewards, up to timeMax in the past 
        timeTmpL = []; timeTmpR = []; nTimeTmpL = []; nTimeTmpR = [];
        while j-k > 0 && behSessionData(responseInds(j)).respondTime - behSessionData(responseInds(j-k)).rewardTime < timeMax
            if behSessionData(responseInds(j-k)).rewardL == 1
                timeTmpL = [timeTmpL (behSessionData(responseInds(j)).respondTime - behSessionData(responseInds(j-k)).rewardTime)];
            end
            if behSessionData(responseInds(j-k)).rewardR == 1
                timeTmpR = [timeTmpR (behSessionData(responseInds(j)).respondTime - behSessionData(responseInds(j-k)).rewardTime)];
            end
            if behSessionData(responseInds(j-k)).rewardL == 0
                nTimeTmpL = [nTimeTmpL (behSessionData(responseInds(j)).respondTime - behSessionData(responseInds(j-k)).rewardTime)];
            end
            if behSessionData(responseInds(j-k)).rewardR == 0
                nTtimeTmpR = [nTimeTmpR (behSessionData(responseInds(j)).respondTime - behSessionData(responseInds(j-k)).rewardTime)];
            end
            k = k + 1;
        end
        %bin outcome times and use to fill matrices
        if ~isempty(timeTmpL)
            binnedRwds = discretize(timeTmpL,timeBinEdges);
            for k = 1:tMax
                if ~isempty(binnedRwds == k)
                    rwdTmpMatx(k,j) = -1*sum(binnedRwds == k);
                else
                    rwdTmpMatx(k,j) = 0;
                end
            end
        end
        if ~isempty(timeTmpR)
            binnedRwds = discretize(timeTmpR,timeBinEdges);
            for k = 1:tMax
                if ~isempty(binnedRwds == k) & isnan(rwdTmpMatx(k,j))
                    rwdTmpMatx(k,j) = sum(binnedRwds == k);
                elseif ~isempty(binnedRwds == k) & ~isnan(rwdTmpMatx(k,j))
                    rwdTmpMatx(k,j) = rwdTmpMatx(k,j) + sum(binnedRwds == k);
                else
                    rwdTmpMatx(k,j) = 0;
                end
            end
        end
        if isempty(timeTmpL) & isempty(timeTmpR)
            rwdTmpMatx(:,j) = 0;
        end
        if ~isempty(nTimeTmpL)
            binnedNoRwds = discretize(nTimeTmpL,timeBinEdges);
            for k = 1:tMax
                if ~isempty(binnedNoRwds == k)
                    noRwdTmpMatx(k,j) = -1*sum(binnedNoRwds == k);
                else
                    noRwdTmpMatx(k,j) = 0;
                end
            end
        end
        if ~isempty(nTimeTmpR)
            binnedNoRwds = discretize(nTimeTmpR,timeBinEdges);
            for k = 1:tMax
                if ~isempty(binnedNoRwds == k) & isnan(noRwdTmpMatx(k,j))
                    noRwdTmpMatx(k,j) = sum(binnedNoRwds == k);
                elseif ~isempty(binnedNoRwds == k) & ~isnan(noRwdTmpMatx(k,j))
                    noRwdTmpMatx(k,j) = noRwdTmpMatx(k,j) + sum(binnedNoRwds == k);
                else
                    noRwdTmpMatx(k,j) = 0;
                end
            end
        end
        if isempty(nTimeTmpL) & isempty(nTimeTmpR)
            noRwdTmpMatx(:,j) = 0;
        end
    end
    
    %concatenate temp matrix with combined matrix
    rwdTmpMatx(:,1) = NaN;
    rwdMatx = [rwdMatx NaN(length(timeBinEdges)-1, 100) rwdTmpMatx];
    noRwdMatx = [noRwdMatx NaN(length(timeBinEdges)-1, 100) noRwdTmpMatx];
    combinedAllChoice_R = [combinedAllChoice_R NaN(1,100) allChoice_R];
    
    seshLength = [seshLength length(rwdTmpMatx)];
end

seshInd = zeros(length(rwdMatx), length(dayList));
tmpInd = 101;
for i = 1:length(dayList)-1
    seshInd(tmpInd:tmpInd+seshLength(i)-1,i) = 1;
    tmpInd = tmpInd + seshLength(i) + 100;
end


%logistic regression models
glm_rwd = fitglm([rwdMatx]', combinedAllChoice_R,'distribution','binomial','link','logit');
glm_noRwd = fitglm([noRwdMatx]', combinedAllChoice_R,'distribution','binomial','link','logit');
glm_all = fitglm([rwdMatx' noRwdMatx'], combinedAllChoice_R,'distribution','binomial','link','logit');
glm_seshInd = fitglm([rwdMatx' noRwdMatx', seshInd], combinedAllChoice_R,'distribution','binomial','link','logit');
glm_auto = fitglm([rwdMatx+noRwdMatx]', combinedAllChoice_R,'distribution','binomial','link','logit');

t = struct;
t.binSize = binSize;
t.timeMax = timeMax;
t.tMax = tMax;
t.timeBinEdges = timeBinEdges;

if p.Results.plotFlag
    figure; hold on;
    relevInds = 2:tMax+1;
    coefVals = glm_all.Coefficients.Estimate(relevInds);
    CIbands = coefCI(glm_all);
    errorL = abs(coefVals - CIbands(relevInds,1));
    errorU = abs(coefVals - CIbands(relevInds,2));
    errorbar(((1:tMax)*binSize/1000),coefVals,errorL,errorU,'Color', [0.7 0 1],'linewidth',2)

    relevInds = tMax+2:2*tMax+1;
    coefVals = glm_all.Coefficients.Estimate(relevInds);
    CIbands = coefCI(glm_all);
    errorL = abs(coefVals - CIbands(relevInds,1));
    errorU = abs(coefVals - CIbands(relevInds,2));
    errorbar(((1:tMax)*binSize/1000),coefVals,errorL,errorU,'b','linewidth',2)
    line([0 tMax*binSize/1000 + binSize/1000], [0 0], 'Color','red','LineStyle','--')

    legend('Reward', 'No Reward')
    xlabel('Outcome n seconds back')
    ylabel('\beta Coefficient')
    xlim([0 tMax*binSize/1000 + binSize/1000])
    title([animal ' ' category])
    set(gca, 'tickdir', 'out')
    set(gcf, 'renderer', 'painters')
end