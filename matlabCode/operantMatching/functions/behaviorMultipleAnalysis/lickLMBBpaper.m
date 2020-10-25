function [rwdMatx, lickMatx, combinedLickLat, sessionID, aniID] = lickLMBBpaper(xlFile,animal, category, varargin)
%task and model parameters
p = inputParser;
% default parameters if none given
p.addParameter('plotFlag', 1);
p.addParameter('maxTrials', 600);
p.addParameter('numTrial', 10);
p.parse(varargin{:});

numTrial = p.Results.numTrial;
maxTrial = p.Results.maxTrials;
colors = cool(5);
[root, ~] = currComputer();

[~, dayList, ~] = xlsread([root xlFile], animal);
col = contains(dayList(1,:),category);
dayList = dayList(2:end,col);
endInd = find(cellfun(@isempty,dayList),1);
if ~isempty(endInd)
    dayList = dayList(1:endInd-1,:);
end

combinedLickLat = [];
lickMatx = [];
rwdMatx = [];
combinedTrialNumber = [];
aniID = [];
sessionLen = zeros(size(dayList));
for i = 1:length(dayList)
    % load data
    [animalName, ~] = strtok(dayList{i}, 'd'); 
    animalName = animalName(2:end);
    behSessionData = loadBehavioralData([dayList{i} '.asc']);
    behavStruct = parseBehavioralData(behSessionData,maxTrial);
    outcome = behavStruct.allRewards(1:min(maxTrial,length(behavStruct.allChoices)));
    outcome = abs(outcome);
    lickLat = behavStruct.lickLatMlog(1:min(maxTrial,length(behavStruct.allChoices)));
    lickMatTemp = NaN(numTrial,length(lickLat));
    rwdMatTemp = NaN(numTrial,length(lickLat));
    for j = 1:numTrial
        lickMatTemp(j,j+1:end) = lickLat(1:end-j);
        rwdMatTemp(j,j+1:end) = outcome(1:end-j);       
    end
    aniTemp = cell(1,length(lickLat));
    aniTemp(:) = {animalName};
    lickMatx = [lickMatx lickMatTemp];
    rwdMatx = [rwdMatx rwdMatTemp];
    combinedLickLat = [combinedLickLat lickLat];
    combinedTrialNumber = [combinedTrialNumber 1:length(lickLat)];
    aniID = [aniID aniTemp];
    sessionLen(i)=length(lickLat);
end
% session ID vector
sessionID = zeros(1,length(combinedLickLat));
for i = 1:length(dayList)
    if i > 1
        sessionID(1,sum(sessionLen(1:i-1))+1:sum(sessionLen(1:i-1))+sessionLen(i))=i;
    else
        sessionID(1,1:sessionLen(i))=1;
    end

end

glm_LickRwd = fitlm([rwdMatx' lickMatx'], combinedLickLat);
glm_Rwd = fitlm([rwdMatx', combinedTrialNumber'], combinedLickLat);
%% plot everything
c = 2;
b = 3;
a = 2;
figure2('position',[0 0 1200 400]); 
subplot(1,3,1); hold on;
CIbands = coefCI(glm_LickRwd);
relevInds = 2:numTrial+1;
coefVals = glm_LickRwd.Coefficients.Estimate(relevInds);
errorL = abs(coefVals - CIbands(relevInds,1));
errorU = abs(coefVals - CIbands(relevInds,2));
errorbar(1:numTrial,coefVals,errorL,errorU,'Color', colors(c,:),'linewidth',2)
line([0 numTrial+1],[0 0], 'Color', [0.5 0.5 0.5], 'LineStyle','--');
xlim([0 numTrial+1])
xlabel('Trials')
ylabel('\beta Coefficient')
title('rewards')

subplot(1,3,2); hold on;
relevInds = numTrial + 2:2*numTrial+1;
coefVals = glm_LickRwd.Coefficients.Estimate(relevInds);
errorL = abs(coefVals - CIbands(relevInds,1));
errorU = abs(coefVals - CIbands(relevInds,2));
errorbar(1:numTrial,coefVals,errorL,errorU,'Color', colors(c,:),'linewidth',2)
line([0 numTrial+1],[0 0], 'Color', [0.5 0.5 0.5], 'LineStyle','--');
xlim([0 numTrial+1])
xlabel('Trials')
ylabel('\beta Coefficient')
title('lickLats')

% subplot(2,2,3); hold on;
% CIbands = coefCI(glm_preLickRwd);
% relevInds = 2:numTrial;
% coefVals = glm_preLickRwd.Coefficients.Estimate(relevInds);
% errorL = abs(coefVals - CIbands(relevInds,1));
% errorU = abs(coefVals - CIbands(relevInds,2));
% errorbar(1:numTrial-1,coefVals,errorL,errorU,'Color', colors(5,:),'linewidth',2)
% line([0 numTrial+1],[0 0], 'Color', [1 0 0], 'LineStyle','--');
% xlim([0 numTrial+1])
% xlabel('Trials')
% ylabel('\beta Coefficient')
% title('rewards')
% suptitle([animal '  ' category])

subplot(1,3,3); hold on;
CIbands = coefCI(glm_LickRwd);
int = glm_LickRwd.Coefficients.Estimate(1);
bar(a, int, 'FaceColor', colors(b,:), 'EdgeColor', colors(c,:))
errorL = abs(int - CIbands(1,1));
errorU = abs(int - CIbands(1,2));
errorbar(a,int,errorL,errorU,'Color', colors(c,:),'linewidth',2)
xlim([0 3])
xlabel('Int')
ylabel('\beta Coefficient')
title('Intercept')
suptitle([animal '  ' category])





























