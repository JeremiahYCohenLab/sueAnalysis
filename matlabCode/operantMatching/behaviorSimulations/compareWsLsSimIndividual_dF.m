function compareWsLsSimIndividual_dF(xlFile, sheet, category, beh, varargin)

p = inputParser;
% default parameters if none given
p.addParameter('rwdProbs', [90 50 10])
p.addParameter('maxTrials', 500)
p.addParameter('runs', 1) %runs per sample
p.addParameter('randomSeed', 318625)
p.addParameter('modelName', [{'fiveParam_bias'} {'sevenParam_absPePeAN_scale_int_bias_ord'}])
p.addParameter('sessionParamsFlag', 0)
p.addParameter('samps', 1000)
p.addParameter('taskType', 'decoupled')
p.addParameter('revForFlag', 0)
p.addParameter('tMax', 5)
p.parse(varargin{:});

[root, sep] = currComputer();

%initialize analysis parameters
rSeed = p.Results.randomSeed;
numMdls = length(p.Results.modelName);
tMax = p.Results.tMax;


%extract session list from excel sheet
[~, dayList, ~] = xlsread(xlFile, sheet);
[~,col] = find(strcmp(dayList, category));
dayList = dayList(2:end, col);
endInd = find(cellfun(@isempty,dayList),1);
if ~isempty(endInd)
    dayList = dayList(1:endInd-1,:);
end

prevAnimal = [];
animals = [];
for currSesh = 1:length(dayList)
    sessionName = dayList{currSesh};
    [animal, ~] = strtok(sessionName, 'd'); 
    animal = animal(2:end);
    if strcmp(animal, prevAnimal) == 0
        animals = [animals; {animal}];
    end
    prevAnimal = animal;
end
numAnimals = length(animals);

%initialize structures for actual ws-ls
wS = cell(numAnimals, 1);
lS = cell(numAnimals, 1);
wsRwdHx = cell(numAnimals, 1);
lsRwdHx = cell(numAnimals, 1);

%initialize structures for simulated ws-ls
wsSim = cell(numAnimals, numMdls);
lsSim = cell(numAnimals, numMdls);
wsRwdHxSim = cell(numAnimals, numMdls);
lsRwdHxSim = cell(numAnimals, numMdls);

prevAnimal = [];
currA = 0;
for currSesh = 1:length(dayList)
    sessionName = dayList{currSesh};
    [animal, ~] = strtok(sessionName, 'd'); 
    animal = animal(2:end);
    if strcmp(animal, prevAnimal) == 0 || p.Results.sessionParamsFlag
        if strcmp(animal, prevAnimal) == 0
            fprintf('Simulating animal %s \n', animal);
            currA = currA + 1;
            %get ws-ls from actual behavior
            [wsRwdHx{currA}, lsRwdHx{currA}, wS{currA}, lS{currA}] = wslsRwdCount_dF(xlFile, animal, category, 'tMax', p.Results.tMax, 'revForFlag', p.Results.revForFlag);
        end
        prevAnimal = animal;
        
        for currM = 1:numMdls
            modelPath = [root animal sep animal 'sorted' sep 'stan' sep 'bernoulli' sep p.Results.modelName{currM} sep animal...
                    beh '_' p.Results.modelName{currM} '.mat'];
            [t, ~] = getStanModelParams_samps(p.Results.modelName{currM}, modelPath, p.Results.samps,...
                            'sessionParamsFlag', p.Results.sessionParamsFlag, 'sessionName', sessionName);
            
            wsTmp = nan(p.Results.runs*p.Results.samps, 1);  
            lsTmp = nan(p.Results.runs*p.Results.samps, 1);  
            wsRwdHxTmp = nan(p.Results.runs*p.Results.samps, tMax+1); 
            lsRwdHxTmp = nan(p.Results.runs*p.Results.samps, tMax+1);
            for currS = 1:p.Results.samps
                sInd = (currS - 1) * p.Results.runs;
                for currR = 1:p.Results.runs
                    rSeed = rSeed + 1;
                    [~,allRewards, allChoices, blockProbs, blockSwitch] = runSim_dF(p.Results.modelName{currM}, t.params(currS,:),...
                            p.Results.maxTrials, rSeed, p.Results.rwdProbs);
                    allRewards = abs(allRewards(1:end-1));
                    changeChoice = [abs(diff(allChoices)) > 0];
                    rwdHx = movsum(allRewards, [tMax 0]);
                    rwdHx = [0 rwdHx(1:end-1)];

                    wsTmp(currR+sInd) = 1 - (sum(changeChoice(allRewards==1))/sum(allRewards==1)); 
                    lsTmp(currR+sInd) = sum(changeChoice(allRewards==0))/sum(allRewards==0);
                    for j = 1:tMax+1
                        tmpInds = logical(rwdHx == j-1);
                        wsRwdHxTmp(currR+sInd, j) = 1 - (sum(changeChoice(allRewards==1 & tmpInds))/sum(allRewards==1 & tmpInds));
                        lsRwdHxTmp(currR+sInd, j) = sum(changeChoice(allRewards==0 & tmpInds))/sum(allRewards==0 & tmpInds);
                    end 
                end
            end
            wsSim{currA, currM} = [wsSim{currA, currM}; mean(wsTmp)];
            lsSim{currA, currM} = [lsSim{currA, currM}; mean(lsTmp)];
            wsRwdHxSim{currA, currM} = [wsRwdHxSim{currA, currM}; nanmean(wsRwdHxTmp,1)];
            lsRwdHxSim{currA, currM} = [lsRwdHxSim{currA, currM}; nanmean(lsRwdHxTmp,1)];
        end
    end 
end

%convert actual data to matrix
wS = cell2mat(wS);
lS = cell2mat(lS);
wsRwdHx = cell2mat(wsRwdHx);
lsRwdHx = cell2mat(lsRwdHx);

%find within animal mean if running by session
if p.Results.sessionParamsFlag
    wrapper = @(x) nanmean(x, 1);
    wsSim = cellfun(wrapper, wsSim);
    lsSim = cellfun(wrapper, wsSim);
    wsRwdHxSim = cellfun(wrapper, wsSim);
    lsRwdHxSim = cellfun(wrapper, wsSim);
end

%organize animals into single matrix for each model
wsSim = cell2mat(wsSim);
lsSim = cell2mat(lsSim);
wsRwdHxSimTmp = cell2mat(wsRwdHxSim);
wsRwdHxSim = cell(1,numMdls);
lsRwdHxSimTmp = cell2mat(lsRwdHxSim);
lsRwdHxSim = cell(1,numMdls);

for currM = 1:numMdls
    wsRwdHxSim{currM} = wsRwdHxSimTmp(:,(tMax+1)*(currM-1)+1 : (tMax+1)*(currM));
    lsRwdHxSim{currM} = lsRwdHxSimTmp(:,(tMax+1)*(currM-1)+1 : (tMax+1)*(currM));
end

colors = cool(numMdls);

figure;
subplot(2,2,1); hold on;
x = [0:tMax];
plotFilled(x, wsRwdHx, 'k');
for currM = 1:numMdls
    plotFilled(x, wsRwdHxSim{currM}, colors(currM,:));
end
xlim([-0.5 tMax+0.5])
xlabel(['number of rewards in last ' num2str(tMax) ' trials'])
ylabel('probability')
title('win-stay')
set(gca, 'tickdir', 'out')

subplot(2,2,3); hold on;
plotFilled(x, lsRwdHx, 'k');
legTmp = [{'actual'} {' '}];
for currM = 1:numMdls
    plotFilled(x, lsRwdHxSim{currM}, colors(currM,:));
    legTmp = [legTmp strrep({p.Results.modelName{currM}}, '_', ' ') ' '];
end
xlim([-0.5 tMax+0.5])
xlabel(['number of rewards in last ' num2str(tMax) ' trials'])
ylabel('probability')
title('lose-shift')
set(gca, 'tickdir', 'out')
legend(legTmp)

subplot(2,2,2); hold on;
errorbar(0, mean(wS), sem(wS), 'o', 'color', [0 0 0], 'markerfacecolor', [0 0 0], 'markersize', 20)
tickLbls = {'actual'};
for currM = 1:numMdls
    errorbar(currM, mean(wsSim(:, currM)), sem(wsSim(:,currM)), 'o', 'color', colors(currM,:),...
        'markerfacecolor', colors(currM,:), 'markersize', 20)
    tickLbls = [tickLbls strrep({p.Results.modelName{currM}}, '_', ' ')];
end
xlim([-0.5 numMdls+0.5])
ylabel('win-stay')
xticks([0:numMdls])
xticklabels(tickLbls)
set(gca, 'tickdir', 'out')

subplot(2,2,4); hold on;
errorbar(0, mean(lS), sem(lS), 'o', 'color', [0 0 0], 'markerfacecolor', [0 0 0], 'markersize', 20)
for currM = 1:numMdls
    errorbar(currM, mean(lsSim(:,currM)), sem(lsSim(:,currM)), 'o', 'color', colors(currM,:),...
        'markerfacecolor', colors(currM,:), 'markersize', 20)
end
xlim([-0.5 numMdls+0.5])
ylabel('lose-shift')
xticks([0:numMdls])
xticklabels(tickLbls)
set(gca, 'tickdir', 'out')

set(gcf, 'renderer', 'painters', 'position', [-1919 41 1920 963])