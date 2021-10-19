function compareBlockSimSamp_dF(xlFile, sheet, category, varargin)

p = inputParser;
% default parameters if none given
p.addParameter('rwdProbs', [90 50 10])
p.addParameter('samps', 1000)
p.addParameter('randomSeed', 98773)
p.addParameter('modelName', 'sevenParam_absPePeAN_scale_int_bias_ord')
p.addParameter('mdlBeh', 'clean')
p.addParameter('bernFlag', 1)
p.addParameter('sessionParamsFlag', 1)
p.parse(varargin{:});

[root, sep] = currComputer();

if p.Results.sessionParamsFlag
    biasFlag = 1;
else
    biasFlag = 0;
end

pHigh = p.Results.rwdProbs(1);
pMed = p.Results.rwdProbs(2);
pLow = p.Results.rwdProbs(3);


%extract session list from excel sheet
[~, dayList, ~] = xlsread(xlFile, sheet);
[~,col] = find(strcmp(dayList, category));
dayList = dayList(2:end, col);
endInd = find(cellfun(@isempty,dayList),1);
if ~isempty(endInd)
    dayList = dayList(1:endInd-1,:);
end

numA = 0;
prevAnimal = [];
for currSesh = 1:length(dayList)
    sessionName = dayList{currSesh};
    [animal, ~] = strtok(sessionName, 'd'); 
    animal = animal(2:end);
    if strcmp(animal, prevAnimal) == 0
        numA = numA + 1;
        prevAnimal = animal;
    end
end

choices = cell(numA, 2);
probChoice = cell(numA, 2);
prevRwds = cell(numA, 2);
changeChoice = cell(numA, 2);
wS = cell(numA, 2);
lS = cell(numA, 2);

rSeed = p.Results.randomSeed;

prevAnimal = [];
aInd = 0;
for currSesh = 1:length(dayList)
    sessionName = dayList{currSesh};
    [animal, ~] = strtok(sessionName, 'd'); 
    animal = animal(2:end);
    if strcmp(animal, prevAnimal) == 0 
        fprintf('Analyzing animal %s \n', animal);
        aInd = aInd + 1;
        prevAnimal = animal;
    end

    %get model variables
    if p.Results.bernFlag
        modelPath = [root animal sep animal 'sorted' sep 'stan' sep 'bernoulli' sep p.Results.modelName sep animal...
        p.Results.mdlBeh '_' p.Results.modelName '.mat'];
    else
        modelPath = [root animal sep animal 'sorted' sep 'stan' sep p.Results.modelName sep animal...
        p.Results.mdlBeh '_' p.Results.modelName '.mat'];
    end
    t = getStanModelParams_samps(p.Results.modelName, modelPath, p.Results.samps, 'sessionParamsFlag', p.Results.sessionParamsFlag,...
        'sessionName', sessionName, 'varFlag', 1, 'biasFlag', biasFlag);
    
    %get session behavior data
    [behSessionData, blockSwitch, ~] = loadBehavioralData([sessionName '.asc'], 0);
    o = parseBehavioralData(behSessionData, blockSwitch);
    
    prevRwdsTmp = [0 abs(o.allRewards(1:end-1))];
    changeChoiceTmp = [0 abs(diff(o.allChoices)) > 0];
    
    rwdProb_R = [behSessionData(o.responseInds).rewardProbR]; 
    rwdProb_L = [behSessionData(o.responseInds).rewardProbL];
    
    blockSwitch = [o.blockSwitch length(o.responseInds)];
    for currB = 1:length(blockSwitch)-1
        sInd = blockSwitch(currB); eInd = blockSwitch(currB+1); 
        if (eInd - sInd) > 99
            eInd = sInd + 99;
            buff = 0;
        else
            buff = 99 - (eInd - sInd);
        end
        if rwdProb_L(sInd) == pHigh & rwdProb_R(sInd) == pLow
            choices{aInd, 2} = [choices{aInd, 2}; o.allChoice_L(sInd:eInd) nan(1, buff)];
            probChoice{aInd, 2} = [probChoice{aInd, 2}; t.probChoice(sInd:eInd, 1)' nan(1, buff)];
            prevRwds{aInd, 2} = [prevRwds{aInd, 2}; prevRwdsTmp(sInd:eInd) nan(1, buff)];
            changeChoice{aInd, 2} = [changeChoice{aInd, 2}; changeChoiceTmp(sInd:eInd) nan(1, buff)];    
        elseif rwdProb_L(sInd) == pMed & rwdProb_R(sInd) == pLow
            choices{aInd, 1} = [choices{aInd, 1}; o.allChoice_L(sInd:eInd) nan(1, buff)];
            probChoice{aInd, 1} = [probChoice{aInd, 1}; t.probChoice(sInd:eInd, 1)' nan(1, buff)];
            prevRwds{aInd, 1} = [prevRwds{aInd, 1}; prevRwdsTmp(sInd:eInd) nan(1, buff)];
            changeChoice{aInd, 1} = [changeChoice{aInd, 1}; changeChoiceTmp(sInd:eInd) nan(1, buff)];
        elseif rwdProb_R(sInd) == pHigh & rwdProb_L(sInd) == pLow
            choices{aInd, 2} = [choices{aInd, 2}; o.allChoice_R(sInd:eInd) nan(1, buff)];
            probChoice{aInd, 2} = [probChoice{aInd, 2}; t.probChoice(sInd:eInd, 2)' nan(1, buff)];
            prevRwds{aInd, 2} = [prevRwds{aInd, 2}; prevRwdsTmp(sInd:eInd) nan(1, buff)];
            changeChoice{aInd, 2} = [changeChoice{aInd, 2}; changeChoiceTmp(sInd:eInd) nan(1, buff)];
        elseif rwdProb_R(sInd) == pMed & rwdProb_L(sInd) == pLow
            choices{aInd, 1} = [choices{aInd, 1}; o.allChoice_R(sInd:eInd) nan(1, buff)];
            probChoice{aInd, 1} = [probChoice{aInd, 1}; t.probChoice(sInd:eInd, 2)' nan(1, buff)];
            prevRwds{aInd, 1} = [prevRwds{aInd, 1}; prevRwdsTmp(sInd:eInd) nan(1, buff)];
            changeChoice{aInd, 1} = [changeChoice{aInd, 1}; changeChoiceTmp(sInd:eInd) nan(1, buff)];
        end
    end
    
 
end

for currA = 1:numA
    for currP = 1:2
        for rInd = 1:size(choices{currA, currP}, 1)
            tmp = [];
            for tInd = 1:99
                if choices{currA, currP}(rInd, tInd:tInd+1) == [1 1] | choices{currA, currP}(rInd, tInd:tInd+1) == [1 0]
                    tmp = [tmp tInd tInd+1];
                end
            end
            tmp = setdiff([1:100], tmp); 
            prevRwds{currA, currP}(rInd, tmp) = NaN;
            changeChoice{currA, currP}(rInd, tmp) = NaN;
        end
        for tInd = 1:100
            wS{currA, currP}(tInd) = 1 - ((sum(changeChoice{currA, currP}(find(prevRwds{currA, currP}(:,tInd)==1), tInd)))/sum(prevRwds{currA, currP}(:,tInd)==1));
            lS{currA, currP}(tInd) = sum(changeChoice{currA, currP}(find(prevRwds{currA, currP}(:,tInd)==0), tInd))/sum(prevRwds{currA, currP}(:,tInd)==0);
        end

        choices{currA, currP} = nanmean(choices{currA, currP}, 1);
        probChoice{currA, currP} = nanmean(probChoice{currA, currP}, 1);
    end
end

choices{1} = cell2mat(choices(:,1));
choices{2} = cell2mat(choices(:,2));
probChoice{1} = cell2mat(probChoice(:,1));
probChoice{2} = cell2mat(probChoice(:,2));
wS{1} = cell2mat(wS(:,1));
wS{2} = cell2mat(wS(:,2));
lS{1} = cell2mat(lS(:,1));
lS{2} = cell2mat(lS(:,2));


colors = cool(4);
figure; hold on;
plotFilled([0:29], probChoice{2}(:, 1:30), colors(4,:));
plotFilled([0:29], choices{2}(:, 1:30), colors(3,:));
plotFilled([0:29], probChoice{1}(:, 1:30), colors(2,:));
plotFilled([0:29], choices{1}(:, 1:30), colors(1,:));
xlabel('trials from block begin')
ylabel('choice average')
set(gca, 'tickdir', 'out')
legend('sim 90/10', '', 'actual 90/10', '', 'sim 50/10', '', 'actual 50/10', '')
ylim([0.5 0.75])

suptitle(strrep(p.Results.modelName, '_', ' '))
set(gcf, 'renderer', 'painters', 'position', [-1578 294 1211 562])

figure;
subplot(2,3,[1:3]); hold on;
plotFilled([0:29], wS{1}(:, 1:30), colors(3,:));
plotFilled([0:29], wS{2}(:, 1:30), colors(4,:));
xlabel('trials from block begin')
ylabel('probability')
set(gca, 'tickdir', 'out')
legend('w-s 50/10', '', 'w-s 90/10', '')

subplot(2,3,[4:6]); hold on;
plotFilled([0:29], lS{1}(:, 1:30), colors(2,:));
plotFilled([0:29], lS{2}(:, 1:30), colors(1,:));
xlabel('trials from block begin')
ylabel('probability')
set(gca, 'tickdir', 'out')
legend('l-s 50/10', '', 'l-s 90/10', '')

suptitle(strrep(p.Results.modelName, '_', ' '))
set(gcf, 'renderer', 'painters', 'position', [-1578 294 1211 562])




