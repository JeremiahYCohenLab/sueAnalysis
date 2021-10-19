function [mdl] = simulationLogReg(xlFile, sheet, category, varargin)


%task and model parameters
p = inputParser;
% default parameters if none given
p.addParameter('beh', 'clean')
p.addParameter('modelNames', {'fiveParam_bias', 'sevenParam_absPePeAN_scale_int_bias_ord'})
p.addParameter('sessionParamsFlag', 0)
p.addParameter('runs', 1)
p.addParameter('samps', 1000)
p.addParameter('taskType', 'decoupled')
p.addParameter('rwdProbs', [90 50 10])
p.addParameter('maxTrials', 300)
p.addParameter('randomSeed', 27)
p.addParameter('revForFlag', 0)
p.addParameter('tMax', 12)
p.parse(varargin{:});

%get params from inputs
rSeed = p.Results.randomSeed;
numMdls = length(p.Results.modelNames);
if p.Results.sessionParamsFlag
    biasFlag = 1;
else
    biasFlag = 0;
end

%get approrpriate root and sep for computer
[root, sep] = currComputer;

%initialze cell arrays for model data
rwdMatx = cell(1,numMdls);
noRwdMatx = cell(1,numMdls);
combinedAllChoice_R = cell(1,numMdls);

if regexp(p.Results.taskType, 'switch')
    [actual, tMax] = combineLogReg_opMD(xlFile, sheet, 'switch', 'revForFlag', p.Results.revForFlag,...
                                        'tMax', p.Results.tMax);
else
    [actual, tMax] = combineLogReg_opMD(xlFile, sheet, category, 'revForFlag', p.Results.revForFlag,...
                                        'tMax', p.Results.tMax);
end

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

rwdMatx = cell(numMdls, numA);
noRwdMatx = cell(numMdls, numA);
combinedAllChoice_R = cell(numMdls, numA);

prevAnimal = [];
animals = [];
currA = 0;
for currSesh = 1:length(dayList)
    sessionName = dayList{currSesh};
    [animal, ~] = strtok(sessionName, 'd'); 
    animal = animal(2:end);
    if strcmp(animal, prevAnimal) == 0 || p.Results.sessionParamsFlag
        if strcmp(animal, prevAnimal) == 0
            fprintf('Simulating animal %s \n', animal);
            animals = [animals; {animal}];
            currA = currA + 1;
        end
        prevAnimal = animal;

        for currM = 1:numMdls
            modelPath = [root animal sep animal 'sorted' sep 'stan' sep 'bernoulli' sep p.Results.modelNames{currM} sep animal...
                    p.Results.beh '_' p.Results.modelNames{currM} '.mat'];
            [t, ~] = getStanModelParams_samps(p.Results.modelNames{currM}, modelPath, p.Results.samps,...
                            'sessionParamsFlag', p.Results.sessionParamsFlag, 'sessionName', sessionName, 'biasFlag', biasFlag);

            noRwdMatxSesh = []; rwdMatxSesh = []; allChoiceRsesh = [];
            for currS = 1:p.Results.samps
                for currR = 1:p.Results.runs
                    rSeed = rSeed + 1;
                    [~,allRewards, allChoices, blockProbs, blockSwitch] = runSim_dF(p.Results.modelNames{currM}, t.params(currS,:),...
                            p.Results.maxTrials, rSeed, p.Results.rwdProbs);

                    allNoRewards = allChoices;
                    allNoRewards(allRewards == 1) = 0;
                    allNoRewards(allRewards == -1) = 0;

                    allChoice_R = allChoices;
                    allChoice_R(allChoice_R == -1) = 0;

                    rwdMatxTmp = [];
                    noRwdMatxTmp = [];
                    for tInd = 1:tMax
                        rwdMatxTmp(tInd,:) = [NaN(1,tInd) allRewards(1:end-tInd)];
                        noRwdMatxTmp(tInd,:) = [NaN(1,tInd) allNoRewards(1:end-tInd)];
                    end
                    rwdMatxSesh = [rwdMatxSesh NaN(tMax,100) rwdMatxTmp];
                    noRwdMatxSesh = [noRwdMatxSesh NaN(tMax,100) noRwdMatxTmp];
                    allChoiceRsesh = [allChoiceRsesh NaN(1,100) allChoice_R];
                end
            end
            rwdMatx{currM, currA} = [rwdMatx{currM, currA} rwdMatxSesh];
            noRwdMatx{currM, currA} = [noRwdMatx{currM, currA} noRwdMatxSesh];
            combinedAllChoice_R{currM, currA} = [combinedAllChoice_R{currM, currA} allChoiceRsesh];
        end
    end
end

for currM = 1:numMdls
    rwdMatxC{currM} = cell2mat(rwdMatx(currM,:));
    noRwdMatxC{currM} = cell2mat(noRwdMatx(currM,:));
    combinedAllChoice_RC{currM} = cell2mat(combinedAllChoice_R(currM,:));
end


%run regression for each model simulation
for currM = 1:numMdls
    mdl.(p.Results.modelNames{currM}) = fitglm([rwdMatxC{currM}' noRwdMatxC{currM}'], combinedAllChoice_RC{currM},'distribution','binomial','link','logit');
end

%plot regression coefficients
colors = cool(numMdls + 1);
regInds = [2:tMax+1; tMax+2:tMax*2+1];
legTxt = {'actual'};
figure;
for currReg = 1:2
    subplot(1,2,currReg); hold on;
    relevInds = regInds(currReg,:);
    coefVals = actual.Coefficients.Estimate(relevInds);
    CIbands = coefCI(actual);
    errorL = abs(coefVals - CIbands(relevInds,1));
    errorU = abs(coefVals - CIbands(relevInds,2));
    errorbar([1:tMax],coefVals,errorL,errorU,'Color',colors(end,:),'linewidth',2)

    for currM = 1:numMdls
        coefVals = mdl.(p.Results.modelNames{currM}).Coefficients.Estimate(relevInds);
        CIbands = coefCI(mdl.(p.Results.modelNames{currM}));
        errorL = abs(coefVals - CIbands(relevInds,1));
        errorU = abs(coefVals - CIbands(relevInds,2));
        errorbar([1:tMax],coefVals,errorL,errorU,'Color',colors(currM,:),'linewidth',2)
        if currReg == 2
            legTxt = [legTxt {strrep(p.Results.modelNames{currM}, '_', ' ')}];
        end
    end
    if currReg == 1
        xlabel('reward n trials back')
        yl = ylim;
    else
        xlabel('no reward n trials back')
        ylim(yl)
        legend(legTxt)
    end
    ylabel('\beta Coefficient')
    xlim([0 tMax+1])
    set(gca, 'tickdir', 'out')
end

set(gcf, 'renderer', 'painters')
