function [mdl] = simulationLogRegIndividual(xlFile, sheet, category, varargin)


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
p.addParameter('randomSeed', 27184)
p.addParameter('revForFlag', 0)
p.addParameter('tMax', 12)
p.parse(varargin{:});

%get params from inputs
rSeed = p.Results.randomSeed;
numMdls = length(p.Results.modelNames);

%get approrpriate root and sep for computer
[root, sep] = currComputer;


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

%initialze cell arrays for model data
rwdMatx = cell(numAnimals, numMdls);
noRwdMatx = cell(numAnimals, numMdls);
combinedAllChoice_R = cell(numAnimals, numMdls);
regCoeffs = cell(1,2);

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
            if regexp(p.Results.taskType, 'switch')
                [mdl, tMax] = combineLogReg_opMD(xlFile, animal, 'switch', 'revForFlag', p.Results.revForFlag,...
                                                    'tMax', p.Results.tMax);
            else
                [mdl, tMax] = combineLogReg_opMD(xlFile, animal, category, 'revForFlag', p.Results.revForFlag,...
                                                    'tMax', p.Results.tMax);
            end
            regCoeffs{1} = [regCoeffs{1}; mdl.Coefficients.Estimate(2:tMax+1)'];
            regCoeffs{2} = [regCoeffs{2}; mdl.Coefficients.Estimate(tMax+2:end)'];
        end
        prevAnimal = animal;

        for currM = 1:numMdls
            modelPath = [root animal sep animal 'sorted' sep 'stan' sep 'bernoulli' sep p.Results.modelNames{currM} sep animal...
                    p.Results.beh '_' p.Results.modelNames{currM} '.mat'];
            [t, ~] = getStanModelParams_samps(p.Results.modelNames{currM}, modelPath, p.Results.samps,...
                            'sessionParamsFlag', p.Results.sessionParamsFlag, 'sessionName', sessionName);

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

                    rwdMatx{currA, currM} = [rwdMatx{currA, currM} NaN(tMax,100) rwdMatxTmp];
                    noRwdMatx{currA, currM} = [noRwdMatx{currA, currM} NaN(tMax,100) noRwdMatxTmp];
                    combinedAllChoice_R{currA, currM} = [combinedAllChoice_R{currA, currM} NaN(1,100) allChoice_R];
                end
            end
        end
    end
end

numAnimals = length(animals);
%run regression for each model simulation
regCoeffsSim = cell(2, numMdls);
for currM = 1:numMdls
    for currA = 1:numAnimals
        mdl = fitglm([rwdMatx{currA, currM}' noRwdMatx{currA, currM}'], combinedAllChoice_R{currA, currM},'distribution','binomial','link','logit');
        regCoeffsSim{1, currM} = [regCoeffsSim{1, currM}; mdl.Coefficients.Estimate(2:tMax+1)'];
        regCoeffsSim{2, currM} = [regCoeffsSim{2, currM}; mdl.Coefficients.Estimate(tMax+2:end)'];
    end
end

%plot regression coefficients
colors = cool(numMdls + 1);
regInds = [2:tMax+1; tMax+2:tMax*2+1];
figure;
for currReg = 1:2
    legTxt = {'actual', ' '};
    subplot(1,2,currReg); hold on;
    plotFilled([1:tMax], regCoeffs{currReg}, colors(1,:));
    for currM = 1:numMdls
        plotFilled([1:tMax], regCoeffsSim{currReg, currM}, colors(currM+1,:));
        [rho, pVal] = corr(nanmean(regCoeffsSim{currReg, currM}, 1)', nanmean(regCoeffs{currReg}, 1)');
        legTxt = [legTxt {strrep(p.Results.modelNames{currM},'_',' ')} {['\rho = ' num2str(rho) ', p = ' num2str(pVal)]}];
    end
    if currReg == 1
        xlabel('reward n trials back')
        yl = ylim;
    else
        xlabel('no reward n trials back')
        ylim(yl)
    end
    legend(legTxt)
    ylabel('\beta Coefficient')
    xlim([0 tMax+1])
    set(gca, 'tickdir', 'out')
end

set(gcf, 'renderer', 'painters')
