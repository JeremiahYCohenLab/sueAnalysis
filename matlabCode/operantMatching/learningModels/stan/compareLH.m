function compareLH(xlFile, sheet, beh, varargin)

%task and model parameters
p = inputParser;
% default parameters if none given
p.addParameter('modelNames', {'fourParam', 'sixParam_absPePeAN_bi', 'eightParam_absPePeAN_bi_bias_k', 'eightParam_absPePeAN_exp_bias_k'})
p.addParameter('compareType', 'same');
p.addParameter('modelBeh', 'clean');
p.addParameter('revForFlag', 0);
p.addParameter('plotFlag', 1);
p.addParameter('bernFlag', []);
p.addParameter('sessionFlag', 1);
p.parse(varargin{:});

if isempty(p.Results.bernFlag)
    bernFlag = ones(1,length(p.Results.modelNames));
else
    bernFlag = p.Results.bernFlag;
end

modelNames = p.Results.modelNames;
[root, sep] = currComputer();

[weights, dayList, ~] = xlsread(xlFile, sheet);
[~,col] = find(~cellfun(@isempty,strfind(dayList, beh)) == 1);
dayList = dayList(2:end,col);
endInd = find(cellfun(@isempty,dayList),1);
if ~isempty(endInd)
    dayList = dayList(1:endInd-1,:);
end

LH = nan(length(dayList), length(p.Results.modelNames));
aic = nan(length(dayList), length(p.Results.modelNames));
bic = nan(length(dayList), length(p.Results.modelNames));

switch p.Results.compareType
    case 'all'
        animalName = 'all';
        for currMod = 1:length(p.Results.modelNames)
            if bernFlag(currMod)
                modelPath = [root animalName sep animalName 'sorted' sep 'stan' sep 'bernoulli' sep modelNames{currMod} sep animalName...
                beh '_' modelNames{currMod} '.mat'];
            else
                modelPath = [root animalName sep animalName 'sorted' sep 'stan' sep modelNames{currMod} sep animalName...
                beh '_' modelNames{currMod} '.mat'];
            end
            for currSesh = 1:length(dayList)
                sessionName = dayList{currSesh};
                t = generateStanModelTerms_opMD(modelNames{currMod}, modelPath, sessionName, p.Results.sesionFlag, p.Results.revForFlag);
                LH(currSesh, currMod) = t.LH;
                [aic(currSesh, currMod), bic(currSesh, currMod)] = aicbic(-1*t.LH, length(t.params), length(t.probChoice));
            end
        end
    case 'across'   %if looking at sessions from different animals
        prevAnimal = [];
        for currSesh = 1:length(dayList)
            sessionName = dayList{currSesh};
            [animalName, date] = strtok(sessionName, 'd'); 
            animalName = animalName(2:end);
            if strcmp(animalName, prevAnimal) == 0
                for currMod = 1:length(p.Results.modelNames)
                    if bernFlag(currMod)
                        modelPath = [root animalName sep animalName 'sorted' sep 'stan' sep 'bernoulli' sep modelNames{currMod} sep animalName...
                        p.Results.modelBeh '_' modelNames{currMod} '.mat'];
                    else
                        modelPath = [root animalName sep animalName 'sorted' sep 'stan' sep modelNames{currMod} sep animalName...
                        p.Results.modelBeh '_' modelNames{currMod} '.mat'];
                    end
                    for currS = 1:length(dayList)
                        if regexp(dayList{currS}, animalName)
                            sessionName = dayList{currS};
                            t = generateStanModelTerms_opMD(modelNames{currMod}, modelPath, sessionName, p.Results.sesionFlag, p.Results.revForFlag);
                            LH(currS, currMod) = t.LH;
                            [aic(currS, currMod), bic(currS, currMod)] = aicbic(-1*t.LH, length(t.params), length(t.probChoice));
                        end
                    end
                end
            end
            prevAnimal = animalName;
        end
    case 'same'  %if only looking at sessions from the same animal
        for currSesh = 1:length(dayList)
            sessionName = dayList{currSesh};
            [animalName, date] = strtok(sessionName, 'd'); 
            animalName = animalName(2:end);

            for currMod = 1:length(p.Results.modelNames)
                if bernFlag(currMod)
                    modelPath = [root animalName sep animalName 'sorted' sep 'stan' sep 'bernoulli' sep modelNames{currMod} sep animalName...
                    beh '_' modelNames{currMod} '.mat'];
                else
                    modelPath = [root animalName sep animalName 'sorted' sep 'stan' sep modelNames{currMod} sep animalName...
                    beh '_' modelNames{currMod} '.mat'];
                end
                t = generateStanModelTerms_opMD(modelNames{currMod}, modelPath, sessionName, p.Results.sessionFlag, p.Results.revForFlag);
                LH(currSesh, currMod) = t.LH;
                meow(currSesh, :) = t.params;
                [aic(currSesh, currMod), bic(currSesh, currMod)] = aicbic(-1*t.LH, length(t.params), length(t.probChoice));
                if currSesh == 1
                    numParams(currMod) = length(t.params);
                end
            end
        end
end


numMdls = length(modelNames);
colors = cool(numMdls);

figure; 
subplot(1,3,1); hold on;
for currMod = 1:numMdls
    histogram(LH(:,currMod), 20, 'FaceColor', colors(currMod,:), 'Normalization', 'Probability')
end
set(gca, 'tickdir', 'out')
ylabel('Probability')
xlabel('-log(likelihood)')
legend(modelNames, 'interpreter', 'none')

subplot(1,3,2); hold on;
for currMod = 1:numMdls
    histogram(bic(:,currMod), 20, 'FaceColor', colors(currMod,:), 'Normalization', 'Probability')
end
set(gca, 'tickdir', 'out')
ylabel('Probability')
xlabel('BIC')

subplot(1,3,3); hold on;
for currMod = 1:numMdls
    histogram(aic(:,currMod), 20, 'FaceColor', colors(currMod,:), 'Normalization', 'Probability')
end
set(gca, 'tickdir', 'out')
ylabel('Probability')
xlabel('AIC')

set(gcf, 'renderer', 'painters')



