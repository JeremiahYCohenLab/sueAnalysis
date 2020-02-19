function compareSuccessSim_opMD(xlFile, animals, category, varargin)

p = inputParser;
% default parameters if none given
p.addParameter('rwdProbs', [90 50 10]);
p.addParameter('params', []);
p.addParameter('maxTrials', 350);
p.addParameter('runs', 400);
p.addParameter('randomSeed', 89465);
p.addParameter('modelName', [{'sixParam_absPePeAN_bi'}, {'fourParam'}]);
p.addParameter('bernFlag', [1 1]);
p.addParameter('revForFlag', 0);
p.parse(varargin{:});

[root, sep] = currComputer();

numMdls = length(p.Results.modelName);
numAnimals = length(animals);

if isempty(p.Results.runs)
    runs = length(dayList);  
else
    runs = p.Results.runs;
end

sumRwds = 0;
sumChoices = 0;
sumRwds_sim = zeros(1,length(p.Results.modelName));
sumChoices_sim = zeros(1,length(p.Results.modelName));
sumRwds_c = 0;
sumChoices_c = 0;
rSeed = p.Results.randomSeed;

for currA = 1:numAnimals

    [~, dayList, ~] = xlsread(xlFile, animals{currA});
    [~,col] = find(~cellfun(@isempty,strfind(dayList, category)) == 1);
    dayList = dayList(2:end,col);
    endInd = find(cellfun(@isempty,dayList),1);
    if ~isempty(endInd)
        dayList = dayList(1:endInd-1,:);
    end
    
    for currS = 1:length(dayList)
        sessionName = dayList{currS};
        filename = [sessionName '.asc'];
        [behSessionData, unCorrectedBlockSwitch, out] = loadBehavioralData(filename, p.Results.revForFlag);
        behavStruct = parseBehavioralData(behSessionData, unCorrectedBlockSwitch);

        sumRwds = sumRwds + sum(abs(behavStruct.allRewards)); 
        sumChoices = sumChoices + length(behavStruct.allRewards);

        if isempty(p.Results.params)
            for currM = 1:length(p.Results.modelName)
                if p.Results.bernFlag(currM)
                    modelPath = [root animals{currA} sep animals{currA} 'sorted' sep 'stan' sep 'bernoulli' sep p.Results.modelName{currM}...
                        sep animals{currA} category '_' p.Results.modelName{currM} '.mat'];
                else
                    modelPath = [root animals{currA} sep animals{currA} 'sorted' sep 'stan' sep p.Results.modelName{currM} sep animals{currM}...
                                category '_' p.Results.modelName{currM} '.mat'];
                end
                t = generateStanModelTerms_opMD(p.Results.modelName{currM}, modelPath, sessionName, 1);
                params{currM} = t.params;
            end
        else
            for currM = 1:length(p.Results.modelName)
                params{currM} = p.Results.params;
            end
        end


        for i = 1:runs
            if rem(i,100) == 0
                fprintf('Running simulation %d of %d for animal %d of %d \n', i, runs, currA, numAnimals);
            end
            rSeed = rSeed + 1;

            for currM = 1:numMdls
                switch p.Results.modelName{currM}
                    case 'fourParam'
                        [allRewards, allChoices, ~, ~] = qLearningModel_simNoPlot('params', params{currM},...
                            'maxTrials', p.Results.maxTrials, 'randomSeed', rSeed, 'rwdProbs', p.Results.rwdProbs);
                    case 'fiveParam_rBeta_scale'
                        [rBar, allRewards, allChoices, ~, ~] = qLearningModel_rBeta_simNoPlot('params', params{currM},...
                            'maxTrials', p.Results.maxTrials, 'randomSeed', rSeed, 'rwdProbs', p.Results.rwdProbs);
                    case 'sixParam_rBeta_scale_min'
                        [rBar, allRewards, allChoices, ~, ~] = qLearningModel_rBeta_min_simNoPlot('params', params{currM},...
                            'maxTrials', p.Results.maxTrials, 'randomSeed', rSeed, 'rwdProbs', p.Results.rwdProbs);
                    case 'sixParam_rBeta_kappa'
                        [rBar, allRewards, allChoices, ~, ~] = qLearningModel_rBetaKappa_simNoPlot('params', params{currM},...
                            'maxTrials', p.Results.maxTrials, 'randomSeed', rSeed, 'rwdProbs', p.Results.rwdProbs);
                    case 'fiveParam_peBeta_avg'
                        [~, allRewards, allChoices, ~, ~] = qLearningModel_peBeta_avg_simNoPlot('params', params{currM},...
                            'maxTrials', p.Results.maxTrials, 'randomSeed', rSeed, 'rwdProbs', p.Results.rwdProbs);
                    case 'sixParam_absPePeAN_bi'
                        [~, allRewards, allChoices, ~, ~] = qLearningModel_absPePeAN_bi_simNoPlot('params', params{currM},...
                            'maxTrials', p.Results.maxTrials, 'randomSeed', rSeed, 'rwdProbs', p.Results.rwdProbs);
                end

                sumRwds_sim(currM) = sumRwds_sim(currM) + sum(abs(allRewards));
                sumChoices_sim(currM) = sumChoices_sim(currM) + length(allChoices);
            end

            %find reward rate of clairvoyant mouse
            [allRewards_c, allChoices_c, ~, ~] = qLearningModel_clairvoyant_simNoPlot('maxTrials', p.Results.maxTrials,...
            'randomSeed', rSeed, 'rwdProbs', p.Results.rwdProbs);

            sumRwds_c = sumRwds_c + sum(abs(allRewards_c));
            sumChoices_c = sumChoices_c + length(allChoices_c);

        end
    end
end



rwdRate = sumRwds/sumChoices;
rwdRate_c = sumRwds_c/sumChoices_c;
for currM = 1:numMdls
    rwdRate_sim(currM) = sumRwds_sim(currM)/sumChoices_sim(currM);
end

colors = cool(numMdls + 2);
figure; hold on;
scatter(0, rwdRate, 2000, colors(1,:),  'filled')
text(0, rwdRate, num2str(rwdRate), 'horizontalalignment', 'center');
scatter(1, rwdRate_c, 2000, colors(2,:),  'filled')
text(1, rwdRate_c, num2str(rwdRate_c), 'horizontalalignment', 'center');
tickLbls = {'actual', 'clairvoyant'};
for currM = 1:numMdls
    scatter(currM+1, rwdRate_sim(currM), 2000, colors(currM+2,:),  'filled')
    text(currM+1, rwdRate_sim(currM), num2str(rwdRate_sim(currM)), 'horizontalalignment', 'center');
    tickLbls = [tickLbls {p.Results.modelName{currM}}];
end
tickLbls = strrep(tickLbls, '_', ' ');
xlim([-0.5 numMdls+1.5])
ylabel('reward rate')
xticks([0:numMdls+1])
xticklabels(tickLbls)
xtickangle(20)
set(gca, 'tickdir', 'out')

