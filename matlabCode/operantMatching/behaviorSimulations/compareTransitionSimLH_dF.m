function compareTransitionSimLH_dF(xlFile, sheet, category, beh, varargin)

p = inputParser;
% default parameters if none given
p.addParameter('rwdProbs', [90 50 10])
p.addParameter('tranWin', 5)
p.addParameter('runs', 10) %max runs per sample per condition
p.addParameter('randomSeed', 698512)
p.addParameter('modelNames', {'fiveParam_bias', 'sevenParam_absPePeAN_scale_int_bias_ord'})
p.addParameter('sessionParamsFlag', 0)
p.addParameter('bernFlag', 1)
p.addParameter('samps', 1000)
p.parse(varargin{:});

[root, sep] = currComputer();

%get task params for finding special transitions
pHigh = p.Results.rwdProbs(1);
pMed = p.Results.rwdProbs(2);
pLow = p.Results.rwdProbs(3);
probDiffH = pHigh - pLow;
tranWin = p.Results.tranWin;

%get transition results from actual data
[wsls, transMed, transHigh, aNames] = transitionAnalysis_opMD(xlFile, sheet, category, p.Results.rwdProbs, p.Results.tranWin);

%initialize cells for choice probs around transitions
numMdls = length(p.Results.modelNames);
cpHigh = cell(1,numMdls);
cpMed = cell(1,numMdls);
cp = cell(1,numMdls);
LHs = cell(1,numMdls);

%set range for number of trials on either side of transition
range = 15;

%extract session list from excel file
[~, seshList, ~] = xlsread(xlFile, sheet);
[~,col] = find(strcmp(seshList, category));
seshList = seshList(2:end, col);
endInd = find(cellfun(@isempty,seshList),1);
if ~isempty(endInd)
    seshList = seshList(1:endInd-1,:);
end

prevAnimal = [];
%get choice probs around transitions for all models
for currM = 1:numMdls
    for currSesh = 1:length(seshList)
        sessionName = seshList{currSesh};
        [animal, ~] = strtok(sessionName, 'd'); 
        animal = animal(2:end);
        if strcmp(animal, prevAnimal) == 0 
            fprintf('Simulating animal %s, model %d of %d \n', animal, currM, numMdls);
            prevAnimal = animal;
        end
        if p.Results.bernFlag
            modelPath = [root animal sep animal 'sorted' sep 'stan' sep 'bernoulli' sep p.Results.modelNames{currM}...
                sep animal beh '_' p.Results.modelNames{currM} '.mat'];
        else
            modelPath = [root animal sep animal 'sorted' sep 'stan' sep p.Results.modelNames{currM}...
                sep animal beh '_' p.Results.modelNames{currM} '.mat'];
        end
        t = getStanModelParams_samps(p.Results.modelNames{currM}, modelPath, p.Results.samps, 'sessionParamsFlag', p.Results.sessionParamsFlag,...
            'sessionName', sessionName, 'varFlag', 1, 'biasFlag', p.Results.sessionParamsFlag);
%         t = getStanModelParams_mode(p.Results.modelNames{currM}, modelPath, sessionName, p.Results.sessionParamsFlag);
        LHs{currM} = [LHs{currM} t.LH];
        
        %get session behavior data
        [behSessionData,o.blockSwitch,~] = loadBehavioralData([sessionName '.asc'], 0);
        o = parseBehavioralData(behSessionData, o.blockSwitch);
        allChoices = o.allChoice_R - o.allChoice_L;
        rwdProb_L = [behSessionData(o.responseInds).rewardProbL]';
        rwdProb_R = [behSessionData(o.responseInds).rewardProbR]';
        
        
        choiceProbs = nan(1, length(allChoices));
        choiceProbs(allChoices==1) = t.probChoice(allChoices==1, 2);
        choiceProbs(allChoices==-1) = t.probChoice(allChoices==-1, 1);

        for j = 2:(length(o.blockSwitch) - 1)
            tmpInd = o.blockSwitch(j);
            if tmpInd-tranWin > 0 & tmpInd+tranWin <= length(allChoices)
                if rwdProb_R(tmpInd-1) == pHigh & rwdProb_R(tmpInd) == pLow & any(diff(rwdProb_L(tmpInd-tranWin:tmpInd)) == probDiffH)
                        if (tmpInd - range - 1) > 0 & length(allChoices) > (tmpInd + range)
                            cpHigh{currM} = [cpHigh{currM}; choiceProbs((tmpInd-range+1):(tmpInd+range))];
                            cp{currM} = [cp{currM}; choiceProbs((tmpInd-range+1):(tmpInd+range))];
                        end
                elseif rwdProb_R(tmpInd-1) == pMed & rwdProb_R(tmpInd) == pLow & any(diff(rwdProb_L(tmpInd-tranWin:tmpInd)) == probDiffH)
                        if (tmpInd - range - 1) > 0 & length(allChoices) > (tmpInd + range)
                            cpMed{currM} = [cpMed{currM}; choiceProbs((tmpInd-range+1):(tmpInd+range))];
                            cp{currM} = [cp{currM}; choiceProbs((tmpInd-range+1):(tmpInd+range))];
                        end
                elseif rwdProb_L(tmpInd-1) == pHigh & rwdProb_L(tmpInd) == pLow & any(diff(rwdProb_R(tmpInd-tranWin:tmpInd)) == probDiffH)
                        if (tmpInd - range - 1) > 0 & length(allChoices) > (tmpInd + range)
                            cpHigh{currM} = [cpHigh{currM}; choiceProbs((tmpInd-range+1):(tmpInd+range))];
                            cp{currM} = [cp{currM}; choiceProbs((tmpInd-range+1):(tmpInd+range))];
                        end 
                elseif rwdProb_L(tmpInd-1) == pMed & rwdProb_L(tmpInd) == pLow & any(diff(rwdProb_R(tmpInd-tranWin:tmpInd)) == probDiffH)
                        if (tmpInd - range - 1) > 0 & length(allChoices) > (tmpInd + range)
                            cpMed{currM} = [cpMed{currM}; choiceProbs((tmpInd-range+1):(tmpInd+range))];
                            cp{currM} = [cp{currM}; choiceProbs((tmpInd-range+1):(tmpInd+range))];
                        end
                end   
            end
        end
    end
end


x = [-range+1:range];

figure; 
subplot(1,4,1); hold on;
plotFilledBern(x, transMed, [0.7 0.7 0.7]);
plotFilledBern(x, transHigh, [0 0 0]);
legend('medium -> low', '', 'high -> low', '')
plot([0 0], [0 1], ':k')
plot([x(1) x(end)], [0.5 0.5], ':k')
ylim([0 1])
set(gca, 'tickdir', 'out')
xlabel('Trials from switch')
ylabel('Choice probability')
title('actual')

colors = cool(numMdls*2);
legTxt = [];
subplot(1,4,2); hold on;
for currM = 1:numMdls
    plotFilled(x, cpMed{currM}, colors((currM-1)*numMdls + 1,:));
    plotFilled(x, cpHigh{currM}, colors((currM-1)*numMdls + 2,:));
    legTxt = [legTxt {[strrep(p.Results.modelNames{currM}, '_', ' ') ' medium']} {' '}...
                {[strrep(p.Results.modelNames{currM}, '_', ' ') ' high']} {' '}];
end
legend(legTxt)
plot([0 0], [0 1], ':k')
ylim([0 1])
set(gca, 'tickdir', 'out')

legTxt = [];
subplot(1,4,3); hold on;
for currM = 1:numMdls
    plotFilled(x, cp{currM}, colors((currM-1)*numMdls + 1,:));
    legTxt = [legTxt {[strrep(p.Results.modelNames{currM}, '_', ' ')]} {' '}];
end
legend(legTxt)
plot([0 0], [0 1], ':k')
ylim([0 1])
set(gca, 'tickdir', 'out')

subplot(1,4,4); hold on;
for currM = 1:numMdls
    histogram(LHs{currM}, 10, 'facecolor', colors((currM-1)*numMdls + 1,:));
end
set(gca, 'tickdir', 'out', 'box', 'off')
ylabel('session -LLHs')

set(gcf, 'renderer', 'painters', 'position', [-1824 344 1799 457])

if numMdls == 2
    figure;
    subplot(1,2,1);
    LHdiffs = LHs{1} - LHs{2};
    histogram(LHdiffs, 10, 'facecolor', [1 1 1])
    set(gca, 'tickdir', 'out', 'box', 'off')
    xlabel('differences in session LHs')
    ylabel('count')
    
    subplot(1,2,2); hold on;
    for currS = 1:length(LHs{1})
        plot([1 2], [LHs{1}(currS) LHs{2}(currS)], '-k')
    end
    set(gca, 'tickdir', 'out')
    xlim([0 3])
    xticks([1 2])
    tickLbls = [{strrep(p.Results.modelNames{1}, '_', ' ')} {strrep(p.Results.modelNames{2}, '_', ' ')}];
    xticklabels(tickLbls)
    xtickangle(15)
    ylabel('-LLHs')
end
set(gcf, 'renderer', 'painters', 'position', [-1544 298 1203 485])

