function [wsls, transChoiceMatxMed, transChoiceMatxHigh, sessionNames, meanTrials, rwdMatxMed, rwdMatxHigh] = ...
    transitionAnalysis_opMD(xlFile, sheet, category, probs, tranWin, figFlag)

if nargin < 6
    figFlag = 1;
end

if nargin < 5
    tranWin = 10;
end

if nargin < 4
    pHigh = 90;
    pMed = 50;
else
    pHigh = probs(1);
    pMed = probs(2);
end

probDiffH = pHigh - 10;

[root, sep] = currComputer();
transChoiceMatxMed = []; 
transChoiceMatxHigh = [];
changeChoiceMatxMed = [];
changeChoiceMatxHigh = [];
prevRwdMatxMed = [];
prevRwdMatxHigh = [];
rwdMatxMed = [];
rwdMatxHigh = [];
range = 15;

[~, dayList, ~] = xlsread([root xlFile], sheet);
[~,col] = find(strcmp(dayList, category));
dayList = dayList(2:end,col);
endInd = find(cellfun(@isempty,dayList),1);
if ~isempty(endInd)
    dayList = dayList(1:endInd-1,:);
end

sessionNames = [];
numTrials = [];

for currS = 1: length(dayList)
    sessionName = dayList{currS};

    [behSessionData, unCorrectedBlockSwitch, out] = loadBehavioralData([sessionName '.asc'], 0);
    s = parseBehavioralData(behSessionData, unCorrectedBlockSwitch);

    numTrials = [numTrials length(s.responseInds)];
    prevRewardsBin = [0 abs(s.allRewards(1:end-1))];
    changeChoice = [0 abs(diff(s.allChoices)) > 0];

    rwdProb_R = [behSessionData(s.responseInds).rewardProbR]; 
    rwdProb_L = [behSessionData(s.responseInds).rewardProbL]; 

    for j = 2:(length(s.blockSwitch) - 1)
        tmpInd = s.blockSwitch(j);
        if tmpInd-tranWin > 0 && tmpInd+tranWin <= length(s.allChoices)
            if rwdProb_R(tmpInd-1) == pHigh && rwdProb_R(tmpInd) == 10 && any(diff(rwdProb_L(tmpInd-tranWin:tmpInd)) == probDiffH)
                    if (tmpInd - range - 1) > 0 && length(s.allChoices) > (tmpInd + range)
                        transChoiceMatxHigh = [transChoiceMatxHigh; s.allChoices((tmpInd-range+1):(tmpInd+range))];
                        changeChoiceMatxHigh = [changeChoiceMatxHigh; changeChoice((tmpInd-range+1):(tmpInd+range))];
                        prevRwdMatxHigh = [prevRwdMatxHigh; prevRewardsBin((tmpInd-range+1):(tmpInd+range))];
                        rwdMatxHigh = [rwdMatxHigh; s.allRewards((tmpInd-range+1):(tmpInd+range))];
                        sessionNames = [sessionNames {sessionName}];
                    end
            elseif rwdProb_R(tmpInd-1) == pMed && rwdProb_R(tmpInd) == 10 && any(diff(rwdProb_L(tmpInd-tranWin:tmpInd)) == probDiffH)
                    if (tmpInd - range - 1) > 0 && length(s.allChoices) > (tmpInd + range)
                        transChoiceMatxMed = [transChoiceMatxMed; s.allChoices((tmpInd-range+1):(tmpInd+range))];
                        changeChoiceMatxMed = [changeChoiceMatxMed; changeChoice((tmpInd-range+1):(tmpInd+range))];
                        prevRwdMatxMed = [prevRwdMatxMed; prevRewardsBin((tmpInd-range+1):(tmpInd+range))];
                        rwdMatxMed = [rwdMatxMed; s.allRewards((tmpInd-range+1):(tmpInd+range))];
                        sessionNames = [sessionNames {sessionName}];
                    end
            elseif rwdProb_L(tmpInd-1) == pHigh && rwdProb_L(tmpInd) == 10 && any(diff(rwdProb_R(tmpInd-tranWin:tmpInd)) == probDiffH)
                    if (tmpInd - range - 1) > 0 && length(s.allChoices) > (tmpInd + range)
                        transChoiceMatxHigh = [transChoiceMatxHigh; (s.allChoices((tmpInd-range+1):(tmpInd+range))*-1)];
                        changeChoiceMatxHigh = [changeChoiceMatxHigh; changeChoice((tmpInd-range+1):(tmpInd+range))];
                        prevRwdMatxHigh = [prevRwdMatxHigh; prevRewardsBin((tmpInd-range+1):(tmpInd+range))];
                        rwdMatxHigh = [rwdMatxHigh; s.allRewards((tmpInd-range+1):(tmpInd+range))*-1];
                        sessionNames = [sessionNames {sessionName}];
                    end 
            elseif rwdProb_L(tmpInd-1) == pMed && rwdProb_L(tmpInd) == 10 && any(diff(rwdProb_R(tmpInd-tranWin:tmpInd)) == probDiffH)
                    if (tmpInd - range - 1) > 0 && length(s.allChoices) > (tmpInd + range)
                        transChoiceMatxMed = [transChoiceMatxMed; (s.allChoices((tmpInd-range+1):(tmpInd+range))*-1)];
                        changeChoiceMatxMed = [changeChoiceMatxMed; changeChoice((tmpInd-range+1):(tmpInd+range))];
                        prevRwdMatxMed = [prevRwdMatxMed; prevRewardsBin((tmpInd-range+1):(tmpInd+range))];
                        rwdMatxMed = [rwdMatxMed; s.allRewards((tmpInd-range+1):(tmpInd+range))*-1];
                        sessionNames = [sessionNames {sessionName}];
                    end
            end   
        end
    end
end   

meanTrials = ceil(mean(numTrials));

transChoiceMatxMed(transChoiceMatxMed == -1) = 0;
transChoiceMatxHigh(transChoiceMatxHigh == -1) = 0;

%find only high-first side choices for ws-ls analysis
for rInd = 1:size(transChoiceMatxHigh,1)
    tmp = [];
    for tInd = 1:range*2-1
        if transChoiceMatxHigh(rInd, tInd:tInd+1) == [1 1] | transChoiceMatxHigh(rInd, tInd:tInd+1) == [1 0]
            tmp = [tmp tInd tInd+1];
        end
    end
    tmp = setdiff([1:range*2], tmp); 
    prevRwdMatxHigh(rInd, tmp) = NaN;
    changeChoiceMatxHigh(rInd, tmp) = NaN;
end
for rInd = 1:size(transChoiceMatxMed,1)
    tmp = [];
    for tInd = 1:range*2-1
        if transChoiceMatxMed(rInd, tInd:tInd+1) == [1 1] | transChoiceMatxMed(rInd, tInd:tInd+1) == [1 0]
            tmp = [tmp tInd tInd+1];
        end
    end
    tmp = setdiff([1:range*2], tmp); 
    prevRwdMatxMed(rInd, tmp) = NaN;
    changeChoiceMatxMed(rInd, tmp) = NaN;
end
%find ws-ls rates and bernoullli sem's
for tInd = 1:range*2
    if ~isempty(transChoiceMatxMed)
        lS_med(tInd) = sum(changeChoiceMatxMed(prevRwdMatxMed(:,tInd)==0, tInd))/sum(prevRwdMatxMed(:,tInd)==0);
        sem_lS_med(tInd) = sem_bernoulli(sum(changeChoiceMatxMed(prevRwdMatxMed(:,tInd)==0, tInd)), sum(prevRwdMatxMed(:,tInd)==0));
        wS_med(tInd) = 1 - ((sum(changeChoiceMatxMed(prevRwdMatxMed(:,tInd)==1, tInd)))/sum(prevRwdMatxMed(:,tInd)==1));
        sem_wS_med(tInd) = sem_bernoulli(sum(~changeChoiceMatxMed(prevRwdMatxMed(:,tInd)==1, tInd)), sum(prevRwdMatxMed(:,tInd)==1));
    end
    if ~isempty(transChoiceMatxHigh)
        lS_high(tInd) = sum(changeChoiceMatxHigh(find(prevRwdMatxHigh(:,tInd)==0), tInd))/sum(prevRwdMatxHigh(:,tInd)==0);
        sem_lS_high(tInd) = sem_bernoulli(sum(changeChoiceMatxHigh(find(prevRwdMatxHigh(:,tInd)==0), tInd)), sum(prevRwdMatxHigh(:,tInd)==0));
        wS_high(tInd) = 1 - ((sum(changeChoiceMatxHigh(find(prevRwdMatxHigh(:,tInd)==1), tInd)))/sum(prevRwdMatxHigh(:,tInd)==1));
        sem_wS_high(tInd) = sem_bernoulli(sum(~changeChoiceMatxHigh(find(prevRwdMatxHigh(:,tInd)==1), tInd)), sum(prevRwdMatxHigh(:,tInd)==1));
    end
end

cWin = 8;
limVal = 1;
transChoiceMatxMed_lim = [];
transChoiceMatxHigh_lim = [];
for currT = 1:size(transChoiceMatxMed,1)
    if sum(transChoiceMatxMed(currT,range-cWin+1:range)) >= limVal*cWin  & sum(prevRwdMatxMed(currT,range-cWin+2:range+1)) <= ceil(pMed/100 * cWin)
        transChoiceMatxMed_lim = [transChoiceMatxMed_lim; transChoiceMatxMed(currT,:)];
    end
end
for currT = 1:size(transChoiceMatxHigh,1)
    if sum(transChoiceMatxHigh(currT,range-cWin+1:range)) >= limVal*cWin  & sum(prevRwdMatxHigh(currT,range-cWin+2:range+1)) >= floor(pHigh/100 * cWin)
        transChoiceMatxHigh_lim = [transChoiceMatxHigh_lim; transChoiceMatxHigh(currT,:)];
    end
end

%run an LME to predict choice average after transitions as a function of trial and transition type
medAvg = mean(transChoiceMatxMed(:,range:range+range));
highAvg = mean(transChoiceMatxHigh(:,range:range+range));

allTrans = [medAvg highAvg];
numTT = length([medAvg]);
transType = [ones(1, numTT) ones(1, numTT)*2];
trial = repmat([1:range+1], 1, length(allTrans)/(range+1));
tt = table(zscore(allTrans)', trial', transType', ...
        'VariableNames', {'choiceProbs', 'trial', 'trans'});
tt.trans = nominal(tt.trans);
mdl = fitlme(tt, 'choiceProbs~trial*trans')
stats = anova(mdl)


%run lme on choice history controlled transitions
stayWin = 8;
transAll = [transChoiceMatxMed; transChoiceMatxHigh];
choiceX = sum(transAll(:,range-stayWin+1:range), 2);
xInds = find(choiceX == stayWin);
transAll = transAll(xInds, :);
rwdHx = [rwdMatxMed; rwdMatxHigh];
rwdHx = sum(rwdHx(:,range-stayWin+1:range), 2);
rwdHx = rwdHx(xInds);
sortIndsI = find(rwdHx <= 0.5*stayWin);
sortIndsII = find(rwdHx > 0.5*stayWin);

allTrans = [mean(transAll(sortIndsI,range:range+range)) mean(transAll(sortIndsII,range:range+range))];
tt = table(zscore(allTrans)', trial', transType', ...
        'VariableNames', {'choiceProbs', 'trial', 'trans'});
tt.trans = nominal(tt.trans);
mdl = fitlme(tt, 'choiceProbs~trial*trans')
stats = anova(mdl)


if figFlag && ~isempty(transChoiceMatxMed) && ~isempty(transChoiceMatxHigh)
    figure; 
    subplot(2,3,1); hold on;
    x = [-range+1:range];
    plotFilledBern(x, transChoiceMatxMed, 'r');
    plotFilledBern(x, transChoiceMatxHigh, 'b');
    ylabel('Choice probability')
    plot([0 0], [0 1], ':k')
    plot([x(1) x(end)], [0.5 0.5], ':k')
    legend({'median-low', '', 'high-low', ''})
    ylim([0 1])
    title('all transitions')
    set(gca, 'tickdir', 'out')
    
    subplot(2,3,2); hold on;
    if ~isempty(transChoiceMatxMed_lim)
        plotFilledBern(x, transChoiceMatxMed_lim, 'r');
    end
    if ~isempty(transChoiceMatxHigh_lim)
        plotFilledBern(x, transChoiceMatxHigh_lim, 'b');
    end
    plot([0 0], [0 1], ':k')
    plot([x(1) x(end)], [0.5 0.5], ':k')
    ylim([0 1])
    title(['P(C=higher) >= ' num2str(limVal) ' in ' num2str(cWin) ' trials before transition'])
    set(gca, 'tickdir', 'out')
    
    subplot(2,3,4); hold on
    errorbar(x, wS_med, sem_wS_med, 'r', 'linewidth', 1.3)
    errorbar(x, wS_high, sem_wS_high, 'b', 'linewidth', 1.3)
    ylabel('Win-stay')
    xlabel('Trials from switch')
    y = ylim;
    plot([0 0], y, ':k')
    set(gca, 'tickdir', 'out')

    subplot(2,3,5); hold on
    errorbar(x, lS_med, sem_lS_med, 'r', 'linewidth', 1.3)
    errorbar(x, lS_high, sem_lS_high, 'b', 'linewidth', 1.3)
    ylabel('Lose-shift')
    xlabel('Trials from switch')
    y = ylim;
    plot([0 0], y, ':k')
    legend('medium -> med', 'high -> med')
    set(gca, 'tickdir', 'out')
    
    subplot(2,3,3); hold on;
    colors2 = cool(2);
    plotFilled(x, transAll(sortIndsI, :), colors2(1,:));
    plotFilled(x, transAll(sortIndsII, :), colors2(2,:));
    plot([0 0], [0 1], ':k')
    plot([x(1) x(end)], [0.5 0.5], ':k')
    legend({'<=5','','>5',''})
    ylim([0 1])
    ylabel('choice average')
    title('transitions by rwd hist P(C=higher)=1)')
    
    subplot(2,3,6); hold on;
    transAll = [transChoiceMatxMed; transChoiceMatxHigh];
    choiceX = sum(transAll(:,range-stayWin+1:range), 2);
    xInds = find(choiceX >= 0.8*stayWin);
    transAll = transAll(xInds, :);
    rwdHx = [rwdMatxMed; rwdMatxHigh];
    rwdHx = sum(rwdHx(:,range-stayWin+1:range), 2);
    rwdHx = rwdHx(xInds);
    sortIndsI = find(rwdHx <= 0.5*stayWin);
    sortIndsII = find(rwdHx > 0.5*stayWin);
    colors2 = cool(2);
    plotFilled(x, transAll(sortIndsI, :), colors2(1,:));
    plotFilled(x, transAll(sortIndsII, :), colors2(2,:));
    plot([0 0], [0 1], ':k')
    plot([x(1) x(end)], [0.5 0.5], ':k')
    legend({'<=5','','>5',''})
    ylim([0 1])
    ylabel('choice average')
    title('transitions by rwd hist P(C=higher)>=0.8')

    set(gcf, 'position', [0 157 917 745], 'renderer', 'painters')
    suptitle([sheet ' ' category])
    
end


if ~isempty(transChoiceMatxMed)
    wsls.wS_med = wS_med;
    wsls.sem_wS_med = sem_wS_med;
    wsls.lS_med = lS_med;
    wsls.sem_lS_med = sem_lS_med;
end
if ~isempty(transChoiceMatxHigh)
    wsls.wS_high = wS_high;
    wsls.sem_wS_high = sem_wS_high;
    wsls.lS_high = lS_high;
    wsls.sem_lS_high = sem_lS_high;
end

end
    