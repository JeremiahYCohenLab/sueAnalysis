session = 'mZS078d20220605';
photometryAnalysisSession(session);
figure2;
hold on;
for j = 1:psthBinNum
    plotFilled(midPointsFilter, LCmatChoiceFiltered(target>=edges(j)&target<edges(j+1),:), colorPSTH(j,:));
end

set(gca, 'Box', 'off' )
set(gca, 'TickDir', 'Out')
xlabel('time from respond (s)','FontSize', 15)
ylabel('dF/F, zscored', 'FontSize', 15)
xlim([-1 2.5])

title(['LC'])
%%
figure2;
hold on;
for j = 1:psthBinNum
    plotFilled(midPointsFilter, mPFCmatChoiceFiltered(target>=edges(j)&target<edges(j+1),:), colorPSTH(j,:));
end

set(gca, 'Box', 'off' )
set(gca, 'TickDir', 'Out')
xlabel('time from respond (s)','FontSize', 15)
ylabel('dF/F, zscored', 'FontSize', 15)
xlim([-1 2.5])

title(['mPFC'])
%%
currWin = focusWins{2};
focusMean = mean(mPFCmatChoice(:, midPoints>=currWin(1)&midPoints<currWin(2)), 2);
focusMean = focusMean(sortedInd);
[sigSorted] = sort(focusMean);
lowLim = sigSorted(round(0.05*length(focusMean)));
highLim = sigSorted(round(0.95*length(focusMean)));
focusMean = focusMean(focusMean > lowLim & focusMean < highLim);
figure2;
imagesc(zscore(focusMean));
myMap = [linspace(0, 1, 200);linspace(0, 1, 200);linspace(0, 1, 200)]';
colormap(myMap);
title('mPFC')
ylabel('1<--- rpe --->-1', 'FontSize',15)
%%
currWin = focusWins{2};
focusMean = mean(LCmatChoice(:, midPoints>=currWin(1)&midPoints<currWin(2)), 2);
focusMean = focusMean(sortedInd);
[sigSorted] = sort(focusMean);
lowLim = sigSorted(round(0.05*length(focusMean)));
highLim = sigSorted(round(0.95*length(focusMean)));
focusMean = focusMean(focusMean > lowLim & focusMean < highLim);
figure2;
imagesc(zscore(focusMean));
myMap = [linspace(0, 1, 200);linspace(0, 1, 200);linspace(0, 1, 200)]';
colormap(myMap);
title('LC')
ylabel('1<--- rpe --->-1', 'FontSize',15)
%%