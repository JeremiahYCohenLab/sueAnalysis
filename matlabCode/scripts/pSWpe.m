xlFile = {'ZS059', 'ZS060', 'ZS061', 'ZS062'};
sheet = {'ZS059', 'ZS060', 'ZS061', 'ZS062'};
col = {'good', 'good','good','good'};
numBins = 6;
modelName = '5params';
numSamps = 2000;
edgesRpe = linspace(-1-0.01, 0, 4);
edgesRpe = [edgesRpe(1:end-1) linspace(0, 1+0.01, 4)];
pSWNext = [];
pSWNextZS = [];
peAll = [];
[root, sep] = currComputer();
for ani = 1:length(xlFile)
    dayList = getDayList(xlFile{ani}, sheet{ani}, col{ani});
    for i = 1:length(dayList)
        session = dayList{i};
        pd = parseSessionString_df(session, root, sep);
        params = getStanModelParams_sampsOnly(pd.animalName, col{ani}, modelName, numSamps, 'sessionName', session);
        t = inferModelVar(session, params, modelName);
        pe = t.pe;
        s = behAnalysisNoPlot_opMD(session, 'simpleFlag', 1);
        svsNext = NaN(size(s.allChoices));
        svsNext(s.changeChoice_Inds-1) = 1;
        svsNext(s.stayChoice_Inds-1) = 0;
        tempSW = NaN(1,numBins);
        tempPe = NaN(1,numBins);
        for j = 1:numBins
            tempSW(j) = mean(svsNext(pe>=edgesRpe(j) & pe<edgesRpe(j+1)), 'omitnan');
            tempPe(j) = mean(pe(pe>=edgesRpe(j) & pe<edgesRpe(j+1)), 'omitnan');
        end
        pSWNext = [pSWNext; tempSW];
        tempSW(~isnan(tempSW)) = zscore(tempSW(~isnan(tempSW)));
        pSWNextZS = [pSWNextZS; tempSW];
        peAll = [peAll; tempPe];
    end
end
%%
meanPe = mean(peAll, 'omitnan');
meanSW = mean(pSWNext, 'omitnan');
semSW = sem(pSWNext);
figure2;
plot(meanPe, meanSW, 'Color', [0.2 0.2 0.2], 'LineWidth', 2);
patch([meanPe, flip(meanPe)], [meanSW-semSW flip(meanSW+semSW)], [0.2 0.2 0.2], 'edgeColor', 'none', 'FaceAlpha', 0.5);

set(gca, 'Box', 'off');
set(gca,'tickdir', 'out');
set(gca, 'XTick', [-1 0 1])
set(gca, 'YTick', [0:0.2:0.5], 'FontSize',14)
xlim([0 0.5])
xlim([-1.1 1.1])
xlabel('Rpe', 'FontSize', 18);
ylabel('P(switchNext)', 'FontSize', 18);
%%