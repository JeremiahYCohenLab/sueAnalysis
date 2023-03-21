xlFile = 'allDBh-cre';
sheet = 'all-DBh';
col = 'allLick';
dayList = getDayList(xlFile, sheet, col);
%%
xlFile = 'inhibitionGt';
sheet = 'allGt';
col = 'cueOnGood';
dayList = getDayList(xlFile, sheet, col);
%% lickLat switch vs stay
combineLickLats = [];
combineSvs = [];
for i = 1:length(dayList)
    s = behAnalysisNoPlot_opMD(dayList{i}, 'simpleFlag', 1);
    combineLickLats = [combineLickLats s.lickLat];
    svs = NaN(size(s.allChoices));
    svs(s.changeChoice_Inds) = 1;
    svs(s.stayChoice_Inds) = 0;
    combineSvs = [combineSvs svs];
end
%%
figure2; hold on;
edges = linspace(min(combineLickLats)-0.01, max(combineLickLats)+0.01, 100);
histogram(combineLickLats(combineSvs==1), edges, 'FaceColor', 'r', 'EdgeColor', 'none', 'Normalization', 'probability');
histogram(combineLickLats(combineSvs==0), edges, 'FaceColor', 'k', 'EdgeColor', 'none', 'Normalization', 'probability');
xlabel('lickLat (ms)')
legend({'switch', 'stay'})
set(gca, 'TickDir', 'out');
%% ILI sorted by Qchosen

%% lickLat sorted by Qchosen


%% lickRate aligned to rwd in absence of reward 
tb = 1;
tf = 2;
stepSize = 20;
binSize = 50;
allMeansNrwd300 = [];
allMeansRwd300 = [];
allMeansNrwd200 = [];
allMeansRwd200 = [];
for i = 1:length(dayList)
    session = dayList{i};
    [cellChoice, matChoice, matChoiceSlide, slideTime] = getLickMatChoice(session, tb, tf, stepSize, binSize);
    s = behAnalysisNoPlot_opMD(dayList{i}, 'simpleFlag', 1);
    [~, sortInd] = sort(abs(s.allRewards));
    if s.rwdDelay>190 && s.rwdDelay<210
        allMeansNrwd200 = [allMeansNrwd200; mean(matChoiceSlide(s.nrwd_Inds,:))];
        allMeansRwd200 = [allMeansRwd200; mean(matChoiceSlide(s.rwd_Inds,:))];
    else
        if s.rwdDelay>290 && s.rwdDelay<310
            allMeansNrwd300 = [allMeansNrwd300; mean(matChoiceSlide(s.nrwd_Inds,:))];
            allMeansRwd300 = [allMeansRwd300; mean(matChoiceSlide(s.rwd_Inds,:))];      
        end
    end
    figure;
    subplot(1,4,1); hold on;
    plotSpikeRaster(cellChoice(sortInd),'PlotType','vertline'); hold on;
    plot([300 300], [0 length(cellChoice)], 'Color', 'r', 'LineStyle', '--')
    xlim([-100 1500])
    title(session);
    subplot(3,4,2); hold on;
    plotFilled(slideTime, matChoiceSlide(abs(s.allRewards)==1, :), 'b');
    plotFilled(slideTime, matChoiceSlide(abs(s.allRewards)==0, :), 'r');
    plot([300 300], [0 10], 'Color', 'r', 'LineStyle', '--')
    subplot(3,4,3); hold on;
    matChoice(:, tb*1000) = 0;
    plot(-tb*1000:tf*1000, mean(matChoice(abs(s.allRewards)==1, :)), 'b');
    plot(-tb*1000:tf*1000, mean(matChoice(abs(s.allRewards)==0, :)), 'r');
    plot([300 300], [0 0.1], 'Color', 'k', 'LineStyle', '--')
    
    screen = get(0,'Screensize');
    screen(4) = screen(4) - 100;
    set(gcf, 'Position', screen)
end
%%
figure2;
subplot(1,2,1);
plotFilled(slideTime, allMeansNrwd300, 'r')
hold on;
plotFilled(slideTime, allMeansRwd300, 'b')
plot([300 300], [0 10], 'Color', 'k', 'LineStyle', '--')
xlabel('time from first lick (ms)')
set(gca, 'TickDir', 'out');
set(gca, 'box', 'off')
title('rwdDelay 300')
subplot(1,2,2);
plotFilled(slideTime, allMeansNrwd200, 'r')
hold on;
plotFilled(slideTime, allMeansRwd200, 'b')
plot([200 200], [0 10], 'Color', 'k', 'LineStyle', '--')
xlabel('time from first lick (ms)')
set(gca, 'TickDir', 'out');
set(gca, 'box', 'off')
title('rwdDelay 200')
%%