% get all sessions
xlFile = {'ZS059', 'ZS060', 'ZS061', 'ZS062'};
sheet = {'ZS059', 'ZS060', 'ZS061', 'ZS062'};
col = 'good';
modelName = '5params_k_bias_LaserNegRPE';
numSamps = 2000;
[root, sep] = currComputer();
dayList = {};
for i = 1:length(xlFile)
    dayListCurr = getDayList(xlFile{i}, sheet{i}, col);
    dayList = [dayList; dayListCurr];
end
%%
numBins = 6;
tf = 2;
minILI = 50;
postRwdTime = 500;
ILIDist = [];
ILIDistCleaned = [];
lickDiffCombined =[];
lickLatCombined = [];
lickNumPostRwdCombined = [];
laserCombined = [];
rwdCombined = [];
rpeCombined = [];
svsCombined = [];
Qcombined = [];
lickLatBinned = NaN(length(dayList), numBins);
lickDiffBinned = NaN(length(dayList), numBins);
lickNumBinned = NaN(length(dayList), numBins);
peBinned = NaN(length(dayList), numBins);
lickNumPeBinned = NaN(length(dayList), numBins);
lickNumPeBinnedLaser = NaN(length(dayList), numBins);
lickNumPeBinnedNoLaser = NaN(length(dayList), numBins);
peBinnedNrwd = NaN(length(dayList), 0.5*numBins);
lickNumPeBinnedNrwd = NaN(length(dayList), 0.5*numBins);
lickNumPeBinnedLaserNrwd = NaN(length(dayList), 0.5*numBins);
lickNumPeBinnedNoLaserNrwd = NaN(length(dayList), 0.5*numBins);
peOriBinnedNrwd = NaN(length(dayList), 0.5*numBins);
lickNumPeOriBinnedNrwd  = NaN(length(dayList), 0.5*numBins);
lickNumPeOriBinnedLaserNrwd  = NaN(length(dayList), 0.5*numBins);
lickNumPeOriBinnedNoLaserNrwd  = NaN(length(dayList), 0.5*numBins);
for sess = 1:length(dayList)
    ILIDist = [];
    ILIDistCleaned = [];
    lickDiff = [];
    lickNumPostRwd = [];
    session = dayList{sess};
    s = behAnalysisNoPlot_opMD(session, 'simpleFlag', 1);
    sessionData = s.behSessionData;
    allTrial_lick = {};
    pd = parseSessionString_df(session, root, sep);
    params = getStanModelParams_sampsOnly(pd.animalName, col, modelName, numSamps, 'sessionName', [], 'sessionParamsFlag', 0);
    t = inferModelVar(session, params, modelName, 'perturb', 1);
    Qchosen = t.Q(:,2);
    Qchosen(s.allChoices<0) = t.Q(s.allChoices<0,1);
    Qchosen = zscore(Qchosen);
    for k = 1:length(s.responseInds)
        if ~isnan(sessionData(s.responseInds(k)).rewardL)
            currTrial_lickInd = [sessionData(s.responseInds(k)).licksL] < (sessionData(s.responseInds(k)).respondTime + tf*1000);
            currTrial_lick = sessionData(s.responseInds(k)).licksL(currTrial_lickInd) - sessionData(s.responseInds(k)).respondTime;
        elseif ~isnan(sessionData(s.responseInds(k)).rewardR)
            currTrial_lickInd = [sessionData(s.responseInds(k)).licksR] < (sessionData(s.responseInds(k)).respondTime + tf*1000);
            currTrial_lick = sessionData(s.responseInds(k)).licksR(currTrial_lickInd) - sessionData(s.responseInds(k)).respondTime;  
        else
            currTrial_lick = 0;
        end
        % raw licks
        allTrial_lick{k} = [currTrial_lick];
        ILIDist = [ILIDist diff(currTrial_lick)];
        % cleaned licks
        ILI = diff(currTrial_lick);
        if ~isnan(ILI)
            while sum(~isnan(ILI))>0 && ILI(1) <= minILI
                currTrial_lick = currTrial_lick(currTrial_lick~=currTrial_lick(2));
                ILI = diff(currTrial_lick);
            end
            if ~isempty(ILI)
                lickDiff(k) = ILI(1);
            else
                lickDiff(k) = NaN;
            end
        else
            lickDiff(k) = NaN;
        end
        if lickDiff(k) > s.rwdDelay
            lickDiff(k) = NaN;
        end
        % get post-rwd licks
        lickNumPostRwd(k) = sum(currTrial_lick>=(s.rwdDelay+50)&currTrial_lick<=(s.rwdDelay+50+postRwdTime));
        allTrial_lickClean{k} = [currTrial_lick];
        ILIDistCleaned = [ILIDistCleaned diff(currTrial_lick)];
    end
    svs = NaN(1,length(s.allChoices));
    svs(s.changeChoice_Inds-1) = 1;
    svs(s.stayChoice_Inds-1) = 0;
    lickNumPostRwd = zscore(lickNumPostRwd);
    lickDiff(~isnan(lickDiff)&s.allChoices>0) = zscore(lickDiff(~isnan(lickDiff)&s.allChoices>0));
    lickDiff(~isnan(lickDiff)&s.allChoices<0) = zscore(lickDiff(~isnan(lickDiff)&s.allChoices<0));
    lickDiffCombined = [lickDiffCombined lickDiff];
    lickLatCombined = [lickLatCombined s.lickLatZ];
    lickNumPostRwdCombined = [lickNumPostRwdCombined lickNumPostRwd];
    rwdCombined = [rwdCombined abs(s.allRewards)];
    svsCombined = [svsCombined svs];
    laserCombined = [laserCombined s.laser];
    Qcombined = [Qcombined Qchosen'];
    rpeCombined = [rpeCombined t.peOri'];

    figure;
    % raw figure
    subplot(4,4,[1 5 9])
    [~, sortInd] = sort(s.lickLatZ);
    plotSpikeRaster(allTrial_lick(sortInd),'PlotType','vertline'); 
    xlim([-200 2000])
    subplot(4,4, 13)
    edges = 0:20:max(ILIDist)+1;
    histogram(ILIDist, edges, 'FaceColor', 'c', 'EdgeColor', 'none')
    set(gca, 'YScale', 'log')

    % cleaned figure
    subplot(4,4,[2 6 10])
    [~, sortInd] = sort(s.lickLatZ);
    plotSpikeRaster(allTrial_lickClean(sortInd),'PlotType','vertline'); 
    xlim([-200 2000])
    subplot(4,4, 14)
    edges = 0:20:max(ILIDistCleaned)+1;
    histogram(ILIDistCleaned, edges, 'FaceColor', 'm', 'EdgeColor', 'none')
    set(gca, 'YScale', 'log')

    subplot(4,4,[3, 7, 11]); hold on;
    [~, sortInd] = sort(t.pe);    
    plot([-200 2000], [length(s.nrwd_Inds) length(s.nrwd_Inds)], 'LineStyle', '--', 'Color', 'r'); hold on;
    plotSpikeRaster(allTrial_lick(sortInd),'PlotType','vertline'); hold on;
    plot(([s.rwdDelay+50 s.rwdDelay+50]), [0 length(s.allChoices)], 'LineStyle', '--', 'Color', 'r');
    plot(([s.rwdDelay+50+500 s.rwdDelay+50+500]), [0 length(s.allChoices)], 'LineStyle', '--', 'Color', 'r');
    xlim([-200 2000])

    sgtitle(session)
    screen = get(0,'Screensize');
    screen(4) = screen(4) - 100;
    set(gcf, 'Position', screen);
    edges = quantile(s.lickLatLogZ, linspace(0, 1, numBins+1));
%     edges = linspace(min(s.lickLatLogZ)-0.01, max(s.lickLatLogZ)+0.01, numBins+1);
    edges(1) = edges(1)-0.01;
    edges(end) = edges(end) +0.01;
    tempLat = NaN(1,numBins);
    tempDiff = NaN(1,numBins);
    tempNum = NaN(1,numBins);
    for j = 1:numBins
        tempLat(j) = mean(s.lickLatLogZ(s.lickLatLogZ>edges(j)&s.lickLatLogZ<=edges(j+1)), 'omitnan');
        tempDiff(j) = mean(lickDiff(s.lickLatLogZ>edges(j)&s.lickLatLogZ<=edges(j+1)), 'omitnan');
        tempNum(j) = mean(lickNumPostRwd(s.lickLatLogZ>edges(j)&s.lickLatLogZ<=edges(j+1)), 'omitnan');
    end
    lickLatBinned(sess, :) = tempLat;
    lickDiffBinned(sess, :) = tempDiff;
    lickNumBinned(sess, :) = tempNum;
    % bin by RPE
%     edges = quantile(s, linspace(0, 1, numBins+1));
    edges = [linspace(min(t.pe)-0.01, 0, 0.5*numBins+1) linspace(0, max(t.pe)+0.01, 0.5*numBins+1)];
    edges = edges([1:0.5*numBins+1 0.5*numBins+3:end]); 
    tempPe = NaN(1,numBins);
    tempNumPe = NaN(1,numBins);
    for j = 1:numBins
        tempPe(j) = mean(t.pe(t.pe>edges(j)&t.pe<=edges(j+1)), 'omitnan');
        tempNumPe(j) = mean(lickNumPostRwd(t.pe>edges(j)&t.pe<=edges(j+1)), 'omitnan');
        tempNumPeLaser(j) = mean(lickNumPostRwd(t.pe>edges(j)&t.pe<=edges(j+1)&s.laser'==1), 'omitnan');
        tempNumPeNoLaser(j) = mean(lickNumPostRwd(t.pe>edges(j)&t.pe<=edges(j+1)&s.laser'==0), 'omitnan');
       
    end
    peBinned(sess,:) = tempPe;
    lickNumPeBinned(sess,:) = tempNumPe;
    lickNumPeBinnedLaser(sess, :) = tempNumPeLaser;
    lickNumPeBinnedNoLaser(sess, :) = tempNumPeNoLaser;
    % only focus on nrwd trials
    edges = linspace(min(t.pe(s.nrwd_Inds))-0.01, max(t.pe(s.nrwd_Inds))+0.01, 0.5*numBins+1);
    tempPe = NaN(1,0.5*numBins);
    tempNumPe = NaN(1,0.5*numBins);
    tempNumPeLaser = NaN(1,0.5*numBins);
    tempNumPeNoLaser = NaN(1,0.5*numBins);
    for j = 1:0.5*numBins
        tempPe(j) = mean(t.pe(t.pe>edges(j)&t.pe<=edges(j+1) & abs(s.allRewards)'==0), 'omitnan');
        tempNumPe(j) = mean(lickNumPostRwd(t.pe>edges(j)&t.pe<=edges(j+1) & abs(s.allRewards)'==0), 'omitnan');
        tempNumPeLaser(j) = mean(lickNumPostRwd(t.pe>edges(j)&t.pe<=edges(j+1)&s.laser'==1 & abs(s.allRewards)'==0), 'omitnan');
        tempNumPeNoLaser(j) = mean(lickNumPostRwd(t.pe>edges(j)&t.pe<=edges(j+1)&s.laser'==0 & abs(s.allRewards)'==0), 'omitnan');     
    end
    
    peBinnedNrwd(sess,:) = tempPe;
    lickNumPeBinnedNrwd(sess,:) = tempNumPe;
    lickNumPeBinnedLaserNrwd(sess, :) = tempNumPeLaser;
    lickNumPeBinnedNoLaserNrwd(sess, :) = tempNumPeNoLaser;
    
%     focus on nrwd and use unupdated rpe
    edges = linspace(min(t.peOri(s.nrwd_Inds))-0.01, max(t.peOri(s.nrwd_Inds))+0.01, 0.5*numBins+1);
    tempPe = NaN(1,0.5*numBins);
    tempNumPe = NaN(1,0.5*numBins);
    tempNumPeLaser = NaN(1,0.5*numBins);
    tempNumPeNoLaser = NaN(1,0.5*numBins);
    for j = 1:0.5*numBins
        tempPe(j) = mean(t.peOri(t.peOri>edges(j)&t.peOri<=edges(j+1) & abs(s.allRewards)'==0), 'omitnan');
        tempNumPe(j) = mean(lickNumPostRwd(t.peOri>edges(j)&t.peOri<=edges(j+1) & abs(s.allRewards)'==0), 'omitnan');
        tempNumPeLaser(j) = mean(lickNumPostRwd(t.peOri>edges(j)&t.peOri<=edges(j+1)&s.laser'==1 & abs(s.allRewards)'==0), 'omitnan');
        tempNumPeNoLaser(j) = mean(lickNumPostRwd(t.peOri>edges(j)&t.peOri<=edges(j+1)&s.laser'==0 & abs(s.allRewards)'==0), 'omitnan');     
    end
    
    peOriBinnedNrwd(sess,:) = tempPe;
    lickNumPeOriBinnedNrwd(sess,:) = tempNumPe;
    lickNumPeOriBinnedLaserNrwd(sess, :) = tempNumPeLaser;
    lickNumPeOriBinnedNoLaserNrwd(sess, :) = tempNumPeNoLaser;
end
%% lick diff - lick lat
figure2; hold on;
% scatter(lickLatCombined, lickDiffCombined, 5, [0.7 0.7 0.7], 'filled', 'MarkerFaceAlpha',0.2)
plot(lickLatBinned', lickDiffBinned', 'Color', [0.8 0.8 0.8])
plotFilled(mean(lickLatBinned,'omitnan'), lickDiffBinned, 'm')
xlabel('lick latency, log, zscored');
ylabel('ILI, zscored')
%% plot lickRate-rpe
figure2;
subplot(1, 2, 1);
plotFilled(mean(peBinned(:,1:0.5*numBins), 'omitnan'), lickNumPeBinned(:,1:0.5*numBins), 'r')
hold on;
subplot(1,2,2)
plotFilled(mean(peBinned(:,0.5*numBins+1:end), 'omitnan'), lickNumPeBinned(:,0.5*numBins+1:end), 'b')
%%
lm = fitlm(rpeCombined(rwdCombined==0)', zscore(lickNumPostRwdCombined(rwdCombined==0))');
%% plot lickRate-rpe sep by laser and no laser
figure2;
subplot(1, 2, 1); hold on;
plotFilled(mean(peBinned(:,1:0.5*numBins), 'omitnan'), lickNumPeBinnedLaser(:,1:0.5*numBins), 'b')
plotFilled(mean(peBinned(:,1:0.5*numBins), 'omitnan'), lickNumPeBinnedNoLaser(:,1:0.5*numBins), [0.5 0.5 0.5])
legend({'laser', '', 'no laser', ''})
hold on;
subplot(1,2,2); hold on;
plotFilled(mean(peBinned(:,0.5*numBins+1:end), 'omitnan'), lickNumPeBinnedLaser(:,0.5*numBins+1:end), 'b')
plotFilled(mean(peBinned(:,0.5*numBins+1:end), 'omitnan'), lickNumPeBinnedNoLaser(:,0.5*numBins+1:end), [0.5 0.5 0.5])
%% focus on nrwd trials
%% plot lickRate-rpe
figure2;
plotFilled(mean(peBinnedNrwd(:,1:0.5*numBins), 'omitnan'), lickNumPeBinnedNrwd(:,1:0.5*numBins), 'r')
%% plot lickRate-rpe sep by laser and no laser
figure2;
hold on;
plotFilled(mean(peBinnedNrwd(:,1:0.5*numBins), 'omitnan'), lickNumPeBinnedLaserNrwd(:,1:0.5*numBins), 'b')
plotFilled(mean(peBinnedNrwd(:,1:0.5*numBins), 'omitnan'), lickNumPeBinnedNoLaserNrwd(:,1:0.5*numBins), [0.5 0.5 0.5])
legend({'laser', '', 'no laser', ''})
%% focus on nrwd trials, use un updated rpe
% plot lickRate-rpe
figure2;
plotFilled(mean(peOriBinnedNrwd(:,1:0.5*numBins), 'omitnan'), lickNumPeOriBinnedNrwd(:,1:0.5*numBins), 'r')
%% plot lickRate-rpe sep by laser and no laser
figure2;
hold on;
plotFilled(mean(peOriBinnedNrwd(:,1:0.5*numBins), 'omitnan'), lickNumPeOriBinnedLaserNrwd(:,1:0.5*numBins), 'b')
plotFilled(mean(peOriBinnedNrwd(:,1:0.5*numBins), 'omitnan'), lickNumPeOriBinnedNoLaserNrwd(:,1:0.5*numBins), [0.5 0.5 0.5])
legend({'laser', '', 'no laser', ''})
%%
indLaser = find(laserCombined == 1 & rwdCombined == 0);
indNoLaser =  find(laserCombined == 0 & rwdCombined == 0);
figure2; hold on;
edges = linspace(min(lickNumPostRwdCombined), max(lickNumPostRwdCombined), 50);
histogram(lickNumPostRwdCombined(indLaser),edges, 'Normalization', 'probability', 'FaceColor', 'm', 'EdgeColor', 'none');
histogram(lickNumPostRwdCombined(indNoLaser),edges, 'Normalization', 'probability', 'FaceColor', 'k', 'EdgeColor', 'none');
%% to analyze no lick window distribution 
%%
numBins = 6;
tf = 10;
tb = 10;
minILI = 50;
postRwdTime = 500;
earlyLickRate = [];
ILIDist = [];
ILIDistCleaned = [];
lickDiffCombined =[];
lickLatCombined = [];
lickNumPostRwdCombined = [];
laserCombined = [];
rwdCombined = [];
rpeCombined = [];
svsCombined = [];
Qcombined = [];
lickLatBinned = NaN(length(dayList), numBins);
lickDiffBinned = NaN(length(dayList), numBins);
lickNumBinned = NaN(length(dayList), numBins);
peBinned = NaN(length(dayList), numBins);
lickNumPeBinned = NaN(length(dayList), numBins);
lickNumPeBinnedLaser = NaN(length(dayList), numBins);
lickNumPeBinnedNoLaser = NaN(length(dayList), numBins);
peBinnedNrwd = NaN(length(dayList), 0.5*numBins);
lickNumPeBinnedNrwd = NaN(length(dayList), 0.5*numBins);
lickNumPeBinnedLaserNrwd = NaN(length(dayList), 0.5*numBins);
lickNumPeBinnedNoLaserNrwd = NaN(length(dayList), 0.5*numBins);
peOriBinnedNrwd = NaN(length(dayList), 0.5*numBins);
lickNumPeOriBinnedNrwd  = NaN(length(dayList), 0.5*numBins);
lickNumPeOriBinnedLaserNrwd  = NaN(length(dayList), 0.5*numBins);
lickNumPeOriBinnedNoLaserNrwd  = NaN(length(dayList), 0.5*numBins);
for sess = 1:length(dayList)
    session = dayList{sess};
    s = behAnalysisNoPlot_opMD(session, 'simpleFlag', 1);
    sessionData = s.behSessionData;
    allLicks = sessionData.allLicks;
    allTrial_lick = {};
    for k = 1:length(s.responseInds)
        currLickInds = allLicks>(sessionData(s.responseInds(k)).CSon - tb*1000) & (allLicks<sessionData(s.responseInds(k)).CSon + tf*1000);
        currLicks = allLicks(currLickInds) - sessionData(s.responseInds(k)).CSon;
        allTrial_lick{k} = currLicks;
    end

    figure2;hold on
    % raw figure licks
    [~, sortInd] = sort(s.lickLatZ);
    plotSpikeRaster(allTrial_lick(sortInd),'PlotType','vertline'); 
    title(session)
    hold on
    xlim([-tb*1000 tf*1000])
    plot([0, 0], [0, length(sessionData)], 'Color', 'r', "LineStyle", '--')
    plot([-1000, -1000], [0, length(sessionData)], 'Color', 'b', "LineStyle", '--')
    % P(lick) - time since ITI start, 1s window
    
    


end
%%
plotSessionLick('m689514d20240109')
plotSessionLick('m689514d20231129') 
plotSessionLick('m689514d20231129') 
plotSessionLick('mZS060d20210426')
plotSessionLick('m699462d20240119')
plotSessionLick('m684890d20231211')
%%
function plotSessionLick(session)
    tb = 3;
    tf = 5;
    binSize = 100;
    stepSize = 50;
    edgesLick = [-tb*1000+binSize*0.5]:stepSize:[tf*1000-binSize*0.5];
    s = behAnalysisNoPlot_opMD(session, 'simpleFlag', 1);
    sessionData = s.behSessionData;
    allLicks = sessionData.allLicks;
    allTrial_lick = {};
    allLick_mat = zeros(length(sessionData), length(edgesLick));
    for k = 1:length(s.responseInds)
        currLickInds = allLicks>(sessionData(s.responseInds(k)).CSon - tb*1000) & (allLicks<sessionData(s.responseInds(k)).CSon + tf*1000);
        currLicks = allLicks(currLickInds) - sessionData(s.responseInds(k)).CSon;
        allTrial_lick{k} = currLicks;
        countsPre = searchsort(currLicks, edgesLick-binSize*0.5);
        countsPost = searchsort(currLicks, edgesLick+binSize*0.5);
        allLick_mat(k,:) = 1000 * (countsPost - countsPre)/binSize;
    end
    figure; 
    subplot(3, 1, 2); hold on
    % raw figure licks
    [~, sortInd] = sort(s.lickLatZ);
    LineFormat = struct;
    LineFormat.LineWidth = 2;
    plotSpikeRaster(allTrial_lick(sortInd),'PlotType','vertline', 'LineFormat', LineFormat); 
    title(session)
    hold on
    xlim([-tb*1000 tf*1000])
    plot([0, 0], [0, length(sessionData)], 'Color', 'r', "LineStyle", '--')
    plot([-1000, -1000], [0, length(sessionData)], 'Color', 'b', "LineStyle", '--')

    subplot(3, 1, 1); hold on
    edges = 0:20:2000;
    small = find(diff(allLicks)<50) + 1;
    allLicks(small) = NaN;
    allLicks = allLicks(~isnan(allLicks));
    histogram(diff(allLicks), edges)
    title('inter lick interval')  

    subplot(3, 1, 3); hold on
    plotFilled(edgesLick, allLick_mat, 'k')
    xlim([-tb*1000 tf*1000])
    title('Mean lick rate')
    xlabel('Time from go cue')

    sgtitle(session)

end
%%
function locs = searchsort(x, edges)
    locs = zeros(size(edges));
    for i = 1:length(locs)
        locs(i) = sum(x<=edges(i));
    end
end


    
