ani = 'allGt';
sheet = 'inhibitionGt';
col = 'cueOnShamLate';
dayList = getDayList(sheet, ani, col);
[root, sep] = currComputer();
%% get rid of days with big bias
modelName = '5params';
bias = cell(length(dayList),1);
for i = 1:length(dayList)
    session = dayList{i};
    pd = parseSessionString_df(session, root, sep);
    sampFile = [pd.animalName col '_' modelName];
    path = [root pd.animalName sep pd.animalName 'sorted' sep 'stan' sep 'bernoulli' sep modelName sep col sep];
    [t,~,noSession] = getStanModelParams_samps(modelName, [path sampFile '.mat'], 4000, 'sessionName', session);
    bias{i} = t.params(:,end);
    figure2;
    subplot(1,2,1)
    histogram(bias{i})
    title('bias');
    subplot(1,2,2)
    histogram(t.params(:,1))
    title('aN');
    sgtitle(session);
    
    if abs(mean(bias{i})) >= 0.75
        dayList{i} = {};
    end
end
dayList = dayList(~cellfun(@isempty, dayList));
%%
switchCurrControl = [];
switchCurrLaser = [];
lcCurrControl = [];
lcCurrLaser = [];
wsCurrControl = [];
wsCurrLaser = [];
lcControl = [];
lcControlAfterNrwd = [];
lcControlAfterRwd = [];
wsControl = [];
lcLaser = [];
lcLaserlAfterNrwd = [];
lcLaserlAfterRwd = [];
wsLaser = [];
missControl = [];
missLaser = [];
crControl = [];
crLaser = [];
meanLatCtrl = [];
meanLatLaser = [];
meanLickRateRwdLaser = [];
meanLickRateRwdCtrl = [];
meanLickRateNoRwdLaser = [];
meanLickRateNoRwdCtrl = [];
lenC = [];
lenL = [];
numNrwdLaser = zeros(size(dayList));
numNrwdControl = zeros(size(dayList));

for i = 1:length(dayList)
    os = behAnalysisNoPlot_opMD(dayList{i},  'simpleFlag', 1);
    svs = zeros(size(os.allChoices));
    svs(os.changeChoice_Inds) = 1;
    laserCue = [os.behSessionData.laser];   
    laserNrwdInd = intersect(os.nrwd_Inds, find(os.laser==1));
    ctrlNrwdInd = intersect(os.nrwd_Inds, find(os.laser==0));
    laserRwdInd = intersect(os.rwd_Inds, find(os.laser==1));
    ctrlRwdInd = intersect(os.rwd_Inds, find(os.laser==0));
    
    lcCtemp = length(mintersect(os.nrwd_Inds+1, find(os.laser==0)+1, os.changeChoice_Inds))/length(ctrlNrwdInd);
    lcLtemp = length(mintersect(os.nrwd_Inds+1, find(os.laser==1)+1, os.changeChoice_Inds))/length(laserNrwdInd);
    
%     lcCtemp = sum(ismember(ctrlNrwdInd+1, os.changeChoice_Inds))/length(ctrlNrwdInd);
%     lcLtemp = sum(ismember(laserNrwdInd+1, os.changeChoice_Inds))/length(laserNrwdInd);
    
    
    lcCtempRwd = length(mintersect(os.rwd_Inds+2, ctrlNrwdInd+1, os.changeChoice_Inds, os.stayChoice_Inds+1))/length(mintersect(os.rwd_Inds+1, ctrlNrwdInd, os.stayChoice_Inds));
    lcLtempRwd = length(mintersect(os.rwd_Inds+2, laserNrwdInd+1, os.changeChoice_Inds, os.stayChoice_Inds+1))/length(mintersect(os.rwd_Inds+1, laserNrwdInd,  os.stayChoice_Inds));
    lcCtempNrwd = length(mintersect(os.nrwd_Inds+2, ctrlNrwdInd+1, os.changeChoice_Inds, os.stayChoice_Inds+1))/length(mintersect(os.nrwd_Inds+1, ctrlNrwdInd,  os.stayChoice_Inds));
    lcLtempNrwd = length(mintersect(os.nrwd_Inds+2, laserNrwdInd+1, os.changeChoice_Inds, os.stayChoice_Inds+1))/length(mintersect(os.nrwd_Inds+1, laserNrwdInd,  os.stayChoice_Inds));
% 
%     lenC(i) = length(mintersect(os.rwd_Inds+1, ctrlNrwdInd, os.stayChoice_Inds));
%     lenL(i) = length(mintersect(os.rwd_Inds+1, laserNrwdInd, os.stayChoice_Inds));
    lenC(i) = length(ctrlNrwdInd);
    lenL(i) = length(laserNrwdInd);
    
    wsCtemp = sum(ismember(intersect(os.rwd_Inds, find(os.laser==0))+1, os.stayChoice_Inds))/length(intersect(os.rwd_Inds, find(os.laser==0)));
    wsLtemp = sum(ismember(intersect(os.rwd_Inds, find(os.laser>0))+1, os.stayChoice_Inds))/length(intersect(os.rwd_Inds, find(os.laser>0)));
    
    missCtemp = 1 - length(intersect(find(laserCue==0 & os.CSplus), os.responseInds))/sum(laserCue==0 & os.CSplus);
    missLtemp = 1 - length(intersect(find(laserCue > 0 & os.CSplus), os.responseInds))/sum(laserCue > 0 & os.CSplus);

    crCtemp = sum([false laserCue(1:end-1)] == 0 & os.CSminus & isnan(os.lickSide))/sum([false laserCue(1:end-1)]==0 & os.CSminus);
    crLtemp = sum([false laserCue(1:end-1)] > 0 & os.CSminus & isnan(os.lickSide))/sum([false laserCue(1:end-1)] > 0 & os.CSminus);
    
    
    numNrwdLaser(i) = length(intersect(os.nrwd_Inds, find(os.laser==1)));
    numNrwdControl(i) = length(intersect(os.nrwd_Inds, find(os.laser==0)));
    
    switchCurrControl(i) = sum(svs(os.laser==0))/sum(os.laser==0);
    switchCurrLaser(i) = sum(svs(os.laser==1))/sum(os.laser==1);
    
    lcCurrControl(i) = length(mintersect(os.nrwd_Inds+1, find(os.laser==0), os.changeChoice_Inds))/... 
                        length(intersect(os.nrwd_Inds+1, find(os.laser==0)));
    lcCurrLaser(i) = length(mintersect(os.nrwd_Inds+1, find(os.laser==1), os.changeChoice_Inds))/... 
                        length(intersect(os.nrwd_Inds+1, find(os.laser==1)));
                    
    wsCurrControl(i) = length(mintersect(os.rwd_Inds+1, find(os.laser==0), os.stayChoice_Inds))/... 
                        length(intersect(os.rwd_Inds+1, find(os.laser==0)));
    wsCurrLaser(i) = length(mintersect(os.rwd_Inds+1, find(os.laser==1), os.stayChoice_Inds))/... 
                        length(intersect(os.rwd_Inds+1, find(os.laser==1)));                   

    lcControl(i) = lcCtemp;
    lcLaser(i) = lcLtemp;
    
    lcControlAfterNrwd(i) = lcCtempNrwd;
    lcLaserlAfterNrwd(i) = lcLtempNrwd;
    
    lcControlAfterRwd(i) = lcCtempRwd;
    lcLaserlAfterRwd(i) = lcLtempRwd;
    
    wsControl(i) = wsCtemp;
    wsLaser(i) = wsLtemp;
     
    missControl(i) = missCtemp;
    missLaser(i) = missLtemp;
    
    crControl(i)= crCtemp;
    crLaser(i)= crLtemp;
    
    
    meanLatCtrl(i) = mean(os.lickLat(os.laser==0), 'omitnan');
    meanLatLaser(i) = mean(os.lickLat(os.laser==1), 'omitnan');
    
    meanLatCtrlNext(i) = mean(os.lickLat(intersect(find(os.laser==0)+1, os.changeChoice_Inds)), 'omitnan');
    meanLatLaserNext(i) = mean(os.lickLat(intersect(find(os.laser==1)+1, os.changeChoice_Inds)), 'omitnan');
    
    meanLatCtrlNrwdNext(i) = mean(os.lickLat([false os.laser(1:end-1)==0]&[false os.allRewards(1:end-1)==0]), 'omitnan');
    meanLatLaserNrwdNext(i) = mean(os.lickLat([false os.laser(1:end-1)==1]&[false os.allRewards(1:end-1)==0]), 'omitnan');

    meanLatCtrlNrwdStay(i) = mean(os.lickLat([false os.laser(1:end-1)==0]&[false os.allRewards(1:end-1)==0] & svs==0), 'omitnan');
    meanLatLaserNrwdStay(i) = mean(os.lickLat([false os.laser(1:end-1)==1]&[false os.allRewards(1:end-1)==0] & svs==0), 'omitnan');

    meanLatCtrlNrwdChange(i) = mean(os.lickLat([false os.laser(1:end-1)==0]&[false os.allRewards(1:end-1)==0] & svs==1), 'omitnan');
    meanLatLaserNrwdChange(i) = mean(os.lickLat([false os.laser(1:end-1)==1]&[false os.allRewards(1:end-1)==0] & svs==1), 'omitnan');
    
    meanLickRateRwdLaser(i) = mean(os.lickNumRwd(os.laser==1 & abs(os.allRewards)==1), 'omitnan');
    meanLickRateRwdCtrl(i) = mean(os.lickNumRwd(os.laser==0 & abs(os.allRewards)==1), 'omitnan');
    
    meanLickRateNoRwdLaser(i) = mean(os.lickNumRwd(os.laser==1 & abs(os.allRewards)==0), 'omitnan');
    meanLickRateNoRwdCtrl(i) = mean(os.lickNumRwd(os.laser==0 & abs(os.allRewards)==0), 'omitnan');
end
%% curr switch
figure; hold on;
scatter(switchCurrControl, switchCurrLaser, 30, 'r', 'filled')
plot([0 0.25], [0 0.25],'LineStyle', '--', 'LineWidth', 2, 'Color', [0.6 0.6 0.6])
p = signrank(switchCurrControl, switchCurrLaser); 
set(gca, 'TickDir', 'Out')
set(gca, 'Box', 'off')
set(gca, 'XTick', [0 0.1 0.2], 'FontSize', 14)
set(gca, 'YTick', [0 0.1 0.2], 'FontSize', 14)
xlabel('Control', 'FontSize', 18)
ylim([0 0.25]); xlim([0 0.25])
title(['pSwitchCurr' num2str(p,3)], 'FontSize', 13)
xlabel('Control', 'FontSize', 18)
ylabel('Laser', 'FontSize', 18)
%% curr nrwd switch
figure; hold on;
scatter(lcCurrControl, lcCurrLaser, 30, 'r', 'filled')
plot([0 0.5], [0 0.5],'LineStyle', '--', 'LineWidth', 2, 'Color', [0.6 0.6 0.6])
p = signrank(lcCurrControl, lcCurrLaser); 
set(gca, 'TickDir', 'Out')
set(gca, 'Box', 'off')
set(gca, 'XTick', 0:0.2:0.6, 'FontSize', 14)
set(gca, 'YTick', 0:0.2:0.6, 'FontSize', 14)
xlabel('Control', 'FontSize', 18)
% ylim([0 0.25]); xlim([0 0.25])
title(['lcCurr' num2str(p,3)], 'FontSize', 13)
xlabel('Control', 'FontSize', 18)
ylabel('Laser', 'FontSize', 18)
%% curr rwd stay
figure; hold on;
scatter(wsCurrControl, wsCurrLaser, 30, 'r', 'filled')
plot([0 1], [0 1],'LineStyle', '--', 'LineWidth', 2, 'Color', [0.6 0.6 0.6])
p = signrank(wsCurrControl, wsCurrLaser); 
set(gca, 'TickDir', 'Out')
set(gca, 'Box', 'off')
set(gca, 'XTick', 0:0.1:1, 'FontSize', 14)
set(gca, 'YTick', 0:0.1:1, 'FontSize', 14)
xlabel('Control', 'FontSize', 18)
ylim([0.8 1]); xlim([0.8 1])
title(['wsCurr' num2str(p,3)], 'FontSize', 13)
xlabel('Control', 'FontSize', 18)
ylabel('Laser', 'FontSize', 18)
%% next trial
figure; hold on;
scatter(lcControl, lcLaser, 30, 'r', 'filled')
plot([0 0.6], [0 0.6],'LineStyle', '--', 'LineWidth', 2, 'Color', [0.6 0.6 0.6])
p = signrank(lcControl, lcLaser); 
set(gca, 'TickDir', 'Out')
set(gca, 'Box', 'off')
set(gca, 'XTick', [0 0.25 0.5], 'FontSize', 14)
set(gca, 'YTick', [0 0.25 0.5], 'FontSize', 14)
xlabel('Control', 'FontSize', 18)
ylim([0 0.65]); xlim([0 0.65])
title(['pSwitch' num2str(p,3)], 'FontSize', 13)
xlabel('Control', 'FontSize', 18)
ylabel('Laser', 'FontSize', 18)

%%  
figure; hold on;
scatter(lcControlAfterNrwd, lcLaserlAfterNrwd, 30, 'r', 'filled')
plot([0 0.6], [0 0.6],'LineStyle', '--', 'LineWidth', 2, 'Color', [0.6 0.6 0.6])
p = signrank(lcControlAfterNrwd, lcLaserlAfterNrwd); 
set(gca, 'TickDir', 'Out')
set(gca, 'Box', 'off')
set(gca, 'XTick', [0 0.25 0.5], 'FontSize', 14)
set(gca, 'YTick', [0 0.25 0.5], 'FontSize', 14)
xlabel('Control', 'FontSize', 18)
ylim([0 0.65]); xlim([0 0.65])
title(['pSwitchAfterNrwd' num2str(p,3)], 'FontSize', 13)
xlabel('Control', 'FontSize', 18)
ylabel('Laser', 'FontSize', 18)
%%
figure; hold on;
scatter(lcControlAfterRwd, lcLaserlAfterRwd, 30, 'r', 'filled')
plot([0 0.6], [0 0.6],'LineStyle', '--', 'LineWidth', 2, 'Color', [0.6 0.6 0.6])
p = signrank(lcControlAfterRwd, lcLaserlAfterRwd); 
set(gca, 'TickDir', 'Out')
set(gca, 'Box', 'off')
set(gca, 'XTick', [0 0.25 0.5], 'FontSize', 14)
set(gca, 'YTick', [0 0.25 0.5], 'FontSize', 14)
xlabel('Control', 'FontSize', 18)
ylim([0 0.65]); xlim([0 0.65])
title(['pSwitchAfterRwd' num2str(p,3)], 'FontSize', 13)
xlabel('Control', 'FontSize', 18)
ylabel('Laser', 'FontSize', 18)
%%
figure; hold on;
scatter(wsControl, wsLaser, 30, 'r', 'filled')
plot([0.8 1.1], [0.8 1.1],'LineStyle', '--', 'LineWidth', 2, 'Color', [0.6 0.6 0.6])
p = signrank(wsControl, wsLaser); 
set(gca, 'TickDir', 'Out')
set(gca, 'Box', 'off')
set(gca, 'XTick', [0.95 1], 'FontSize', 14)
set(gca, 'YTick', [0.95 1], 'FontSize', 14)
xlabel('Control', 'FontSize', 18)
ylim([0.9 1.05]); xlim([0.9 1.05])
title(['pStay' num2str(p,3)], 'FontSize', 13)
xlabel('Control', 'FontSize', 18)
ylabel('Laser', 'FontSize', 18)
%% focus on licks
figure; hold on;
scatter(meanLatCtrl, meanLatLaser, 30, 'r', 'filled');
p = signrank(meanLatCtrl, meanLatLaser); 
xlim([150 350])
ylim([150 350])
plot([0 350], [0 350],'LineStyle', '--', 'LineWidth', 2, 'Color', [0.6 0.6 0.6])
title(['lick lat ' ani ' ' num2str(p)])
%%
figure; hold on;
scatter(meanLickRateRwdCtrl, meanLickRateRwdLaser, 30, 'r', 'filled');
p = signrank(meanLickRateRwdCtrl, meanLickRateRwdLaser); 

plot([3 8], [3 8],'LineStyle', '--', 'LineWidth', 2, 'Color', [0.6 0.6 0.6])
title(['lickRateRwd ' ani ' ' num2str(p)])

%%
figure; hold on;
scatter(meanLickRateNoRwdCtrl, meanLickRateNoRwdLaser, 30, 'r', 'filled');
p = signrank(meanLickRateNoRwdCtrl, meanLickRateNoRwdLaser); 

plot([1 4], [1 4],'LineStyle', '--', 'LineWidth', 2, 'Color', [0.6 0.6 0.6])
title(['lickRateNoRwd ' ani ' ' num2str(p)])

%%
figure; hold on;
scatter(meanLatCtrlNext, meanLatLaserNext, 'r', 'filled');
plot([0 350], [0 350],'LineStyle', '--', 'LineWidth', 2, 'Color', [0.6 0.6 0.6])
p = signrank(meanLatCtrlNext, meanLatLaserNext); 
xlim([150 350])
xlim([150 350])
title(['lick lat next (stay) ' ani ' ' num2str(p)])
xlabel('pre no laser')
ylabel('pre laser')
%%
figure; hold on;
scatter(meanLatCtrlNrwdNext, meanLatLaserNrwdNext, 'r', 'filled');
plot([0 500], [0 500],'LineStyle', '--', 'LineWidth', 2, 'Color', [0.6 0.6 0.6])
title(['lick lat Nrwd next ' ani])
%% laser effect on lick latencies
ani = 'allGt';
dayList = getDayList('inhibitionGt', ani, 'cueOnGood');
latControl = [];
latControlAfterNrwd = [];
latControlAfterRwd = [];
wsControl = [];
latLaser = [];
latLaserlAfterNrwd = [];
latLaserlAfterRwd = [];
latLaserlWSAfterNrwd = [];
latLaserlWStayAfterRwd = [];
wsLaser = [];
missControl = [];
missLaser = [];
crControl = [];
crLaser = [];
meanLatCtrl = [];
meanLatLaser = [];
numNrwdLaser = zeros(size(dayList));
numNrwdControl = zeros(size(dayList));

for i = 1:length(dayList)
    os = behAnalysisNoPlot_opMD(dayList{i},  'simpleFlag', 1);
    svs = zeros(size(os.allChoices));
    svs(os.changeChoice_Inds) = 1;
    laserCue = [os.behSessionData.laser];
    
    laserNrwdInd = intersect(os.nrwd_Inds, find(os.laser==1));
    ctrlNrwdInd = intersect(os.nrwd_Inds, find(os.laser==0));
    laserRwdInd = intersect(os.rwd_Inds, find(os.laser==1));
    ctrlRwdInd = intersect(os.rwd_Inds, find(os.laser==0));
    

    latCtemp = mean(os.lickLat(mintersect(os.nrwd_Inds+1, find(os.laser==0)+1, os.stayChoice_Inds)));
    latLtemp = mean(os.lickLat(mintersect(os.nrwd_Inds+1, find(os.laser==1)+1, os.stayChoice_Inds)));

    
    latCtempRwd = mean(os.lickLat(mintersect(os.rwd_Inds+2, ctrlNrwdInd+1, os.stayChoice_Inds, os.stayChoice_Inds+1)));
    latLtempRwd = mean(os.lickLat(mintersect(os.rwd_Inds+2, laserNrwdInd+1, os.stayChoice_Inds, os.stayChoice_Inds+1)));
    latCtempNrwd = mean(os.lickLat(mintersect(os.nrwd_Inds+2, ctrlNrwdInd+1, os.stayChoice_Inds, os.stayChoice_Inds+1)));
    latLtempNrwd = mean(os.lickLat(mintersect(os.nrwd_Inds+2, laserNrwdInd+1, os.stayChoice_Inds, os.stayChoice_Inds+1)));

    wsCtemp = sum(ismember(intersect(os.rwd_Inds, find(os.laser==0))+1, os.stayChoice_Inds))/length(intersect(os.rwd_Inds, find(os.laser==0)));
    wsLtemp = sum(ismember(intersect(os.rwd_Inds, find(os.laser>0))+1, os.stayChoice_Inds))/length(intersect(os.rwd_Inds, find(os.laser>0)));
    
    missCtemp = 1 - length(intersect(find(laserCue==0 & os.CSplus), os.responseInds))/sum(laserCue==0 & os.CSplus);
    missLtemp = 1 - length(intersect(find(laserCue > 0 & os.CSplus), os.responseInds))/sum(laserCue > 0 & os.CSplus);

    crCtemp = sum([false laserCue(1:end-1)] == 0 & os.CSminus & isnan(os.lickSide))/sum([false laserCue(1:end-1)]==0 & os.CSminus);
    crLtemp = sum([false laserCue(1:end-1)] > 0 & os.CSminus & isnan(os.lickSide))/sum([false laserCue(1:end-1)] > 0 & os.CSminus);
    
    
    numNrwdLaser(i) = length(intersect(os.nrwd_Inds, find(os.laser==1)));
    numNrwdControl(i) = length(intersect(os.nrwd_Inds, find(os.laser==0)));
    
    latControl(i) = latCtemp;
    latLaser(i) = latLtemp;
    
    latControlAfterNrwd(i) = latCtempNrwd;
    latLaserlAfterNrwd(i) = latLtempNrwd;
    
    latControlAfterRwd(i) = latCtempRwd;
    latLaserlAfterRwd(i) = latLtempRwd;
    
    wsControl(i) = wsCtemp;
    wsLaser(i) = wsLtemp;
     
    missControl(i) = missCtemp;
    missLaser(i) = missLtemp;
    
    crControl(i)= crCtemp;
    crLaser(i)= crLtemp;
    
    meanLatCtrl(i) = mean(os.lickLat(os.laser==0), 'omitnan');
    meanLatLaser(i) = mean(os.lickLat(os.laser==1), 'omitnan');
    
    meanLatCtrlNext(i) = mean(os.lickLat([false os.laser(1:end-1)==0]), 'omitnan');
    meanLatLaserNext(i) = mean(os.lickLat([false os.laser(1:end-1)==1]), 'omitnan');
    
    meanLatCtrlNrwdNext(i) = mean(os.lickLat([false os.laser(1:end-1)==0]&[false os.allRewards(1:end-1)==0]), 'omitnan');
    meanLatLaserNrwdNext(i) = mean(os.lickLat([false os.laser(1:end-1)==1]&[false os.allRewards(1:end-1)==0]), 'omitnan');

    meanLatCtrlNrwdStay(i) = mean(os.lickLat([false os.laser(1:end-1)==0]&[false os.allRewards(1:end-1)==0] & svs==0), 'omitnan');
    meanLatLaserNrwdStay(i) = mean(os.lickLat([false os.laser(1:end-1)==1]&[false os.allRewards(1:end-1)==0] & svs==0), 'omitnan');

    meanLatCtrlNrwdChange(i) = mean(os.lickLat([false os.laser(1:end-1)==0]&[false os.allRewards(1:end-1)==0] & svs==1), 'omitnan');
    meanLatLaserNrwdChange(i) = mean(os.lickLat([false os.laser(1:end-1)==1]&[false os.allRewards(1:end-1)==0] & svs==1), 'omitnan');

end
%% figure
figure2;
scatter(latControlAfterRwd(lenC>=50), latLaserlAfterRwd(lenC>=50), 20, 'b', 'filled');
hold on; plot([250 600], [250 600], 'LineStyle', '--', 'LineWidth', 2, 'Color', [0.7 0.7 0.7]); title('after no rwd')
p = signrank(latControlAfterNrwd(lenC>=50), latLaserlAfterNrwd(lenC>=50));
title(['lat on switch trial on N+1, nrwd on N' num2str(p,3)], 'FontSize', 13)
%% consecutive rwd + laser vs control nrwd
len = 2:6;
lickLats = cell(1,length(len));
lickRates = cell(1,length(len));
lickRatesRwd = cell(1,length(len));
lickRatesNoRwd = cell(1,length(len));
lcChoiceCtrl = cell(1,length(len));
lcChoiceLaser = cell(1,length(len));
lcChoiceCtrlNrwd = cell(1,length(len));
lcChoiceLaserNrwd = cell(1,length(len));
wsChoice = cell(1,length(len));
lickLatsAfterNoRwd = cell(1,length(len));
numTrials = zeros(1, length(len));

for sess = 1:length(dayList)
    os = behAnalysisNoPlot_opMD(dayList{sess},'simpleFlag',1);
    switches = zeros(1, length(os.allChoices));
    switches(os.changeChoice_Inds) = 1;
    lickInds = cell(1,length(len));
    lickRwdInds = cell(1, length(len));
    lickNoRwdInds = cell(1, length(len));
    lickIndsNrwd = cell(1,length(len));
    lickNoRwdsNrwdInds = cell(1,length(len));
    lickRwdsRNrwdInds = cell(1,length(len));
    for i = 1:length(len)
        % generate choice and reward history
        choicesHis = conv(os.allChoices, ones(1,len(i)));
        choicesHis = choicesHis(1:end-(len(i)-1)); % number of consecutive choices to current trial
        rewardsHis = conv(abs(os.allRewards), ones(1,len(i)-1));
        rewardsHis = [0 rewardsHis(1:end-(len(i)-1))]; % number of consecutive rewards before current trial
        noRewardsHis = conv(abs(os.allNoRewards), ones(1,len(i)-1));
        noRewardsHis = [0 noRewardsHis(1:end-(len(i)-1))]; % number of consecutive no rewards before current trial
        % detect len consecutive choices and len-1 consecutive rewards
        conChoicesInds = find(abs(choicesHis)>=len(i));
        conPreRewardInds = find(rewardsHis>=(len(i)-1));
        conPreNoRwdInds = find(noRewardsHis>=(len(i)-1));
        lickInds{i} = intersect(conChoicesInds, conPreRewardInds); 
        lickNoRwdInds{i} = intersect(lickInds{i}, os.nrwd_Inds); 
        lickRwdInds{i} = intersect(lickInds{i}, os.rwd_Inds); 
        lickIndsNrwd{i} = intersect(conChoicesInds, conPreNoRwdInds);
        lickNoRwdsNrwdInds{i} = intersect(lickIndsNrwd{i}, os.nrwd_Inds);
        lickRwdsRNrwdInds{i} = intersect(lickIndsNrwd{i}, os.rwd_Inds);
    end
    for i = 1:(length(len)-1)
        lickInds{i} = setxor(lickInds{i},lickInds{i+1}); % take longer ones away from shorter groups
        lickNoRwdInds{i} = setxor(lickNoRwdInds{i},lickNoRwdInds{i+1});
        lickRwdInds{i} = setxor(lickRwdInds{i},lickRwdInds{i+1});
        
        lickIndsNrwd{i} = setxor(lickIndsNrwd{i},lickIndsNrwd{i+1});
        lickNoRwdsNrwdInds{i} = setxor(lickNoRwdsNrwdInds{i},lickNoRwdsNrwdInds{i+1});
        lickRwdsRNrwdInds{i} = setxor(lickRwdsRNrwdInds{i},lickRwdsRNrwdInds{i+1});
    end
    for i = 1:length(len)
        lickLats{i} = [lickLats{i} os.lickLatLogZ(lickInds{i})]; % take lickLatz of certain groups
        lickRates{i} = [lickRates{i} os.lickRateZ(lickInds{i})]; % take lick rate of certain groups
        lickRatesNoRwd{i} = [lickRatesNoRwd{i} os.lickRateRwdZ(lickNoRwdInds{i})]; % take cons lick Rate of certain groups
        lickRatesRwd{i} = [lickRatesRwd{i} os.lickRateRwdZ(lickRwdInds{i})]; % take cons lick Rate of certain groups
        tempInds = lickNoRwdInds{i};
        tempInds = tempInds(tempInds<length(os.allChoices)); % take away the last trial
        lcChoiceCtrl{i} = [lcChoiceCtrl{i} switches(intersect(find(os.laser==0), tempInds)+1)]; % take next choice no laser
        lcChoiceLaser{i} = [lcChoiceLaser{i} switches(intersect(find(os.laser==1), tempInds)+1)]; % take next choice
        lickLatsAfterNoRwd{i} = [lickLatsAfterNoRwd{i} os.lickLatLogZ(tempInds+1)];
        tempInds = lickRwdInds{i};
        tempInds = tempInds(tempInds<length(os.allChoices)); % take away the last trial
        wsChoice{i} = [wsChoice{i} 1-switches(tempInds+1)]; % take next choice
        tempInds = lickNoRwdsNrwdInds{i};
        numTrials(i) = numTrials(i) + length(tempInds);
        tempInds = tempInds(tempInds<length(os.allChoices)); % take away the last trial
        lcChoiceCtrlNrwd{i} = [lcChoiceCtrlNrwd{i} switches(intersect(find(os.laser==0), tempInds)+1)];
        lcChoiceLaserNrwd{i} = [lcChoiceLaserNrwd{i} switches(intersect(find(os.laser==1), tempInds)+1)];
    end
end
%% plot p(switch|nrwd) after rwds
figure2; hold on
tempProb = zeros(1, length(len)-1);
tempProbSem = zeros(1, length(len)-1);

for i = 1:length(len)-1
    tempProb(i) = mean(lcChoiceCtrl{i});
    tempProbSem(i) = sem_bern(lcChoiceCtrl{i});
end
plot(len(1:end-1), tempProb, 'Color', [0.7 0.7 0.7], 'LineWidth', 2)
fill([len(1:end-1), flip(len(1:end-1))], [tempProb-tempProbSem flip(tempProb+tempProbSem)], [0.7 0.7 0.7], 'FaceAlpha', 0.2, 'EdgeColor', 'none');

tempProb = zeros(1, length(len)-1);
tempProbSem = zeros(1, length(len)-1);

for i = 1:length(len)-1
    tempProb(i) = mean(lcChoiceLaser{i});
    tempProbSem(i) = sem_bern(lcChoiceLaser{i});
end
plot(len(1:end-1), tempProb, 'Color', [0 0 1], 'LineWidth', 2)
fill([len(1:end-1), flip(len(1:end-1))], [tempProb-tempProbSem flip(tempProb+tempProbSem)], [0 0 1], 'FaceAlpha', 0.2, 'EdgeColor', 'none');
%% plot p(switch|nrwd) after nrwds
figure2; hold on
tempProb = zeros(1, length(len));
tempProbSem = zeros(1, length(len));

for i = 1:length(len)
    tempProb(i) = mean(lcChoiceCtrlNrwd{i});
    tempProbSem(i) = sem_bern(lcChoiceCtrlNrwd{i});
end
plot(1:length(len), tempProb, 'Color', [0.7 0.7 0.7], 'LineWidth', 2)
fill([1:length(len), flip(1:length(len))], [tempProb-tempProbSem flip(tempProb+tempProbSem)], [0.7 0.7 0.7], 'FaceAlpha', 0.2, 'EdgeColor', 'none');

tempProb = zeros(1, length(len));
tempProbSem = zeros(1, length(len));

for i = 1:length(len)
    tempProb(i) = mean(lcChoiceLaserNrwd{i});
    tempProbSem(i) = sem_bern(lcChoiceLaserNrwd{i});
end
plot(1:length(len), tempProb, 'Color', [0 0 1], 'LineWidth', 2)
fill([1:length(len), flip(1:length(len))], [tempProb-tempProbSem flip(tempProb+tempProbSem)], [0 0 1], 'FaceAlpha', 0.2, 'EdgeColor', 'none');

%% plot P switch given Q chosen
modelName = '5params';
ani = 'allGt';
col = 'cueOnGood';
dayList = getDayList('inhibitionGt', ani, col);
numBins = 3;
QmeansC = zeros(length(dayList), numBins);
pSwitchMeanC = zeros(length(dayList), numBins);
QmeansL = zeros(length(dayList), numBins);
pSwitchMeanL = zeros(length(dayList), numBins);
[root, sep] = currComputer();
%%
for i = 1:length(dayList)
    session = dayList{i};
    pd = parseSessionString_df(session, root, sep);
    sampFile = [pd.animalName col '_' modelName];
    path = [root pd.animalName sep pd.animalName 'sorted' sep 'stan' sep 'bernoulli' sep modelName sep col sep];
    os = behAnalysisNoPlot_opMD(session, 'simpleFlag', 1);
    choice = os.allChoices';
    choice(choice<0) = 0;
    outcome = abs(os.allRewards)';
    responseInds = os.responseInds; 
    ITI = os.timeBtwn(1:length(choice)); 
    preRwd = [NaN abs(os.allRewards(1:end-1))]';
    % behavior
    % switch
    svsTemp = find(choice(2:end) ~= choice(1:end-1)) + 1;
    svs = zeros(1,length(responseInds));
    svs(svsTemp) = 1;
    % model 
    % generate best estimates of parameters
    [t,~,badSessionFlag] = getStanModelParams_samps(modelName, [path sampFile '.mat'], 2000, 'sessionName', session);
    % diff value
    Qdiff = abs(t.Q(:,2)-t.Q(:,1));
    % total value
    Qsum = sum(t.Q,2);
    % chosen valie
    Qchosen  = zeros(length(choice),1);
    Qunchosen  = zeros(length(choice),1);
    for j = 1:length(choice)
        if choice(j)>0
            Qchosen(j) = t.Q(j,2);
            Qunchosen(j) = t.Q(j,1);
        else
            Qchosen(j) = t.Q(j,1);
            Qunchosen(j) = t.Q(j,2);
        end
    end
    % bias
    paramNames = getParamNames_dF(modelName,1);
    biasSide = zeros(size(responseInds))';
    biasInd = contains(paramNames,'bias');
    if mean(t.params(:,biasInd))>0
        biasSide(os.lickR_Inds)=1;
    else
        biasSide(os.lickL_Inds)=1;
    end
    rightSide = zeros(length(choice),1);
    rightSide(os.lickR_Inds) = 1;
    choiceConf = 2*t.probChoice - 1;
    pe = t.pe;
    dawExp = double(t.probChoice <= 0.5);
    % laser
    laser = os.laser;
    % sham
%     laser = rand(1,length(laser));
%     laser(laser>=0.7) = 1;
%     laser(laser<0.7) = 0;
    % bin trials by Qchosen with laser or not
    sessIndRange = 1:length(os.responseInds);
    Qchosen = zscore(Qchosen);
    Qsum = zscore(Qsum);
    Qdiff = zscore(Qdiff);
    Qunchosen = zscore(Qunchosen);
    choiceConf = zscore(choiceConf);
%     figure; 
%     histogram(Qchosen);
%     title(session);
    target = Qdiff;
    edges = linspace(min(target)-0.01, max(target)+0.01, numBins+1);
%     edges = binEqualSize(target, numBins+1);
    % lasered trials
    tempQ = target(mintersect(os.nrwd_Inds, find(laser==1), os.stayChoice_Inds, sessIndRange(1:end-1)));
    tempSw = svs(mintersect(os.nrwd_Inds+1, find(laser==1)+1,os.stayChoice_Inds+1, sessIndRange(2:end)));
%     edges = linspace(min(tempQ)-0.01, max(tempQ)+0.01, numBins+1);
    for j = 1:numBins
        QmeansL(i,j) = mean(tempQ(tempQ>=edges(j) & tempQ<edges(j+1)));
        pSwitchMeanL(i,j) = mean(tempSw(tempQ>=edges(j) & tempQ<edges(j+1)));
    end
    % control trials
    tempQ = target(mintersect(os.nrwd_Inds, find(laser==0), os.stayChoice_Inds, sessIndRange(1:end-1)));
    tempSw = svs(mintersect(os.nrwd_Inds+1, find(laser==0)+1, os.stayChoice_Inds+1, sessIndRange(2:end)));
%     edges = linspace(min(tempQ)-0.01, max(tempQ)+0.01, numBins+1);
    for j = 1:numBins
        QmeansC(i,j) = mean(tempQ(tempQ>=edges(j) & tempQ<edges(j+1)));
        pSwitchMeanC(i,j) = mean(tempSw(tempQ>=edges(j) & tempQ<edges(j+1)));
    end    
    
    
end
%% plot pSwitch against Qchosen, split by laser or not
semSWL = sem(pSwitchMeanL);
semSWC = sem(pSwitchMeanC);
meanSWL = mean(pSwitchMeanL, 'omitnan');
meanSWC = mean(pSwitchMeanC, 'omitnan');
meanQL = mean(QmeansL, 'omitnan');
meanQC = mean(QmeansC, 'omitnan');

figure2; hold on;
plot(meanQC, meanSWC, 'LineWidth',2, 'Color', [0.6 0.6 0.6]);
fill([meanQC flip(meanQC)], [meanSWC-semSWC flip(meanSWC+semSWC)], [0.6 0.6 0.6], 'FaceAlpha', 0.25, 'EdgeColor', 'none');

plot(meanQL, meanSWL, 'LineWidth',2, 'Color', [0 0.5 1]);
fill([meanQL flip(meanQL)], [meanSWL-semSWL flip(meanSWL+semSWL)], [0 0.5 1], 'FaceAlpha', 0.25, 'EdgeColor', 'none');
%% plot all session
for i = 1:length(dayList)
    session = dayList{i};
    s = behAnalysisNoPlot_opMD(session, 'simpleFlag', 1);
    laser = find(s.laser==1);
    
    figure; 
    subplot(3,1,1); hold on;
    plot([1:length(s.responseInds); 1:length(s.responseInds)], 0.5*[zeros(1, length(s.responseInds)); s.allChoices + s.allRewards], 'k');
    scatter(laser, s.allChoices(laser), 20, 'b', 'filled')
    title(session)
    screen = get(0,'Screensize');
    screen(4) = screen(4) - 100;
    set(gcf, 'Position', screen)   
end
%%





%%
ani = 'allGt';
sheet = 'inhibitionGt';
col = 'cueOnShamLate';
dayList = getDayList(sheet, ani, col);
[root, sep] = currComputer();


%% pupil validation
outcomeCombine = [];
laserCombine = [];
outcomeSig = [];
laserSig = []; 
for i = 1:length(dayList)
    session = dayList{i};
    pd = parseSessionString_df(session, root, sep);
    pupilAlignPath = [pd.sortedFolder session '_pupil.mat'];
    if exist(pupilAlignPath, 'file')
       load(pupilAlignPath);
    else
        continue
    end
    stepSize = round(1*FR); % in frames
    binSize = round(0.25*FR);  % in frames
    s = behAnalysisNoPlot_opMD(session, 'simpleFlag', 1);
    svs = NaN(size(s.allChoices));
    svs(s.changeChoice_Inds) = 1;
    svs(s.stayChoice_Inds) = 0;
    midpoints = round(0.5*binSize+1):stepSize:size(sessionPupilCue,2)-round(0.5*binSize)+1;
    pupilSlide = zeros(size(sessionPupilCue,1), length(midpoints));
    combMat = [abs(s.allRewards)', s.laser'];
    for j = 1:length(midpoints)
        pupilSlideTemp = mean(sessionPupilCueZ(:,midpoints(j)-round(0.5*binSize):midpoints(j)+round(0.5*binSize)-1),2,'omitnan');
        lm = fitlm(combMat(qualInd(s.responseInds),:),pupilSlideTemp(s.responseInds(ismember(s.responseInds, find(qualInd>0)))));
        outcomeCombine(i, j) = lm.Coefficients.tStat(2);
        laserCombine(i,j) = lm.Coefficients.tStat(3);
        outcomeSig(i, j) = lm.Coefficients.pValue(2)<0.05;
        laserSig(i,j) = lm.Coefficients.pValue(3)<0.05;    
        pupilSlide(:,j) = pupilSlideTemp;
    end
    pupilSlideAll{i} = pupilSlide;
end
%%
slideTime = linspace(-2+0.5*binSize/FR, 10-0.5*binSize/FR, length(midpoints));
figure;
subplot(1,2,1);
plot(slideTime, mean(outcomeSig));
title('outcome')
subplot(1,2,2)
plot(slideTime, mean(laserSig));
title('laser')
%%
[~, focusWin] = max(abs(laserCombine),[],2, 'linear');
maxCoeff = laserCombine(focusWin);
[~, focusWin] = max(abs(laserCombine),[],2);
figure;
edges = linspace(min(maxCoeff)-0.01, max(maxCoeff)+0.01, 15);
histogram(maxCoeff, edges, 'FaceColor', 'm', 'EdgeColor', 'none');
set(gca, 'TickDir', 'out');
set(gca, 'box', 'off')
xlabel('laser Coeff')
title('all gt sessions')
%%
laserMean = [];
laserSem = [];
noLaserMean = [];
noLaserSem = [];
for i = 1:length(dayList)
    session = dayList{i};
    s = behAnalysisNoPlot_opMD(dayList{i}, 'simpleFlag', 1);
    pd = parseSessionString_df(session, root, sep);
    pupilAlignPath = [pd.sortedFolder session '_pupil.mat'];
    if exist(pupilAlignPath, 'file')
       load(pupilAlignPath);
    else
        continue
    end
    pupilTemp = pupilSlideAll{i};
    pupilTemp = pupilTemp(:, focusWin(i));
    pupilTemp = pupilTemp(s.responseInds);
    laserMean(i) = mean(pupilTemp(s.laser==1 & qualInd(s.responseInds)), 'omitnan');
    laserSem(i) = sem(pupilTemp(s.laser==1 & qualInd(s.responseInds)));
    noLaserMean(i) = mean(pupilTemp(s.laser==0 & qualInd(s.responseInds)), 'omitnan');
    noLaserSem(i) = sem(pupilTemp(s.laser==1 & qualInd(s.responseInds)));
    combMat = [abs(s.allRewards)', s.laser'];
%     lm = fitlm(combMat(qualInd(s.responseInds),:),pupilTemp(s.responseInds(ismember(s.responseInds, find(qualInd>0)))));
end
%%
figure; hold on;
scatter(noLaserMean, laserMean, 25, 'm', 'filled');
plot([-1 0.8], [-1 0.8], 'Color', [0.6 0.6 0.6], 'LineStyle', '--', 'LineWidth', 2);
set(gca, 'TickDir', 'out');
set(gca, 'box', 'off')
xlabel('mean dia without light (zscored)')
ylabel('mean dia with light (zscored)')
%%



















