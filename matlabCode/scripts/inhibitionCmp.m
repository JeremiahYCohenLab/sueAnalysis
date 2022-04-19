%% load session lists

xlFile = 'dynaVol.xlsx';
sheet = 'ZS076';
col1 = 'dmso';
col2 =  'CNO';

control = getDayList(xlFile, sheet, col1);
inhibition = getDayList(xlFile, sheet, col2);
%%
% control
cSummary = struct();
for i = 1:length(control)
    os = behAnalysisNoPlot_opMD(control{i},'simpleFlag',0);
    cSummary(i).wc = length(intersect(os.rwd_Inds(os.rwd_Inds<length(os.responseInds))+1, os.changeChoice_Inds))/sum(os.rwd_Inds<length(os.responseInds));
    cSummary(i).ws = length(intersect(os.rwd_Inds(os.rwd_Inds<length(os.responseInds))+1, os.stayChoice_Inds))/sum(os.rwd_Inds<length(os.responseInds));
    cSummary(i).lc = length(intersect(os.nrwd_Inds(os.nrwd_Inds<length(os.responseInds))+1, os.changeChoice_Inds))/sum(os.nrwd_Inds<length(os.responseInds));
    cSummary(i).ls = length(intersect(os.nrwd_Inds(os.nrwd_Inds<length(os.responseInds))+1, os.stayChoice_Inds))/sum(os.nrwd_Inds<length(os.responseInds));
    cSummary(i).hmmStates = sum(os.hmmStates==1)/length(os.hmmStates);
    cSummary(i).pStaySw = length(intersect(os.changeChoice_Inds(os.changeChoice_Inds<length(os.responseInds))+1, os.stayChoice_Inds))/sum(os.changeChoice_Inds<length(os.responseInds));
    cSummary(i).pSwSw = length(intersect(os.changeChoice_Inds(os.changeChoice_Inds<length(os.responseInds))+1, os.changeChoice_Inds))/sum(os.changeChoice_Inds<length(os.responseInds));
    cSummary(i).pSwSwNrwd =  length(intersect(intersect(os.changeChoice_Inds(os.changeChoice_Inds<length(os.responseInds)), os.nrwd_Inds)+1, os.changeChoice_Inds))/length(intersect(os.changeChoice_Inds(os.changeChoice_Inds<length(os.responseInds)), os.nrwd_Inds));
    cSummary(i).pStSwRwd =  length(intersect(intersect(os.changeChoice_Inds(os.changeChoice_Inds<length(os.responseInds)), os.rwd_Inds)+1, os.stayChoice_Inds))/length(intersect(os.changeChoice_Inds(os.changeChoice_Inds<length(os.responseInds)), os.rwd_Inds));
end
%% inhibition
iSummary = struct();
for i = 1:length(inhibition)
    os = behAnalysisNoPlot_opMD(inhibition{i}, 'simpleFlag', 0);
    iSummary(i).wc = length(intersect(os.rwd_Inds(os.rwd_Inds<length(os.responseInds))+1, os.changeChoice_Inds))/sum(os.rwd_Inds<length(os.responseInds));
    iSummary(i).ws = length(intersect(os.rwd_Inds(os.rwd_Inds<length(os.responseInds))+1, os.stayChoice_Inds))/sum(os.rwd_Inds<length(os.responseInds));
    iSummary(i).lc = length(intersect(os.nrwd_Inds(os.nrwd_Inds<length(os.responseInds))+1, os.changeChoice_Inds))/sum(os.nrwd_Inds<length(os.responseInds));
    iSummary(i).ls = length(intersect(os.nrwd_Inds(os.nrwd_Inds<length(os.responseInds))+1, os.stayChoice_Inds))/sum(os.nrwd_Inds<length(os.responseInds));
    iSummary(i).hmmStates = sum(os.hmmStates==1)/length(os.hmmStates);
    iSummary(i).pStaySw = length(intersect(os.changeChoice_Inds(os.changeChoice_Inds<length(os.responseInds))+1, os.stayChoice_Inds))/sum(os.changeChoice_Inds<length(os.responseInds));
    iSummary(i).pSwSw = length(intersect(os.changeChoice_Inds(os.changeChoice_Inds<length(os.responseInds))+1, os.changeChoice_Inds))/sum(os.changeChoice_Inds<length(os.responseInds));
    iSummary(i).pSwSwNrwd =  length(intersect(intersect(os.changeChoice_Inds(os.changeChoice_Inds<length(os.responseInds)), os.nrwd_Inds)+1, os.changeChoice_Inds))/length(intersect(os.changeChoice_Inds(os.changeChoice_Inds<length(os.responseInds)), os.nrwd_Inds));
    iSummary(i).pStSwRwd =  length(intersect(intersect(os.changeChoice_Inds(os.changeChoice_Inds<length(os.responseInds)), os.rwd_Inds)+1, os.stayChoice_Inds))/length(intersect(os.changeChoice_Inds(os.changeChoice_Inds<length(os.responseInds)), os.rwd_Inds));
end
%% day inds
dayList = [control; inhibition];
dates = zeros(size(dayList));
for i = 1:length(dates)
    tmp = dayList{i};
    tmp = str2double(tmp(end-3:end));
    if isnan(tmp)
        tmp = dayList{i};
        tmp = str2double(tmp(end-4:end-1));
    end
    dates(i) = tmp;
end

[~, ind] = sort(dates);
ind(ind) = 1:length(dates);
indControl = ind(1:length(control));
indInhibition =  ind(length(control)+1:end);
%%
figure; 
subplot(3,4,1); hold on; 
xbins = 0.7:0.02:1.0;
histogram([cSummary.ws], xbins, 'Normalization','probability', 'FaceColor', 'c');
histogram([iSummary.ws], xbins, 'Normalization','probability', 'FaceColor', 'm');
legend({'control' ,'inhibition'})
title('win-stay')

subplot(3,4,2); hold on; 
xbins = 0:0.05:0.5;
histogram([cSummary.lc], xbins, 'Normalization','probability', 'FaceColor', 'c');
histogram([iSummary.lc], xbins, 'Normalization','probability', 'FaceColor', 'm');
title('loss-change')

subplot(3,4,3); hold on; 
xbins = 0:0.05:1;
histogram([cSummary.hmmStates], xbins, 'Normalization','probability', 'FaceColor', 'c');
histogram([iSummary.hmmStates], xbins, 'Normalization','probability', 'FaceColor', 'm');
title('hmm explore')

subplot(3,4,4); hold on; 
xbins = 0:0.05:1;
histogram([cSummary.ws] - [cSummary.lc], xbins, 'Normalization','probability', 'FaceColor', 'c');
histogram([iSummary.ws] - [iSummary.lc], xbins, 'Normalization','probability', 'FaceColor', 'm');
title('ws - lc')

subplot(3,4,5); hold on; 
xbins = 0:0.05:1;
histogram([cSummary.pSwSw], xbins, 'Normalization','probability', 'FaceColor', 'c');
histogram([iSummary.pSwSw], xbins, 'Normalization','probability', 'FaceColor', 'm');
title('P(switch|switch)')

subplot(3,4,6); hold on;
xbins = 0:0.05:1;
histogram([cSummary.pSwSwNrwd], xbins, 'Normalization','probability', 'FaceColor', 'c');
histogram([iSummary.pSwSwNrwd], xbins, 'Normalization','probability', 'FaceColor', 'm');
title('P(switch|switch,nrwd)')

subplot(3,4,7); hold on;
xbins = -1:0.05:1;
histogram([cSummary.pSwSwNrwd]-[cSummary.pSwSw], xbins, 'Normalization','probability', 'FaceColor', 'c');
histogram([iSummary.pSwSwNrwd]-[iSummary.pSwSw], xbins, 'Normalization','probability', 'FaceColor', 'm');
title('P(switch|switch,nrwd)-P(switch|switch)')

subplot(3,4,8); hold on;
xbins = -1:0.05:1;
histogram([cSummary.pSwSwNrwd]-[cSummary.lc], xbins, 'Normalization','probability', 'FaceColor', 'c');
histogram([iSummary.pSwSwNrwd]-[iSummary.lc], xbins, 'Normalization','probability', 'FaceColor', 'm');
title('P(switch|switch,nrwd)-P(switch|nrwd)')
suptitle(sheet)

subplot(3,4,9); hold on; 
xbins = 0:0.05:1;
histogram([cSummary.pStaySw], xbins, 'Normalization','probability', 'FaceColor', 'c');
histogram([iSummary.pStaySw], xbins, 'Normalization','probability', 'FaceColor', 'm');
title('P(stay|switch)')

subplot(3,4,10); hold on;
xbins = 0:0.05:1;
histogram([cSummary.pStSwRwd], xbins, 'Normalization','probability', 'FaceColor', 'c');
histogram([iSummary.pStSwRwd], xbins, 'Normalization','probability', 'FaceColor', 'm');
title('P(stay|switch,rwd)')

subplot(3,4,11); hold on;
xbins = -1:0.05:1;
histogram([cSummary.pStSwRwd]-[cSummary.pStaySw], xbins, 'Normalization','probability', 'FaceColor', 'c');
histogram([iSummary.pStSwRwd]-[iSummary.pStaySw], xbins, 'Normalization','probability', 'FaceColor', 'm');
title('P(stay|switch,rwd)-P(stay|switch)')

subplot(3,4,12); hold on;
xbins = -1:0.05:1;
histogram([cSummary.pStSwRwd]-[cSummary.ws], xbins, 'Normalization','probability', 'FaceColor', 'c');
histogram([iSummary.pStSwRwd]-[iSummary.ws], xbins, 'Normalization','probability', 'FaceColor', 'm');
title('P(stay|switch,rwd)-P(stay|rwd)')
suptitle(sheet)

screen = get(0,'Screensize');
screen(4) = screen(4) - 100;
set(gcf, 'Position', screen)
%%
figure; % only for single animal
subplot(2,4,1); hold on; 
scatter(indControl', [cSummary.ws], 12, 'c', 'filled');
scatter(indInhibition, [iSummary.ws], 12, 'm', 'filled');
legend({'control' ,'inhibition'})
title('win-stay')

subplot(2,4,2); hold on; 
scatter(indControl, [cSummary.lc], 12, 'c', 'filled');
scatter(indInhibition, [iSummary.lc], 12, 'm', 'filled');
title('loss-change')

subplot(2,4,3); hold on; 
scatter(indControl, [cSummary.hmmStates], 12, 'c', 'filled');
scatter(indInhibition, [iSummary.hmmStates], 12, 'm', 'filled');
title('hmm explore')

subplot(2,4,4); hold on; 
scatter(indControl, [cSummary.ws] - [cSummary.lc], 12, 'c', 'filled');
scatter(indInhibition, [iSummary.ws] - [iSummary.lc], 12, 'm', 'filled');
title('ws - lc')

subplot(2,4,5); hold on; 
scatter(indControl, [cSummary.pSwSw], 12, 'c', 'filled');
scatter(indInhibition, [iSummary.pSwSw], 12, 'm', 'filled');
title('P(switch|switch)')

subplot(2,4,6); hold on;
scatter(indControl, [cSummary.pSwSwNrwd], 12, 'c', 'filled');
scatter(indInhibition, [iSummary.pSwSwNrwd], 12, 'm', 'filled');
title('P(switch|switch,nrwd)')

subplot(2,4,7); hold on;
scatter(indControl, [cSummary.pSwSwNrwd]-[cSummary.pSwSw], 12, 'c', 'filled');
scatter(indInhibition, [iSummary.pSwSwNrwd]-[iSummary.pSwSw], 12, 'm', 'filled');
title('P(switch|switch,nrwd)-P(switch|switch)')

subplot(2,4,8); hold on;
scatter(indControl, [cSummary.pSwSwNrwd]-[cSummary.lc], 12, 'c', 'filled');
scatter(indInhibition, [iSummary.pSwSwNrwd]-[iSummary.lc], 12, 'm', 'filled');
title('P(switch|switch,nrwd)-P(switch|nrwd)')
suptitle(sheet)
%% simulation of control and inhibition
cSummarySim = struct();
for i = 1:length(control)
    sessionName = control{i};
    [animalName, date] = strtok(sessionName, 'd'); 
    animalName = animalName(2:end);
    [params, ~, ~, noSession] = getStanModelParams_sampsOnly(animalName, col1, '5params', 200, 'sessionParamsFlag', 1, 'sessionName', control{i});
    if noSession
        fprintf([dayList{i} ' no good behavior \n'])
        continue
    end
    params = mean(params);
    [allRewards, allChoices] = qLearningModel_simNoPlot('maxTrials', 400, 'params', params,'plotFlag',0);
    os = struct;
    os.rwd_Inds = find(abs(allRewards)>0);
    os.nrwd_Inds = find(abs(allRewards)==0);
    os.changeChoice_Inds = find(allChoices(2:end)~=allChoices(1:end-1))+1;
    os.stayChoice_Inds = find(allChoices(2:end)==allChoices(1:end-1))+1;
    cSummarySim(i).wc = length(intersect(os.rwd_Inds(os.rwd_Inds<length(allChoices))+1, os.changeChoice_Inds))/sum(os.rwd_Inds<length(allChoices));
    cSummarySim(i).ws = length(intersect(os.rwd_Inds(os.rwd_Inds<length(allChoices))+1, os.stayChoice_Inds))/sum(os.rwd_Inds<length(allChoices));
    cSummarySim(i).lc = length(intersect(os.nrwd_Inds(os.nrwd_Inds<length(allChoices))+1, os.changeChoice_Inds))/sum(os.nrwd_Inds<length(allChoices));
    cSummarySim(i).ls = length(intersect(os.nrwd_Inds(os.nrwd_Inds<length(allChoices))+1, os.stayChoice_Inds))/sum(os.nrwd_Inds<length(allChoices));
%     cSummarySim(i).hmmStates = sum(os.hmmStates==1)/length(os.hmmStates);
    cSummarySim(i).pStaySw = length(intersect(os.changeChoice_Inds(os.changeChoice_Inds<length(allChoices))+1, os.stayChoice_Inds))/sum(os.changeChoice_Inds<length(allChoices));
    cSummarySim(i).pSwSw = length(intersect(os.changeChoice_Inds(os.changeChoice_Inds<length(allChoices))+1, os.changeChoice_Inds))/sum(os.changeChoice_Inds<length(allChoices));
    cSummarySim(i).pSwSwNrwd =  length(intersect(intersect(os.changeChoice_Inds(os.changeChoice_Inds<length(allChoices)), os.nrwd_Inds)+1, os.changeChoice_Inds))/length(intersect(os.changeChoice_Inds(os.changeChoice_Inds<length(allChoices)), os.nrwd_Inds));
    cSummarySim(i).pStSwRwd =  length(intersect(intersect(os.changeChoice_Inds(os.changeChoice_Inds<length(allChoices)), os.rwd_Inds)+1, os.stayChoice_Inds))/length(intersect(os.changeChoice_Inds(os.changeChoice_Inds<length(allChoices)), os.rwd_Inds));
end
%% compare simulation and real
figure; 
subplot(3,4,1); hold on; 
xbins = 0.7:0.02:1.0;
histogram([cSummary.ws], xbins, 'Normalization','probability', 'FaceColor', 'c');
histogram([cSummarySim.ws], xbins, 'Normalization','probability', 'FaceColor', 'm');
legend({'control' ,'controlSim'})
title('win-stay')

subplot(3,4,2); hold on; 
xbins = 0:0.05:0.5;
histogram([cSummary.lc], xbins, 'Normalization','probability', 'FaceColor', 'c');
histogram([cSummarySim.lc], xbins, 'Normalization','probability', 'FaceColor', 'm');
title('loss-change')

% subplot(3,4,3); hold on; 
% xbins = 0:0.05:1;
% histogram([cSummary.hmmStates], xbins, 'Normalization','probability', 'FaceColor', 'c');
% histogram([cSummarySim.hmmStates], xbins, 'Normalization','probability', 'FaceColor', 'm');
% title('hmm explore')

subplot(3,4,4); hold on; 
xbins = 0:0.05:1;
histogram([cSummary.ws] - [cSummary.lc], xbins, 'Normalization','probability', 'FaceColor', 'c');
histogram([cSummarySim.ws] - [iSummary.lc], xbins, 'Normalization','probability', 'FaceColor', 'm');
title('ws - lc')

subplot(3,4,5); hold on; 
xbins = 0:0.05:1;
histogram([cSummary.pSwSw], xbins, 'Normalization','probability', 'FaceColor', 'c');
histogram([cSummarySim.pSwSw], xbins, 'Normalization','probability', 'FaceColor', 'm');
title('P(switch|switch)')

subplot(3,4,6); hold on;
xbins = 0:0.05:1;
histogram([cSummary.pSwSwNrwd], xbins, 'Normalization','probability', 'FaceColor', 'c');
histogram([cSummarySim.pSwSwNrwd], xbins, 'Normalization','probability', 'FaceColor', 'm');
title('P(switch|switch,nrwd)')

subplot(3,4,7); hold on;
xbins = -1:0.05:1;
histogram([cSummary.pSwSwNrwd]-[cSummary.pSwSw], xbins, 'Normalization','probability', 'FaceColor', 'c');
histogram([cSummarySim.pSwSwNrwd]-[iSummary.pSwSw], xbins, 'Normalization','probability', 'FaceColor', 'm');
title('P(switch|switch,nrwd)-P(switch|switch)')

subplot(3,4,8); hold on;
xbins = -1:0.05:1;
histogram([cSummary.pSwSwNrwd]-[cSummary.lc], xbins, 'Normalization','probability', 'FaceColor', 'c');
histogram([cSummarySim.pSwSwNrwd]-[iSummary.lc], xbins, 'Normalization','probability', 'FaceColor', 'm');
title('P(switch|switch,nrwd)-P(switch|nrwd)')
suptitle(sheet)

subplot(3,4,9); hold on; 
xbins = 0:0.05:1;
histogram([cSummary.pStaySw], xbins, 'Normalization','probability', 'FaceColor', 'c');
histogram([cSummarySim.pStaySw], xbins, 'Normalization','probability', 'FaceColor', 'm');
title('P(stay|switch)')

subplot(3,4,10); hold on;
xbins = 0:0.05:1;
histogram([cSummary.pStSwRwd], xbins, 'Normalization','probability', 'FaceColor', 'c');
histogram([cSummarySim.pStSwRwd], xbins, 'Normalization','probability', 'FaceColor', 'm');
title('P(stay|switch,rwd)')

subplot(3,4,11); hold on;
xbins = -1:0.05:1;
histogram([cSummary.pStSwRwd]-[cSummary.pStaySw], xbins, 'Normalization','probability', 'FaceColor', 'c');
histogram([cSummarySim.pStSwRwd]-[iSummary.pStaySw], xbins, 'Normalization','probability', 'FaceColor', 'm');
title('P(stay|switch,rwd)-P(stay|switch)')

subplot(3,4,12); hold on;
xbins = -1:0.05:1;
histogram([cSummary.pStSwRwd]-[cSummary.ws], xbins, 'Normalization','probability', 'FaceColor', 'c');
histogram([cSummarySim.pStSwRwd]-[iSummary.ws], xbins, 'Normalization','probability', 'FaceColor', 'm');
title('P(stay|switch,rwd)-P(stay|rwd)')
suptitle(sheet)

screen = get(0,'Screensize');
screen(4) = screen(4) - 100;
set(gcf, 'Position', screen)
%%
figure; % only for single animal
subplot(2,4,1); hold on; 
scatter(indControl', [cSummary.ws], 12, 'c', 'filled');
scatter(indInhibition, [cSummarySim.ws], 12, 'm', 'filled');
legend({'control' ,'controlSim'})
title('win-stay')

subplot(2,4,2); hold on; 
scatter(indControl, [cSummary.lc], 12, 'c', 'filled');
scatter(indControl, [cSummarySim.lc], 12, 'm', 'filled');
title('loss-change')

% subplot(2,4,3); hold on; 
% scatter(indControl, [cSummary.hmmStates], 12, 'c', 'filled');
% scatter(indControl, [cSummarySim.hmmStates], 12, 'm', 'filled');
% title('hmm explore')

subplot(2,4,4); hold on; 
scatter(indControl, [cSummary.ws] - [cSummary.lc], 12, 'c', 'filled');
scatter(indControl, [cSummarySim.ws] - [iSummary.lc], 12, 'm', 'filled');
title('ws - lc')

subplot(2,4,5); hold on; 
scatter(indControl, [cSummary.pSwSw], 12, 'c', 'filled');
scatter(indControl, [cSummarySim.pSwSw], 12, 'm', 'filled');
title('P(switch|switch)')

subplot(2,4,6); hold on;
scatter(indControl, [cSummary.pSwSwNrwd], 12, 'c', 'filled');
scatter(indControl, [cSummarySim.pSwSwNrwd], 12, 'm', 'filled');
title('P(switch|switch,nrwd)')

subplot(2,4,7); hold on;
scatter(indControl, [cSummary.pSwSwNrwd]-[cSummary.pSwSw], 12, 'c', 'filled');
scatter(indControl, [cSummarySim.pSwSwNrwd]-[iSummary.pSwSw], 12, 'm', 'filled');
title('P(switch|switch,nrwd)-P(switch|switch)')

subplot(2,4,8); hold on;
scatter(indControl, [cSummary.pSwSwNrwd]-[cSummary.lc], 12, 'c', 'filled');
scatter(indControl, [cSummarySim.pSwSwNrwd]-[iSummary.lc], 12, 'm', 'filled');
title('P(switch|switch,nrwd)-P(switch|nrwd)')
suptitle(sheet)


% for i = 1:length(control)
%     session = control{i};
%     fitHmmOpt(session,1);
%     
% end
%%