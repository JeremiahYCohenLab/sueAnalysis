%% load session lists

xlFile = 'inhibitionAll';
sheet = 'ZS071';
col1 = 'control1';
col2 =  'inhibition';

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
end
%%
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
subplot(2,4,1); hold on; 
xbins = 0.7:0.02:1.0;
histogram([cSummary.ws], xbins, 'Normalization','probability', 'FaceColor', 'c');
histogram([iSummary.ws], xbins, 'Normalization','probability', 'FaceColor', 'm');
legend({'control' ,'inhibition'})
title('win-stay')

subplot(2,4,2); hold on; 
xbins = 0:0.05:0.5;
histogram([cSummary.lc], xbins, 'Normalization','probability', 'FaceColor', 'c');
histogram([iSummary.lc], xbins, 'Normalization','probability', 'FaceColor', 'm');
title('loss-change')

subplot(2,4,3); hold on; 
xbins = 0:0.05:1;
histogram([cSummary.hmmStates], xbins, 'Normalization','probability', 'FaceColor', 'c');
histogram([iSummary.hmmStates], xbins, 'Normalization','probability', 'FaceColor', 'm');
title('hmm explore')

subplot(2,4,4); hold on; 
xbins = 0:0.05:1;
histogram([cSummary.ws] - [cSummary.lc], xbins, 'Normalization','probability', 'FaceColor', 'c');
histogram([iSummary.ws] - [iSummary.lc], xbins, 'Normalization','probability', 'FaceColor', 'm');
title('ws - lc')

subplot(2,4,5); hold on; 
xbins = 0:0.05:1;
histogram([cSummary.pSwSw], xbins, 'Normalization','probability', 'FaceColor', 'c');
histogram([iSummary.pSwSw], xbins, 'Normalization','probability', 'FaceColor', 'm');
title('P(switch|switch)')

subplot(2,4,6); hold on;
xbins = 0:0.05:1;
histogram([cSummary.pSwSwNrwd], xbins, 'Normalization','probability', 'FaceColor', 'c');
histogram([iSummary.pSwSwNrwd], xbins, 'Normalization','probability', 'FaceColor', 'm');
title('P(switch|switch,nrwd)')

subplot(2,4,7); hold on;
xbins = -1:0.05:1;
histogram([cSummary.pSwSwNrwd]-[cSummary.pSwSw], xbins, 'Normalization','probability', 'FaceColor', 'c');
histogram([iSummary.pSwSwNrwd]-[iSummary.pSwSw], xbins, 'Normalization','probability', 'FaceColor', 'm');
title('P(switch|switch,nrwd)-P(switch|switch)')

subplot(2,4,8); hold on;
xbins = -1:0.05:1;
histogram([cSummary.pSwSwNrwd]-[cSummary.lc], xbins, 'Normalization','probability', 'FaceColor', 'c');
histogram([iSummary.pSwSwNrwd]-[iSummary.lc], xbins, 'Normalization','probability', 'FaceColor', 'm');
title('P(switch|switch,nrwd)-P(switch|nrwd)')
suptitle(sheet)
%%
figure; 
subplot(2,4,1); hold on; 
plot(indControl', [cSummary.ws], 'c', 'LineWidth', 1.5);
plot(indInhibition, [iSummary.ws], 'm', 'lineWidth', 1.5);
legend({'control' ,'inhibition'})
title('win-stay')

subplot(2,4,2); hold on; 
plot(indControl, [cSummary.lc], 'c', 'LineWidth', 1.5);
plot(indInhibition, [iSummary.lc],  'm', 'LineWidth', 1.5);
title('loss-change')

subplot(2,4,3); hold on; 
plot(indControl, [cSummary.hmmStates],  'c', 'LineWidth', 1.5);
plot(indInhibition, [iSummary.hmmStates], 'm', 'LineWidth', 1.5);
title('hmm explore')

subplot(2,4,4); hold on; 
plot(indControl, [cSummary.ws] - [cSummary.lc],  'c', 'LineWidth', 1.5);
plot(indInhibition, [iSummary.ws] - [iSummary.lc], 'm', 'LineWidth', 1.5);
title('ws - lc')

subplot(2,4,5); hold on; 
plot(indControl, [cSummary.pSwSw],   'c', 'LineWidth', 1.5);
plot(indInhibition, [iSummary.pSwSw], 'm', 'LineWidth', 1.5);
title('P(switch|switch)')

subplot(2,4,6); hold on;
plot(indControl, [cSummary.pSwSwNrwd],  'c', 'LineWidth', 1.5);
plot(indInhibition, [iSummary.pSwSwNrwd], 'm', 'LineWidth', 1.5);
title('P(switch|switch,nrwd)')

subplot(2,4,7); hold on;
plot(indControl, [cSummary.pSwSwNrwd]-[cSummary.pSwSw],  'c', 'LineWidth', 1.5);
plot(indInhibition, [iSummary.pSwSwNrwd]-[iSummary.pSwSw], 'm', 'LineWidth', 1.5);
title('P(switch|switch,nrwd)-P(switch|switch)')

subplot(2,4,8); hold on;
plot(indControl, [cSummary.pSwSwNrwd]-[cSummary.lc],   'c', 'LineWidth', 1.5);
plot(indInhibition, [iSummary.pSwSwNrwd]-[iSummary.lc], 'm', 'LineWidth', 1.5);
title('P(switch|switch,nrwd)-P(switch|nrwd)')
suptitle(sheet)
%%

% for i = 1:length(control)
%     session = control{i};
%     fitHmmOpt(session,1);
%     
% end
%%