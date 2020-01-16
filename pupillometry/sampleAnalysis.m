session = 'mZS034d20191120';
trainingSession = 'mZS034d20191116';
trainingdate = 'Nov16';
itr = '55000';
[animalName, date] = strtok(session, 'd'); 
animalName = animalName(2:end);
date = date(1:9);
sessionFolder = ['m' animalName date];
[root,sep] = currComputer_df;
videopath = ['C:\Users\zhixiao\Documents\pupil\' animalName sep session '\pupil'];
savePath = [videopath sep 'figures' sep];
skeleton = [videopath sep session 'DLC_resnet50_' trainingSession trainingdate 'shuffle1_' itr '_skeleton.csv'];
position = [videopath sep session 'DLC_resnet50_' trainingSession trainingdate 'shuffle1_' itr '.csv'];
%%
diaRaw = csvread(skeleton,2,0);
positionRaw = csvread(position,3,0);
ll = positionRaw(:,4).*positionRaw(:,7);
ind = ll<0.95; ind = ~ind;
diameter = interp1(diaRaw(ind,1),diaRaw(ind,2),diaRaw(:,1));
c = interp1(positionRaw(ind,1),positionRaw(ind,[2,3,5,6]),positionRaw(:,1));
x = (c(:,1)+c(:,3))/2;
y = (c(:,2)+c(:,4))/2;
figure;hold on;
plot(diaRaw(:,2));
plot(diameter);
%%
pupilStruct = pupil_aviToMat(session);

%%
%load pupil

pathData = parseSessionString_df(session, root, '\');
matFile = [pathData.baseFolder 'pupil\' erase(session, '\') '.mat'];
load(matFile)
clc
fprintf('Loaded %s\n', matFile)
%%
fits = struct;
fits = generatePupilTime_wrapper(session, root, pupilStruct, fits);
appendFits_pupil(session, root, fits);
%%
%load behavior
load(['C:\Users\zhixiao\Documents\pupil\' animalName  sep session '\sorted\session\' session '_sessionData_behav.mat'])
fR = fits.cleaned.timing.(session).frameRate*fits.cleaned.timing.(session).scaling_factor;
fits.cleaned.timing.(session).rwdTiming_inFrames = fR*([behSessionData.rewardTime]-behSessionData(1).CSon)/1000+fits.cleaned.timing.(session).trialTiming_inFrames(1);
%%
len = length(behSessionData);
postlen = round(fR*3);
prelen = round(fR*1);
pupilDiamatrix = NaN(len,postlen+prelen+1);
pupilXmatrix = NaN(len,postlen+prelen+1);
pupilYmatrix = NaN(len,postlen+prelen+1);
tt = round(fits.cleaned.timing.(session).trialTiming_inFrames);
tr = round(fits.cleaned.timing.(session).rwdTiming_inFrames);
for i = 1: len
    if ~(isnan(tt(i)))
     ftemp = tt(i);
     behSessionData(i).pDia = diameter(ftemp-prelen:ftemp+postlen)';
     behSessionData(i).pDiapre = mean(diameter(ftemp-round(fR/2):ftemp));
     behSessionData(i).pX = x(ftemp-prelen:ftemp+postlen)';
     behSessionData(i).pY = y(ftemp-prelen:ftemp+postlen)';
     pupilDiamatrix(i,:) = behSessionData(i).pDia;
     pupilXmatrix(i,:) = behSessionData(i).pX;
     pupilYmatrix(i,:) = behSessionData(i).pY;
    end
end

%%
postlenr = round(fR*2.5);
prelenr = round(fR*1.5);
pupilDiarwdmatrix = NaN(len,postlenr+prelenr+1);
pupilXrwdmatrix = NaN(len,postlenr+prelenr+1);
pupilYrwdmatrix = NaN(len,postlenr+prelenr+1);
for i = 1: len
    if ~isnan(behSessionData(i).rewardTime) && ~isnan(tt(i))
        ftemp = round(tt(i) + fR*((behSessionData(i).rewardTime - behSessionData(i).CSon)/1000));
        behSessionData(i).pDiarwd = diameter(ftemp-prelenr:ftemp+postlenr)';
        behSessionData(i).pDiaprer = mean(diameter(ftemp-round(fR/2):ftemp));
        behSessionData(i).pXrwd = x(ftemp-prelenr:ftemp+postlenr)';
        behSessionData(i).pYrwd = y(ftemp-prelenr:ftemp+postlenr)';
        pupilDiarwdmatrix(i,:) = behSessionData(i).pDiarwd;
        pupilXrwdmatrix(i,:) = behSessionData(i).pXrwd;
        pupilYrwdmatrix(i,:) = behSessionData(i).pYrwd;
    end
end


%%
% matrix
len = length(behSessionData);
postlen = round(fR*3);
prelen = round(fR*1);
pupilDiamatrix = NaN(len,postlen+prelen+1);
pupilXmatrix = NaN(len,postlen+prelen+1);
pupilYmatrix = NaN(len,postlen+prelen+1);
postlenr = round(fR*2.5);
prelenr = round(fR*1.5);
pupilDiarwdmatrix = NaN(len,postlenr+prelenr+1);
pupilXrwdmatrix = NaN(len,postlenr+prelenr+1);
pupilYrwdmatrix = NaN(len,postlenr+prelenr+1);
pre = NaN(len,1);
prer = NaN(len,1);
for i = 1: len
    if ~(isnan(tt(i)))
     pupilDiamatrix(i,:) = behSessionData(i).pDia;
     pupilXmatrix(i,:) = behSessionData(i).pX;
     pupilYmatrix(i,:) = behSessionData(i).pY;
     pre(i) = behSessionData(i).pDiapre;
    end
end

for i = 1: len
    if ~isnan(behSessionData(i).rewardTime) && ~isnan(tt(i))
        pupilDiarwdmatrix(i,:) = behSessionData(i).pDiarwd;
        pupilXrwdmatrix(i,:) = behSessionData(i).pXrwd;
        pupilYrwdmatrix(i,:) = behSessionData(i).pYrwd;
        prer(i) = behSessionData(i).pDiaprer;
    end
end
savepath = [root animalName sep sessionFolder sep 'sorted' sep 'session' sep];
save([savepath session '_sessionData_behav.mat'], 'behSessionData','blockSwitch','blockSwitchL','blockSwitchR');

%%
%model

mdIdx = ~isnan([behSessionData.rewardTime]);
mdIdx = find(mdIdx > 0);

for i = 1:length(mdIdx)
    behSessionData(mdIdx(i)).rpe = ms.Q3_5.modOutput.rpe(i);
    behSessionData(mdIdx(i)).qR = ms.Q3_5.modOutput.Q(i,1);
    behSessionData(mdIdx(i)).qL = ms.Q3_5.modOutput.Q(i,2);
end

[~,rpeRank] = sort(ms.Q3_5.modOutput.rpe);
rpeRank = mdIdx(rpeRank);
[~,rpeAbsRank] = sort(abs(ms.Q3_5.modOutput.rpe));
rpeAbsRank = mdIdx(rpeAbsRank);

pDiarwdmd = pupilDiarwdmatrix(rpeRank,:);
zpDiarwdmd = zscore(pDiarwdmd')';
prer = mean(pDiarwdmd(:,round(prelenr - 0.5*fR):prelenr)');
npDiarwdmd = pDiarwdmd./prer';


pDiarwdmdabs = pupilDiarwdmatrix(rpeAbsRank,:);
zpDiarwdmdabs = zscore(pDiarwdmdabs')';
prer = mean(pDiarwdmdabs(:,round(prelenr - 0.5*fR):prelenr)');
npDiarwdmdabs = pDiarwdmdabs./prer';
%%
%idx = ~ (pupilDiarwdmatrix(:,1) == 0);
%pupilDiamatrix = pupilDiamatrix(idx,:);
%pupilXmatrix = pupilXmatrix(idx,:);
%pupilYmatrix = pupilYmatrix(idx,:);
pDiamatrix = pupilDiamatrix(~isnan(pupilDiamatrix(:,1)),:);
normpupilDiamatrix = pupilDiamatrix(~isnan(pupilDiamatrix(:,1)),:) ./ pre(~isnan(pre));
zpupilDiamatrix = zscore(pupilDiamatrix(~isnan(pupilDiamatrix(:,1)),:)')';
zmin = min(zpupilDiamatrix,[],'all');
zmax = max(zpupilDiamatrix,[],'all');
nmin = min(normpupilDiamatrix,[],'all');
nmax = max(normpupilDiamatrix,[],'all');

%%
pDiarwdmatrix = pupilDiarwdmatrix(~isnan(pupilDiarwdmatrix(:,1)),:);
pupilXrwdmatrix = pupilXrwdmatrix(~isnan(pupilDiarwdmatrix(:,1)),:);
pupilYrwdmatrix = pupilYrwdmatrix(~isnan(pupilDiarwdmatrix(:,1)),:);
normpupilDiarwdmatrix = pDiarwdmatrix ./ prer(~isnan(pupilDiarwdmatrix(:,1)));
zpupilDiarwdmatrix = zscore(pupilDiarwdmatrix(~isnan(pupilDiarwdmatrix(:,1)),:)')';
zrmin = min(zpupilDiarwdmatrix,[],'all');
zrmax = max(zpupilDiarwdmatrix,[],'all');
nrmin = min(normpupilDiarwdmatrix,[],'all');
nrmax = max(normpupilDiarwdmatrix,[],'all');
lenr = size(pupilXrwdmatrix,1);
%%
figure;hold on
suptitle(session);
subplot(6,2,[1,3]);
plotmin = min([zmin,zrmin]);
plotmax = max([zmax,zrmax]);
imagesc([-prelen/fR postlen/fR],[1 len],zpupilDiamatrix,[plotmin plotmax])
line([0 0],[0 len],'Color','k');
title('zscore CS alignment')
% line([1.4 1.4],[0 len],'Color','b');

subplot(6,2,5);hold on;
plotFilled(linspace(-prelen/fR, postlen/fR, prelen + postlen + 1),zpupilDiamatrix,[0 0 1]);hold on;
plot(linspace(-prelen/fR, postlen/fR, prelen + postlen + 1),mean(zpupilDiamatrix,1),'Color','b');
xlim([-prelen/fR postlen/fR]);
ylim([plotmin plotmax]);

subplot(6,2,[2,4]);
imagesc([-prelenr/fR postlenr/fR],[1 lenr],zpupilDiarwdmatrix,[plotmin plotmax])
line([0 0],[0 lenr],'Color','k');
title('zscore lick alignment')
% line([1.4 1.4],[0 lenr],'Color','b');

subplot(6,2,6);hold on;
plotFilled(linspace(-prelenr/fR, postlenr/fR, prelenr + postlenr + 1),zpupilDiarwdmatrix,[0 0 1]);hold on;
plot(linspace(-prelenr/fR, postlenr/fR, prelenr + postlenr + 1),mean(zpupilDiarwdmatrix,1),'Color','b');
xlim([-prelenr/fR postlenr/fR]);
ylim([plotmin plotmax]);

subplot(6,2,[7,9]);
nplotmin = min([nmin,nrmin]);
nplotmax = max([nmax,nrmax])-0.1;
imagesc([-prelen/fR postlen/fR],[1 len],normpupilDiamatrix,[nplotmin nplotmax])
line([0 0],[0 len],'Color','k');
title('dilation CS alignment')
% line([1.4 1.4],[0 len],'Color','b');

subplot(6,2,11);hold on;
plotFilled(linspace(-prelen/fR, postlen/fR, prelen + postlen + 1),normpupilDiamatrix,[0 0 1]);hold on;
plot(linspace(-prelen/fR, postlen/fR, prelen + postlen + 1),mean(normpupilDiamatrix,1),'Color','b');
xlim([-prelen/fR postlen/fR]);
ylim([nplotmin nplotmax]);

subplot(6,2,[8,10]);
imagesc([-prelenr/fR postlenr/fR],[1 lenr],normpupilDiarwdmatrix,[nplotmin nplotmax])
line([0 0],[0 lenr],'Color','k');
title('dilation lick alignment')
% line([1.4 1.4],[0 lenr],'Color','b');

subplot(6,2,12);hold on;
plotFilled(linspace(-prelenr/fR, postlenr/fR, prelenr + postlenr + 1),normpupilDiarwdmatrix,[0 0 1]); hold on;
plot(linspace(-prelenr/fR, postlenr/fR, prelenr + postlenr + 1),mean(normpupilDiarwdmatrix,1),'Color','b');
xlim([-prelenr/fR postlenr/fR]);
ylim([nplotmin nplotmax]);
%%

if isempty(dir(savePath))
    mkdir(savePath)
end
saveFigurePDF(gcf,[savePath session '_pupil'])


%%
%plot against model value
figure;hold on
suptitle([session 'rpeRank']);

plotmin = min([zmin,zrmin]);
plotmax = max([zmax,zrmax]);


subplot(2,2,1);

imagesc([-prelenr/fR postlenr/fR],[1 lenr],zpDiarwdmd(~isnan(zpDiarwdmd(:,1)),:),[plotmin plotmax])
line([0 0],[0 lenr],'Color','k');
neg = length(find(ms.Q3_5.modOutput.rpe < 0));
line([-prelenr/fR postlenr/fR],[neg neg])
title('zscore lick alignment, rpeRank')
% line([1.4 1.4],[0 lenr],'Color','b');


nplotmin = min([nmin,nrmin]);
nplotmax = max([nmax,nrmax])-0.1;


subplot(2,2,2);
imagesc([-prelenr/fR postlenr/fR],[1 lenr],npDiarwdmd(~isnan(npDiarwdmd(:,1)),:),[nplotmin nplotmax])
line([0 0],[0 lenr],'Color','k');
line([-prelenr/fR postlenr/fR],[neg neg])
title('dilation lick alignment, rpeRank')
% line([1.4 1.4],[0 lenr],'Color','b');




subplot(2,2,3);
imagesc([-prelenr/fR postlenr/fR],[1 lenr],zpDiarwdmdabs(~isnan(zpDiarwdmdabs(:,1)),:),[plotmin plotmax])
line([0 0],[0 lenr],'Color','k');
title('zscore lick alignment, rpeAbsRank')
% line([1.4 1.4],[0 lenr],'Color','b');

subplot(2,2,4);
imagesc([-prelenr/fR postlenr/fR],[1 lenr],npDiarwdmdabs(~isnan(npDiarwdmdabs(:,1)),:),[nplotmin nplotmax])
line([0 0],[0 lenr],'Color','k');
title('dilation lick alignment, rpeAbsRank')
% line([1.4 1.4],[0 lenr],'Color','b');

%%

if isempty(dir(savePath))
    mkdir(savePath)
end
saveFigurePDF(gcf,[savePath session '_pupilrpe'])

%%
%rwd nrwd
mdIdx = ~isnan([behSessionData.rewardTime]);
mdIdx = find(mdIdx > 0);

responseInds = find(~isnan([behSessionData.rewardTime])); % find CS+ trials with a response in the lick window
allReward_R = [behSessionData(responseInds).rewardR]; 
allReward_L = [behSessionData(responseInds).rewardL]; 
allReward_R = [behSessionData(responseInds).rewardR]; 
allReward_L = [behSessionData(responseInds).rewardL];  
allChoices = NaN(1,length(behSessionData(responseInds)));
allChoices(~isnan(allReward_R)) = 1;
allChoices(~isnan(allReward_L)) = -1;

allReward_R(isnan(allReward_R)) = 0;
allReward_L(isnan(allReward_L)) = 0;
allChoice_R = double(allChoices == 1);
allChoice_L = double(allChoices == -1);

allRewards = zeros(1,length(allChoices));
allRewards(logical(allReward_R)) = 1;
allRewards(logical(allReward_L)) = -1;

allITIs = [behSessionData(responseInds).trialEnd] - [behSessionData(responseInds).CSon];
%%
% rwd vs non rwd stay vs switch raw trace
figure; hold on;
suptitle([session ' raw trace' ' rwd vs nrwd stay vs switch']);
rwdInd = logical(allReward_L + allReward_R);
nrwdInd = ~rwdInd;
% subplot(3,2,2);hold on;
% rwdInd = logical(allReward_L + allReward_R);
% rnormpD = pupilDiamatrix(mdIdx,:);
% rpre = pre(mdIdx);
% rpre = rpre(rwdInd);
% rnormpD = rnormpD(rwdInd,:)./rpre;
% plotFilled(linspace(-prelen/fR, postlen/fR, prelen + postlen + 1),rnormpD,[1 0 0]);

% nrnormpD = pupilDiamatrix(mdIdx,:);
% rpre = pre(mdIdx);
% rpre = rpre(nrwdInd);
% nrnormpD = nrnormpD(nrwdInd,:)./rpre;
% plotFilled(linspace(-prelen/fR, postlen/fR, prelen + postlen + 1),nrnormpD,[0 0 1]);
% 
% xlim([-prelen/fR postlen/fR]);
% title('Dilation');
% legend(['rwdTrial'],[''],['nrwdTrial'],['']);

subplot(4,2,1);hold on;
M1 = pupilDiamatrix(mdIdx,:);
M1 = M1(rwdInd,:);
plotFilled(linspace(-prelen/fR, postlen/fR, prelen + postlen + 1),M1,[1 0 0]);

M2 = pupilDiamatrix(mdIdx,:);
M2 = M2(nrwdInd,:);
plotFilled(linspace(-prelen/fR, postlen/fR, prelen + postlen + 1),M2,[0 0 1]);

xlim([-prelen/fR postlen/fR]);
title('Rwd vs nRwd');
legend(['rwdTrial'],[''],['nrwdTrial'],[''],'location','northwest','fontsize',12);
% subplot(2,2,2);hold on;
% rnormpD = pupilDiarwdmatrix(mdIdx,:);
% rpre = prer(mdIdx);
% rpre = rpre(rwdInd);
% rnormpD = rnormpD(rwdInd,:)./rpre;
% plotFilled(linspace(-prelenr/fR, postlenr/fR, prelenr + postlenr + 1),rnormpD,[0 0 1]);
% 
% nrnormpD = pupilDiarwdmatrix(mdIdx,:);
% rpre = prer(mdIdx);
% rpre = rpre(nrwdInd);
% nrnormpD = nrnormpD(nrwdInd,:)./rpre;
% plotFilled(linspace(-prelenr/fR, postlenr/fR, prelenr + postlenr + 1),nrnormpD,[1 0 0]);

% rwd vs non rwd

subplot(4,2,2);hold on;

preRwdInd = logical([0,rwdInd(1:end-1)]);
rpDpre = pupilDiamatrix(mdIdx,:);
rpDpre = rpDpre(preRwdInd,:);
plotFilled(linspace(-prelen/fR, postlen/fR, prelen + postlen + 1),rpDpre,[1 0 0]);

nrwdInd = ~preRwdInd;
nrpDpre = pupilDiamatrix(mdIdx,:);
nrpDpre = nrpDpre(nrwdInd,:);
plotFilled(linspace(-prelen/fR, postlen/fR, prelen + postlen + 1),nrpDpre,[0 0 1]);

xlim([-prelen/fR postlen/fR]);
title('preRwd vs prenRwd');
legend(['prerwdTrial'],[''],['prenrwdTrial'],[''],'location','northwest','fontsize',12);
% subplot(2,2,4);hold on;
% rpD = pupilDiarwdmatrix(mdIdx,:);
% rpD = rpD(preRwdInd,:);
% plotFilled(linspace(-prelenr/fR, postlenr/fR, prelenr + postlenr + 1),rpD,[0 0 1]);
% 
% nrpD = pupilDiarwdmatrix(mdIdx,:);
% nrpD = nrpD(nrwdInd,:);
% plotFilled(linspace(-prelenr/fR, postlenr/fR, prelenr + postlenr + 1),nrpD,[1 0 0]);

%stay trial
subplot(4,2,3);hold on;

r_nr_stay = find((allChoices(2:end).* allChoices(1:end-1)) > 0 & allRewards(2:end) == 0 & allRewards(1:end-1) ~= 0)+1;
nr_nr_stay = find((allChoices(2:end).* allChoices(1:end-1)) > 0 & allRewards(2:end) == 0 & allRewards(1:end-1) == 0)+1;

M1 = pupilDiamatrix(mdIdx,:);
M1 = M1(r_nr_stay,:);
plotFilled(linspace(-prelen/fR, postlen/fR, prelen + postlen + 1),M1,[1 0 0]);

M2 = pupilDiamatrix(mdIdx,:);
M2 = M2(nr_nr_stay,:);
plotFilled(linspace(-prelen/fR, postlen/fR, prelen + postlen + 1),M2,[0 0 1]);

xlim([-prelen/fR postlen/fR]);
title('rwd-nrwd vs nrwd-nrwd');
legend(['rwd-nrwd'],[''],['nrwd-nrwd'],[''],'location','northwest','fontsize',12);

subplot(4,2,4);hold on;

r_r_stay = find((allChoices(2:end).* allChoices(1:end-1)) > 0 & allRewards(2:end) ~= 0 & allRewards(1:end-1) ~= 0)+1;
nr_r_stay = find((allChoices(2:end).* allChoices(1:end-1)) > 0 & allRewards(2:end) ~= 0 & allRewards(1:end-1) == 0)+1;

M1 = pupilDiamatrix(mdIdx,:);
M1 = M1(r_r_stay,:);
plotFilled(linspace(-prelen/fR, postlen/fR, prelen + postlen + 1),M1,[1 0 0]);

M2 = pupilDiamatrix(mdIdx,:);
M2 = M2(nr_r_stay,:);
plotFilled(linspace(-prelen/fR, postlen/fR, prelen + postlen + 1),M2,[0 0 1]);

xlim([-prelen/fR postlen/fR]);
title('rwd-rwd vs nrwd-rwd');
legend(['rwd-rwd'],[''],['nrwd-rwd'],[''],'location','northwest','fontsize',12);

%switch trial
subplot(4,2,5);hold on;

r_nr_switch = find((allChoices(2:end).* allChoices(1:end-1)) < 0 & allRewards(2:end) == 0 & allRewards(1:end-1) ~= 0)+1;
nr_nr_switch = find((allChoices(2:end).* allChoices(1:end-1)) < 0 & allRewards(2:end) == 0 & allRewards(1:end-1) == 0)+1;

M1 = pupilDiamatrix(mdIdx,:);
M1 = M1(r_nr_switch,:);
plotFilled(linspace(-prelen/fR, postlen/fR, prelen + postlen + 1),M1,[1 0 0]);

M2 = pupilDiamatrix(mdIdx,:);
M2 = M2(nr_nr_switch,:);
plotFilled(linspace(-prelen/fR, postlen/fR, prelen + postlen + 1),M2,[0 0 1]);

xlim([-prelen/fR postlen/fR]);
title('rwd-switch-nrwd vs nrwd-switch-nrwd');
legend(['rwd-switch-nrwd'],[''],['nrwd-switch-nrwd'],[''],'location','northwest','fontsize',12);

subplot(4,2,6);hold on;

r_r_switch = find((allChoices(2:end).* allChoices(1:end-1)) < 0 & allRewards(2:end) ~= 0 & allRewards(1:end-1) ~= 0)+1;
nr_r_switch = find((allChoices(2:end).* allChoices(1:end-1)) < 0 & allRewards(2:end) ~= 0 & allRewards(1:end-1) == 0)+1;

M1 = pupilDiamatrix(mdIdx,:);
M1 = M1(r_r_switch,:);
plotFilled(linspace(-prelen/fR, postlen/fR, prelen + postlen + 1),M1,[1 0 0]);

M2 = pupilDiamatrix(mdIdx,:);
M2 = M2(nr_r_switch,:);
plotFilled(linspace(-prelen/fR, postlen/fR, prelen + postlen + 1),M2,[0 0 1]);

xlim([-prelen/fR postlen/fR]);
title('rwd-switch-rwd vs nrwd-switch-rwd');
legend(['rwd-switch-rwd'],[''],['nrwd-switch-rwd'],[''],'location','northwest','fontsize',12);

subplot(4,2,7);hold on;

r_switch = find((allChoices(2:end).* allChoices(1:end-1)) < 0 & allRewards(1:end-1)~=0)+1;
nr_switch = find((allChoices(2:end).* allChoices(1:end-1)) < 0 & allRewards(1:end-1) == 0)+1;
M1 = pupilDiamatrix(mdIdx,:);
M1 = M1(r_switch,:);
plotFilled(linspace(-prelen/fR, postlen/fR, prelen + postlen + 1),M1,[1 0 0]); hold on; 

M2 = pupilDiamatrix(mdIdx,:);
M2 = M2(nr_switch,:);
plotFilled(linspace(-prelen/fR, postlen/fR, prelen + postlen + 1),M2,[0 0 1]);

xlim([-prelen/fR postlen/fR]);
title('rwd-switch vs nrwd-switch');
legend(['rwd-switch'],[''],['nrwd-switch'],[''],'location','northwest','fontsize',12);

subplot(4,2,8);hold on;

switches = find((allChoices(2:end).* allChoices(1:end-1)) < 0)+1;
stay = find((allChoices(2:end).* allChoices(1:end-1)) > 0)+1;

M1 = pupilDiamatrix(mdIdx,:);
M1 = M1(switches,:);
plotFilled(linspace(-prelen/fR, postlen/fR, prelen + postlen + 1),M1,[1 0 0]);

M2 = pupilDiamatrix(mdIdx,:);
M2 = M2(stay,:);
plotFilled(linspace(-prelen/fR, postlen/fR, prelen + postlen + 1),M2,[0 0 1]);

xlim([-prelen/fR postlen/fR]);
title('switch vs stay');
legend(['switch'],[''],['stay'],[''],'location','northwest','fontsize',12);
%%
if isempty(dir(savePath))
    mkdir(savePath)
end
saveFigurePDF(gcf,[savePath session '_pupilswitchstay'])
%%
% rwd vs non rwd stay vs switch normalized
figure; hold on;
suptitle([session ' normalized' ' rwd vs nrwd stay vs switch']);
rwdInd = logical(allReward_L + allReward_R);
nrwdInd = ~rwdInd;

npupilDiamatrix = pupilDiamatrix./pre - 1;

subplot(4,2,1);hold on;
M1 = npupilDiamatrix(mdIdx,:);
M1 = M1(rwdInd,:);
plotFilled(linspace(-prelen/fR, postlen/fR, prelen + postlen + 1),M1,[1 0 0]);

M2 = npupilDiamatrix(mdIdx,:);
M2 = M2(nrwdInd,:);
plotFilled(linspace(-prelen/fR, postlen/fR, prelen + postlen + 1),M2,[0 0 1]);

xlim([-prelen/fR postlen/fR]);
title('Rwd vs nRwd');
legend(['rwdTrial'],[''],['nrwdTrial'],[''],'location','northwest','fontsize',12);
% subplot(2,2,2);hold on;
% rnormpD = pupilDiarwdmatrix(mdIdx,:);
% rpre = prer(mdIdx);
% rpre = rpre(rwdInd);
% rnormpD = rnormpD(rwdInd,:)./rpre;
% plotFilled(linspace(-prelenr/fR, postlenr/fR, prelenr + postlenr + 1),rnormpD,[0 0 1]);
% 
% nrnormpD = pupilDiarwdmatrix(mdIdx,:);
% rpre = prer(mdIdx);
% rpre = rpre(nrwdInd);
% nrnormpD = nrnormpD(nrwdInd,:)./rpre;
% plotFilled(linspace(-prelenr/fR, postlenr/fR, prelenr + postlenr + 1),nrnormpD,[1 0 0]);

% rwd vs non rwd

subplot(4,2,2);hold on;

preRwdInd = logical([0,rwdInd(1:end-1)]);
rpDpre = npupilDiamatrix(mdIdx,:);
rpDpre = rpDpre(preRwdInd,:);
plotFilled(linspace(-prelen/fR, postlen/fR, prelen + postlen + 1),rpDpre,[1 0 0]);

nrwdInd = ~preRwdInd;
nrpDpre = npupilDiamatrix(mdIdx,:);
nrpDpre = nrpDpre(nrwdInd,:);
plotFilled(linspace(-prelen/fR, postlen/fR, prelen + postlen + 1),nrpDpre,[0 0 1]);

xlim([-prelen/fR postlen/fR]);
title('preRwd vs prenRwd');
legend(['prerwdTrial'],[''],['prenrwdTrial'],[''],'location','northwest','fontsize',12);

%stay trial
subplot(4,2,3);hold on;

r_nr_stay = find((allChoices(2:end).* allChoices(1:end-1)) > 0 & allRewards(2:end) == 0 & allRewards(1:end-1) ~= 0)+1;
nr_nr_stay = find((allChoices(2:end).* allChoices(1:end-1)) > 0 & allRewards(2:end) == 0 & allRewards(1:end-1) == 0)+1;

M1 = npupilDiamatrix(mdIdx,:);
M1 = M1(r_nr_stay,:);
plotFilled(linspace(-prelen/fR, postlen/fR, prelen + postlen + 1),M1,[1 0 0]);

M2 = npupilDiamatrix(mdIdx,:);
M2 = M2(nr_nr_stay,:);
plotFilled(linspace(-prelen/fR, postlen/fR, prelen + postlen + 1),M2,[0 0 1]);

xlim([-prelen/fR postlen/fR]);
title('rwd-nrwd vs nrwd-nrwd');
legend(['rwd-nrwd'],[''],['nrwd-nrwd'],[''],'location','northwest','fontsize',12);

subplot(4,2,4);hold on;

r_r_stay = find((allChoices(2:end).* allChoices(1:end-1)) > 0 & allRewards(2:end) ~= 0 & allRewards(1:end-1) ~= 0)+1;
nr_r_stay = find((allChoices(2:end).* allChoices(1:end-1)) > 0 & allRewards(2:end) ~= 0 & allRewards(1:end-1) == 0)+1;

M1 = npupilDiamatrix(mdIdx,:);
M1 = M1(r_r_stay,:);
plotFilled(linspace(-prelen/fR, postlen/fR, prelen + postlen + 1),M1,[1 0 0]);

M2 = npupilDiamatrix(mdIdx,:);
M2 = M2(nr_r_stay,:);
plotFilled(linspace(-prelen/fR, postlen/fR, prelen + postlen + 1),M2,[0 0 1]);

xlim([-prelen/fR postlen/fR]);
title('rwd-rwd vs nrwd-rwd');
legend(['rwd-rwd'],[''],['nrwd-rwd'],[''],'location','northwest','fontsize',12);

%switch trial
subplot(4,2,5);hold on;

r_nr_switch = find((allChoices(2:end).* allChoices(1:end-1)) < 0 & allRewards(2:end) == 0 & allRewards(1:end-1) ~= 0)+1;
nr_nr_switch = find((allChoices(2:end).* allChoices(1:end-1)) < 0 & allRewards(2:end) == 0 & allRewards(1:end-1) == 0)+1;

M1 = npupilDiamatrix(mdIdx,:);
M1 = M1(r_nr_switch,:);
plotFilled(linspace(-prelen/fR, postlen/fR, prelen + postlen + 1),M1,[1 0 0]);

M2 = npupilDiamatrix(mdIdx,:);
M2 = M2(nr_nr_switch,:);
plotFilled(linspace(-prelen/fR, postlen/fR, prelen + postlen + 1),M2,[0 0 1]);

xlim([-prelen/fR postlen/fR]);
title('rwd-switch-nrwd vs nrwd-switch-nrwd');
legend(['rwd-switch-nrwd'],[''],['nrwd-switch-nrwd'],[''],'location','northwest','fontsize',12);

subplot(4,2,6);hold on;

r_r_switch = find((allChoices(2:end).* allChoices(1:end-1)) < 0 & allRewards(2:end) ~= 0 & allRewards(1:end-1) ~= 0)+1;
nr_r_switch = find((allChoices(2:end).* allChoices(1:end-1)) < 0 & allRewards(2:end) ~= 0 & allRewards(1:end-1) == 0)+1;

M1 = npupilDiamatrix(mdIdx,:);
M1 = M1(r_r_switch,:);
plotFilled(linspace(-prelen/fR, postlen/fR, prelen + postlen + 1),M1,[1 0 0]);

%%
M2 = npupilDiamatrix(mdIdx,:);
M2 = M2(nr_r_switch,:);
plotFilled(linspace(-prelen/fR, postlen/fR, prelen + postlen + 1),M2,[0 0 1]);

xlim([-prelen/fR postlen/fR]);
title('rwd-switch-rwd vs nrwd-switch-rwd');
legend(['rwd-switch-rwd'],[''],['nrwd-switch-rwd'],[''],'location','northwest','fontsize',12);

subplot(4,2,7);hold on;

r_switch = find((allChoices(2:end).* allChoices(1:end-1)) < 0 & allRewards(1:end-1)~=0)+1;
nr_switch = find((allChoices(2:end).* allChoices(1:end-1)) < 0 & allRewards(1:end-1) == 0)+1;
M1 = npupilDiamatrix(mdIdx,:);
M1 = M1(r_switch,:);
plotFilled(linspace(-prelen/fR, postlen/fR, prelen + postlen + 1),M1,[1 0 0]); hold on; 

M2 = npupilDiamatrix(mdIdx,:);
M2 = M2(nr_switch,:);
plotFilled(linspace(-prelen/fR, postlen/fR, prelen + postlen + 1),M2,[0 0 1]);

xlim([-prelen/fR postlen/fR]);
title('rwd-switch vs nrwd-switch');
legend(['rwd-switch'],[''],['nrwd-switch'],[''],'location','northwest','fontsize',12);

subplot(4,2,8);hold on;

switches = find((allChoices(2:end).* allChoices(1:end-1)) < 0)+1;
stay = find((allChoices(2:end).* allChoices(1:end-1)) > 0)+1;

M1 = npupilDiamatrix(mdIdx,:);
M1 = M1(switches,:);
plotFilled(linspace(-prelen/fR, postlen/fR, prelen + postlen + 1),M1,[1 0 0]);

M2 = npupilDiamatrix(mdIdx,:);
M2 = M2(stay,:);
plotFilled(linspace(-prelen/fR, postlen/fR, prelen + postlen + 1),M2,[0 0 1]);

xlim([-prelen/fR postlen/fR]);
title('switch vs stay');
legend(['switch'],[''],['stay'],[''],'location','northwest','fontsize',12);
%%
if isempty(dir(savePath))
    mkdir(savePath)
end
saveFigurePDF(gcf,[savePath session '_pupilswitchstaynorm'])
%%
%lick histogram
lickLat = [];       lickRate = [];
lickLat_L = [];     lickRate_L = [];
lickLat_R = [];     lickRate_R = [];
for i = 1:length(behSessionData)
    if ~isempty(behSessionData(i).rewardTime)
        lickLat = [lickLat behSessionData(i).rewardTime - behSessionData(i).CSon];
        if ~isnan(behSessionData(i).rewardL)
            lickLat_L = [lickLat_L behSessionData(i).rewardTime - behSessionData(i).CSon];
            if behSessionData(i).rewardL == 1
                if length(behSessionData(i).licksL) > 1
                    lickRateTemp = 1000/(min(diff(behSessionData(i).licksL)));
                    lickRate = [lickRate lickRateTemp];
                    lickRate_L = [lickRate_L lickRateTemp];
                else
                   lickRate = [lickRate 0];
                   lickRate_L = [lickRate_L 0]; 
                end
            end
        elseif ~isnan(behSessionData(i).rewardR)
            lickLat_R = [lickLat_R behSessionData(i).rewardTime - behSessionData(i).CSon];      %make single licks zeros for easier indexing
            if behSessionData(i).rewardR == 1
                if length(behSessionData(i).licksR) > 1
                    lickRateTemp = 1000/(min(diff(behSessionData(i).licksR)));
                    lickRate = [lickRate lickRateTemp];
                    lickRate_R = [lickRate_R lickRateTemp];
                else
                    lickRate = [lickRate 0];
                    lickRate_R = [lickRate_R 0];
                end
            end
        end
    end
end
%lickLat = lickLat(mdIdx);

%%
stay = find((allChoices(2:end).* allChoices(1:end-1)) > 0)+1;
switches = find((allChoices(2:end).* allChoices(1:end-1)) < 0 & allRewards(1:end-1) ~= 0)+1;
r_stay = find((allChoices(2:end).* allChoices(1:end-1)) > 0 & allRewards(1:end-1) ~= 0)+1;
nr_stay = find((allChoices(2:end).* allChoices(1:end-1)) > 0 & allRewards(1:end-1) == 0)+1;
r_nr_stay_stay = find((allChoices(3:end).* allChoices(2:end-1)) > 0 & (allChoices(3:end).* allChoices(1:end-2)) > 0 & allRewards(1:end-2) ~= 0 & allRewards(2:end-1) == 0)+2;
r_nr_stay_nr = find((allChoices(3:end).* allChoices(2:end-1)) > 0 & (allChoices(3:end).* allChoices(1:end-2)) > 0 & allRewards(1:end-2) ~= 0 & allRewards(2:end-1) == 0)+1;
nr_nr_stay_stay = find((allChoices(3:end).* allChoices(2:end-1)) > 0 & (allChoices(3:end).* allChoices(1:end-2)) > 0 & allRewards(1:end-2) == 0 & allRewards(2:end-1) == 0)+2;

figure; hold on;
subplot(3,2,1); hold on;
histogram(lickLat(stay),0:50:1000,'Normalization','probability', 'FaceColor', 'm');
histogram(lickLat(switches),0:50:1000,'Normalization','probability', 'FaceColor', 'c');
legend(['stay'],['switch'],'location','northwest');

subplot(3,2,2); hold on;
histogram(lickLat(r_stay),0:50:1000,'Normalization','probability', 'FaceColor', 'm');
histogram(lickLat(nr_stay),0:50:1000,'Normalization','probability', 'FaceColor', 'c');
legend(['r-stay'],['nr-stay'],'location','northwest');

subplot(3,2,3); hold on;
histogram(lickLat(r_nr_stay_stay),0:50:1000,'Normalization','probability', 'FaceColor', 'm');
histogram(lickLat(r_nr_stay_nr),0:50:1000,'Normalization','probability', 'FaceColor', 'c');
legend(['r-nr-stay-stay'],['r-nr-stay-nr'],'location','northwest');

subplot(3,2,4); hold on;
histogram(lickLat(r_nr_stay_stay),0:50:1000,'Normalization','probability', 'FaceColor', 'm');
histogram(lickLat(nr_nr_stay_stay),0:50:1000,'Normalization','probability', 'FaceColor', 'c');
legend(['r-nr-stay-stay'],['nr-nr-stay-stay'],'location','northwest');

subplot(3,2,5);hold on;
scatter(lickLat(r_nr_stay_stay),lickLat(r_nr_stay_nr))
xlim([0 max(lickLat)])
ylim([0 max(lickLat)])
line([0 max(lickLat)],[0 max(lickLat)])
xlabel('r-nr-stay-stay');
ylabel('r-nr-stay-nr');

%%
if isempty(dir(savePath))
    mkdir(savePath)
end
saveFigurePDF(gcf,[savePath session '_lick'])
%%
switch_r_stay = find((allChoices(2:end-1).* allChoices(1:end-2)) < 0 & allRewards(2:end-1)~= 0 & (allChoices(2:end-1).* allChoices(3:end)) > 0)+1;
switch_r = find((allChoices(2:end-1).* allChoices(1:end-2)) < 0 & allRewards(2:end-1)~= 0)+1;
p = size(switch_r_stay)/size(switch_r)