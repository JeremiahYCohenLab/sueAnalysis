%%
mZS034d20191116 = load('Z:\ZS034\mZS034d20191116\sorted\session\mZS034d20191116_sessionData_behav.mat');
mZS034d20191120 = load('Z:\ZS034\mZS034d20191120\sorted\session\mZS034d20191120_sessionData_behav.mat');
mZS022d20191118 = load('Z:\ZS022\mZS022d20191118\sorted\session\mZS022d20191118_sessionData_behav.mat');
mZS022d20191120 = load('Z:\ZS022\mZS022d20191120\sorted\session\mZS022d20191120_sessionData_behav.mat');

mZS034d20191119 = load('Z:\ZS034\mZS034d20191119\sorted\session\mZS034d20191119_sessionData_behav.mat');
mZS034d20191121 = load('Z:\ZS034\mZS034d20191121\sorted\session\mZS034d20191121_sessionData_behav.mat');
mZS022d20191119 = load('Z:\ZS022\mZS022d20191119\sorted\session\mZS022d20191119_sessionData_behav.mat');
mZS022d20191121 = load('Z:\ZS022\mZS022d20191121\sorted\session\mZS022d20191121_sessionData_behav.mat');
%%
fields = fieldnames(mZS034d20191116.behSessionData)';
fields{2,8} = NaN;
gap = struct(fields{:});
gap(2:10) =gap(1);
behSessionData = [mZS022d20191119.behSessionData,gap,mZS022d20191121.behSessionData, gap,mZS034d20191119.behSessionData,gap,mZS034d20191121.behSessionData];    
%behSessionData = [mZS022d20191118.behSessionData,gap,mZS022d20191120.behSessionData, gap,mZS034d20191116.behSessionData,gap,mZS034d20191120.behSessionData];
session = 'ZS022Ctrl';
%% List load control

workbookFile = 'Z:\combineAnimals';
ani = 'combine';
lst = 'ctrlp';
[~, dayList, ~] = xlsread(workbookFile, ani);
col = contains(dayList(1,:),lst);
dayList = dayList(2:end,col);
endInd = find(cellfun(@isempty,dayList),1);
if ~isempty(endInd)
    dayList = dayList(1:endInd-1,:);
end

[root,sep] = currComputer_df;

for i = 1:length(dayList)
    session = dayList{i};
    [animalName, ~] = strtok(session, 'd'); 
    animalName = animalName(2:end);
    tempBeh = load([root animalName sep session sep 'sorted' sep 'session' sep session '_sessionData_behav_pupil.mat']);
    fields = fieldnames(tempBeh.behSessionData);
    if ~isempty(find(contains(fields,'hmm') == 1, 1))
        break
    end
end
    fields = fields';
    fields{2,8} = NaN;
    gap = struct(fields{:});
    gap(2:10) =gap(1);

    session = dayList{1};
    [animalName, ~] = strtok(session, 'd'); 
    animalName = animalName(2:end);
    tempBeh = load([root animalName sep session sep 'sorted' sep 'session' sep session '_sessionData_behav_pupil.mat']);
    fields = fieldnames(tempBeh.behSessionData);
    if isempty(find(contains(fields,'hmm') == 1, 1))
        [~, states] = fitHmmOpt(dayList{i},1);
        tempBeh = load([root animalName sep session sep 'sorted' sep 'session' sep session '_sessionData_behav_pupil.mat']);
    end
    behaSessionDataCnb = tempBeh.behSessionData;
    responseInds = find(~isnan([tempBeh.behSessionData.rewardTime]));
    allITIs = [tempBeh.behSessionData(responseInds(1:end-1) + 1).CSon] - [tempBeh.behSessionData(responseInds(1:end-1)).CSon];
    allITIs = [allITIs 20000];
    
for i = 2:length(dayList)
    session = dayList{i};
    [animalName, ~] = strtok(session, 'd'); 
    animalName = animalName(2:end);
    tempBeh = load([root animalName sep session sep 'sorted' sep 'session' sep session '_sessionData_behav_pupil.mat']);
    fields = fieldnames(tempBeh.behSessionData);
    if isempty(find(contains(fields,'hmm') == 1, 1))
        [~, states] = fitHmmOpt(dayList{i},1);
        tempBeh = load([root animalName sep session sep 'sorted' sep 'session' sep session '_sessionData_behav_pupil.mat']);
    end
    behaSessionDataCnb = [behaSessionDataCnb, gap, tempBeh.behSessionData];
    responseInds = find(~isnan([tempBeh.behSessionData.rewardTime]));
    allITItemp = [tempBeh.behSessionData(responseInds(1:end-1) + 1).CSon] - [tempBeh.behSessionData(responseInds(1:end-1)).CSon];
    allITItemp = [allITItemp  20000];
    allITIs = [allITIs allITItemp];
end

%%
% matrix
fR = 21;
len = length(behaSessionDataCnb);
postlen = round(fR*6);
prelen = round(fR*1);
tracelen = postlen + prelen;
pupilDiamatrix = NaN(len,tracelen);
pupilXmatrix = NaN(len,tracelen);
pre = NaN(len,1);
prer = NaN(len,1);
for i = 1: len
    if ~(isnan(behaSessionDataCnb(i).pDiapre))
     pupilDiamatrix(i,:) = behaSessionDataCnb(i).pDia(1:tracelen);
     pupilXmatrix(i,:) = behaSessionDataCnb(i).pX(1:tracelen);
     pre(i) = behaSessionDataCnb(i).pDiapre;
    end
end
npupilDiamatrix = pupilDiamatrix./pre - 1;
cpupilDiamatrix = pupilDiamatrix - pre ;


%%
%rwd nrwd
mdIdx = ~isnan([behaSessionDataCnb.rewardTime]);
mdIdx = find(mdIdx > 0);

responseInds = find(~isnan([behaSessionDataCnb.rewardTime])); % find CS+ trials with a response in the lick window
allReward_R = [behaSessionDataCnb(responseInds).rewardR]; 
allReward_L = [behaSessionDataCnb(responseInds).rewardL]; 
allReward_R = [behaSessionDataCnb(responseInds).rewardR]; 
allReward_L = [behaSessionDataCnb(responseInds).rewardL];  
allChoices = NaN(1,length(behaSessionDataCnb(responseInds)));
allChoices(~isnan(allReward_R)) = 1;
allChoices(~isnan(allReward_L)) = -1;

allReward_R(isnan(allReward_R)) = 0;
allReward_L(isnan(allReward_L)) = 0;
allChoice_R = double(allChoices == 1);
allChoice_L = double(allChoices == -1);

allRewards = zeros(1,length(allChoices));
allRewards(logical(allReward_R)) = 1;
allRewards(logical(allReward_L)) = -1;
%%
% rwd vs non rwd stay vs switch raw trace
figure; hold on;
suptitle([session ' raw trace' ' rwd vs nrwd stay vs switch']);
M = pupilDiamatrix(mdIdx,:);
shortITIIdx = allITIs < 1000*(postlen/fR) ;
M(shortITIIdx, :) = NaN;

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

M1 = M(rwdInd,:);
plotFilled(linspace(-prelen/fR, postlen/fR, tracelen),M1,[1 0 0]);

M2 = M(nrwdInd,:);
plotFilled(linspace(-prelen/fR, postlen/fR, tracelen),M2,[0 0 1]);

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
rpDpre = M(preRwdInd,:);
plotFilled(linspace(-prelen/fR, postlen/fR, tracelen),rpDpre,[1 0 0]);

nrwdInd = ~preRwdInd;
nrpDpre = M(nrwdInd,:);
plotFilled(linspace(-prelen/fR, postlen/fR,tracelen),nrpDpre,[0 0 1]);

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

M1 = M(r_nr_stay,:);
plotFilled(linspace(-prelen/fR, postlen/fR, tracelen),M1,[1 0 0]);

M2 = M(nr_nr_stay,:);
plotFilled(linspace(-prelen/fR, postlen/fR, tracelen),M2,[0 0 1]);

xlim([-prelen/fR postlen/fR]);
title('rwd-nrwd vs nrwd-nrwd');
legend(['rwd-nrwd'],[''],['nrwd-nrwd'],[''],'location','northwest','fontsize',12);

subplot(4,2,4);hold on;

r_r_stay = find((allChoices(2:end).* allChoices(1:end-1)) > 0 & allRewards(2:end) ~= 0 & allRewards(1:end-1) ~= 0)+1;
nr_r_stay = find((allChoices(2:end).* allChoices(1:end-1)) > 0 & allRewards(2:end) ~= 0 & allRewards(1:end-1) == 0)+1;

M1 = M(r_r_stay,:);
plotFilled(linspace(-prelen/fR, postlen/fR, tracelen),M1,[1 0 0]);

M2 = M(nr_r_stay,:);
plotFilled(linspace(-prelen/fR, postlen/fR, tracelen),M2,[0 0 1]);

xlim([-prelen/fR postlen/fR]);
title('rwd-rwd vs nrwd-rwd');
legend(['rwd-rwd'],[''],['nrwd-rwd'],[''],'location','northwest','fontsize',12);

%switch trial
subplot(4,2,5);hold on;

r_nr_switch = find((allChoices(2:end).* allChoices(1:end-1)) < 0 & allRewards(2:end) == 0 & allRewards(1:end-1) ~= 0)+1;
nr_nr_switch = find((allChoices(2:end).* allChoices(1:end-1)) < 0 & allRewards(2:end) == 0 & allRewards(1:end-1) == 0)+1;

M1 = M(r_nr_switch,:);
plotFilled(linspace(-prelen/fR, postlen/fR, tracelen),M1,[1 0 0]);

M2 = M(nr_nr_switch,:);
plotFilled(linspace(-prelen/fR, postlen/fR, tracelen),M2,[0 0 1]);

xlim([-prelen/fR postlen/fR]);
title('rwd-switch-nrwd vs nrwd-switch-nrwd');
legend(['rwd-switch-nrwd'],[''],['nrwd-switch-nrwd'],[''],'location','northwest','fontsize',12);

subplot(4,2,6);hold on;

r_r_switch = find((allChoices(2:end).* allChoices(1:end-1)) < 0 & allRewards(2:end) ~= 0 & allRewards(1:end-1) ~= 0)+1;
nr_r_switch = find((allChoices(2:end).* allChoices(1:end-1)) < 0 & allRewards(2:end) ~= 0 & allRewards(1:end-1) == 0)+1;

M1 = M(r_r_switch,:);
plotFilled(linspace(-prelen/fR, postlen/fR, tracelen),M1,[1 0 0]);

M2 = M(nr_r_switch,:);
plotFilled(linspace(-prelen/fR, postlen/fR, tracelen),M2,[0 0 1]);

xlim([-prelen/fR postlen/fR]);
title('rwd-switch-rwd vs nrwd-switch-rwd');
legend(['rwd-switch-rwd'],[''],['nrwd-switch-rwd'],[''],'location','northwest','fontsize',12);

subplot(4,2,7);hold on;

r_switch = find((allChoices(2:end).* allChoices(1:end-1)) < 0 & allRewards(1:end-1)~=0)+1;
nr_switch = find((allChoices(2:end).* allChoices(1:end-1)) < 0 & allRewards(1:end-1) == 0)+1;

M1 = M(r_switch,:);
plotFilled(linspace(-prelen/fR, postlen/fR, tracelen),M1,[1 0 0]); hold on; 

M2 = M(nr_switch,:);
plotFilled(linspace(-prelen/fR, postlen/fR, tracelen),M2,[0 0 1]);

xlim([-prelen/fR postlen/fR]);
title('rwd-switch vs nrwd-switch');
legend(['rwd-switch'],[''],['nrwd-switch'],[''],'location','northwest','fontsize',12);

subplot(4,2,8);hold on;

switches = find((allChoices(2:end).* allChoices(1:end-1)) < 0)+1;
stay = find((allChoices(2:end).* allChoices(1:end-1)) > 0)+1;

M1 = M(switches,:);
plotFilled(linspace(-prelen/fR, postlen/fR, tracelen),M1,[1 0 0]);

M2 = M(stay,:);
plotFilled(linspace(-prelen/fR, postlen/fR, tracelen),M2,[0 0 1]);

xlim([-prelen/fR postlen/fR]);
title('switch vs stay');
legend(['switch'],[''],['stay'],[''],'location','northwest','fontsize',12);
%%
if isempty(dir(savePath))
    mkdir(savePath)
end
saveFigurePDF(gcf,[savePath session '_pupilswitchstay'])
%%
figure; hold on; 
suptitle(session);
M = pupilDiamatrix(mdIdx,:);
shortITIIdx = allITIs < 1000*(postlen/fR) ;
M(shortITIIdx, :) = NaN;
subplot(1,2,1);hold on;
switch_r = find((allChoices(2:end).* allChoices(1:end-1)) < 0 & allRewards(2:end) ~= 0)+1;
switch_nr = find((allChoices(2:end).* allChoices(1:end-1)) < 0 & allRewards(2:end) == 0)+1;

M1 = M(r_nr_switch,:);
plotFilled(linspace(-prelen/fR, postlen/fR, tracelen),M1,[1 0 0]);

M2 = M(nr_nr_switch,:);
plotFilled(linspace(-prelen/fR, postlen/fR, tracelen),M2,[0 0 1]);

xlim([-prelen/fR postlen/fR]);
title('switch-rwd vs switch-nrwd raw');
legend(['switch-rwd'],[''],['switch-nrwd'],[''],'location','northwest','fontsize',12);

M = npupilDiamatrix(mdIdx,:);
shortITIIdx = allITIs < 1000*(postlen/fR) ;
M(shortITIIdx, :) = NaN;
subplot(1,2,2);hold on;
switch_r = find((allChoices(2:end).* allChoices(1:end-1)) < 0 & allRewards(2:end) ~= 0)+1;
switch_nr = find((allChoices(2:end).* allChoices(1:end-1)) < 0 & allRewards(2:end) == 0)+1;

M1 = M(r_nr_switch,:);
plotFilled(linspace(-prelen/fR, postlen/fR, tracelen),M1,[1 0 0]);

M2 = M(nr_nr_switch,:);
plotFilled(linspace(-prelen/fR, postlen/fR, tracelen),M2,[0 0 1]);

xlim([-prelen/fR postlen/fR]);
title('switch-rwd vs switch-nrwd norm');
legend(['switch-rwd'],[''],['switch-nrwd'],[''],'location','northwest','fontsize',12);

%%
% rwd vs non rwd stay vs switch change
figure; hold on;
suptitle([session ' norm' ' rwd vs nrwd stay vs switch']);
rwdInd = logical(allReward_L + allReward_R);
nrwdInd = ~rwdInd;

M = npupilDiamatrix(mdIdx,:);
shortITIIdx = allITIs < 1000*(postlen/fR) ;
M(shortITIIdx, :) = NaN;

subplot(4,2,1);hold on;
M1 = M(rwdInd,:);
plotFilled(linspace(-prelen/fR, postlen/fR, tracelen),M1,[1 0 0]);

M2 = M(nrwdInd,:);
plotFilled(linspace(-prelen/fR, postlen/fR, tracelen),M2,[0 0 1]);

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
rpDpre = M(preRwdInd,:);
plotFilled(linspace(-prelen/fR, postlen/fR, tracelen),rpDpre,[1 0 0]);

nrwdInd = ~preRwdInd;
nrpDpre = M(nrwdInd,:);
plotFilled(linspace(-prelen/fR, postlen/fR, tracelen),nrpDpre,[0 0 1]);

xlim([-prelen/fR postlen/fR]);
title('preRwd vs prenRwd');
legend(['prerwdTrial'],[''],['prenrwdTrial'],[''],'location','northwest','fontsize',12);

%stay trial
subplot(4,2,3);hold on;

r_nr_stay = find((allChoices(2:end).* allChoices(1:end-1)) > 0 & allRewards(2:end) == 0 & allRewards(1:end-1) ~= 0)+1;
nr_nr_stay = find((allChoices(2:end).* allChoices(1:end-1)) > 0 & allRewards(2:end) == 0 & allRewards(1:end-1) == 0)+1;

M1 = M(r_nr_stay,:);
plotFilled(linspace(-prelen/fR, postlen/fR, tracelen),M1,[1 0 0]);

M2 = M(nr_nr_stay,:);
plotFilled(linspace(-prelen/fR, postlen/fR, tracelen),M2,[0 0 1]);

xlim([-prelen/fR postlen/fR]);
title('rwd-nrwd vs nrwd-nrwd');
legend(['rwd-nrwd'],[''],['nrwd-nrwd'],[''],'location','northwest','fontsize',12);

subplot(4,2,4);hold on;

r_r_stay = find((allChoices(2:end).* allChoices(1:end-1)) > 0 & allRewards(2:end) ~= 0 & allRewards(1:end-1) ~= 0)+1;
nr_r_stay = find((allChoices(2:end).* allChoices(1:end-1)) > 0 & allRewards(2:end) ~= 0 & allRewards(1:end-1) == 0)+1;

M1 = M(r_r_stay,:);
plotFilled(linspace(-prelen/fR, postlen/fR, tracelen),M1,[1 0 0]);

M2 = M(nr_r_stay,:);
plotFilled(linspace(-prelen/fR, postlen/fR, tracelen),M2,[0 0 1]);

xlim([-prelen/fR postlen/fR]);
title('rwd-rwd vs nrwd-rwd');
legend(['rwd-rwd'],[''],['nrwd-rwd'],[''],'location','northwest','fontsize',12);

%switch trial
subplot(4,2,5);hold on;

r_nr_switch = find((allChoices(2:end).* allChoices(1:end-1)) < 0 & allRewards(2:end) == 0 & allRewards(1:end-1) ~= 0)+1;
nr_nr_switch = find((allChoices(2:end).* allChoices(1:end-1)) < 0 & allRewards(2:end) == 0 & allRewards(1:end-1) == 0)+1;

M1 = M(r_nr_switch,:);
plotFilled(linspace(-prelen/fR, postlen/fR, tracelen),M1,[1 0 0]);

M2 = M(nr_nr_switch,:);
plotFilled(linspace(-prelen/fR, postlen/fR, tracelen),M2,[0 0 1]);

xlim([-prelen/fR postlen/fR]);
title('rwd-switch-nrwd vs nrwd-switch-nrwd');
legend(['rwd-switch-nrwd'],[''],['nrwd-switch-nrwd'],[''],'location','northwest','fontsize',12);

subplot(4,2,6);hold on;

r_r_switch = find((allChoices(2:end).* allChoices(1:end-1)) < 0 & allRewards(2:end) ~= 0 & allRewards(1:end-1) ~= 0)+1;
nr_r_switch = find((allChoices(2:end).* allChoices(1:end-1)) < 0 & allRewards(2:end) ~= 0 & allRewards(1:end-1) == 0)+1;

M1 = M(r_r_switch,:);
plotFilled(linspace(-prelen/fR, postlen/fR, tracelen),M1,[1 0 0]);

M2 = M(nr_r_switch,:);
plotFilled(linspace(-prelen/fR, postlen/fR, tracelen),M2,[0 0 1]);

xlim([-prelen/fR postlen/fR]);
title('rwd-switch-rwd vs nrwd-switch-rwd');
legend(['rwd-switch-rwd'],[''],['nrwd-switch-rwd'],[''],'location','northwest','fontsize',12);

subplot(4,2,7);hold on;

r_switch = find((allChoices(2:end).* allChoices(1:end-1)) < 0 & allRewards(1:end-1)~=0)+1;
nr_switch = find((allChoices(2:end).* allChoices(1:end-1)) < 0 & allRewards(1:end-1) == 0)+1;

M1 = M(r_switch,:);
plotFilled(linspace(-prelen/fR, postlen/fR, tracelen),M1,[1 0 0]); hold on; 

M2 = M(nr_switch,:);
plotFilled(linspace(-prelen/fR, postlen/fR, tracelen),M2,[0 0 1]);

xlim([-prelen/fR postlen/fR]);
title('rwd-switch vs nrwd-switch');
legend(['rwd-switch'],[''],['nrwd-switch'],[''],'location','northwest','fontsize',12);

subplot(4,2,8);hold on;

switches = find((allChoices(2:end).* allChoices(1:end-1)) < 0)+1;
stay = find((allChoices(2:end).* allChoices(1:end-1)) > 0)+1;

M1 = M(switches,:);
plotFilled(linspace(-prelen/fR, postlen/fR, tracelen),M1,[1 0 0]);

M2 = M(stay,:);
plotFilled(linspace(-prelen/fR, postlen/fR, tracelen),M2,[0 0 1]);

xlim([-prelen/fR postlen/fR]);
title('switch vs stay');
legend(['switch'],[''],['stay'],[''],'location','northwest','fontsize',12);
%%
figure; hold on; 
suptitle(session);
subplot(1,2,1);hold on;
switch_r = find((allChoices(2:end).* allChoices(1:end-1)) < 0 & allRewards(2:end) ~= 0)+1;
switch_nr = find((allChoices(2:end).* allChoices(1:end-1)) < 0 & allRewards(2:end) == 0)+1;
M = pupilDiamatrix(mdIdx,:);
shortITIIdx = allITIs < 1000*(postlen/fR) ;
M(shortITIIdx, :) = NaN;

M1 = M(switch_r,:);
plotFilled(linspace(-prelen/fR, postlen/fR, tracelen),M1,[1 0 0]);

M2 = M(switch_nr,:);
plotFilled(linspace(-prelen/fR, postlen/fR, tracelen),M2,[0 0 1]);

xlim([-prelen/fR postlen/fR]);
title('switch-rwd vs switch-nrwd raw');
legend(['switch-rwd'],[''],['switch-nrwd'],[''],'location','northwest','fontsize',12);

subplot(1,2,2);hold on;
M = cpupilDiamatrix(mdIdx,:);
shortITIIdx = allITIs < 1000*(postlen/fR) ;
M(shortITIIdx, :) = NaN; 

M1 = M(switch_r,:);
plotFilled(linspace(-prelen/fR, postlen/fR, tracelen),M1,[1 0 0]);

M2 = M(switch_nr,:);
plotFilled(linspace(-prelen/fR, postlen/fR, tracelen),M2,[0 0 1]);

xlim([-prelen/fR postlen/fR]);
title('switch-rwd vs switch-nrwd norm');
legend(['switch-rwd'],[''],['switch-nrwd'],[''],'location','northwest','fontsize',12);

%%
% eye position
Xpre = mean(pupilXmatrix(:,10:17)');
Xpre = Xpre(responseInds);

tMax = 10;
rwdMatx = [];
for i = 1:tMax
    rwdMatx(i,:) = [NaN(1,i) allRewards(1:end-i)];
end

allNoRewards = allChoices;
allNoRewards(allRewards~=0) = 0;
noRwdMatx = [];
for i = 1:tMax
    noRwdMatx(i,:) = [NaN(1,i) allNoRewards(1:end-i)];
end

glm_rwdX = fitlm([rwdMatx(:,~isnan(Xpre));noRwdMatx(:,~isnan(Xpre))]', -Xpre(~isnan(Xpre)));                 

%plot beta coefficients from lrm's
 figure; hold on;
 suptitle(session);
 subplot(2,2,1); hold on;
 relevInds = 2:tMax+1;
 coefVals = glm_rwdX.Coefficients.Estimate(relevInds);
 CIbands = coefCI(glm_rwdX);
 errorL = abs(coefVals - CIbands(relevInds,1));
 errorU = abs(coefVals - CIbands(relevInds,2));
 errorbar(1:10,coefVals,errorL,errorU,'Color', [0.7 0 1],'linewidth',2)
 
 relevInds = tMax+2:2*tMax+1;
 coefVals = glm_rwdX.Coefficients.Estimate(relevInds);
 CIbands = coefCI(glm_rwdX);
 errorL = abs(coefVals - CIbands(relevInds,1));
 errorU = abs(coefVals - CIbands(relevInds,2));
 errorbar(1:10,coefVals,errorL,errorU,'Color',[30,144,255]/255,'linewidth',2)

 xlabel('Right reward n trials back')
 ylabel('\beta Coefficient')
 title('LRM: rewards in trial on eye position')
 line([0 10], [0 0], 'Color','k','LineStyle','--')
 legend('Rwd','nRwd');
 
xl = Xpre(logical(allChoice_L));
xr = Xpre(logical(allChoice_R)); 

[h,p,ci,stats] = ttest2(xl, xr);
[hpre,ppre,cipre,statspre] = ttest2(Xpre(logical([0,allReward_L(1:end-1)])),Xpre(logical([0,allReward_R(1:end-1)])));

seml = std(xl(~isnan(xl)))/sqrt(length(xl(~isnan(xl)))); 
semr = std(xr(~isnan(xr)))/sqrt(length(xr(~isnan(xr)))); 

subplot(2,2,4); hold on; 
bar(1:2,[nanmean(xl) nanmean(xr)]); hold on; 
er = errorbar(1:2,[nanmean(xl) nanmean(xr)],[seml semr]);                                      
er.LineStyle = 'none';  
hold off
xlabel('L lick   R lick')
ylabel('eye position')
ylim([135 140])
title('next pretrial eye position')
xl = Xpre(logical([0,allReward_L(1:end-1)]));
xr = Xpre(logical([0,allReward_R(1:end-1)]));
seml = std(xl(~isnan(xl)))/sqrt(length(xl(~isnan(xl)))); 
semr = std(xr(~isnan(xr)))/sqrt(length(xr(~isnan(xr)))); 
legend(sprintf('p = %s',ppre),'location','northeast');

subplot(2,2,2); hold on; 
bar(1:2,[nanmean(xl) nanmean(xr)]); hold on; 
er = errorbar(1:2,[nanmean(xl) nanmean(xr)],[seml semr]);                                      
er.LineStyle = 'none';  
hold off
xlabel('L Reward   R Reward')
ylabel('eye position')
ylim([135 140])
title('pre trial eye position')
legend(sprintf('p = %s',p),'location','northeast');

subplot(2,2,3)
tMax = 10;
xMatx = [];
for i = 1:tMax
    xMatx(i,:) = [NaN(1,i-1) -Xpre(1:end-i+1)];
end

glm_rwdX = fitlm(xMatx', allChoices);                 

%plot beta coefficients from lrm's
 subplot(2,2,3); hold on;
 relevInds = 2:tMax+1;
 coefVals = glm_rwdX.Coefficients.Estimate(relevInds);
 CIbands = coefCI(glm_rwdX);
 errorL = abs(coefVals - CIbands(relevInds,1));
 errorU = abs(coefVals - CIbands(relevInds,2));
 errorbar(1:10,coefVals,errorL,errorU,'Color', [0.7 0 1],'linewidth',2)
 xlabel('pretrial eye position n trials back')
 ylabel('\beta Coefficient')
 title('LRM: eye position on choices')
 line([0 10], [0 0], 'Color','k','LineStyle','--')
%%
%CSplus vs CS minus
CSplusIn = strcmp({behaSessionDataCnb.trialType},'CSplus');
CSminusIn = strcmp({behaSessionDataCnb.trialType},'CSminus');
figure; hold on; 
suptitle([session 'CSplus vs CS minus'])

subplot(1,2,1);hold on;
M1 = pupilDiamatrix(CSplusIn,:);
plotFilled(linspace(-prelen/fR, postlen/fR, tracelen),M1,[1 0 0]);

M2 = pupilDiamatrix(CSminusIn,:);
plotFilled(linspace(-prelen/fR, postlen/fR, tracelen),M2,[0 0 1]);

xlim([-prelen/fR postlen/fR]);
title('CS+ vs CSminus raw');
legend(['CS+'],[''],['CS-'],[''],'location','northwest','fontsize',12);

npupilDiamatrix = pupilDiamatrix./pre - 1;
cpupilDiamatrix = pupilDiamatrix - pre ;
subplot(1,2,2);hold on;

M1 = cpupilDiamatrix(CSplusIn,:);
plotFilled(linspace(-prelen/fR, postlen/fR, tracelen),M1,[1 0 0]);

M2 = cpupilDiamatrix(CSminusIn,:);
plotFilled(linspace(-prelen/fR, postlen/fR, tracelen),M2,[0 0 1]);

xlim([-prelen/fR postlen/fR]);
title('CS+ vs CSminus change');
legend(['CS+'],[''],['CS-'],[''],'location','northwest','fontsize',12);

%%
%lick latency
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
lickLat = lickLat(~isnan(lickLat));
%%
%lick lat
figure;hold on
suptitle(session);
subplot(1,2,1);
scatter(lickLat, pre(mdIdx),10,'filled');
xlabel('lickLat');
ylabel('baselineDia')
o = pre(mdIdx)';
o = o(~isnan(o));
[R,P,RL,RU] = corrcoef([lickLat(~isnan(pre(mdIdx))); o]');
legend(sprintf('r = %s \n p =%s',R(1,2),P(1,2)),'location','northeast')

subplot(1,2,2);
m = max(npupilDiamatrix(:,22:end)');
scatter(lickLat, m(mdIdx),10,'filled');
xlabel('lickLat');
ylabel('dilation')
o = m(mdIdx)';
o = o(~isnan(o));
[R,P,RL,RU] = corrcoef([lickLat(~isnan(pre(mdIdx))); o']');
legend([sprintf('r = %s \n p =%s',R(1,2),P(1,2))],'location','northeast')


%%
%model compare
mod = 'Q10_2';
figure;
suptitle(mod);
for i = 1:length(modStruct.(mod).params)
    subplot(3,2,i);hold on;
    [p,tbl,stats] = kruskalwallis([x1(:,i),x2(:,i)],[],'off');
    boxplot([x1(:,i),x2(:,i)]);
    [h,pt,ci,stats] = ttest2(x1(:,i),x2(:,i));
    xlabel('Ctrl (8)                                             Inhi(8)')
    title([modStruct.(mod).params{1,i} sprintf('\n K-W: p = %d T: p = %d',round(p,3),round(pt,3))])
end

figure;
subplot(1,2,1);hold on;
boxplot([x1(:,6),x2(:,6)]);
[h,pt,ci,stats] = ttest2(x1(:,6),x2(:,6));
xlabel('Ctrl (8)                                             Inhi(8)')
title(['pstay/rwd' sprintf('\nT: p = %d',round(pt,3))])

subplot(1,2,2);hold on;
boxplot([x1(:,7),x2(:,7)]);
[h,pt,ci,stats] = ttest2(x1(:,7),x2(:,7));
xlabel('Ctrl (8)                                             Inhi(8)')
title(['pswitch/nrwd' sprintf('\nT: p = %d',round(pt,3))])

%%
%model fitting
[root,sep] = currComputer_df;
sessionFolder = 'session';
mod = 'Q10_2';
workbookFile = 'Z:\combineAnimals'; 
[A, B, C] = xlsread(workbookFile, mod);
model = struct;
%modStruct = importLearningModels_df('Q');
%%
for i = [7,9,15,16]
    model(i).session = [B{i,1}];
    session = model(i).session;
    [animalName, date] = strtok(model(i).session, 'd'); 
    animalName = animalName(2:end);
% load behavioral data
    behSessionDataPath = [root animalName sep session sep 'sorted' sep 'session' sep session '_sessionData_behav.mat'];
    load(behSessionDataPath);
    fullStruct.s = behSessionData;
    os = parseBehavioralData_df(fullStruct, blockSwitch, blockSwitch);
% fit behavioral model
    ms = fitLearningMods_df(os, 'modStruct', modStruct, 'particularModel', {mod});
    mdIdx = ~isnan([behSessionData.rewardTime]);
    mdIdx = find(mdIdx > 0);
    for j = 1:length(mdIdx)
        behSessionData(mdIdx(j)).rpe = ms.(mod).modOutput.rpe(j);
        behSessionData(mdIdx(j)).qL = ms.(mod).modOutput.Q(j,1);
        behSessionData(mdIdx(j)).qR = ms.(mod).modOutput.Q(j,2);
        behSessionData(mdIdx(j)).pchoice = ms.(mod).modOutput.probChosenChoice(j);
    end
    save(behSessionDataPath, 'behSessionData', 'blockSwitch', 'blockSwitchL', 'blockSwitchR')
end
%%
%dilation vs rpe
%rpe
rpe = [behSessionData(mdIdx).rpe];
figure;hold on
suptitle(session);

subplot(3,2,1);
m = max(npupilDiamatrix(:,22:end)');
m = m(mdIdx);
rpe = rpe(~isnan(m));
o = m(~isnan(m));
absrpe = abs(rpe);
scatter(abs(rpe), o,10,'filled');
xlabel('rpe');
ylabel('change');
[R,P,RL,RU] = corrcoef([absrpe; o]');
legend(sprintf('r = %s \n p =%s',R(1,2),P(1,2)),'location','northeast')

tMax = 5;
arpeMatx = [];
for i = 1:tMax
    arpeMatx(i,:) = [NaN(1,i-1) absrpe(1:end-i+1)];
end

glm_rpe = fitlm(arpeMatx',o);
relevInds = 2:tMax+1;
subplot(3,2,2);
coefVals = glm_rpe.Coefficients.Estimate(relevInds);
CIbands = coefCI(glm_rpe);
errorL = abs(coefVals - CIbands(relevInds,1));
errorU = abs(coefVals - CIbands(relevInds,2));
errorbar(1:tMax,coefVals,errorL,errorU,'Color', [0.7 0 1],'linewidth',2)
xlim([0,tMax+1])
xlabel('absrpe n-1 trials back')
ylabel('\beta Coefficient')
title('LRM: absrpe on dilation')
line([0 10], [0 0], 'Color','k','LineStyle','--')

subplot(3,2,3);
m = max(npupilDiamatrix(:,22:end)');
m = m(mdIdx);
o = m(~isnan(m));
prpe = rpe(rpe>0);
scatter(prpe, o((rpe > 0)),10,'filled');
xlabel('prpe');
ylabel('change')

[R,P,RL,RU] = corrcoef([prpe; o(rpe>0)]');
legend(sprintf('r = %s \n p =%s',R(1,2),P(1,2)),'location','northeast')

tMax = 5;
arpeMatx = [];
for i = 1:tMax
    arpeMatx(i,:) = [NaN(1,i-1) prpe(1:end-i+1)];
end

glm_rpe = fitlm(arpeMatx',o(rpe>0));
relevInds = 2:tMax+1;
subplot(3,2,4);
coefVals = glm_rpe.Coefficients.Estimate(relevInds);
CIbands = coefCI(glm_rpe);
errorL = abs(coefVals - CIbands(relevInds,1));
errorU = abs(coefVals - CIbands(relevInds,2));
errorbar(1:tMax,coefVals,errorL,errorU,'Color', [0.7 0 1],'linewidth',2)
xlim([0,tMax+1])
xlabel('absrpe n-1 trials back')
ylabel('\beta Coefficient')
title('LRM: prpe on dilation')
line([0 10], [0 0], 'Color','k','LineStyle','--')

subplot(3,2,5);
m = max(npupilDiamatrix(:,22:end)');
m = m(mdIdx);
o = m(~isnan(m));
nrpe = rpe(rpe<0);
scatter(nrpe, o((rpe < 0)),10,'filled');
xlabel('nrpe');
ylabel('dilation')

[R,P,RL,RU] = corrcoef([nrpe; o(rpe<0)]');
legend(sprintf('r = %s \n p =%s',R(1,2),P(1,2)),'location','northeast')

tMax = 5;
arpeMatx = [];
for i = 1:tMax
    arpeMatx(i,:) = [NaN(1,i-1) nrpe(1:end-i+1)];
end

glm_rpe = fitlm(arpeMatx',o(rpe<0));
relevInds = 2:tMax+1;
subplot(3,2,6);
coefVals = glm_rpe.Coefficients.Estimate(relevInds);
CIbands = coefCI(glm_rpe);
errorL = abs(coefVals - CIbands(relevInds,1));
errorU = abs(coefVals - CIbands(relevInds,2));
errorbar(1:tMax,coefVals,errorL,errorU,'Color', [0.7 0 1],'linewidth',2)
xlim([0,tMax+1])
xlabel('absrpe n-1 trials back')
ylabel('\beta Coefficient')
title('LRM: nrpe on dilation')
line([0 10], [0 0], 'Color','k','LineStyle','--')

%%
%baseline vs exp
mdIdx = ~isnan([behSessionData.rewardTime]);
mdIdx = find(mdIdx > 0);
eve = [behSessionData.pchoice]/0.5 - 1;
m = max(pupilDiamatrix(:,22:end)');

figure;hold on;
suptitle(session);
subplot(1,2,1);
m = max(cpupilDiamatrix(:,22:end)');
scatter(eve, m(mdIdx),10,'filled');
xlabel('explore -- exploit');
ylabel('change')
o = m(mdIdx)';
o = o(~isnan(o));
[R,P,RL,RU] = corrcoef([eve(~isnan(m(mdIdx)))', o]);
glm_eve = fitlm(eve(~isnan(m(mdIdx)))',o);
refline([glm_eve.Coefficients.Estimate(2),glm_eve.Coefficients.Estimate(1)])
legend(sprintf('r = %s p =%s \n', R(1,2),P(1,2)), sprintf('b Coefficient = %s \n p =%s', glm_eve.Coefficients.Estimate(2),glm_eve.Coefficients.pValue(2)),'location','northeast');

subplot(1,2,2);
scatter(eve, pre(mdIdx),10,'filled');
xlabel('explore -- exploit');
ylabel('baseline')
o = pre(mdIdx)';
o = o(~isnan(o));
[R,P,RL,RU] = corrcoef([eve(~isnan(m(mdIdx)))', o']);
glm_eve = fitlm(eve(~isnan(m(mdIdx)))',o);
refline([glm_eve.Coefficients.Estimate(2),glm_eve.Coefficients.Estimate(1)])
legend(sprintf('r = %s p =%s \n', R(1,2),P(1,2)), sprintf('b Coefficient = %s \n p =%s', glm_eve.Coefficients.Estimate(2),glm_eve.Coefficients.pValue(2)),'location','northeast')
%%
%rwd-nrwd - nrwd-nrwd ctrl vs inhi
figure;
suptitle('ctrl')
subplot(1,2,1);hold on;
%plot(linspace(-prelen/fR, postlen/fR, tracelen),A,'LineWidth',4)
plot(linspace(-prelen/fR, postlen/fR, tracelen),D,'LineWidth',4)
% legend('Inhi','Ctrl')
line([-1 3], [0 0], 'Color','k','LineStyle','--')
title('rwd-nrwd - nrwd-nrwd')
ylim([-0.025 0.05])

subplot(1,2,2);hold on;
%plot(linspace(-prelen/fR, postlen/fR, tracelen),B,'LineWidth',4)
plot(linspace(-prelen/fR, postlen/fR, tracelen),E,'LineWidth',4)
legend('Inhi','Ctrl')
line([-1 3], [0 0], 'Color','k','LineStyle','--')
% title('nrwd-rwd - rwd-rwd')
ylim([-0.025 0.05])
%%
%raw trace
figure;hold on; 
subplot(1,2,1);hold on;
plotFilled(linspace(-prelen/fR, postlen/fR, tracelen),F,[1 0 0]);
plotFilled(linspace(-prelen/fR, postlen/fR, tracelen),C,[0 0 1]);
legend('ctrl',[''],'inhi',[''])
title('raw trace')

subplot(1,2,2);hold on;
plotFilled(linspace(-prelen/fR, postlen/fR, tracelen),Fn,[1 0 0]);
plotFilled(linspace(-prelen/fR, postlen/fR, tracelen),Cn,[0 0 1]);
legend('ctrl',[''],'inhi',[''])
title('Dilation')
%% stay-nrwd vs switch-nrwd baseline opens plasticity
% stay switch dilation
figure;hold on;
M = cpupilDiamatrix(mdIdx,:);
shortITIIdx = allITIs < 1000*(postlen/fR) ;
M(shortITIIdx, :) = NaN;

subplot(2,2,1); hold on;
stay_r = find((allChoices(2:end).* allChoices(1:end-1)) > 0 & allRewards(2:end) ~= 0)+1;
switch_r = find((allChoices(2:end).* allChoices(1:end-1)) <0 & allRewards(2:end) ~= 0)+1;

M1 = M(stay_r,:);
plotFilled(linspace(-prelen/fR, postlen/fR, tracelen),M1,[1 0 0]);

M2 = M(switch_r,:);
plotFilled(linspace(-prelen/fR, postlen/fR, tracelen),M2,[0 0 1]);

xlim([-prelen/fR postlen/fR]);
title('stay-rwd vs switch-rwd');
legend(['stay-rwd'],[''],['switch-rwd'],[''],'location','northwest','fontsize',12);

subplot(2,2,2); hold on;
stay_nr = find((allChoices(2:end).* allChoices(1:end-1)) > 0 & allRewards(2:end) == 0)+1;
switch_nr = find((allChoices(2:end).* allChoices(1:end-1)) <0 & allRewards(2:end) == 0)+1;

M1 = M(stay_nr,:);
plotFilled(linspace(-prelen/fR, postlen/fR, tracelen),M1,[1 0 0]);

M2 = M(switch_nr,:);
plotFilled(linspace(-prelen/fR, postlen/fR, tracelen),M2,[0 0 1]);

xlim([-prelen/fR postlen/fR]);
title('stay-nrwd vs switch-nrwd');
legend(['stay-nrwd'],[''],['switch-nrwd'],[''],'location','northwest','fontsize',12);

M = pupilDiamatrix(mdIdx,:);
shortITIIdx = allITIs < 1000*(postlen/fR) ;
M(shortITIIdx, :) = NaN;

subplot(2,2,3); hold on;
stay_r = find((allChoices(2:end).* allChoices(1:end-1)) > 0 & allRewards(2:end) ~= 0)+1;
switch_r = find((allChoices(2:end).* allChoices(1:end-1)) <0 & allRewards(2:end) ~= 0)+1;

M1 = M(stay_r,:);
plotFilled(linspace(-prelen/fR, postlen/fR, tracelen),M1,[1 0 0]);

M2 = M(switch_r,:);
plotFilled(linspace(-prelen/fR, postlen/fR, tracelen),M2,[0 0 1]);

xlim([-prelen/fR postlen/fR]);
title('stay-rwd vs switch-rwd');
legend(['stay-rwd'],[''],['switch-rwd'],[''],'location','northwest','fontsize',12);

subplot(2,2,4); hold on;
stay_nr = find((allChoices(2:end).* allChoices(1:end-1)) > 0 & allRewards(2:end) == 0)+1;
switch_nr = find((allChoices(2:end).* allChoices(1:end-1)) <0 & allRewards(2:end) == 0)+1;

M1 = M(stay_nr,:);
plotFilled(linspace(-prelen/fR, postlen/fR, tracelen),M1,[1 0 0]);

M2 = M(switch_nr,:);
plotFilled(linspace(-prelen/fR, postlen/fR, tracelen),M2,[0 0 1]);

xlim([-prelen/fR postlen/fR]);
title('stay-nrwd vs switch-nrwd');
legend(['stay-nrwd'],[''],['switch-nrwd'],[''],'location','northwest','fontsize',12);
%% nr-stay-r vs nr-switch-r
% stay switch dilation
figure;hold on;
M = npupilDiamatrix(mdIdx,:);
shortITIIdx = allITIs < 1000*(postlen/fR) ;
M(shortITIIdx, :) = NaN;

subplot(2,2,1); hold on;
nr_stay_r = find((allChoices(2:end).* allChoices(1:end-1)) > 0 & allRewards(2:end) ~= 0 & allRewards(1:end-1) == 0)+1;
nr_switch_r = find((allChoices(2:end).* allChoices(1:end-1)) <0 & allRewards(2:end) ~= 0 & allRewards(1:end-1) == 0)+1;

M1 = M(nr_stay_r,:);
plotFilled(linspace(-prelen/fR, postlen/fR, tracelen),M1,[1 0 0]);

M2 = M(nr_switch_r,:);
plotFilled(linspace(-prelen/fR, postlen/fR, tracelen),M2,[0 0 1]);
                 
xlim([-prelen/fR postlen/fR]);
title('nr-stay-r vs nr-switch-r');
legend(['nr-stay-r'],[''],['nr-switch-r'],[''],'location','northwest','fontsize',12);

subplot(2,2,2); hold on;
switch_r = find((allChoices(2:end).* allChoices(1:end-1)) < 0 & allRewards(2:end) ~= 0)+1;
switch_nr = find((allChoices(2:end).* allChoices(1:end-1)) < 0 & allRewards(2:end) == 0)+1;

M1 = M(switch_r,:);
plotFilled(linspace(-prelen/fR, postlen/fR, tracelen),M1,[1 0 0]);

M2 = M(switch_nr,:);
plotFilled(linspace(-prelen/fR, postlen/fR, tracelen),M2,[0 0 1]);
                 
xlim([-prelen/fR postlen/fR]);
title('switch-r vs switch-nr');
legend(['switch-r'],[''],['switch-nr'],[''],'location','northwest','fontsize',12);

% stay switch raw

M = pupilDiamatrix(mdIdx,:);
shortITIIdx = allITIs < 1000*(postlen/fR) ;
M(shortITIIdx, :) = NaN;

subplot(2,2,3); hold on;

M1 = M(nr_stay_r,:);
plotFilled(linspace(-prelen/fR, postlen/fR, tracelen),M1,[1 0 0]);

M2 = M(nr_switch_r,:);
plotFilled(linspace(-prelen/fR, postlen/fR, tracelen),M2,[0 0 1]);
                 
xlim([-prelen/fR postlen/fR]);
title('nr-stay-r vs nr-switch-r');
legend(['nr-stay-r'],[''],['nr-switch-r'],[''],'location','northwest','fontsize',12);

subplot(2,2,4); hold on;

M1 = M(switch_r,:);
plotFilled(linspace(-prelen/fR, postlen/fR, tracelen),M1,[1 0 0]);

M2 = M(switch_nr,:);
plotFilled(linspace(-prelen/fR, postlen/fR, tracelen),M2,[0 0 1]);
                 
xlim([-prelen/fR postlen/fR]);
title('switch-r vs switch-nr');
legend(['switch-r'],[''],['switch-nr'],[''],'location','northwest','fontsize',12);


%%
%p switch - stay
[root,sep] = currComputer_df;
sessionFolder = 'session';
mod = 'Q10_2';
workbookFile = 'Z:\combineAnimals'; 
[A, B, C] = xlsread(workbookFile, mod);
%%
p = struct;
for i = 1:16
    session = [B{i+1,1}];
    p(i).session = session; 
    [animalName, date] = strtok(session, 'd'); 
    animalName = animalName(2:end);
 % load behavioral data
    behSessionDataPath = [root animalName sep session sep 'sorted' sep 'session' sep session '_sessionData_behav.mat'];
    load(behSessionDataPath);
% fit behavioral model
responseInds = find(~isnan([behSessionData.rewardTime])); % find CS+ trials with a response in the lick window
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
    switch_r_stay = find((allChoices(2:end-1).* allChoices(1:end-2)) < 0 & allRewards(2:end-1)~= 0 & (allChoices(2:end-1).* allChoices(3:end)) > 0)+1;
    switch_r = find((allChoices(2:end-1).* allChoices(1:end-2)) < 0 & allRewards(2:end-1)~= 0)+1;
    ptemp = size(switch_r_stay)/size(switch_r);
    p(i).pstayswitchrwd = ptemp; 
    nr_switch = find((allChoices(2:end).* allChoices(1:end-1)) < 0 & allRewards(1:end-1)== 0);
    nr= find(allRewards(1:end-1)== 0);
    ptemp = size(nr_switch)/size(nr);
    p(i).pswitchnrwd = ptemp; 
    
end

compare = [p.pstayswitchrwd];
x1 = compare(1:8);x2 = compare(9:16);
[h,P,ci,stats] = ttest2(x1,x2);

figure;                                    
boxplot([x1;x2]');
xlabel('Contrl                                Inhibition')
title(['P(stay/switched&rewarded)' sprintf('p = %d',P)])
%%
%switch_r_stay switch_r_switch
M = npupilDiamatrix(mdIdx,:);
shortITIIdx = allITIs < 1000*(postlen/fR) ;
M(shortITIIdx, :) = NaN;

figure; hold on;
subplot(1,2,1);hold on;
switch_r_stay = find((allChoices(2:end-1).* allChoices(1:end-2)) < 0 & allRewards(2:end-1)~= 0 & (allChoices(2:end-1).* allChoices(3:end)) > 0)+1;
switch_r_switch = find((allChoices(2:end-1).* allChoices(1:end-2)) < 0 & allRewards(2:end-1)~= 0 & (allChoices(3:end).* allChoices(2:end-1)) < 0)+1;

M1 = M(switch_r_stay,:);
plotFilled(linspace(-prelen/fR, postlen/fR, tracelen),M1,[1 0 0]);


M2 = M(switch_r_switch,:);
plotFilled(linspace(-prelen/fR, postlen/fR, tracelen),M2,[0 0 1]);

xlim([-prelen/fR postlen/fR]);
title('Ctrl norm switch-r-stay vs switch-r-switch');
legend(['switch-r-stay'],[''],['switch-r-switch'],[''],'location','northwest','fontsize',12);

subplot(1,2,2);hold on;
M = pupilDiamatrix(mdIdx,:);
shortITIIdx = allITIs < 1000*(postlen/fR) ;
M(shortITIIdx, :) = NaN;

M1 = M(switch_r_stay,:);
plotFilled(linspace(-prelen/fR, postlen/fR, tracelen),M1,[1 0 0]);


M2 = M(switch_r_switch,:);
plotFilled(linspace(-prelen/fR, postlen/fR, tracelen),M2,[0 0 1]);

xlim([-prelen/fR postlen/fR]);
title('Ctrl norm switch-r-stay vs switch-r-switch');
legend(['switch-r-stay'],[''],['switch-r-switch'],[''],'location','northwest','fontsize',12);

%%
%regression on baseline and maxDia
%baseline
premd = pre(mdIdx);

%diff for stay trials
diffstay = abs(diff(abs(allRewards)));
diffstay(find((allChoices(2:end).* allChoices(1:end-1)) < 0)) = 0;
diffstay = [NaN,diffstay];
%switch or not
svs = 0.5*[NaN,abs(diff(allChoices))];
%lick lat
lickLat = NaN(1,len);
for i = 1:len
    if ~isnan(behSessionData(i).rewardTime)
        lickLat(i) = behSessionData(i).rewardTime - behSessionData(i).CSon;
    end
end
lickLat = lickLat(mdIdx);
lickLat = zscore(lickLat);
% maxD
m = max(pupilDiamatrix(:,21:end)');
m = m(mdIdx);
%dilation
d = max(npupilDiamatrix(:,21:end)');
d = d(mdIdx);
%change
c = m - premd'; 
%rpe
rpe = [behSessionData(mdIdx).rpe];
absrpe = abs(rpe);
%rpe(find((allChoices(2:end).* allChoices(1:end-1)) < 0)+1) = NaN;
regMat = [[NaN,premd(1:end-1)']; [NaN,diffstay(1:end-1)];  svs;  lickLat];
regLick = [[NaN, lickLat(1:end-1)];[NaN, abs(allRewards(1:end-1))];premd';svs];
%%
regMat = regMat(:,~isnan(premd));
regLick = regLick(:,~isnan(premd));
m = m(~isnan(premd));
d = d(~isnan(premd));
b = premd(~isnan(premd));
c = c(~isnan(premd));
lickLat = lickLat(~isnan(premd));
glm_all = fitlm(regMat', b');
glm_lick = fitlm(regLick', lickLat');
%%
coefVals = glm_lick.Coefficients.Estimate(2:5);
CIbands = coefCI(glm_lick);
errorL = abs(coefVals - CIbands(2:5,1));
errorU = abs(coefVals - CIbands(2:5,2));
figure;
bar(1:4,coefVals); hold on; 
er = errorbar(1:4,coefVals,errorL,errorU,'Color', [0 0 0],'linewidth',1.5);
er.LineStyle = 'none'; 
title('lrm: on licklat')
xlabel('pre Licklat  prereward baseline switch vs stay')
ylabel('\beta Coefficient')
text(0.5,0.35,sprintf('R^2 = %d',glm_lick.Rsquared.Adjusted))
hold off
%% prepDia


