%% change pupil trace length

%% list
workbookFile = 'Z:\combineAnimals';
ani = 'ZS036';
lst = 'Inhi';
[~, dayList, ~] = xlsread(workbookFile, ani);
col = contains(dayList(1,:),lst);
dayList = dayList(2:end,col);
endInd = find(cellfun(@isempty,dayList),1);
if ~isempty(endInd)
    dayList = dayList(1:endInd-1,:);
end

%% training set
trainingSession = 'mZS036d20191207';
trainingdate = 'Dec7';
itr = '90000';

[root,sep] = currComputer_df;
   
for i = 1:length(dayList)
    %% behavior
    session = dayList{i};
    [animalName, ~] = strtok(session, 'd'); 
    animalName = animalName(2:end);
    load([root animalName sep session sep 'sorted' sep 'session' sep session '_sessionData_behav.mat']);
    load([root animalName sep session sep 'pupil' sep session '.mat'], 'fits')
    fR = fits.cleaned.timing.(session).frameRate*fits.cleaned.timing.(session).scaling_factor;
    fits.cleaned.timing.(session).rwdTiming_inFrames = fR*([behSessionData.rewardTime]-behSessionData(1).CSon)/1000+fits.cleaned.timing.(session).trialTiming_inFrames(1);
    fields = fieldnames(behSessionData);
    if isempty(find(contains(fields,'hmm') == 1, 1))
        [~, states] = fitHmmOpt(dayList{i});
        load([root animalName sep session sep 'sorted' sep 'session' sep session '_sessionData_behav.mat']);
    end
    
    %% diameter
    videopath = ['C:\Users\zhixiao\Documents\pupil\' animalName sep session '\pupil'];
    skeleton = [videopath sep session 'DLC_resnet50_' trainingSession trainingdate 'shuffle1_' itr '_skeleton.csv'];
    position = [videopath sep session 'DLC_resnet50_' trainingSession trainingdate 'shuffle1_' itr '.csv'];    
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
    title(session);
    
    %% append
    len = length(behSessionData);
    postlen = round(fR*6);
    prelen = round(fR*1);
    postlenr = round(fR*5.5);
    prelenr = round(fR*1.5);
    tt = round(fits.cleaned.timing.(session).trialTiming_inFrames);
    tr = round(fits.cleaned.timing.(session).rwdTiming_inFrames);
    for j = 1: len
        if ~(isnan(tt(j)))
            ftemp = tt(j);
            behSessionData(j).pDia = diameter(ftemp-prelen:ftemp+postlen)';
            behSessionData(j).pDiapre = mean(diameter(ftemp-round(fR/2):ftemp));
            behSessionData(j).pX = x(ftemp-prelen:ftemp+postlen)';
            behSessionData(j).pY = y(ftemp-prelen:ftemp+postlen)';
        end
        
        if ~isnan(behSessionData(j).rewardTime) && ~isnan(tt(j))
            ftemp = round(tt(j) + fR*((behSessionData(j).rewardTime - behSessionData(j).CSon)/1000));
            behSessionData(j).pDiarwd = diameter(ftemp-prelenr:ftemp+postlenr)';
            behSessionData(j).pDiaprer = mean(diameter(ftemp-round(fR/2):ftemp));
            behSessionData(j).pXrwd = x(ftemp-prelenr:ftemp+postlenr)';
            behSessionData(j).pYrwd = y(ftemp-prelenr:ftemp+postlenr)';
        end
    end
    savepath = [root animalName sep session sep 'sorted' sep 'session' sep];
    save([savepath session '_sessionData_behav_pupil.mat'], 'behSessionData','blockSwitch','blockSwitchL','blockSwitchR');
end