%% List load control

workbookFile = 'Z:\combineAnimals';
ani = 'combine';
lst = 'Ctrl';
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
    tempBeh = load([root animalName sep session sep 'sorted' sep 'session' sep session '_sessionData_behav.mat']);
    fields = fieldnames(tempBeh.behSessionData);
    if ~isempty(find(contains(fields,'hmm') == 1, 1))
        break
    end
end
    fields = fields';
    fields{2,8} = NaN;
    gap = struct(fields{:});
    gap(2:10) =gap(1);
%% start wtih first day
    session = dayList{1};
    [animalName, ~] = strtok(session, 'd'); 
    animalName = animalName(2:end);
    tempBeh = load([root animalName sep session sep 'sorted' sep 'session' sep session '_sessionData_behav.mat']);
    fields = fieldnames(tempBeh.behSessionData);
    if isempty(find(contains(fields,'hmm') == 1, 1))
        [~, states] = fitHmmOpt(dayList{i},1);
        tempBeh = load([root animalName sep session sep 'sorted' sep 'session' sep session '_sessionData_behav.mat']);
    end
    behaSessionDataCnb = tempBeh.behSessionData;
    responseInds = find(~isnan([tempBeh.behSessionData.rewardTime]));
    allITIs = [tempBeh.behSessionData(responseInds(1:end-1) + 1).CSon] - [tempBeh.behSessionData(responseInds(1:end-1)).CSon];
    allITIs = [allITIs 20000];
    evetemp = [tempBeh.behSessionData(responseInds).hmm];
    evetemp(evetemp ~= 1) = 0;
    hmm = [tempBeh.behSessionData(responseInds).hmm];
    swtemp = hmm(2:end) == 1 & hmm(1:end-1) ~=1;
    swtemp = [0 swtemp];
    eve = evetemp;
    sw = swtemp;
    %%
    
for i = 2:length(dayList)
    session = dayList{i};
    [animalName, ~] = strtok(session, 'd'); 
    animalName = animalName(2:end);
    tempBeh = load([root animalName sep session sep 'sorted' sep 'session' sep session '_sessionData_behav.mat']);
    fields = fieldnames(tempBeh.behSessionData);
    if isempty(find(contains(fields,'hmm') == 1, 1))
        [~, states] = fitHmmOpt(dayList{i},1);
        tempBeh = load([root animalName sep session sep 'sorted' sep 'session' sep session '_sessionData_behav.mat']);
    end
    behaSessionDataCnb = [behaSessionDataCnb, gap, tempBeh.behSessionData];
    responseInds = find(~isnan([tempBeh.behSessionData.rewardTime]));
    allITItemp = [tempBeh.behSessionData(responseInds(1:end-1) + 1).CSon] - [tempBeh.behSessionData(responseInds(1:end-1)).CSon];
    allITItemp = [allITItemp  20000];
    allITIs = [allITIs allITItemp];
    evetemp = [tempBeh.behSessionData(responseInds).hmm];
    evetemp(evetemp ~= 1) = 0;
    eve = [eve evetemp];
    hmm = [tempBeh.behSessionData(responseInds).hmm];
    swtemp = hmm(2:end) == 1 & hmm(1:end-1) ~=1;
    swtemp = [0 swtemp];
    sw = [sw swtemp];
end
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

%% logistic regression to switch/not
%pre rwd
rwd = abs(allRewards)-0.5;
prerwd = [NaN,rwd(1:end-1)];
preprerwd = [NaN,prerwd(1:end-1)];
prepreprerwd = [NaN,preprerwd(1:end-1)];
preprepreprerwd = [NaN,prepreprerwd(1:end-1)];
inter = prerwd .* [NaN, allITIs(1:end-1)]/1000; 
%switch
svs = 0.5 * [NaN,abs(diff(allChoices))];
regMat = [prerwd; [NaN, allITIs(1:end-1)]/10000]';
%regMat = [preprerwd; prepreprerwd; [NaN, allITIs(1:end-1)]/1000]';
%%
nprerwdIdx = find(prerwd < 0);
regMat = regMat(nprerwdIdx,:);
svs = svs(nprerwdIdx);
%%
glm_all = fitglm(regMat, sw, 'distribution','binomial','link','logit');
%%
figure;hold on;

coefVals = glm_all_inhi.Coefficients.Estimate(2:end);
CIbands = coefCI(glm_all_inhi);
errorL = abs(coefVals - CIbands(2:end,1));
errorU = abs(coefVals - CIbands(2:end,2));

bar(5,coefVals(end)); hold on;
errorbar(1:4,coefVals(1:end-1),errorL(1:end-1),errorU(1:end-1),'Color', [0 0 0],'linewidth',1.5); 
errorbar(5,coefVals(end),errorL(end),errorU(end),'Color', [0 0 0],'linewidth',1.5);

title('lrm: on licklat')
xlabel('pre Licklat  prereward baseline switch vs stay')
ylabel('\beta Coefficient')
text(0.5,0.35,sprintf('R^2 = %d',glm_lick.Rsquared.Adjusted))
hold off
%% 
