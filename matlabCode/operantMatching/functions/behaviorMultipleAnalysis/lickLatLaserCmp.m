function lickLatLaserCmp(xlFile, sheet, col, varargin)
p = inputParser;
p.addParameter('modelName', '5params');
p.addParameter('target', 'Qchosen')
p.addParameter('numBins', 2);
p.addParameter('shuffle', 0);
p.parse(varargin{:});

dayList = getDayList(xlFile, sheet, col);
numBins = p.Results.numBins;
QmeansC = zeros(length(dayList), numBins);
latMeansC = zeros(length(dayList), numBins);
QmeansL = zeros(length(dayList), numBins);
latMeansL = zeros(length(dayList), numBins);
[root, sep] = currComputer(); 
svsCombined = [];
targetCombined = [];
laserCombined = [];
outcomeCombined = [];
%%
for i = 1:length(dayList)
    session = dayList{i};
    pd = parseSessionString_df(session, root, sep);
    sampFile = [pd.animalName col '_' p.Results.modelName];
    path = [root pd.animalName sep pd.animalName 'sorted' sep 'stan' sep 'bernoulli' sep p.Results.modelName sep col sep];
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
    [t,~,badSessionFlag] = getStanModelParams_samps(p.Results.modelName, [path sampFile '.mat'], 2000, 'sessionName', session);
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
    paramNames = getParamNames_dF(p.Results.modelName,1);
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
    % shuffle
    if p.Results.shuffle
        laser = rand(1,length(laser));
        laser(laser>=0.7) = 1;
        laser(laser<0.7) = 0;
    end
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
    laserNrwdInd = intersect(os.nrwd_Inds, find(os.laser==1));
    ctrlNrwdInd = intersect(os.nrwd_Inds, find(os.laser==0));
    laserRwdInd = intersect(os.rwd_Inds, find(os.laser==1));
    ctrlRwdInd = intersect(os.rwd_Inds, find(os.laser==0));
    
    latCtempNrwd = mean(os.lickLat(mintersect(os.nrwd_Inds+2, ctrlNrwdInd+1, os.changeChoice_Inds, os.stayChoice_Inds+1)));
    latLtempNrwd = mean(os.lickLat(mintersect(os.nrwd_Inds+2, laserNrwdInd+1, os.changeChoice_Inds, os.stayChoice_Inds+1)));
    latCtempRwd = mean(os.lickLat(mintersect(os.rwd_Inds+2, ctrlNrwdInd+1, os.changeChoice_Inds, os.stayChoice_Inds+1)));
    latLtempRwd = mean(os.lickLat(mintersect(os.rwd_Inds+2, laserNrwdInd+1, os.changeChoice_Inds, os.stayChoice_Inds+1)));

    latCNrwd(i) = latCtempNrwd;
    latLNrwd(i) = latLtempNrwd;
    latCRwd(i) = latCtempRwd;
    latLRwd(i) = latLtempRwd;    
    
    eval(['target = ' p.Results.target ';']);
    focusLaserTrials = mintersect(os.nrwd_Inds, os.stayChoice_Inds, os.changeChoice_Inds-1, sessIndRange(1:end-1));
    temp = target(focusLaserTrials);
    edgesLat = linspace(min(temp)-0.01, max(temp)+0.01, numBins+1);
%     
%     figure2;
%     subplot(2,1,1)
%     hold on;
%     edgesTemp = linspace(min(target)-0.01, max(target)+0.01, 10);
%     histogram(target(mintersect(os.nrwd_Inds+1, os.nrwd_Inds, os.changeChoice_Inds-1, os.stayChoice_Inds)),  edgesTemp);
%     histogram(target(mintersect(os.rwd_Inds+1, os.nrwd_Inds, os.changeChoice_Inds-1, os.stayChoice_Inds)),  edgesTemp);
%     subplot(2,1,2)
%     hold on;
%     histogram(target(focusLaserTrials),  edgesTemp);
%     plot([edgesLat(2), edgesLat(2)], [0 10]);
%     sgtitle(session);
%     edges = binEqualSize(target, numBins+1);
    % lasered trials
    tempQLat = target(mintersect(find(laser==1), focusLaserTrials));
    tempLat = [os.lickLatZ];
    tempLat = tempLat(mintersect(find(laser==1)+1, focusLaserTrials+1));
%     edges = linspace(min(tempQ)-0.01, max(tempQ)+0.01, numBins+1);
    for j = 1:numBins
        QLatmeansL(i,j) = mean(tempQLat(tempQLat>=edgesLat(j) & tempQLat<edgesLat(j+1)));
        latMeansL(i,j) = mean(tempLat(tempQLat>=edgesLat(j) & tempQLat<edgesLat(j+1)));
    end
    % control trials
    tempQLat = target(mintersect(find(laser==0), focusLaserTrials));
    tempLat = [os.lickLatZ];
    tempLat = tempLat(mintersect(find(laser==0)+1, focusLaserTrials+1));
%     edges = linspace(min(tempQ)-0.01, max(tempQ)+0.01, numBins+1);
    for j = 1:numBins
        QLatmeansC(i,j) = mean(tempQLat(tempQLat>=edgesLat(j) & tempQLat<edgesLat(j+1)));
        latMeansC(i,j) = mean(tempLat(tempQLat>=edgesLat(j) & tempQLat<edgesLat(j+1)));
    end    
    
    targetCombined = [targetCombined; Qchosen];
    outcomeCombined = [outcomeCombined; abs(os.allRewards)'];
    laserCombined = [laserCombined; laser'];
end

%% plot pSwitch against Qchosen, split by laser or not
semLatL = sem(latMeansL);
semLatC = sem(latMeansC);
meanLatL = mean(latMeansL,'omitnan');
meanLatC = mean(latMeansC,'omitnan');
meanQLatL = mean(QLatmeansL, 'omitnan');
meanQLatC = mean(QLatmeansC, 'omitnan');

figure2; hold on;
sig = NaN(1,numBins);
for i = 1:numBins
    sig(i) = signrank(latMeansC(:,i), latMeansL(:,i));
end

plot(meanQLatC, meanLatC, 'LineWidth',2, 'Color', [0.6 0.6 0.6]);
fill([meanQLatC flip(meanQLatC)], [meanLatC-semLatC flip(meanLatC+semLatC)], [0.6 0.6 0.6], 'FaceAlpha', 0.25, 'EdgeColor', 'none');

plot(meanQLatL, meanLatL, 'LineWidth',2, 'Color', [0 0.5 1]);
fill([meanQLatL flip(meanQLatL)], [meanLatL-semLatL flip(meanLatL+semLatL)], [0 0.5 1], 'FaceAlpha', 0.25, 'EdgeColor', 'none');
xlabel(p.Results.target);

scatter(meanQLatL(sig<0.05), meanLatL(sig<0.05)+semLatL(sig<0.05)+0.02, 15, 'r', 'filled');

sigCell = mat2cell(sig, 1, ones(1,numBins));
text(meanQLatL', meanLatL'+semLatL'+0.05, cellfun(@(x)num2str(x,3), sigCell, 'UniformOutput', false), 'FontSize', 10);

sgtitle('lickLat')

end