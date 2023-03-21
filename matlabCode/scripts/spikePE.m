load('/Volumes/Elements/computerDataBackUp/tmpData/catWithOutcomePe.mat')
modelName = '5params';
tb = 2;
tf = 4;
focusWin = [300 1300];
[root, sep] = currComputer();
prevSession = [];
category = 'good';
time = -1000*tb:1000*tf;
maxTrial = 1000;
paramNames = getParamNames_dF(modelName, 1);
%%
color1 = [0 0.8 0.8];
color2 = [1 0.2 1];
[coeff,scores,latent, ~, explained, mu] = pca(zscore(waveformsSession, [], 1));
indAll = {};
dis = {};
for a = 1:10
    [indAll{a}, ~, dis{a}] = kmeans([scores(:, 1:5)], 2);
end
[~,optiInds] = min(cellfun(@mean, dis));
cats = indAll{optiInds};
if sum(cats==1)<sum(cats==2)
    cats = 3-cats;
end

figure2;
subplot(1,2,1)
scatter(scores(cats==1, 1), scores(cats==1, 2), 12, color1, 'filled')
hold on;
scatter(scores(cats==2, 1), scores(cats==2, 2), 12, color2, 'filled')

subplot(1,2,2); hold on
plotFilled(1:size(waveformsSession, 2), waveformsSession(cats==1,:), color1);
plotFilled(1:size(waveformsSession, 2), waveformsSession(cats==2,:), color2);
%%
for i = 1:length(allSessions)
    session = allSessions{i};
    unit = allUnits{i};
    pd = parseSessionString_df(session, root, sep);
    animalName = pd.animalName;
    sampFile = [animalName category '_', modelName];
    path = [root animalName sep animalName 'sorted' sep 'stan' sep 'bernoulli' sep modelName sep category sep];
    fprintf([session unit '\n']);
    % paths
    neuralynxDataPath = [pd.sortedFolder session '_sessionData_nL.mat'];
    % load behavior and neurons
    load(neuralynxDataPath)

    if ~strcmp(session,prevSession) % avoiding recomputing model variables for different unit from same session
        %% behavior preparation 
        % parse behavior
        os = behAnalysisNoPlot_opMD(session, 'simpleFlag', 1);
        choice = os.allChoices';
        choice(choice<0) = 0;
        outcome = abs(os.allRewards)';
        choice = choice(1:min(length(choice), maxTrial));
        outcome = outcome(1:length(choice));
        outcomeL = outcome;
        outcomeL(outcomeL==0) = -1;
        outcomeR = outcome;
        outcomeR(outcomeR==0) = -1;
        outcomeL(choice==1) = 0;
        outcomeR(choice==0) = 0;
        responseInds = os.responseInds(1:min(length(choice), maxTrial)); 
        preRwd = [NaN abs(os.allRewards(1:end-1))]';
        %% behavior
        % switch
        svs = zeros(length(os.responseInds),1);
        svs(os.changeChoice_Inds) = 1;
        svsNext = [svs(2:end); NaN];
        svsWhenNrwd = svsNext;
        svsWhenNrwd(os.rwd_Inds) = NaN;
        [t,~,noSession] = getStanModelParams_samps(modelName, [path sampFile '.mat'], 2000, 'sessionName', session);
        
        % diff value
        Qdiff = abs(t.Q(:,2)-t.Q(:,1));
        % total value
        Qsum = sum(t.Q,2);
        % prepe
        prePe = [NaN; t.pe(1:end-1)];
        % pe
        pe = t.pe;
        % dawExp
        dawExp = double(t.probChoice <= 0.5);
        % confidence
        choiceConf = 2.*t.probChoice - 1;
        %time in session
        timeInSession = [sessionData(responseInds).CSon]' - sessionData(responseInds(1)).CSon;
        % chosen valie
        Qchosen  = zeros(length(choice),1);
        Qunchosen  = zeros(length(choice),1);
        QchosenUpdate = NaN(length(choice),1);
        for j = 1:length(choice)
            if j < length(choice)
                if choice(j)>0
                    Qchosen(j) = t.Q(j,2);
                    Qunchosen(j) = t.Q(j,1);
                    QchosenUpdate(j) = t.Q(j+2);
                else
                    Qchosen(j) = t.Q(j,1);
                    Qunchosen(j) = t.Q(j,2);
                    QchosenUpdate(j) = t.Q(j+1);
                end
            else                
                if choice(j)>0
                    Qchosen(j) = t.Q(j,2);
                    Qunchosen(j) = t.Q(j,1);
                else
                    Qchosen(j) = t.Q(j,1);
                    Qunchosen(j) = t.Q(j,2);
                end
                
            end
        end
        % bias side
        biasSide = zeros(size(responseInds))';
        biasInd = contains(paramNames, 'bias');
        if mean(t.params(:,biasInd))>0
            biasSide(os.lickR_Inds)=1;
        else
            biasSide(os.lickL_Inds)=1;
        end
        lickLat = os.lickLatLogZ';
        rightSide = zeros(size(pe));
        rightSide(os.allChoices>0)=1;
        rightSide(os.allChoices<=0)=-1;
        preITI = os.timeBtwn';
        % consecutive no rewards
        conNrwds = zeros(size(pe));
        for j = 1:length(choice)
            if outcome(j) == 0
                k = 1;
                while j-k>0 
                    if outcome(j-k)==0
                        k = k+1;
                    else
                        break
                    end
                end                
                conNrwds(j) = k;
            end
        end
        
        
        prevSession = session;
    end
    
    %% get neuron activity
    spikeFields = fields(sessionData);
    clust = find(contains(spikeFields,unit));
    allTrial_spike_choice = {};
    for k = 1:length(os.responseInds)
        if os.responseInds(k) == 1
            prevTrial_spike = [];
        else
            prevTrial_spikeInd = [sessionData(os.responseInds(k)-1).(spikeFields{clust})] > (sessionData(os.responseInds(k)).respondTime-tb*1000);
            prevTrial_spike = sessionData(os.responseInds(k)-1).(spikeFields{clust})(prevTrial_spikeInd) - sessionData(os.responseInds(k)).respondTime;
        end

        currTrial_spikeInd = sessionData(os.responseInds(k)).(spikeFields{clust}) < sessionData(os.responseInds(k)).respondTime+tf*1000 ... 
            & sessionData(os.responseInds(k)).(spikeFields{clust}) > sessionData(os.responseInds(k)).respondTime-tb*1000;
        currTrial_spike = sessionData(os.responseInds(k)).(spikeFields{clust})(currTrial_spikeInd) - sessionData(os.responseInds(k)).respondTime;

        allTrial_spike_choice{k} = [prevTrial_spike currTrial_spike];
    end
    allTrial_spike_choice(cellfun(@isempty,allTrial_spike_choice)) = {zeros(1,0)}; 
    %% initialize lick matrices
    for j = 1:length(os.responseInds)
        trialDurDiff(j) = (sessionData(os.responseInds(j)).trialEnd - sessionData(os.responseInds(j)).CSon)- tf*1000;
    end
    trialDurDiff(end) = 0; 

    % spike matric for GLM
    allTrial_spikeMatx_choice = zeros(length(os.responseInds),length(time));         
   for j = 1:length(allTrial_spike_choice)
        tempSpike = allTrial_spike_choice{j};
        tempSpike = tempSpike + tb*1000; % add this to pad time for SDF
        if any(tempSpike == 0)
            tempSpike(tempSpike == 0) = 1;
        end
        allTrial_spikeMatx_choice(j,tempSpike) = 1;
        if trialDurDiff(j) < 0
            allTrial_spikeMatx_choice(j, isnan(allTrial_spikeMatx_choice(j, 1:end+trialDurDiff(j)))) = 0;  %converts within trial duration NaNs to 0's
        else
            allTrial_spikeMatx_choice(j, isnan(allTrial_spikeMatx_choice(j,:))) = 0;
        end
    end
    % focus window
    focusSpikes = nansum(allTrial_spikeMatx_choice(:, (focusWin(1)+tb*1000):(focusWin(2)+tb*1000)),2);
    
    figure;
    subplot(2,1,1);
    yyaxis left
    plot(1:length(os.responseInds), pe, 'LineWidth', 2, 'Color', [0.5 0.5 0.5]);
    yyaxis right
    plot(1:length(os.responseInds), focusSpikes, 'LineWidth', 2, 'Color', [1 0 0]);
    sgtitle([session ' ' unit ' ' num2str(cats(i))]);
    subplot(2,1,2);
    autocorr(focusSpikes);
end