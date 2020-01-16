function [behSessionData, states] = fitHmm(sessionName,plotFlag)    
    %%
    if nargin < 2
    plotFlag = 0;
    end
    %load session
    % Determine if computer is PC or Mac and set roots and separators appropriately
    [root, sep] = currComputer();

    % Generate the correct file path and see what files are available
    [animalName, date] = strtok(sessionName, 'd'); 
    animalName = animalName(2:end);
    date = date(1:9);
    sessionFolder = ['m' animalName date];
    filepath = [root animalName sep sessionFolder sep 'behavior' sep sessionName '.asc'];

    if isstrprop(sessionName(end), 'alpha')
        behSessionDataPath = [root animalName sep sessionFolder sep 'sorted' sep 'session ' sessionName(end) sep sessionName '_behSessionData_behav.mat'];
    else
        behSessionDataPath = [root animalName sep sessionFolder sep 'sorted' sep 'session' sep sessionName '_sessionData_behav.mat'];
    end

    if exist(behSessionDataPath,'file')
        load(behSessionDataPath)
    else
        [behSessionData, blockSwitch, blockSwitchL, blockSwitchR] = generateSessionData_operantMatchingDecoupled(sessionName);
    end

    responseInds = find(~isnan([behSessionData.rewardTime])); % find CS+ trials with a response in the lick window
    allReward_R = [behSessionData(responseInds).rewardR]; 
    allReward_L = [behSessionData(responseInds).rewardL]; 
    allChoices = NaN(1,length(behSessionData(responseInds)));
    allChoices(~isnan(allReward_R)) = 2;
    allChoices(~isnan(allReward_L)) = 1;

    %%
    %fit hmm OIT start
    lastwarn('');
    trans_guess = [0.9, 0, 0.1;...
                 0, 0.9, 0.1;...
                 0.4, 0.4, 0.2];
    emis_guess = [1, 0;...
                0, 1;
                0.5, 0.5];
    if allChoices(1) == 2
        allChoices = 3 - allChoices;
        [trans_est,emis_est] = hmmtrain(allChoices,trans_guess,emis_guess,'Tolerance',1e-7,'Maxiterations', 1e5); % should be done separately for each session because of bias.
        bias = abs(emis_est(3,1)-0.5);
        emis_fit =  2 * bias *emis_guess + (1 - 2 * bias)*emis_est;
        tbias = abs(trans_est(3,1) - trans_est(3,2));
        trans_fit = trans_est;
        trans_fit(3,1) = tbias * 0.5*(1 - trans_est(3,3)) + (1 - tbias)*trans_est(3,1);
        trans_fit(3,2) = tbias * 0.5*(1 - trans_est(3,3)) + (1 - tbias)*trans_est(3,2);
        statestemp = hmmviterbi(allChoices,trans_fit,emis_fit);
        allChoices = 3 - allChoices;
        states(statestemp == 1) = 2;
        states(statestemp == 2) = 1;
        states(statestemp == 3) = 3;
    else
        [trans_est,emis_est] = hmmtrain(allChoices,trans_guess,emis_guess,'Tolerance',1e-7,'Maxiterations', 1e5); % should be done separately for each session because of bias.
        bias = abs(emis_est(3,1)-0.5);
        emis_fit =  2 * bias *emis_guess + (1 - 2 * bias)*emis_est;
        tbias = abs(trans_est(3,1) - trans_est(3,2));
        trans_fit = trans_est;
        trans_fit(3,1) = tbias * 0.5*(1 - trans_est(3,3)) + (1 - tbias)*trans_est(3,1);
        trans_fit(3,2) = tbias * 0.5*(1 - trans_est(3,3)) + (1 - tbias)*trans_est(3,2);
        statestemp = hmmviterbi(allChoices,trans_fit,emis_fit);
        states = statestemp;
    end

    %append to behavior
    for i = 1:length(responseInds)
        behSessionData(responseInds(i)).hmm = states(i);
    end

    save(behSessionDataPath, 'behSessionData', 'blockSwitch', 'blockSwitchL', 'blockSwitchR')
    %%
    % Plot Raw Data
if plotFlag == 1
    rMag = 1;
    nrMag = rMag/2;

    % trial plot
    figure;
    suptitle(sessionName);
    subplot(3,8,[1:8]); hold on
    for i = 1:length(states)
        if states(i) == 1
           fill([i-0.5, i+0.5, i+0.5, i-0.5],[-1*rMag, -1*rMag, rMag, rMag],[1 0.7 1],'LineStyle','none');
        else
            if states(i) == 2
                fill([i-0.5, i+0.5, i+0.5, i-0.5],[-1*rMag, -1*rMag, rMag, rMag],[0.7 1 1],'LineStyle','none');
            else
                fill([i-0.5, i+0.5, i+0.5, i-0.5],[-1*rMag, -1*rMag, rMag, rMag],[0.7 0.7 0.7],'LineStyle','none');
            end
        end
    end

    j = 0;
    for i = 1:length(behSessionData)
        if strcmp(behSessionData(i).trialType,'CSplus')
            if ~isnan(behSessionData(i).rewardR)
                j = j+1;
                if behSessionData(i).rewardR == 1 % R side rewarded
                    plot([j j],[0 rMag],'k')
                else
                    plot([j j],[0 nrMag],'k') % R side not rewarded
                end
            elseif ~isnan(behSessionData(i).rewardL)
                j = j+1;
                if behSessionData(i).rewardL == 1 % L side rewarded
                    plot([j j],[-1*rMag 0],'k')
                else
                    plot([j j],[-1*nrMag 0],'k')
                end
            end
        end
    end

    text(0,1.5,'L/R');
    xlim([0.5 j+0.5]);
    ylabel('<-- L       R  -->')
    title('OIT start');
    text(1,1.2,lastwarn);
%% fit hmm ORE start
    lastwarn('');
    trans_guess = [0.2, 0.4, 0.4;...
                  0.1, 0.9, 0;...
                  0.1, 0, 0.9];
    emis_guess = [0.5, 0.5;
                 1, 0;...
                 0, 1];

    [trans_est,emis_est] = hmmtrain(allChoices,trans_guess,emis_guess,'Tolerance',1e-7,'Maxiterations', 1e5); % should be done separately for each session because of bias.
    % correction
    bias = abs(emis_est(1,1)-0.5);
    emis_fit =  2 * bias *emis_guess + (1 - 2 * bias)*emis_est;
    tbias = abs(trans_est(1,2) - trans_est(1,3));
    trans_fit = trans_est;
    trans_fit(1,2) = tbias * 0.5*(1 - trans_est(1,1)) + (1 - tbias)*trans_est(1,2);
    trans_fit(1,3) = tbias * 0.5*(1 - trans_est(1,1)) + (1 - tbias)*trans_est(1,3);
    
    
    % fit
    states = hmmviterbi(allChoices,trans_fit,emis_fit);

    %append to behavior
    for i = 1:length(responseInds)
        behSessionData(responseInds(i)).hmm = states(i);
    end

    save(behSessionDataPath, 'behSessionData', 'blockSwitch', 'blockSwitchL', 'blockSwitchR')
    %%
    % Plot Raw Data
if plotFlag == 1
    rMag = 1;
    nrMag = rMag/2;

    % trial plot
    subplot(3,8,[9:16]); hold on
    for i = 1:length(states)
        if states(i) == 2
           fill([i-0.5, i+0.5, i+0.5, i-0.5],[-1*rMag, -1*rMag, rMag, rMag],[1 0.7 1],'LineStyle','none');
        else
            if states(i) == 3
                fill([i-0.5, i+0.5, i+0.5, i-0.5],[-1*rMag, -1*rMag, rMag, rMag],[0.7 1 1],'LineStyle','none');
            else
                fill([i-0.5, i+0.5, i+0.5, i-0.5],[-1*rMag, -1*rMag, rMag, rMag],[0.7 0.7 0.7],'LineStyle','none');
            end
        end
    end

    j = 0;
    for i = 1:length(behSessionData)
        if strcmp(behSessionData(i).trialType,'CSplus')
            if ~isnan(behSessionData(i).rewardR)
                j = j+1;
                if behSessionData(i).rewardR == 1 % R side rewarded
                    plot([j j],[0 rMag],'k')
                else
                    plot([j j],[0 nrMag],'k') % R side not rewarded
                end
            elseif ~isnan(behSessionData(i).rewardL)
                j = j+1;
                if behSessionData(i).rewardL == 1 % L side rewarded
                    plot([j j],[-1*rMag 0],'k')
                else
                    plot([j j],[-1*nrMag 0],'k')
                end
            end
        end
    end

    text(0,1.5,'L/R');
    xlim([0.5 j+0.5]);
    ylabel('<-- L       R  -->')
    title('ORE start');
    text(1,1.2,lastwarn);
end
end