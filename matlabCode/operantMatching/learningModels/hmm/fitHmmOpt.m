function [behSessionData, states, trans_fit, emis_fit] = fitHmmOpt(sessionName,plotFlag)    
    %%
    if nargin < 2
    plotFlag = 0;
    end
    %%
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
        behSessionDataPath = [root animalName sep sessionFolder sep 'sorted' sep 'session ' sessionName(end) sep sessionName '_sessionData_behav.mat'];
    else
        behSessionDataPath = [root animalName sep sessionFolder sep 'sorted' sep 'session' sep sessionName '_sessionData_behav.mat'];
    end

    if exist(behSessionDataPath,'file')
        load(behSessionDataPath)
    else
        [behSessionData, blockSwitch, blockSwitchL, blockSwitchR] = generateSessionData_operantMatchingDecoupledRwdDelay(sessionName);
    end
    
    if ~exist('behSessionData', 'var')
        behSessionData = sessionData;
    end
    
    responseInds = find(~isnan([behSessionData.rewardTime])); % find CS+ trials with a response in the lick window
    allReward_R = [behSessionData(responseInds).rewardR]; 
    allReward_L = [behSessionData(responseInds).rewardL]; 
    allChoices = NaN(1,length(behSessionData(responseInds)));
    allChoices(~isnan(allReward_R)) = 2;
    allChoices(~isnan(allReward_L)) = 1;

    %% fit hmm
    trans_guess = [0.4, 0.3, 0.3;...
                   0.2, 0.8, 0;...
                   0.2, 0, 0.8];
    emis_guess = [0.5, 0.5;...
                0, 1;
                1, 0];

    [trans_est,emis_est] = hmmtrainPrior(allChoices,trans_guess,emis_guess); % should be done separately for each session because of bias. 
    %% compare
    % fit
    states = hmmviterbi(allChoices,trans_est,emis_est);
   

    %append to behavior
    for i = 1:length(responseInds)
        behSessionData(responseInds(i)).hmm = states(i);
    end
    if exist('blockSwitchL', 'var')
        save(behSessionDataPath, 'behSessionData', 'blockSwitch', 'blockSwitchL', 'blockSwitchR')
    else
        save(behSessionDataPath, 'behSessionData', 'blockSwitch')
    end
    
    %% switches
    stateswitch = [0 states(2:end) ~= states(1:end-1)];
    choiceswitch = [0 allChoices(2:end) ~= allChoices(1:end-1)];
    blockedges = find([allChoices(2:end) ~= allChoices(1:end-1)] ~= 0)+1;
    isiLen = diff([1 blockedges length(allChoices)+1]);
    [lenDistri,~] = histcounts(isiLen, [1:max(isiLen)+1]-0.5); 
    expfit = singleExpFit1para(lenDistri/sum(lenDistri), 1:max(isiLen));
    b = expfit.b;
    fitcurve1 = 1/expfit.b*exp(-(1/expfit.b)*(linspace(0.5,max(isiLen)+0.5,500)));
    expfit = twoExpFit(lenDistri/sum(lenDistri), 1:max(isiLen));
    bl = expfit.a;
    bs = expfit.b;
    fitcurve2 = expfit.c*(1/expfit.b*exp(-(1/expfit.b)*(linspace(0.5,max(isiLen)+0.5,500)))) + ...
        (1 - expfit.c)*(1/expfit.a*exp(-(1/expfit.a)*(linspace(0.5,max(isiLen)+0.5,500))));
    %% Plot Raw Data
    lastwarn('');
if plotFlag == 1
    figure;
    suptitle(sessionName);
    rMag = 1;
    nrMag = rMag/2;

    % trial plot
    subplot(3,4,[1:4]); hold on
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
    title('Opt');
    text(1,1.2,lastwarn);
%% plot block length distribution
    subplot(3,4,5);hold on; 
    histogram(isiLen, 0.5:1:max(isiLen)+0.5, 'Normalization', 'probability');
    plot(linspace(0.5, max(isiLen) + 0.5,500),fitcurve1, 'LineWidth',2);
    plot(linspace(0.5, max(isiLen) + 0.5,500),fitcurve2, 'LineWidth',2);
    text(0.25*max(isiLen), 0.7*max(lenDistri/sum(lenDistri)), sprintf(' b = %s \n bs = %s \n bl = %s', num2str(round(b,2)), num2str(round(bs,2)), num2str(round(bl,2))));
    
    screen = get(0,'Screensize');
    screen(4) = screen(4) - 100;
    set(gcf, 'Position', screen)    
end 