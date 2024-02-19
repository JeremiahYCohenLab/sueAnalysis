function s = unitAutoCorrNp(session, unit1, unit2, varargin)
%task and model parameters
p = inputParser;
% default parameters if none given
p.addParameter('binSize', 2)% in ms  
p.addParameter('saveFigFlag', 0);
p.addParameter('lag', 100);% in ms
p.addParameter('jitter', 5); % in ms
p.addParameter('numShuffle', 500); % pairs of surrogates to generate
p.addParameter('alpha', 0.01);
p.addParameter('plot', 0)
p.parse(varargin{:});

binSize = p.Results.binSize;
stepSize = binSize;
jitter = 5*binSize;
% basic info
[root, sep] = currComputer();

% load unit 1
fileRoot = 'F:\npOptoRecordings\withOpto\';
unitFile1 = [fileRoot session sep 'sorted' sep 'timestamps' sep unit1 '_allSpikes.npy'];
spikeTime{1} = readNPY(unitFile1)*1000;
spikeTime{1} = spikeTime{1}(spikeTime{1}<=spikeTime{1}(end));

% load unit 2
unitFile2 = [fileRoot session sep 'sorted' sep 'timestamps' sep unit2 '_allSpikes.npy'];
spikeTime{2} = readNPY(unitFile2)*1000;
spikeTime{2} = spikeTime{2}(spikeTime{2}<=spikeTime{2}(end));

% load laser time
laserFile = [fileRoot session sep 'sorted' sep 'timestamps' sep 'laserTimes.npy'];
laserTimes = readNPY(laserFile)*1000;

% spike counts
edgesCorr = (min([spikeTime{1}; spikeTime{2}])-0.01):binSize:(min(laserTimes));
xAllCorr = 1000*histcounts(spikeTime{1}, edgesCorr)/binSize;
yAllCorr = 1000*histcounts(spikeTime{2}, edgesCorr)/binSize;

edges = (min([spikeTime{1}; spikeTime{2}])-0.01):binSize: (max([spikeTime{1}; spikeTime{2}])+0.01);
xAll = 1000*histcounts(spikeTime{1}, edges)/binSize;
yAll = 1000*histcounts(spikeTime{2}, edges)/binSize;

timeCorr = 0.5 * (edgesCorr(1:end-1) + edgesCorr(2:end))/60000;
time = 0.5 * (edges(1:end-1) + edges(2:end))/60000;
%% calculate autoCorr

    numLag = round(p.Results.lag/p.Results.binSize);
 
  
    corrCoAll = crosscorr(xAllCorr,yAllCorr, NumLags=numLag);
    
    % corrCoAll(numLag + 1) = NaN;

    lagTime = (-p.Results.lag):binSize:p.Results.lag;
 
    
    s.corrCoAll = corrCoAll;
    s.fr = 1000*[sum(spikeTime{1}<laserTimes(1))/(laserTimes(1) - spikeTime{1}(1));...
        sum(spikeTime{2}<laserTimes(1))/(laserTimes(1) - spikeTime{2}(1))];
    % corrCoAll(0.5*(length(corrCoAll)+1)) = NaN;
    % 
    % 
    % % consolidate into struct
    % s.corrCo = mean(corrCo, 'omitnan');
    % s.corrCoPre = mean(corrCoPre, 'omitnan');
    % s.corrCoAll = corrCoAll;
    % s.corrCoPreSess = corrCoPreSess;
    % s.lag = (-lag:lag)*stepSize;
    
    
    % generating surrogate
    corrSurr = zeros(p.Results.numShuffle, 2*numLag+1);
    spike2 = spikeTime{2};
    tempSpike = cell(p.Results.numShuffle,1);
    tempCnt = cell(p.Results.numShuffle,1);
    parfor a = 1:p.Results.numShuffle
        % surrogate

        tempSpike{a} = spike2 + (2*jitter+1)*rand(size(spike2,1), size(spike2,2)) - jitter -1;    

        tempCnt{a} = 1000 * histcounts(tempSpike{a}, edgesCorr)/binSize;

        corrSurr(a,:) = crosscorr(xAllCorr, tempCnt{a}, NumLags=numLag);

    end
    % find alpha boundary

    for i = 1:(2*numLag+1)
            tmp = sort(corrSurr(:,i));
            CI(1,i) = tmp(floor(0.5*p.Results.alpha*p.Results.numShuffle));
            CI(2,i) = tmp(ceil((1-0.5*p.Results.alpha)*p.Results.numShuffle));
    end
    s.CI = CI;
    s.mean = mean(corrSurr, 1);
%%
  if p.Results.plot
        figure2;
        subplot(1, 5, 1)
        hold on
        plot(lagTime, corrCoAll, 'LineWidth',2, 'Color','k');
        plot(lagTime, zeros(1, length(lagTime)));
        ylim([min(-0.05, min(corrCoAll, [], 'omitmissing')) max(corrCoAll, [], 'omitmissing')]);
        patch([lagTime flip(lagTime)], [CI(1,:), flip(CI(2,:))], [0.5 0.5 0.5], 'FaceAlpha', 0.25, 'edgecolor', 'none');
        xlabel('ms')
        subplot(1, 5, 2:5)
        hold on
        yyaxis left
        plot(time, xAll);
        yyaxis right
        plot(time, yAll);
       
        
        yyaxis left
        plot([laserTimes(1) laserTimes(1)]/(60 * 1000), [0 10], 'LineStyle','--', 'Color','k', 'LineWidth',2)
        ylabel('Spikes/s')

        sgtitle([session ' ' unit1 ' ' unit2], 'interpreter', 'none')
  end
  %%
%     s.CIPreTrial = CIPreTrial;
%     s.CIPreSess = CIPreSess;
%     s.meanPreTrial = mean(allPreTrial, 1, 'omitnan');
%     s.menPreSession = mean(allPreSession, 1, 'omitnan');
% %% plotting
% if p.Results.plotFlag
%     figure2;
%     subplot(1,3,1); hold on;
%     bar((-lag:lag)*stepSize, corrCoAll,'k');
%     legend('allSpikes');
% 
%     subplot(1,3,2); hold on;
%     bar((-lag:lag)*stepSize,mean(corrCo, 'omitnan'),'m', 'FaceAlpha', 0.7);
%     bar((-lag:lag)*stepSize,mean(corrCoPre, 'omitnan'),'c', 'FaceAlpha', 0.7);
%     plot((-lag:lag)*stepSize, mean(allPreTrial, 1, 'omitnan'), 'Color', 'c', 'LineWidth', 2);
%     plot(repmat((-lag:lag)'*stepSize, 1, 2), CIPreTrial', 'Color', 'c', 'LineWidth', 2, 'LineStyle', '--');
%     legend({'in trial', 'pre trial'})
%     title([unit1 unit2 ' crossCorr'], 'Interpreter','none');
% 
%     subplot(1,3,3); hold on;
%     bar((-lag:lag)*stepSize,corrCoPreSess,'k');
%     plot((-lag:lag)*stepSize, mean(allPreSession, 1, 'omitnan'), 'Color', [0.7 0.7 0.7], 'LineWidth', 2);
%     plot(repmat((-lag:lag)'*stepSize', 1, 2), CIPreSess', 'Color', [0.7 0.7 0.7], 'LineWidth', 2, 'LineStyle', '--');
%     legend('preSession');
%     % subplot(2, length(clust)+1, 1); hold on;
%     % spikeTimes = sessionData(clust(1)-min(allClustsInds)+1).allSpikes;
%     % histogram(diff(spikeTimes), 'Normalization', 'probability');
%     % title('ISI')
%     % xlabel('ms')
%     % set(gca, 'XScale', 'log')
%     % plot([2 2], [0 0.2], 'LineStyle', '--', 'Color', 'r');
%     % xlim([0 max(diff(sessionData(clust(1)-min(allClustsInds)+1).allSpikes))])
% 
%     sgtitle([session ' ' unit1 unit2], 'Interpreter', 'none');
%     screen = get(0,'Screensize');
%     screen(4) = screen(4) - 100;
%     set(gcf, 'Position', screen)
end
   
