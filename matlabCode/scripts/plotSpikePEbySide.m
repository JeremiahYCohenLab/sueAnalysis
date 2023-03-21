load('F:\tmpData\popGLMSepSide.mat')
%% scatter Tstats
figure2;
hold on;
line([-15 10], [0 0], 'color', [0.7 0.7 0.7], 'LineStyle', '--')
line([0 0], [-15 20], 'color', [0.7 0.7 0.7], 'LineStyle', '--')
scatter(populationTStatsL(5,1,:), populationTStatsR(5,1,:), 12, [0.5 0.5 0.5], 'filled');
scatter(populationTStatsL(5,1,abs(populationSigL(5,1,:))>0), populationTStatsR(5,1,abs(populationSigL(5,1,:))>0), 12, 'c', 'filled');
scatter(populationTStatsL(5,1,abs(populationSigR(5,1,:))>0), populationTStatsR(5,1,abs(populationSigR(5,1,:))>0), 12, 'm', 'filled');
scatter(populationTStatsL(5,1,abs(populationSigR(5,1,:))>0&abs(populationSigL(5,1,:))>0), populationTStatsR(5,1,abs(populationSigR(5,1,:))>0&abs(populationSigL(5,1,:))>0), 12, 'b', 'filled');
xlabel('L')
ylabel('R')
title('Tstats')
%% scatter coeff
figure2;
hold on;
line([-1.5 1.5], [0 0], 'color', [0.7 0.7 0.7], 'LineStyle', '--')
line([0 0], [-1.5 1.5], 'color', [0.7 0.7 0.7], 'LineStyle', '--')
scatter(populationCoeffsL(5,1,:), populationCoeffsR(5,1,:), 12, [0.5 0.5 0.5], 'filled');
scatter(populationCoeffsL(5,1,abs(populationSigL(5,1,:))>0), populationCoeffsR(5,1,abs(populationSigL(5,1,:))>0), 12, 'c', 'filled');
scatter(populationCoeffsL(5,1,abs(populationSigR(5,1,:))>0), populationCoeffsR(5,1,abs(populationSigR(5,1,:))>0), 12, 'm', 'filled');
scatter(populationCoeffsL(5,1,abs(populationSigR(5,1,:))>0&abs(populationSigL(5,1,:))>0), populationCoeffsR(5,1,abs(populationSigR(5,1,:))>0&abs(populationSigL(5,1,:))>0), 12, 'b', 'filled');
xlabel('L')
ylabel('R')
title('Coeff')
%% plot spike/s-rpe separate by side
figure2;
% pe bins
binNum = 6;
% decide unit Inds by significants
% unitInd = find((populationSigL(5,1,:).*populationSigR(5,1,:)>=0)&(populationSigL(5,1,:)<0 | populationSigR(5,1,:)<0));
% colors
colorL = [1 0 1];
colorR = [0 0.8 0.8];
% decide bins
binEdges = linspace(-1-0.001, 1+0.001, binNum+1);
% signs 
sigs = [-1 0 1];
for i = 1:3
    for j = 1:3
        % find target population
        unitInd = find(populationSigL(5,1,:) == sigs(i) & populationSigR(5,1,:) == sigs(j));
        % zscore and center by each side
        peMeanL = NaN(length(unitInd), binNum);
        peMeanR = NaN(length(unitInd), binNum);
        spikeMeanL = NaN(length(unitInd), binNum);
        spikeMeanR = NaN(length(unitInd), binNum);

        for u = 1:length(unitInd)
            pe = allPe{unitInd(u)};
            choice = allChoices{unitInd(u)};
            spike = zscore(allSpikes{unitInd(u)});
            spikeL = spike(choice == 0);
            spikeL = spikeL - mean(spikeL); % center to get rid of bias
            spikeR = spike(choice == 1);
            spikeR = spikeR - mean(spikeR); % center to get rid of bias
            peL = pe(choice==0);
            peR = pe(choice==1);
            for b = 1:binNum
                % left
                peMeanL(u, b) = mean(peL(peL>binEdges(b) & peL<=binEdges(b+1)));
                spikeMeanL(u,b) = mean(spikeL(peL>binEdges(b) & peL<binEdges(b+1)));
                % right
                peMeanR(u, b) = mean(peR(peR>binEdges(b) & peR<=binEdges(b+1)));
                spikeMeanR(u, b) = mean(spikeR(peR>binEdges(b) & peR<binEdges(b+1)));
            end
        end

        %% plot
        subplot(3,3,sub2ind([3,3], i, 4-(j)))
        hold on;
        peLMeanAll = mean(peMeanL, 1, 'omitnan');
        peRMeanAll = mean(peMeanR, 1, 'omitnan');

        plotFilled(peLMeanAll, spikeMeanL, colorL, 2);
        plotFilled(peRMeanAll, spikeMeanR, colorR, 2);

        if i == 1 && j ==1
            legend({'L', '', 'R', ''})
        end
               
        title([num2str(sign(sigs(i))) '/' num2str(sign(sigs(j))) ' ' num2str(length(unitInd))]);
    end
end

%%
load('F:\tmpData\popGLMInter.mat')
%% scatter Tstats
binInd = 9;
figure2;
hold on;
line([-15 10], [0 0], 'color', [0.7 0.7 0.7], 'LineStyle', '--')
line([0 0], [-8 10], 'color', [0.7 0.7 0.7], 'LineStyle', '--')
scatter(populationTStats(binInd,strcmp(regressors, 'outcome'),:), populationTStats(binInd,strcmp(regressors, 'outcome:rightSide'),:), 12, [0.5 0.5 0.5], 'filled');
scatter(populationTStats(binInd,strcmp(regressors, 'outcome'),abs(populationSig(binInd,strcmp(regressors, 'outcome'),:))>0), populationTStats(binInd,strcmp(regressors, 'outcome:rightSide'),abs(populationSig(binInd,strcmp(regressors, 'outcome'),:))>0), 12, 'c', 'filled');
scatter(populationTStats(binInd,strcmp(regressors, 'outcome'),abs(populationSig(binInd,strcmp(regressors, 'outcome:rightSide'),:))>0), populationTStats(binInd,strcmp(regressors, 'outcome:rightSide'),abs(populationSig(binInd,strcmp(regressors, 'outcome:rightSide'),:))>0), 12, 'm', 'filled');
scatter(populationTStats(binInd,strcmp(regressors, 'outcome'),abs(populationSig(binInd,strcmp(regressors, 'outcome'),:))>0&abs(populationSig(binInd,strcmp(regressors, 'outcome:rightSide'),:))>0), populationTStats(binInd,strcmp(regressors, 'outcome:rightSide'),abs(populationSig(binInd,strcmp(regressors, 'outcome'),:))>0&abs(populationSig(binInd,strcmp(regressors, 'outcome:rightSide'),:))>0), 12, 'b', 'filled');
line([-15 10], [15 -10], 'color', [0.7 0.7 0.7], 'LineStyle', '--')
xlabel('outcome')
ylabel('outcome:rightSide')
title('Tstats')
%% scatter coeff
figure2;
hold on;
line([-2 1.5], [0 0], 'color', [0.7 0.7 0.7], 'LineStyle', '--')
line([0 0], [-1.5 1.5], 'color', [0.7 0.7 0.7], 'LineStyle', '--')
scatter(populationCoeffs(binInd,strcmp(regressors, 'outcome'),:), populationCoeffs(binInd,strcmp(regressors, 'outcome:rightSide'),:), 12, [0.5 0.5 0.5], 'filled');
scatter(populationCoeffs(binInd,strcmp(regressors, 'outcome'),abs(populationSig(binInd,strcmp(regressors, 'outcome'),:))>0), populationCoeffs(binInd,strcmp(regressors, 'outcome:rightSide'),abs(populationSig(binInd,strcmp(regressors, 'outcome'),:))>0), 12, 'c', 'filled');
scatter(populationCoeffs(binInd,strcmp(regressors, 'outcome'),abs(populationSig(binInd,strcmp(regressors, 'outcome:rightSide'),:))>0), populationCoeffs(binInd,strcmp(regressors, 'outcome:rightSide'),abs(populationSig(binInd,strcmp(regressors, 'outcome:rightSide'),:))>0), 12, 'm', 'filled');
scatter(populationCoeffs(binInd,strcmp(regressors, 'outcome'),abs(populationSig(binInd,strcmp(regressors, 'outcome'),:))>0&abs(populationSig(binInd,strcmp(regressors, 'outcome:rightSide'),:))>0), populationCoeffs(binInd,strcmp(regressors, 'outcome:rightSide'),abs(populationSig(binInd,strcmp(regressors, 'outcome'),:))>0&abs(populationSig(binInd,strcmp(regressors, 'outcome:rightSide'),:))>0), 12, 'b', 'filled');
line([-1.5 1.5], [1.5 -1.5], 'color', [0.7 0.7 0.7], 'LineStyle', '--')
xlabel('outcome')
ylabel('outcome:rightSide')
title('coeff')
%% plot spike/s-rpe separate by side
figure2;
% pe bins
binNum = 6;
% decide unit Inds by significants
% unitInd = find((populationSigL(5,1,:).*populationSigR(5,1,:)>=0)&(populationSigL(5,1,:)<0 | populationSigR(5,1,:)<0));
% colors
colorL = [1 0 1];
colorR = [0 0.8 0.8];
% decide bins
binEdges = linspace(-1-0.001, 1+0.001, binNum+1);
% signs 
sigs = [-1 0 1];
for i = 1:3
    for j = 1:3
        % find target population
        unitInd = find(populationSig(binInd,1,:) == sigs(i) & populationSig(binInd,3,:) == sigs(j));
        % zscore and center by each side
        peMeanL = NaN(length(unitInd), binNum);
        peMeanR = NaN(length(unitInd), binNum);
        spikeMeanL = NaN(length(unitInd), binNum);
        spikeMeanR = NaN(length(unitInd), binNum);

        for u = 1:length(unitInd)
            pe = allPe{unitInd(u)};
            choice = allChoices{unitInd(u)};
            spike = zscore(allSpikes{unitInd(u)});
            spikeL = spike(choice == 0);
            spikeL = spikeL - mean(spikeL); % center to get rid of bias
            spikeR = spike(choice == 1);
            spikeR = spikeR - mean(spikeR); % center to get rid of bias
            peL = pe(choice==0);
            peR = pe(choice==1);
            for b = 1:binNum
                % left
                peMeanL(u, b) = mean(peL(peL>binEdges(b) & peL<=binEdges(b+1)));
                spikeMeanL(u,b) = mean(spikeL(peL>binEdges(b) & peL<binEdges(b+1)));
                % right
                peMeanR(u, b) = mean(peR(peR>binEdges(b) & peR<=binEdges(b+1)));
                spikeMeanR(u, b) = mean(spikeR(peR>binEdges(b) & peR<binEdges(b+1)));
            end
        end

        %% plot
        subplot(3,3,sub2ind([3,3], i, 4-(j)))
        hold on;
        peLMeanAll = mean(peMeanL, 1, 'omitnan');
        peRMeanAll = mean(peMeanR, 1, 'omitnan');

        plotFilled(peLMeanAll, spikeMeanL, colorL, 2);
        plotFilled(peRMeanAll, spikeMeanR, colorR, 2);

        if i == 1 && j ==1
            legend({'L', '', 'R', ''})
        end
               
        title([num2str(sign(sigs(i))) '/' num2str(sign(sigs(j))) ' ' num2str(length(unitInd))]);
    end
end














