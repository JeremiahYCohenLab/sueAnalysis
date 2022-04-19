%% task simualtion and evaluation
% steady state matching
p1 = 0.1:0.2:0.9;
p2 = 0.1:0.2:0.9;
sim = 100;
tests = cell(length(p1), length(p2),sim);
%%
for i = 1:length(p1)
    for j = 1:length(p2)
        
            for s = 1:sim
                rewardProbs = taskSimulate(p1(i), p2(j), 100, [200,200],0);
                beh = taskPerform(rewardProbs,0,'sigma',0.2);
                tests{i,j,s} = beh;
            end
       
    end
end
%% undermatching plot
rewardFr = [];
choiceFr = [];
for i = 1:length(p1)
    for j = 1:length(p2)
        if i ~= j
            for s = 1:sim
                currS = tests{i,j,s};
                rewardFr = [rewardFr sum(currS.choices>0 & currS.rewards>0)/sum(currS.rewards>0) sum(currS.choices<1 & currS.rewards>0)/sum(currS.rewards>0)];
                choiceFr = [choiceFr mean(currS.choices) 1-mean(currS.choices)];
            end
        end
    end
end

figure; hold on;
scatter(rewardFr, choiceFr, 12, 'Filled');
line([0 1], [0 1], 'linestyle', '--', 'color', [0.7 0.7 0.7])
%% proportion choice imagesc
proChoice2 = zeros(length(p1),length(p2));
for i = 1:length(p1)
    for j = 1:length(p2)

            allChoices = [];
            for s = 1:sim
                currS = tests{i,j,s};
                allChoices = [allChoices currS.choices];
            end
            proChoice2(i,j) = mean(allChoices);

    end
end
figure;
imagesc(proChoice2);
colorbar;
%% glm for learning with same p, but different block lengths
lenTest = [20, 40, 60, 80];
p1 = 0.2;
p2 = 0.9;
LR = 0.06;
aConf = 0.25;
aFConf = 0.9;
aSup = 1;
sigma = 0.05;
sim = 20;
test = cell(length(lenTest),length(sigma),sim);

%%
    for j = 1:length(lenTest)
        for i = 1:length(sigma)
            currSigma = sigma(i);
            for s = 1:sim
                rewardProbs = taskSimulate(p1, p2, 500, [lenTest(j), lenTest(j)],0);
                beh = taskPerformDynSigma(rewardProbs,0,'aP', LR, 'aN', LR, 'aSup', aSup, 'sigma', currSigma, 'aConf', aConf, 'aFConf', aFConf);
                test{j,i,s} = beh;
            end
        end
    end
%% mean sigma
sigmaMean = zeros(length(lenTest), length(sigma));
sigmaSem = zeros(length(lenTest), length(sigma));
for j = 1:length(lenTest)
    for i = 1:length(sigma)
        currSigma = [];
        for s = 1:sim
            beh = test{j,i,s};
            currSigma = [currSigma, mean(beh.sigmaCurr)];
        end
        sigmaMean(j,i) = mean(currSigma);
        sigmaSem(j,i) = std(currSigma);
    end
end

figure;
errorbar(lenTest, sigmaMean, sigmaSem, 'color', [0 0.5 1], 'linewidth', 4);
line([min(lenTest), max(lenTest)],[sigma sigma], 'color', [0.5 0.5 0.5], 'linewidth', 4)
xlabel('blockLength', 'FontSize', 24)
ylabel('mean(\sigma)','FontSize', 24)
set(gcf,'color','w');
set(gca,'TickDir','out');
set(gca,'XTick', 20:20:80,'FontSize', 18)
set(gca,'ytick',0:0.05:0.2,'FontSize', 18)
ylim([0 0.215]);
%%
sigmaTest = sigma;
colors = cool(length(sigmaTest));
figure; hold on;
for h = 1:length(lenTest)
    for i = 1:length(sigmaTest)
        allChoices = [];
        allRewards = [];
        for s = 1:sim
            currS = test{h,i,s}; 
            allChoices = [allChoices currS.choices];
            allRewards = [allRewards currS.rewards];
        end

        allRewards(allChoices==0) = -allRewards(allChoices==0);
        rewardsMat = NaN(10,length(allChoices));
        for j = 1:10
            rewardsMat(j,j+1:end) = allRewards(1:end-j);
        end
        nrwdMat = NaN(10,length(allChoices));
        allNoRewards = zeros(size(allRewards));
        allNoRewards(allRewards == 0 & allChoices == 0) = -1;
        allNoRewards(allRewards == 0 & allChoices == 1) = 1;
        for j = 1:10
            nrwdMat(j,j+1:end) = allNoRewards(1:end-j);
        end    

        glm = fitglm([rewardsMat;nrwdMat]', allChoices', 'Distribution', 'binomial');
        subplot(1,length(lenTest),h); hold on;
        idx = 2:11;
        ests = glm.Coefficients.Estimate;
        ci = coefCI(glm);
        ci = ci(:,2)-ci(:,1);
        errorbar(1:10, ests(idx), ci(idx), 'color', colors(i,:), 'linewidth', 2);
    end
    subplot(1,length(lenTest),h); line([0 11], [0 0], 'linestyle', '--', 'color', [0.7 0.7 0.7]);
    legend(cellfun(@num2str,num2cell(sigmaTest), 'UniformOutput', false))
end

%% compare choice block switch with differenet levels of exploration
colors = cool(length(sigmaTest));
preSw = 10;
postSw = 20;
figure;hold on;

for i = 1:length(lenTest)
    for j = 1:length(sigmaTest)
        swChoices = [];
        for s = 1:sim
            currBeh = test{i,j,s};
            allSwitches = find(diff(currBeh.rewardsP(:,2)) ~= 0) + 1; % first switchTrials
            for sw = 1:length(allSwitches)
                start = max([1 allSwitches(sw)-preSw]);
                last = min([length(currBeh.choices) allSwitches(sw)+postSw]);
                currSwChoices = [NaN(1,start-allSwitches(sw)+preSw) currBeh.choices(start:last) NaN(1,allSwitches(sw)+postSw-last)];

                if currBeh.rewardsP(allSwitches(sw),2) < currBeh.rewardsP(allSwitches(sw)-1,2)
                    swChoices = [swChoices; currSwChoices];
                else
                    swChoices = [swChoices; 1-currSwChoices];
                end
            end
        end
        nonnanT = sum(~isnan(swChoices), 1);
        subplot(1,length(lenTest),i); hold on;
        plot(-preSw:postSw, mean(swChoices,1,'omitnan'), 'color', 'm', 'linewidth', 2);
        fill([-preSw:postSw, postSw:-1:-preSw], [mean(swChoices,1,'omitnan') + std(swChoices,1,'omitnan')./sqrt(nonnanT-1), flip(mean(swChoices,1,'omitnan') - std(swChoices,1,'omitnan')./sqrt(nonnanT-1))], 'm','FaceAlpha', 0.2,'EdgeColor', 'none' )
        %errorbar(-preSw:postSw, mean(swChoices,1,'omitnan'), std(swChoices,1,'omitnan')./sqrt(nonnanT-1), 'color', colors(j,:), 'linewidth', 2);
        ylim([0 1]);
    end
    legend(cellfun(@num2str,num2cell(sigmaTest), 'UniformOutput', false))
    title(['blockLength = ' num2str(lenTest(i))])
end
%% compare confidence and surprise block switch with differenet levels of exploration
sigmaTest = sigma;
colors = cool(length(sigmaTest));
preSw = 10;
postSw = 25;
figure;hold on;

for i = 1:length(lenTest)
    for j = 1:length(sigmaTest)
        swCs = [];
        swCB = [];
        swSp = [];
        swSigma = [];
        for s = 1:sim
            currBeh = test{i,j,s};
            allSwitches = find(diff(currBeh.rewardsP(:,2)) ~= 0) + 1; % first switchTrials
            for sw = 1:length(allSwitches)
                start = max([1 allSwitches(sw)-preSw]);
                last = min([length(currBeh.choices) allSwitches(sw)+postSw]);
                currSwCs = [NaN(1,start-allSwitches(sw)+preSw) currBeh.c(start:last,2)'-currBeh.c(start:last,1)' NaN(1,allSwitches(sw)+postSw-last)];
                currConfBar = [NaN(1,start-allSwitches(sw)+preSw) currBeh.confBar(start:last) NaN(1,allSwitches(sw)+postSw-last)];
                currSup = [NaN(1,start-allSwitches(sw)+preSw) currBeh.confSup(start:last) NaN(1,allSwitches(sw)+postSw-last)]; 
                currSigma = [NaN(1,start-allSwitches(sw)+preSw) currBeh.sigmaCurr(start:last) NaN(1,allSwitches(sw)+postSw-last)];
                
                if currBeh.rewardsP(allSwitches(sw),2) < currBeh.rewardsP(allSwitches(sw)-1,2)
                    swCs = [swCs; currSwCs];
                    swCB = [swCB; 1-currConfBar];
                    swSp = [swSp; currSup];
                    swSigma = [swSigma; currSigma];
                else
                    swCs = [swCs; -currSwCs];
                    swCB = [swCB; 1-currConfBar];
                    swSp = [swSp; currSup];
                    swSigma = [swSigma; currSigma];
                end
            end
        end
        nonnanT = sum(~isnan(swCs), 1);
        subplot(4,length(lenTest),i); hold on;
        plot(-preSw:postSw, mean(swCs,1,'omitnan'), 'color', 'm', 'linewidth', 2);
        fill([-preSw:postSw, postSw:-1:-preSw], [mean(swCs,1,'omitnan') + std(swCs,1,'omitnan')./sqrt(nonnanT-1), flip(mean(swCs,1,'omitnan') - std(swCs,1,'omitnan')./sqrt(nonnanT-1))], 'm','FaceAlpha', 0.2,'EdgeColor', 'none' )
        title(['blockLen = ' num2str(lenTest(i))])
        %errorbar(-preSw:postSw, mean(swCs,1,'omitnan'), std(swCs,1,'omitnan')./sqrt(nonnanT-1), 'color', colors(j,:), 'linewidth', 2);
         subplot(4,length(lenTest),length(lenTest)+i); hold on;
        plot(-preSw:postSw, mean(swCB,1,'omitnan'), 'color', 'm', 'linewidth', 2);
        fill([-preSw:postSw, postSw:-1:-preSw], [mean(swCB,1,'omitnan') + std(swCB,1,'omitnan')./sqrt(nonnanT-1), flip(mean(swCB,1,'omitnan') - std(swCB,1,'omitnan')./sqrt(nonnanT-1))], 'm','FaceAlpha', 0.2,'EdgeColor', 'none' )
        subplot(4,length(lenTest),2*length(lenTest)+i); hold on;
        plot(-preSw:postSw, mean(swSp,1,'omitnan'), 'color', 'm', 'linewidth', 2);
        fill([-preSw:postSw, postSw:-1:-preSw], [mean(swSp,1,'omitnan') + std(swSp,1,'omitnan')./sqrt(nonnanT-1), flip(mean(swSp,1,'omitnan') - std(swSp,1,'omitnan')./sqrt(nonnanT-1))], 'm','FaceAlpha', 0.2,'EdgeColor', 'none' )
        subplot(4,length(lenTest),3*length(lenTest)+i); hold on;
        plot(-preSw:postSw, mean(swSigma,1,'omitnan'), 'color', 'm', 'linewidth', 2);
        fill([-preSw:postSw, postSw:-1:-preSw], [mean(swSigma,1,'omitnan') + std(swSigma,1,'omitnan')./sqrt(nonnanT-1), flip(mean(swSigma,1,'omitnan') - std(swSigma,1,'omitnan')./sqrt(nonnanT-1))], 'm','FaceAlpha', 0.2,'EdgeColor', 'none' )

    end
end
%% compare reward rate with different sigmas in different block lengths
sigmaTest = sigma;
rwdRates = zeros(length(lenTest), length(sigmaTest), sim);
colors = cool(length(lenTest));
for i = 1:length(sigmaTest)
    for j = 1:length(lenTest)
        for s = 1:sim
            currBeh = test{j,i,s};
            rwdRates(j,i,s) = mean(currBeh.rewards);
        end
    end
end

figure;
%%
hold on;
for j = 1:length(lenTest)
    bar =  std(squeeze(rwdRates(j,:,:)), 0, 1)/sqrt(sim-1);
    m = mean(squeeze(rwdRates(j,:,:)),1);
    plot([0, 0.3], [m m], 'color', colors(j,:), 'linewidth',2);

    fill([0, 0.3, 0.3, 0], [m-bar, m-bar, m+bar, m+bar], colors(j,:),'FaceAlpha', 0.2,'EdgeColor', 'none' )
   
end
legend(cellfun(@num2str,num2cell(lenTest), 'UniformOutput', false))
%%
xlabel('I_2 - I_1', 'FontSize', 24)
ylabel('P(choice=2)','FontSize', 24)
set(gcf,'color','w');
set(gca,'TickDir','out');
set(gca,'XTick',[])
set(gca,'ytick',[0 1],'FontSize', 18)
xlim([-0.08, 0.08])
%%
plot(1-confBar, 'color', [0.3 1 0.3], 'linewidth', 2.5);
plot(0.3*confSup, 'color', [1 0.3 0.3], 'linewidth', 2.5);
fill([1:length(choices), flip(1:length(choices))], [sigmaCurr, -flip(sigmaCurr)], 'FaceColor', [0 1 1], 'FaceAlpha', 0.5);
%%