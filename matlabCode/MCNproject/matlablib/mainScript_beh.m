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
scatter(rewardFr, choiceFr, 12, 'Filled', 'k');
line([0 1], [0 1], 'linestyle', '--', 'color', [0.7 0.7 0.7])
xlabel('rewardFraction','FontSize', 18)
ylabel('choiceFraction','FontSize', 18)
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
imagesc(p1,p2,proChoice2);
xlabel('P(R|choice = 1)','FontSize', 18)
ylabel('P(R|choice = 2)','FontSize', 18)
set(gcf,'color','w');
colormap cool
colorbar;
%% glm for learning with same p, but different block lengths
lenTest = [20, 40, 60, 80];
sigmaTest = [0.025,0.05,0.075, 0.10, 0.125, 0.15, 0.20, 0.25, 0.3];
p1 = 0.2;
p2 = 0.9;
LR = 0.06;
sim = 150;
test = cell(length(sigmaTest),length(lenTest),sim);
%%
for i = 1:length(sigmaTest)
    for j = 1:length(lenTest)
        for s = 1:sim
            rewardProbs = taskSimulate(p1, p2, 500, [lenTest(j), lenTest(j)],0);
            beh = taskPerform(rewardProbs,0,'sigma',sigmaTest(i),'aP', LR, 'aN',LR);
            test{i,j,s} = beh;
        end
    end
end
%%
colors = cool(length(sigmaTest));
figure; hold on;
for h = 1:length(lenTest)
    for i = 1:length(sigmaTest)
        allChoices = [];
        allRewards = [];
        for s = 1:sim
            currS = test{i,h,s}; 
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

%% compare block switch with differenet levels of exploration
colors = cool(length(sigmaTest));
preSw = 10;
postSw = 20;
figure;hold on;

for i = 1:length(lenTest)
    for j = 1:length(sigmaTest)
        swChoices = [];
        for s = 1:sim
            currBeh = test{j,i,s};
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
        plot(-preSw:postSw, mean(swChoices,1,'omitnan'), 'color', colors(j,:), 'linewidth', 2);
        fill([-preSw:postSw, postSw:-1:-preSw], [mean(swChoices,1,'omitnan') + std(swChoices,1,'omitnan')./sqrt(nonnanT-1), flip(mean(swChoices,1,'omitnan') - std(swChoices,1,'omitnan')./sqrt(nonnanT-1))], colors(j,:),'FaceAlpha', 0.2,'EdgeColor', 'none' )
        %errorbar(-preSw:postSw, mean(swChoices,1,'omitnan'), std(swChoices,1,'omitnan')./sqrt(nonnanT-1), 'color', colors(j,:), 'linewidth', 2);
        ylim([0 1])
    end
    legend(cellfun(@num2str,num2cell(sigmaTest), 'UniformOutput', false))
    title(['blockLength = ' num2str(lenTest(i))])
end
%% compare block switch with differenet levels of exploration
colors = cool(length(sigmaTest));
preSw = 10;
postSw = 25;
figure;hold on;

for i = 1:length(lenTest)
    for j = 6
        swCs = [];
        for s = 1:sim
            currBeh = test{j,i,s};
            allSwitches = find(diff(currBeh.rewardsP(:,2)) ~= 0) + 1; % first switchTrials
            for sw = 1:length(allSwitches)
                start = max([1 allSwitches(sw)-preSw]);
                last = min([length(currBeh.choices) allSwitches(sw)+postSw]);
                currSwCs = [NaN(1,start-allSwitches(sw)+preSw) currBeh.c(start:last,2)'-currBeh.c(start:last,1)' NaN(1,allSwitches(sw)+postSw-last)];

                if currBeh.rewardsP(allSwitches(sw),2) < currBeh.rewardsP(allSwitches(sw)-1,2)
                    swCs = [swCs; currSwCs];
                else
                    swCs = [swCs; -currSwCs];
                end
            end
        end
        nonnanT = sum(~isnan(swCs), 1);
        subplot(1,length(lenTest),i); hold on;
        plot(-preSw:postSw, mean(swCs,1,'omitnan'), 'color', 'c', 'linewidth', 2);
        fill([-preSw:postSw, postSw:-1:-preSw], [mean(swCs,1,'omitnan') + std(swCs,1,'omitnan')./sqrt(nonnanT-1), flip(mean(swCs,1,'omitnan') - std(swCs,1,'omitnan')./sqrt(nonnanT-1))], 'c','FaceAlpha', 0.2,'EdgeColor', 'none' )
    end
    legend(cellfun(@num2str,num2cell(sigmaTest), 'UniformOutput', false))
    title(['sigma = ' num2str(lenTest(i))])

end
%% compare reward rate with different sigmas in different block lengths
rwdRates = zeros(length(sigmaTest), length(lenTest), sim);
colors = cool(length(lenTest));
for i = 1:length(sigmaTest)
    for j = 1:length(lenTest)
        for s = 1:sim
            currBeh = test{i,j,s};
            rwdRates(i,j,s) = mean(currBeh.rewards);
        end
    end
end

figure; hold on;
for j = 1:length(lenTest)
    errorbar(sigmaTest, mean(squeeze(rwdRates(:,j,:)),2), std(squeeze(rwdRates(:,j,:)), 0, 2)/sqrt(sim-1), 'color', colors(j,:), 'linewidth',2);
end
legend(cellfun(@num2str,num2cell(lenTest), 'UniformOutput', false))
%%
tickLen = 0.8;
for t = 1:length(choices)
    line([t t], [0 sign(choices(t)-0.5)*(0.5*tickLen + 0.5*rewards(t)*tickLen)], 'linewidth', 1.5);
end

%%