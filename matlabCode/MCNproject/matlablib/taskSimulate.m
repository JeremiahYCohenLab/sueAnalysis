function [rewardProbs, allBlockLen]= taskSimulate(pLow, pHigh, trialNum, blockLenDist,plotFlag)
currblockLen = round(rand(1,1)*range(blockLenDist))+blockLenDist(1);
allBlockLen = currblockLen;
tInB = 1;
rewardProbs = zeros(trialNum,2);
startChoice = round(rand(1,1))+1;

rwdProb = [pLow, pHigh];


for t = 1:trialNum
    if tInB <= currblockLen
        rewardProbs(t,:) = rwdProb;
        tInB = tInB + 1;
    else
        currblockLen = round(rand(1,1)*range(blockLenDist))+blockLenDist(1);
        allBlockLen = [allBlockLen, currblockLen];
        tInB = 1;
        rwdProb = flip(rwdProb);
        rewardProbs(t,:) = rwdProb;
        tInB = tInB + 1;
    end
end

if plotFlag
    figure; hold on; subplot(1,4,1:3); hold on;
    plot(rewardProbs(:,1),'color','r','linewidth',1.5);
    plot(rewardProbs(:,2),'color','g','linewidth',1.5);
    legend({'red', 'green'})
    subplot(1,4,4)
    histogram(allBlockLen);
end