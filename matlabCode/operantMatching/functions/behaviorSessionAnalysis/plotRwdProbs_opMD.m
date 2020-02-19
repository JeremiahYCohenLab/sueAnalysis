function plotRwdProbs_opMD(sessionName)

s = behAnalysisNoPlot_opMD(sessionName);

lProb = [];
bL = diff(s.blockSwitchL);
for i = 1:length(s.blockSwitchL)
    if i == length(s.blockSwitchL)
        lProb = [lProb s.behSessionData(s.blockSwitchL(i) + 1).rewardProbL * ones(1,(length(s.behSessionData) - s.blockSwitchL(end)))];
    else
        lProb = [lProb s.behSessionData(s.blockSwitchL(i) + 1).rewardProbL * ones(1, bL(i))];
    end
end

rProb = [];
bR = diff(s.blockSwitchR);
for i = 1:length(s.blockSwitchR)
    if i == length(s.blockSwitchR)
        rProb = [rProb s.behSessionData(s.blockSwitchR(i) + 1).rewardProbR * ones(1,(length(s.behSessionData) - s.blockSwitchR(end)))];
    else
        rProb = [rProb s.behSessionData(s.blockSwitchR(i) + 1).rewardProbR * ones(1, bR(i))];
    end
end


figure; hold on;
plot(lProb, '-k', 'linewidth', 3)
plot(rProb, '-', 'linewidth', 3, 'color', [0.65 0.65 0.65])
ylim([0 100])
yticks([0 10 40 70 100])
yticklabels([{'0'} {'0.1'} {'0.4'} {'0.7'} {'1'}])
legend('left', 'right')
set(gca, 'tickdir', 'out')
set(gcf, 'renderer', 'painters')

end
