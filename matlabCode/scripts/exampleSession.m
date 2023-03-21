session = 'mZS061d20210408';
s = behAnalysisNoPlot_opMD('mZS061d20210408', 'simpleFlag', 1);
allProbsL = [s.behSessionData(s.responseInds).rewardProbL];
allProbsR = [s.behSessionData(s.responseInds).rewardProbR];
swL = find(allProbsL(1:end-1) ~= allProbsL(2:end))+1;
swR = find(allProbsR(1:end-1) ~= allProbsR(2:end))+1;
swAll = unique([swL, swR]);
rwds = s.allRewards;
choices = s.allChoices;
kernel = ones(1,5);
kernel = kernel/sum(kernel);
smoothedChoice = conv(kernel, choices);
smoothedChoice = smoothedChoice(3:end-2);
smoothedChoice(1:2) = NaN;
smoothedChoice(end) = NaN;
smoothedChoice(end-1) = NaN;
smoothedRwd = conv(kernel, rwds);
smoothedRwd = smoothedRwd(3:end-2);
smoothedRwd(1:2) = NaN;
smoothedRwd(end-1:end) = NaN;
%%
figure2;hold on;
plot([s.rwd_Inds; s.rwd_Inds], [zeros(1, length(s.rwd_Inds)); rwds(s.rwd_Inds)], 'k')
plot([s.nrwd_Inds; s.nrwd_Inds], 0.5 * [zeros(1, length(s.nrwd_Inds)); choices(s.nrwd_Inds)], 'Color', [0.6 0.6 0.6])

plot([swL; swL]-0.5, [zeros(size(swL)); -ones(size(swL))], 'Color', [0.7 0.7 0.7], 'LineStyle', '--')
plot([swR; swR]-0.5, [zeros(size(swR)); ones(size(swR))], 'Color', [0.7 0.7 0.7], 'LineStyle', '--')

set(gca, 'TickDir', 'Out')
set(gca, 'Box', 'off')
xlabel('Trials')
ylabel('L <--------------> R')
%%
figure2; hold on;
plot(smoothedRwd, 'b', 'LineWidth', 2)
plot(smoothedChoice, 'k', 'LineWidth', 2)

set(gca, 'TickDir', 'Out')
set(gca, 'Box', 'off')
xlabel('Trials')
ylabel('L <--------------> R')
%% model fitting example
[root, sep] = currComputer();
modelName = '5params';
numSamps = 2000;
category = 'good';
pd = parseSessionString_df(session, root, sep);
params = getStanModelParams_sampsOnly(pd.animalName, category, modelName, numSamps, 'sessionName', session);
t = inferModelVar(session, params, modelName);
choicePredict = t.probChoice;
choicePredict(s.allChoices<0) = 1-choicePredict(s.allChoices<0);
choicesBin = choices;
choicesBin(choices<0) = 0;
smoothedChoiceB = conv(kernel, choicesBin);
smoothedChoiceB = smoothedChoiceB(3:end-2);
smoothedChoiceB(1:2) = NaN;
smoothedChoiceB(end) = NaN;
smoothedChoiceB(end-1) = NaN;
figure2; hold on;
plot(1:length(smoothedChoiceB), smoothedChoiceB, 'k', 'LineWidth',2);
plot(1:length(smoothedChoiceB), choicePredict, 'b', 'LineWidth',2);
legend({'smoothed choices', 'RL prediction'})
ylabel('L <--------------> R')
%%

%%