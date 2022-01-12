function cmpDVwithParams(xlFile, sheet, col, modelName, paramCtrl, paramInhi)
dayList = getDayList(xlFile, sheet, col);
%%
[paramsCtrl] = getStanModelParams_sampsOnly(sheet, paramCtrl, modelName, 200, 'sessionParamsFlag', 0);
[paramsInhi] = getStanModelParams_sampsOnly(sheet, paramInhi, modelName, 200, 'sessionParamsFlag', 0);
llCtrl = zeros(1,length(dayList));
llInhi = zeros(1,length(dayList));
for dayInd = 1:length(dayList)
    [tCtrl] = inferModelVar(dayList{dayInd}, paramsCtrl, modelName);
    [tInhi] = inferModelVar(dayList{dayInd}, paramsInhi, modelName);
    llCtrl(dayInd) = tCtrl.LH;
    llInhi(dayInd) = tInhi.LH;
    %% cmpSession
    figure;
    rMag = 1;
    nrMag = rMag/2;
    % trial plot
    subplot(3,1,1);
    hold on;
    [behSessionData] = loadBehavioralData([dayList{dayInd} '.asc']);
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
    plot(1:length(tCtrl.probChoice), tCtrl.Q(:,2)-tCtrl.Q(:,1), 'lineWidth', 2, 'color', [0.3 0.3 0.3]);
    plot(1:length(tInhi.probChoice), tInhi.Q(:,2)-tInhi.Q(:,1), 'lineWidth', 2, 'color', [0.3 1 0.3]);
    title('valueDiff')
    % choiceProb
    subplot(3,1,2);
    hold on;
    [behSessionData] = loadBehavioralData([dayList{dayInd} '.asc']);
    os = behAnalysisNoPlot_opMD(dayList{dayInd}, 'simpleFlag',1);
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
    choiceProbCtrl = tCtrl.probChoice;
    choiceProbCtrl(os.allChoices<0) = 1 - choiceProbCtrl(os.allChoices<0);
    choiceProbInhi = tInhi.probChoice;
    choiceProbInhi(os.allChoices<0) = 1 - choiceProbInhi(os.allChoices<0);

    plot(1:length(tCtrl.probChoice), 2*choiceProbCtrl-1, 'lineWidth', 2, 'color', [0.3 0.3 0.3]);
    plot(1:length(tInhi.probChoice), 2*choiceProbInhi-1, 'lineWidth', 2, 'color', [0.3 1 0.3]);
    title('pronChoice')
    
    % chosenProb
    subplot(3,1,3);
    hold on;
    [behSessionData] = loadBehavioralData([dayList{dayInd} '.asc']);
    os = behAnalysisNoPlot_opMD(dayList{dayInd}, 'simpleFlag',1);
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
    choiceProbCtrl = tCtrl.probChoice;
    choiceProbInhi = tInhi.probChoice;
    patch([1:length(tCtrl.probChoice), flip(1:length(tCtrl.probChoice))], [5*(log(choiceProbInhi./choiceProbCtrl)); zeros(length(tCtrl.probChoice),1)], 'r', 'FaceAlpha', 0.3, 'EdgeColor', 'none')
%     patch([1:length(tCtrl.probChoice), flip(1:length(tCtrl.probChoice))], [choiceProbCtrl; flip(-1*choiceProbCtrl)], [1 0 0], 'FaceAlpha', 0.3, 'EdgeColor', 'none');
%     patch([1:length(tInhi.probChoice), flip(1:length(tInhi.probChoice))], [choiceProbInhi; flip(-1*choiceProbInhi)], [0 1 0], 'FaceAlpha', 0.3, 'EdgeColor', 'none');
%     
    title('probChosen')
    
    
    suptitle(dayList{dayInd});
    screenSize = get(0,'Screensize');
    screenSize(4) = screenSize(4) - 100;
    set(gcf, 'renderer', 'painters', 'position', screenSize)
end
%% plot
paramNames = getParamNames_dF(modelName,0);
colors = cool(length(paramNames));
figure;
for i = 1:length(paramNames)
    subplot(1,length(paramNames)+1,i); hold on;
    low = min([paramsCtrl(:,i); paramsInhi(:,i)]');
    high = max([paramsCtrl(:,i); paramsInhi(:,i)]');
    edges = linspace(low, high, 25);
    histogram(paramsCtrl(:,i), edges, 'Normalization', 'probability','FaceColor', colors(i,:));
    histogram(paramsInhi(:,i), edges, 'Normalization', 'probability','FaceColor', [0 1 0]);
    title(paramNames{i});
end
subplot(1,length(paramNames)+1,length(paramNames)+1); hold on;
scatter(llCtrl,llInhi, 12, [1 0 0], 'filled');
line(minmax([llCtrl; llInhi]'),minmax([llCtrl; llInhi]'), 'lineWidth', 1, 'LineStyle', '--', 'Color', [0.7 0.7 0.7]);
xlabel(paramCtrl);
ylabel(paramInhi);
suptitle([sheet ' ' modelName])
%%