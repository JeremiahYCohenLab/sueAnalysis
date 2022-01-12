function comparePredictionParamsGrad(xlFile, sheet, col, modelName, cmpParam, cmpRange, varargin);
p = inputParser;
% default parameters if none given
p.addParameter('numBins', 20);
p.parse(varargin{:});
dayList = getDayList(xlFile, sheet, col);
paramNames = getParamNames_dF(modelName, 1);
cmpInds = strcmp(paramNames, cmpParam);
paramGrad = linspace(cmpRange(1), cmpRange(2), p.Results.numBins);
%%
ll = zeros(length(dayList),1);
llGrad = zeros(length(dayList),p.Results.numBins);
for dayInd = 1:length(dayList)
    [params,~,l] = getStanModelParams_sampsOnly(sheet, col, modelName, 200, 'sessionParamsFlag', 1, 'sessionName', dayList{dayInd});
    paramsCmp = repmat(params,1,1,p.Results.numBins);
    paramsCmp(:,cmpInds,:) = repmat(paramGrad, 200, 1);
    [t] = inferModelVar(dayList{dayInd}, params, modelName);
    ll(dayInd) = t.LH;
    tGrad = {};
    for i = 1:p.Results.numBins
        [tGrad{i}] = inferModelVar(dayList{dayInd}, squeeze(paramsCmp(:,:,i)), modelName);
        llGrad(dayInd,i) = tGrad{i}.LH;
    end

    %% cmpSession
    figure;
    rMag = 1;
    nrMag = rMag/2;
    colors = cool(p.Results.numBins); 
    % trial plot
    subplot(3,5,1:4);
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
    
    for i = 1:p.Results.numBins
        plot(1:length(tGrad{i}.probChoice), tGrad{i}.Q(:,2)-tGrad{i}.Q(:,1), 'lineWidth', 2, 'color', colors(i,:));
    end
    plot(1:length(t.probChoice), t.Q(:,2)-t.Q(:,1), 'lineWidth', 2, 'color', [0.3 0.3 0.3]);
    title('valueDiff')
    
    % choiceProb
    subplot(3,5,6:9);
    hold on;
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

    
    for i = 1:p.Results.numBins
        choiceProb = tGrad{i}.probChoice;
        choiceProb(os.allChoices<0) = 1 - choiceProb(os.allChoices<0);
        plot(1:length(tGrad{i}.probChoice), 2*choiceProb-1, 'lineWidth', 2, 'color', colors(i,:));
    end
    choiceProb = t.probChoice;
    choiceProb(os.allChoices<0) = 1 - choiceProb(os.allChoices<0);
    plot(1:length(t.probChoice), 2*choiceProb-1, 'lineWidth', 2, 'color', [0.3 0.3 0.3]);
    
    title('probChoice')
    
    % chosenProbLL
    subplot(3,5,11:14);
    hold on;
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
    choiceProb = t.probChoice;
    for i = 1:p.Results.numBins
        choiceProbGrad = tGrad{i}.probChoice;
        plot(1:length(t.probChoice), 4*(log(choiceProbGrad)-log(choiceProb)), "LineWidth", 1, 'color', colors(i,:));
    end
    title('probChosenLL')
    
    subplot(3,5,5); hold on
    histogram(params(:,cmpInds), 20, 'Normalization', 'Probability', 'FaceColor', [0.6 0.6 0.6]);
    for i = 1:p.Results.numBins
        plot([paramGrad(i) paramGrad(i)], [0 0.2], 'LineWidth', 2, 'Color', colors(i,:));
    end
    plot(paramGrad, llGrad(dayInd, :)/200, 'LineWidth', 2, 'Color', colors(i,:));
    plot(minmax(paramGrad), [ll(dayInd), ll(dayInd)]/200, 'LineWidth', 2, 'Color', [0.6 0.6 0.6])
    
    suptitle(dayList{dayInd});
    screenSize = get(0,'Screensize');
    screenSize(4) = screenSize(4) - 100;
    set(gcf, 'renderer', 'painters', 'position', screenSize)
end
%% plot

%%