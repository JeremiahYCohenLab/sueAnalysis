function plotChoiceCurve(xlFile, sheet, col, model, varargin)
p = inputParser;
% default parameters if none given
p.addParameter('plotFlag', 0);
p.addParameter('numBins', 8);
p.addParameter('maxTrials', 1000);
p.parse(varargin{:});

dayList = getDayList(xlFile, sheet, col);

numSamps = 200;
combineChoices = [];
combineOutcomes = [];
combineHmmExp = [];
combineBias = [];
combineV = [];
combinePe = [];
combineLickLat = [];
combineLickRates = [];
combineLickRatesRwd = [];
combineSw = [];

for sess = 1:length(dayList)
%     fprintf([dayList{sess} '\n'])
    os = behAnalysisNoPlot_opMD(dayList{sess},'simpleFlag',1);
    switches = zeros(1, length(os.allChoices));
    switches(os.changeChoice_Inds) = 1;
    [animalName, date] = strtok(dayList{sess}, 'd'); 
    animalName = animalName(2:end);
    date = date(1:9);
    [params, model, ~, noSession] = getStanModelParams_sampsOnly(animalName, col, model, numSamps, 'sessionParamsFlag', 1, 'sessionName', dayList{sess});
    paramNames = getParamNames_dF(model,1);
    betaInd = contains(paramNames, 'beta');
    if noSession
        fprintf([dayList{sess} ' no good behavior \n'])
        continue
    end
    t = inferModelVar(dayList{sess}, params, model);
    combinePe = [combinePe; t.pe];
    combineBias = [combineBias; (mean(params(:,end))/mean(params(:,betaInd)))*ones(length(os.responseInds),1)];
    combineV = [combineV; t.Q];
    combineLickLat = [combineLickLat, os.lickLatZ];
    combineLickRates = [combineLickRates, os.lickRateZ];
    combineLickRatesRwd = [combineLickRatesRwd, os.lickRateRwdZ];
    combineSw = [combineSw, switches];
    combineHmmExp = [combineHmmExp, os.hmmStates==1];
    combineChoices = [combineChoices, 0.5*(os.allChoices+1)];
    combineOutcomes = [combineOutcomes, abs(os.allRewards)];
    
end
    
    % plot explore
    vORE = combineV(combineHmmExp>0,2) - combineV(combineHmmExp>0,1); % + combineBias(combineHmmExp>0,1);
    choicesORE = combineChoices(combineHmmExp>0);
    edges = linspace(min(vORE), max(vORE), p.Results.numBins+1);
%     edges = binEqualSize(vORE,p.Results.numBins);
    choiceOREMean = zeros(p.Results.numBins,1);
    choiceORESem = zeros(p.Results.numBins,1);
    vOREMean = zeros(p.Results.numBins,1);
    for k = 1:(length(edges)-1)
        if k < (length(edges)-1)
            vOREMean(k) = mean(vORE(vORE>=edges(k) & vORE < edges(k+1)));
            choiceOREMean(k) = mean(choicesORE(vORE>=edges(k) & vORE < edges(k+1)));
            choiceORESem(k) = sem_bern(choicesORE(vORE>=edges(k) & vORE < edges(k+1)));
        else
            vOREMean(k) = mean(vORE(vORE>=edges(k) & vORE <= edges(k+1)));
            choiceOREMean(k) = mean(choicesORE(vORE>=edges(k) & vORE <= edges(k+1)));
            choiceORESem(k) = sem_bern(choicesORE(vORE>=edges(k) & vORE <= edges(k+1)));            
        end
    end
    
    % plot exploit
    vOIT = combineV(combineHmmExp<1,2) - combineV(combineHmmExp<1,1)% + combineBias(combineHmmExp<1,1);
    choicesOIT = combineChoices(combineHmmExp<1);
    edges = linspace(min(vOIT), max(vOIT), p.Results.numBins+1);
%     edges = binEqualSize(vOIT,p.Results.numBins);
    choiceOITMean = zeros(p.Results.numBins,1);
    choiceOITSem = zeros(p.Results.numBins,1);
    vOITMean = zeros(p.Results.numBins,1);
    for k = 1:(length(edges)-1)
        if k < (length(edges)-1)
            vOITMean(k) = mean(vOIT(vOIT>=edges(k) & vOIT < edges(k+1)));
            choiceOITMean(k) = mean(choicesOIT(vOIT>=edges(k) & vOIT < edges(k+1)));
            choiceOITSem(k) = sem_bern(choicesOIT(vOIT>=edges(k) & vOIT < edges(k+1)));
        else
            vOITMean(k) = mean(vOIT(vOIT>=edges(k) & vOIT <= edges(k+1)));
            choiceOITMean(k) = mean(choicesOIT(vOIT>=edges(k) & vOIT <= edges(k+1)));
            choiceOITSem(k) = sem_bern(choicesOIT(vOIT>=edges(k) & vOIT <= edges(k+1)));            
        end
    end
    
    figure2;
    hold on;
    errorbar(vOITMean, choiceOITMean, choiceOITSem, 'color', 'b', 'LineWidth', 2);
    errorbar(vOREMean, choiceOREMean, choiceORESem, 'color', 'm', 'LineWidth', 2);
    title([sheet ' ' col ' ' model]);
    
end