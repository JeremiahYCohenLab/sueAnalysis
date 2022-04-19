    function pupilGLM_combined(xlFile, sheet, category, varargin)
%task and model parameters
p = inputParser;
% default parameters if none given
p.addParameter('maxTrial', 1000);
% p.addParameter('modelName','7params_absPePeAN_scale_int_bias_ord')
p.addParameter('modelName','5params')
p.addParameter('regressors', '1+pe+biasSide+pe*biasSide')
p.addParameter('binSize', 1000)% in ms
p.addParameter('stepSize', 1000)
% default time window is -2s to 10s to cue time, preLen+postLen+1 frames
p.addParameter('tb', 2)% in s 1.7
p.addParameter('tf', 10)% in s
p.addParameter('saveFigFlag', 1);
p.parse(varargin{:});
populationSig = []; % the matrix with 1 for positive beta, -1 for negative beta
populationTStats = []; % t statistics for each parameter
populationCoeffs =[]; % coeffs for each regressor
paramNames = getParamNames_dF(p.Results.modelName, 1);
maxTrial = p.Results.maxTrial;
% basic info
[root, sep] = currComputer();
% time window
time = -1000*p.Results.tb:1000*p.Results.tf;
midPoints = (0.5*p.Results.binSize + 1):p.Results.stepSize:(length(time)-0.5*p.Results.binSize);
slideTime = midPoints - p.Results.tb*1000;
dayList = getDayList(xlFile, sheet, category);
%% animal loop 
    % load model fitting results
    

%% session and unit loop
    for ses = 1:length(dayList)
        session = dayList{ses};
    % paths
        [animalName, date] = strtok(session, 'd'); 
        animalName = animalName(2:end);
        sampFile = [animalName category '_', p.Results.modelName];
        path = [root animalName sep animalName 'sorted' sep 'stan' sep 'bernoulli' sep p.Results.modelName sep category sep];
        pd = parseSessionString_df(session, root, sep);
        sortedFolderLocation = [pd.sortedFolder 'session' sep];
        pupilFile = [sortedFolderLocation session '_pupil.mat'];
    % load pupil
    if exist(pupilFile,'file')
        load(pupilFile)
    else
        fprintf(['no pupi in ' session '\n'])
        continue
    end
    % check alignment
    if errorProp > 0.4
        fprintf(['pupil no well aligned in ' session '\n']);
        continue
    end
    %% behavior preparation 
    % parse behavior
    os = behAnalysisNoPlot_opMD(session, 'simpleFlag', 0);
    choice = os.allChoices';
    choice(choice<0) = 0;
    outcome = abs(os.allRewards)';
    choice = choice(1:min(length(choice), maxTrial));
    outcome = outcome(1:length(choice));
    responseInds = os.responseInds(1:min(length(choice), maxTrial)); 
    preRwd = [NaN abs(os.allRewards(1:end-1))]';
    %% behavior
    % switch
    svs = zeros(length(os.responseInds),1);
    svs(os.changeChoice_Inds) = 1;
    svsNext = [svs(2:end); NaN];
    svsWhenNrwd = svsNext;
    svsWhenNrwd(os.rwd_Inds) = NaN;
    [t,~,noSession] = getStanModelParams_samps(p.Results.modelName, [path sampFile '.mat'], 2000, 'sessionName', session);
    if noSession
        fprintf(['no good behavior in ' session '\n']);
        continue
    end

    % diff value
    Qdiff = abs(t.Q(:,2)-t.Q(:,1));
    % total value
    Qsum = sum(t.Q,2);
    % prepe
    prePe = [NaN; t.pe(1:end-1)];
    % pe
    pe = t.pe;
    % dawExp
    dawExp = double(t.probChoice <= 0.5);
    % confidence
    choiceConf = 2.*t.probChoice - 1;
    %time in session
    timeInSession = ([os.behSessionData(responseInds).CSon]' - os.behSessionData(responseInds(1)).CSon)./1000;
    % chosen value
    Qchosen  = zeros(length(choice),1);
    Qunchosen  = zeros(length(choice),1);
    QchosenUpdate = NaN(length(choice),1);
    for j = 1:length(choice)
        if j < length(choice)
            if choice(j)>0
                Qchosen(j) = t.Q(j,2);
                Qunchosen(j) = t.Q(j,1);
                QchosenUpdate(j) = t.Q(j+2);
            else
                Qchosen(j) = t.Q(j,1);
                Qunchosen(j) = t.Q(j,2);
                QchosenUpdate(j) = t.Q(j+1);
            end
        else                
            if choice(j)>0
                Qchosen(j) = t.Q(j,2);
                Qunchosen(j) = t.Q(j,1);
            else
                Qchosen(j) = t.Q(j,1);
                Qunchosen(j) = t.Q(j,2);
            end

        end
    end
    % bias side
    biasSide = zeros(size(responseInds))';
    biasInd = contains(paramNames, 'bias');
    if mean(t.params(:,biasInd))>0
        biasSide(os.lickR_Inds)=1;
    else
        biasSide(os.lickL_Inds)=1;
    end
    hmm = double(os.hmmStates==1)';
    lickLat = os.lickLatLogZ';
    rightSide = zeros(size(pe));
    rightSide(os.allChoices>0)=1;
    preITI = os.timeBtwn';
    % consecutive no rewards
    conNrwds = zeros(size(pe));
    for j = 1:length(choice)
        if outcome(j) == 0
            k = 1;
            while j-k>0 
                if outcome(j-k)==0
                    k = k+1;
                else
                    break
                end
            end                
            conNrwds(j) = k;
        end
    end

    if contains(p.Results.modelName, '7params_absPePeAN_scale_int_bias_ord')
        aN = t.aN;
        peBar = t.peBar;
        pePe = t.pePe;
        scPe = pe.*(1-peBar);
        tbl = table(outcome, pe, prePe, preRwd, Qsum, Qdiff, choiceConf, biasSide, rightSide, timeInSession, lickLat, hmm, Qchosen, Qunchosen, QchosenUpdate, preITI, dawExp, svs, svsNext, svsWhenNrwd, conNrwds, scPe, aN, peBar, pePe);
    else
        tbl = table(outcome, pe, prePe, preRwd, Qsum, Qdiff, choiceConf, biasSide, rightSide, timeInSession, lickLat, hmm, Qchosen, Qunchosen, QchosenUpdate, preITI, dawExp, svs, svsNext, svsWhenNrwd, conNrwds);
    end
    names = tbl.Properties.VariableNames;
    % zscore all regressors
    for cols = 1:length(names)
        tmp = tbl.(names{cols});
        tmp(~isnan(tmp)) = zscore(tmp(~isnan(tmp)));
        tbl.(names{cols}) = tmp;
    end
        


    %% generate pupil matrix from choice matrix
    preLen = round(2*FR);
    postLen = round(10*FR);
    slideFrame = round(midPoints*FR/1000);
    binSizeHalf = floor(0.5 * p.Results.binSize*FR/1000);
    
    % create mean pupil size matrix for regression
    pupilMatSlide = zeros(size(sessionPupilChoice,1), length(slideFrame));
    for t = 1:length(slideFrame)
        pupilMatSlide(:,t) = sum(sessionPupilChoice(:,max(slideFrame(t)-binSizeHalf,1):min(slideFrame(t)+binSizeHalf, size(sessionPupilChoice,2))), 2,'omitnan');
    end
    clear sessionPupilChoice
    
    
    %% glm
    fprintf(['glm of ' session '\n']);
    sigs = [];
    tStats = [];
    coeffs = [];
    for k = 1:length(midPoints)
        currPupil = pupilMatSlide(:,k);
        currPupil(~qualInd(os.responseInds)) = NaN;
        currTbl = addvars(tbl, currPupil);
        lm = fitlm(currTbl, ['currPupil~' p.Results.regressors]);
        sigTmp = zeros(1,length(lm.CoefficientNames)-1);
        tStatsTmp = zeros(1,length(lm.CoefficientNames)-1);
        coeffsTmp = zeros(1,length(lm.CoefficientNames)-1);
        for j = 1:(length(lm.CoefficientNames)-1)
            if lm.Coefficients.pValue(j+1)<0.05
               sigTmp(j) = sign(lm.Coefficients.Estimate(j+1));
            else
               sigTmp(j) = 0;
            end
            tStatsTmp(j) = lm.Coefficients.tStat(j+1);
            coeffsTmp(j) = lm.Coefficients.Estimate(j+1);
        end
        sigs = [sigs; sigTmp];
        tStats = [tStats; tStatsTmp];
        coeffs = [coeffs; coeffsTmp];
        
    end

    populationSig = cat(3, populationSig, sigs);
    populationTStats  = cat(3,populationTStats, tStats);
    populationCoeffs  = cat(3,populationCoeffs, coeffs);
    end 


%% plot everything
regressors = lm.CoefficientNames(2:end);
tFig = figure;
screen = get(0,'Screensize');
screen(4) = screen(4) - 100;
set(tFig, 'Position', screen)
suptitle('tStats distribution')
colors = cool(length(regressors));
subplot(length(regressors)+1,1,1); hold on;
allSig = abs(populationSig);
allSig = mean(allSig,3);
for i = 1:length(regressors)
    plot(slideTime,allSig(:,i),'Color', colors(i,:), 'LineStyle','-', 'Marker','none', 'linewidth', 2);
end
edges = [slideTime - 0.5*p.Results.stepSize slideTime(end)+0.5*p.Results.stepSize];
for i = 1:length(edges)
    line([edges(i) edges(i)], [0 1.2*max(allSig,[],'all')], 'color', [0.7 0.7 0.7], 'LineStyle','--')
end
line([300 300], [0 1.2*max(allSig,[],'all')], 'color', 'r', 'LineStyle','--');
legend(regressors)
ylim([0 1.2*max(allSig,[],'all')])
xlim(minmax(edges));
ylabel('ratio of sig untis')
xlabel('time from respond')


minT = min(populationTStats,[],'all');
maxT = max(populationTStats,[],'all');
bins = linspace(minT,maxT,20);
for k = 1:length(midPoints)
    for j = 1:length(regressors)
        subplot(length(regressors)+1, length(midPoints), length(midPoints)*j+k); hold on;
        tmpTStats = populationTStats(k,j,:);
        tmpSig = populationSig(k,j,:);
        nonSig = tmpTStats(tmpSig == 0);
        Sig = tmpTStats(tmpSig ~= 0);
        histogram(nonSig, bins, 'FaceColor', [0.5 0.5 0.5]);
        histogram(Sig, bins, 'FaceColor', colors(j,:));
        if k == 1
            ylabel(regressors{j})
        end
    end
end
suptitle([sheet '-' category])

% for i = 1:length(midPoints)
%     titleStr = sprintf('From %d To %d', edges(i), edges(i+1)); 
%     figure;
%     scatterAll(squeeze(populationCoeffs(i,:,:))', regressors,7,'m');
%     suptitle(titleStr)
% end  












