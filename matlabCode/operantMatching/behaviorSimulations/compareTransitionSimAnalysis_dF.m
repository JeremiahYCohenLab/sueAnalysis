function [t] = compareTransitionSimAnalysis_dF(xlFile, sheet, pre, post, varargin)

p = inputParser;
% default parameters if none given
p.addParameter('runs', 4)
p.addParameter('samps', 1000)
p.addParameter('modelName', {'sevenParam_absPePeAN_scale_int_bias_ord'})
p.addParameter('bernFlag', 1)
p.addParameter('sessionParamsFlag', 0)
p.addParameter('lesionInd', [])
p.addParameter('lesionVal', 0)
p.parse(varargin{:});

[root, sep] = currComputer();
thresh = 0.4;

if length(p.Results.modelName) == 1
    preMdl = p.Results.modelName{1};
    postMdl = p.Results.modelName{1};
else
    preMdl = p.Results.modelName{1};
    postMdl = p.Results.modelName{2};
end

tPre = analyzeTransitionSimSamp_dF(xlFile, sheet, pre, pre, 'modelName', preMdl, 'runs', p.Results.runs,...
                                    'samps', p.Results.samps, 'sessionParamsFlag', p.Results.sessionParamsFlag);
if ~isempty(p.Results.lesionInd)
    %get simulated data from pre + lesion
    tPost = analyzeTransitionSimSamp_dF(xlFile, sheet, pre, pre, 'modelName', postMdl, 'runs',...
                                    p.Results.runs, 'samps', p.Results.samps, 'lesionInd', p.Results.lesionInd,...
                                    'sessionParamsFlag', p.Results.sessionParamsFlag);

    %get actual data from post behavior
    [~, transMed, transHigh, ~, ~, ~, ~] = transitionAnalysis_opMD(xlFile, sheet, post);
    range = 15;
    postTranWin = 15;
    medAvg = mean(transMed(:,range:range+postTranWin));
    highAvg = mean(transHigh(:,range:range+postTranWin));
    tPost.medAvg = medAvg;
    tPost.highAvg = highAvg;
    
    %find zero crossing
    tmpDiff = [highAvg - medAvg];
    zeroInd = find(tmpDiff < 0, 1);
    if zeroInd == 1
        zeroInd = find(tmpDiff < 0, 2);
        zeroInd = zeroInd(2);
    end
    x = [zeroInd-1 zeroInd];                 
    y = tmpDiff(x);
    c = [[1; 1]  x(:)]\y(:);
    zeroCross = (0 - c(1)) / c(2);
    tPost.zeroCross = zeroCross;

    %find thresh cross
    threshInd = find(highAvg < thresh, 1);
    if isempty(threshInd)
        threshCrossHigh = NaN;
    else
        x = [threshInd-1 threshInd];                        %find line between points on sides of zero
        y = highAvg(x);
        c = [[1; 1]  x(:)]\y(:);                            %calculate parameter vector
        threshCrossHigh = (thresh - c(1)) / c(2);    %x = (y-b)/ m
    end
    threshInd = find(medAvg < thresh, 1);
    if isempty(threshInd)
        threshCrossLow = NaN;
    else
        x = [threshInd-1 threshInd];                       %find line between points on sides of zero
        y = medAvg(x);
        c = [[1; 1]  x(:)]\y(:);                           %calculate parameter vector
        threshCrossLow = (thresh - c(1)) / c(2);    %x = (y-b)/ m
    end
    tPost.threshCrossHigh = threshCrossHigh;
    tPost.threshCrossLow = threshCrossLow;

    %get tau for choice probability curves
    ft = fittype('a*exp((-1/b)*x)');
    sp = [0.5 7];
    x = [0:postTranWin];
    tau = nan(1,2); tauCI = nan(1,2);
    medMdl = fit(x', medAvg', ft, 'start', sp);
    tau(1) = medMdl.b;
    tmp = confint(medMdl); 
    tauCI(1) = medMdl.b - tmp(1,2);
    highMdl = fit(x', highAvg', ft, 'start', sp);
    tau(2) = highMdl.b;
    tmp = confint(highMdl); 
    tauCI(2) = highMdl.b - tmp(1,2);
    tPost.tau = tau;
    tPost.tauCI = tauCI;
else
    tPost = analyzeTransitionSimSamp_dF(xlFile, sheet, post, post, 'modelName', postMdl, 'runs',...
                                    p.Results.runs, 'samps', p.Results.samps, ...
                                    'sessionParamsFlag', p.Results.sessionParamsFlag);
end

if p.Results.sessionParamsFlag
    names = {'zeroCrossSim', 'threshCrossLowSim', 'threshCrossHighSim', 'tauSim'};
    behNames = {pre, post};
  
    for currB = 1:2
        [~, dayList, ~] = xlsread(xlFile, sheet);
        [~,col] = find(strcmp(dayList, behNames{currB}));
        dayList = dayList(2:end, col);
        endInd = find(cellfun(@isempty,dayList),1);
        if ~isempty(endInd)
            dayList = dayList(1:endInd-1,:);
        end
        
        prevAnimal = []; aInd = 0; aInds = [];
        for currS = 1:length(dayList)
            [animal, ~] = strtok(dayList{currS}, 'd'); 
            animal = animal(2:end);
            if strcmp(animal, prevAnimal) == 0
                aInd = aInd+1;
            end
            aInds = [aInds aInd];
            prevAnimal = animal;
        end
        
        for currV = 1:length(names)
            tmp = [];
            for currA = 1:max(aInds)
                if currB == 1
                    if regexp(names{currV}, 'tauSim')
                        tmp(currA,:) = nanmean(tPre.tauSim(aInds == currA, :));
                    else
                        tmp(currA) = nanmean(tPre.(names{currV})(aInds == currA));
                    end
                else
                    if regexp(names{currV}, 'tauSim')
                        tmp(currA,:) = nanmean(tPost.tauSim(aInds == currA, :));
                    else
                        tmp(currA) = nanmean(tPost.(names{currV})(aInds == currA));
                    end
                end
            end
            if currB == 1
                tPre.(names{currV}) = tmp;
            else
                tPost.(names{currV}) = tmp;
            end
        end
    end
end


colors = cool(2);
figure; 
subplot(1,4,1); hold on;
binEdges = [1:10];
histogram(tPre.zeroCrossSim, binEdges, 'FaceColor', colors(1,:), 'Normalization', 'probability')
histogram(tPost.zeroCrossSim, binEdges, 'FaceColor', colors(2,:), 'Normalization', 'probability')
ylabel('probability')
xlabel('trial at 0 cross')
numNoCross(1) = sum(isnan(tPre.zeroCrossSim))/length(tPre.zeroCrossSim);
numNoCross(2) = sum(isnan(tPost.zeroCrossSim))/length(tPost.zeroCrossSim);
legend([{strcat(num2str(numNoCross(1)*100), '% never cross')} {strcat(num2str(numNoCross(2)*100), '% never cross')}])

yl = ylim;
plot([tPre.zeroCross tPre.zeroCross], [0 yl(2)], '--', 'color', colors(1,:), 'linewidth', 2)
plot([tPost.zeroCross tPost.zeroCross], [0 yl(2)], '--', 'color', colors(2,:), 'linewidth', 2)
set(gca, 'tickdir', 'out', 'box', 'off')

%fix this to deal with session parameters
subplot(1,4,2); hold on;
numA = length(tPre.threshCrossLowSim);
for currA = 1:numA
    plot([tPre.threshCrossLowSim(currA) tPost.threshCrossLowSim(currA)], ...
        [tPre.threshCrossHighSim(currA) tPost.threshCrossHighSim(currA)], 'k');
end 
scatter(tPre.threshCrossLowSim, tPre.threshCrossHighSim, [], colors(1,:), 'filled')
scatter(tPost.threshCrossLowSim, tPost.threshCrossHighSim, [], colors(2,:), 'filled')
plot([tPre.threshCrossLow tPost.threshCrossLow], [tPre.threshCrossHigh tPost.threshCrossHigh], '-k')
scatter(tPre.threshCrossLow, tPre.threshCrossHigh, 100, 'markerfacecolor', colors(1,:), 'markeredgecolor', 'k');
scatter(tPost.threshCrossLow, tPost.threshCrossHigh, 100, 'markerfacecolor', colors(2,:), 'markeredgecolor', 'k');
xl = xlim; yl = ylim;
plot([0 max([xl yl])], [0 max([xl yl])], '--k') 
ylim([0 max([xl yl])]);
xlim([0 max([xl yl])]);
xlabel('thresh cross medium')
ylabel('thresh cross high')
set(gca, 'tickdir', 'out')

numNoCross(1) = sum(isnan(tPre.threshCrossLowSim) | isnan(tPre.threshCrossHighSim)) / length(tPre.threshCrossLowSim);
numNoCross(2) = sum(isnan(tPost.threshCrossLowSim) | isnan(tPost.threshCrossHighSim)) / length(tPost.threshCrossLowSim);
legTxt = [{['thresh = ' num2str(thresh)]}, {[num2str(numNoCross(1)*100) '% never cross - pre']}, ...
                {[num2str(numNoCross(2)*100) '% never cross - post']}];
legend(legTxt);

subplot(1,4,3); hold on;
for currA = 1:numA
    plot([tPre.tauSim(currA,1) tPost.tauSim(currA,1)], [tPre.tauSim(currA,2) tPost.tauSim(currA,2)], 'k');
end
scatter(tPre.tauSim(:,1), tPre.tauSim(:,2), [], colors(1,:), 'filled')
scatter(tPost.tauSim(:,1), tPost.tauSim(:,2), [], colors(2,:), 'filled')

plot([tPre.tau(1) tPost.tau(1)], [tPre.tau(2) tPost.tau(2)], 'k');
errorbar(tPre.tau(1), tPre.tau(2), tPre.tauCI(2), tPre.tauCI(2), tPre.tauCI(1), tPre.tauCI(1), 'o',...
    'color','k', 'markerfacecolor', colors(1,:), 'markersize', 10)
errorbar(tPost.tau(1), tPost.tau(2), tPost.tauCI(2), tPost.tauCI(2), tPost.tauCI(1), tPost.tauCI(1), 'o',...
    'color','k', 'markerfacecolor', colors(2,:), 'markersize', 10)

xl = xlim; yl = ylim;
plot([0 max([xl yl])], [0 max([xl yl])], '--k') 
ylim([0 max([xl yl])]);
xlim([0 max([xl yl])]);
xlabel('medium \tau')
ylabel('high \tau')
set(gca, 'tickdir', 'out', 'box', 'off')

subplot(1,4,4); hold on;
plot([tPre.tauSimAvg(1) tPost.tauSimAvg(1)], [tPre.tauSimAvg(2) tPost.tauSimAvg(2)], 'k');
errorbar(tPre.tauSimAvg(1), tPre.tauSimAvg(2), tPre.tauSimAvgCI(2), tPre.tauSimAvgCI(2), tPre.tauSimAvgCI(1), tPre.tauSimAvgCI(1), 'o',...
    'color','k', 'markerfacecolor', colors(1,:), 'markersize', 10)
errorbar(tPost.tauSimAvg(1), tPost.tauSimAvg(2), tPost.tauSimAvgCI(2), tPost.tauSimAvgCI(2), tPost.tauSimAvgCI(1), tPost.tauSimAvgCI(1), 'o',...
    'color','k', 'markerfacecolor', colors(2,:), 'markersize', 10)
plot([tPre.tau(1) tPost.tau(1)], [tPre.tau(2) tPost.tau(2)], 'k');
errorbar(tPre.tau(1), tPre.tau(2), tPre.tauCI(2), tPre.tauCI(2), tPre.tauCI(1), tPre.tauCI(1), 'o',...
    'color','k', 'markerfacecolor', colors(1,:), 'markersize', 10)
errorbar(tPost.tau(1), tPost.tau(2), tPost.tauCI(2), tPost.tauCI(2), tPost.tauCI(1), tPost.tauCI(1), 'o',...
    'color','k', 'markerfacecolor', colors(2,:), 'markersize', 10)
set(gca, 'tickdir', 'out', 'box', 'off')


xl = xlim; yl = ylim;
plot([0 max([xl yl])], [0 max([xl yl])], '--k') 
ylim([0 max([xl yl])]);
xlim([0 max([xl yl])]);
xlabel('medium \tau')
ylabel('high \tau')

set(gcf, 'renderer', 'painters', 'position', [-1720 377 1648 420])


%% glm
% allTrans = [tPre.medAvg  reshape(tPre.medAvgSim', 1, size(tPre.medAvgSim,1)*size(tPre.medAvgSim,2))...
%             tPre.highAvg reshape(tPre.highAvgSim', 1, size(tPre.highAvgSim,1)*size(tPre.highAvgSim,2))...
%             tPost.medAvg  reshape(tPost.medAvgSim', 1, size(tPost.medAvgSim,1)*size(tPost.medAvgSim,2))...
%             tPost.highAvg reshape(tPost.highAvgSim', 1, size(tPost.highAvgSim,1)*size(tPost.highAvgSim,2))];
% 
% numTT = length([tPre.medAvg  reshape(tPre.medAvgSim', 1, size(tPre.medAvgSim,1)*size(tPre.medAvgSim,2))]);
% transType = repmat([ones(1, numTT) ones(1, numTT)*2], 1, 2);
% 
% numC = numTT*2;
% cond = [ones(1,numC) ones(1,numC)*2];
% 
% type = repmat([ones(1,range+1) ones(1, length(reshape(tPre.medAvgSim', 1, size(tPre.medAvgSim,1)*size(tPre.medAvgSim,2))))*2], 1, 4);
% 
% trial = repmat([1:range+1], 1, length(allTrans)/(range+1));
% mouse = repmat([nan(1,range+1) ones(1,range+1) ones(1,range+1)*2 ones(1,range+1)*3], 1,4);
% 
% mdl = fitglm([trial' transType' cond' type' mouse'], allTrans', 'distribution', 'gamma', 'categoricalvars', [2 3 4 5]);


%% lme by taus

% t = table(zscore([tPre.tauSim(:,1); tPre.tauSim(:,2); tPost.tauSim(:,1); tPost.tauSim(:,2)]),  repmat([1:numA]', 4, 1), ...
%             repmat([ones(numA, 1); ones(numA, 1)*2], 2, 1), [ones(numA*2,1); ones(numA*2,1)*2], ...
%             'VariableNames', {'taus', 'mouse', 'trans', 'cond'});
% t = table(zscore([tPre.tauSim(:,1); tPre.tauSim(:,2); tPost.tauSim(:,1); tPost.tauSim(:,2)]),  repmat([1:numA]', 4, 1), ...
%         repmat([repmat('m', numA, 1); repmat('h', numA, 1)], 2, 1), [repmat('vh', numA*2,1); repmat('ag', numA*2,1)], ...
%         'VariableNames', {'taus', 'mouse', 'trans', 'cond'});
% t.mouse = nominal(t.mouse);
% t.trans = nominal(t.trans);
% t.cond = nominal(t.cond);
% mdl = fitlme(t, 'taus~trans*cond+(trans*cond|mouse)')
% stats = anova(mdl)


% t = table(zscore([tPre.tauSim(:,1); tPre.tau(1,1); tPre.tauSim(:,2); tPre.tau(1,2);...
%     tPost.tauSim(:,1); tPost.tauSim(:,2); tPost.tau(1,1); tPost.tau(1,2);]),  repmat([1:numA NaN]', 4, 1), ...
%         repmat([repmat('m', numA+1, 1); repmat('h', numA+1, 1)], 2, 1), [repmat('vh', (numA+1)*2,1); repmat('ag', (numA+1)*2,1)], ...
%         repmat([repmat('s', numA, 1); repmat('a', 1, 1)], 4,1), ...
%         'VariableNames', {'taus', 'mouse', 'trans', 'cond', 'type'});
% t = table(zscore([tPre.tauSim(:,1); tPre.tau(1,1); tPre.tauSim(:,2); tPre.tau(1,2);...
%     tPost.tauSim(:,1); tPost.tauSim(:,2); tPost.tau(1,1); tPost.tau(1,2);]),  repmat([1:numA+1]', 4, 1), ...
%         repmat([repmat('m', numA+1, 1); repmat('h', numA+1, 1)], 2, 1), [repmat('vh', (numA+1)*2,1); repmat('ag', (numA+1)*2,1)], ...
%         'VariableNames', {'taus', 'mouse', 'trans', 'cond'});
% t.mouse = nominal(t.mouse);
% t.trans = nominal(t.trans);
% t.cond = nominal(t.cond);
% t.type = nominal(t.type);

% mdl = fitlme(t, 'taus~trans*cond+(trans*cond|mouse)');
% stats = anova(mdl);
            
         

%% lme on trial-by-trial choice probabilities

allTrans = [reshape(tPre.medAvgSim', 1, size(tPre.medAvgSim,1)*size(tPre.medAvgSim,2))...
            reshape(tPre.highAvgSim', 1, size(tPre.highAvgSim,1)*size(tPre.highAvgSim,2))...
            reshape(tPost.medAvgSim', 1, size(tPost.medAvgSim,1)*size(tPost.medAvgSim,2))...
            reshape(tPost.highAvgSim', 1, size(tPost.highAvgSim,1)*size(tPost.highAvgSim,2))];

numTT = length([reshape(tPre.medAvgSim', 1, size(tPre.medAvgSim,1)*size(tPre.medAvgSim,2))]);
transType = repmat([ones(1, numTT) ones(1, numTT)*2], 1, 2);
numC = numTT*2;
cond = [ones(1,numC) ones(1,numC)*2];
trial = repmat([1:range+1], 1, length(allTrans)/(range+1));
mouse = [];
for currA = 1:numA          %fix this for session params
    mouse = [mouse ones(1,range+1) * currA];
end
mouse = repmat(mouse, 1,4);

tt = table(zscore(allTrans)', trial', transType', cond', mouse', ...
        'VariableNames', {'choiceProbs', 'trial', 'trans', 'cond', 'mouse'});
tt.mouse = nominal(tt.mouse);
tt.trans = nominal(tt.trans);
tt.cond = nominal(tt.cond);

%mdl = fitlme(tt, 'choiceProbs~trial*trans*cond+(trial*trans*cond|mouse)');
mdl = fitlme(tt, 'choiceProbs~trial*trans + trans*cond + (trial*trans|mouse) + (trans*cond|mouse)');
stats = anova(mdl)



allTrans = [tPre.medAvg tPre.highAvg tPost.medAvg tPost.highAvg];

numTT = length([tPre.medAvg]);
transType = repmat([ones(1, numTT) ones(1, numTT)*2], 1, 2);
numC = numTT*2;
cond = [ones(1,numC) ones(1,numC)*2];
trial = repmat([1:range+1], 1, length(allTrans)/(range+1));

tt = table(zscore(allTrans)', trial', transType', cond', ...
        'VariableNames', {'choiceProbs', 'trial', 'trans', 'cond'});
tt.trans = nominal(tt.trans);
tt.cond = nominal(tt.cond);

%mdl = fitlme(tt, 'choiceProbs~trial*trans*cond')
mdl = fitlme(tt, 'choiceProbs~trial*trans + trans*cond');
stats = anova(mdl)


