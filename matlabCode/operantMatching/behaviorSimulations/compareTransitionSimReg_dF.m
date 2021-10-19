function compareTransitionSimReg_dF(animal, category, beh, varargin)

p = inputParser;
% default parameters if none given
p.addParameter('modelNames', {'fiveParam_bias', 'sevenParam_absPePeAN_scale_int_bias_ord'})
p.parse(varargin{:});

[root, sep] = currComputer();

%get info about models
modelNames = p.Results.modelNames;
numMdls = length(modelNames);

%set params for number of trials around transition
range = 15;
x = [-range+1:range];

%load transition data structures
for currM = 1:numMdls
    dataPath = [root 'transitionData' sep modelNames{currM} sep animal sep...
                        'transitionMdl' sep animal beh '_transitionSim.mat'];
    dS.(modelNames{currM}) = load(dataPath);
end
        
%extract choice probability curves for simulated data
for currM = 1:numMdls
    medMdl = fitlm(dS.(modelNames{currM}).t.medAvgSim, dS.(modelNames{currM}).t.medAvg);
    medR(currM) = medMdl.Rsquared.Ordinary;
    highMdl = fitlm(dS.(modelNames{currM}).t.highAvgSim, dS.(modelNames{currM}).t.highAvg);
    highR(currM) = highMdl.Rsquared.Ordinary;
end


%plot simulated curves against actual
colors = cool(numMdls+1);
figure;
subplot(1,2,1); hold on;
plotFilledBern(x, dS.(modelNames{1}).t.transMed, colors(end,:));
for currM = 1:numMdls
    plot(x, medCurves(currM, :), '-', 'color', colors(currM,:), 'linewidth', 1.5);
    text('Position', [-range 0.6-(0.1*currM)  0], 'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom', 'string', ...
        sprintf([[strrep(modelNames{currM}, '_', ' ') ' R^2 = ' num2str(medR(currM))]));
end
set(gca, 'tickdir', 'out')
title('medium to low')

subplot(1,2,2); hold on;
legTxt = {'actual', ' '};
plotFilledBern(x, dS.(modelNames{1}).t.transHigh, colors(end,:));
for currM = 1:numMdls
    plot(x, highCurves(currM, :), '-', 'color', colors(currM,:), 'linewidth', 1.5);
    text('Position', [-range 0.6-(0.1*currM)  0], 'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom', 'string', ...
        sprintf([strrep(modelNames{currM}, '_', ' ') ' R^2 = ' num2str(highR(currM))]));
    legTxt = [legTxt {strrep(modelNames{currM}, '_', ' ')}];
end
legend(legTxt);
set(gca, 'tickdir', 'out')
title('high to low')


set(gcf, 'renderer', 'painters', 'position', [-1743 392 1551 420])
    
    
    
    
    
    
    
    