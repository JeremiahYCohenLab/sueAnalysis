function [fit, samples] = stan_qLearningFit_session(sessionName, varargin)

p = inputParser;
% default parameters if none given
p.addParameter('revForFlag', 0);
p.addParameter('fixedParams', 0);
p.addParameter('params', ['fiveParamOfourStart_preS']);
p.addParameter('paramInds', [1 2 3 4 5]);
p.addParameter('paramNames', {'aN', 'aP', 'aF', 'beta', 'v'});
p.addParameter('animalName', []);
p.addParameter('modelType', ['fourParam']);
p.parse(varargin{:});


[root, sep] = currComputer();
filename = [sessionName '.asc'];
[behSessionData, unCorrectedBlockSwitch, out] = loadBehavioralData(filename, p.Results.revForFlag);
behavStruct = parseBehavioralData(behSessionData, unCorrectedBlockSwitch);

choice = behavStruct.allChoices;
choice(choice == -1) = 2;
choice = choice';
outcome = abs(behavStruct.allRewards); 
outcome = outcome';
T = length(outcome); 


if p.Results.fixedParams
    [allParams, modelNames, ~] = xlsread('stanParams.xlsx', animal);
    [~,col] = find(~cellfun(@isempty,strfind(modelNames, p.Results.params)) == 1);
    params = allParams(p.Results.paramInds,col);
    session_dat = struct('T',T,'choice', choice', 'outcome', outcome', 'params', params);
else
    session_dat = struct('T',T,'choice', choice', 'outcome', outcome');
end

filePath = ['C:\Users\cooper_PC\Desktop\githubRepositories\cooperAnalysis\matlabCode\operantMatching\learningModels\stan\'];
switch p.Results.modelType
    case 'fiveParamO_session'
        fit = stan('file',[filePath 'stan_qLearning_5params_opponency_session.stan'],'data',session_dat,'verbose', true);
    otherwise
        fprintf([p.Results.modelType ' model does not exist'])
end

disp('pres any key to continue once the stan processing has been completed')
pause;


samples = fit.extract('permuted',true);

paramEsts = [];
if p.Results.fixedParams
    [binNum, edges]  = discretize(eval(['samples.' p.Results.paramNames{setdiff([1:5], p.Results.paramInds)}]),100);
    binInd = mode(binNum);
    paramEsts = mean([edges(binInd) edges(binInd + 1)]);
else
    for i = 1:length(p.Results.paramInds)
        [binNum, edges]  = discretize(eval(['samples.' p.Results.paramNames{p.Results.paramInds(i)}]),100);
        binInd = mode(binNum);
        paramEsts(i) = mean([edges(binInd) edges(binInd + 1)]);
    end
end


% if p.Results.fixedParams
%     sampFile = [animal category, '_', p.Results.modelType, '_', p.Results.paramNames{setdiff([1:5], p.Results.paramInds)}];
%     saveFile = [sampFile '.mat'];
%     eval([sampFile,  ' = samples;']);
% else
%     sampFile = [animal category '_', p.Results.modelType];
%     saveFile = [sampFile '.mat'];
%     eval([sampFile,  ' = samples;']);
% end
% 
% savePath = [root animal sep animal 'sorted' sep 'stan' sep];
% if ~exist(savePath)
%     mkdir(savePath);
% end
% save([savePath saveFile], sampFile, 'paramEsts', 'dayList')

figure; 
%title([animal ' - ' p.Results.modelType], 'Interpreter', 'none');
if p.Results.fixedParams
    histogram(eval(['samples.' p.Results.paramNames{setdiff([1:5], p.Results.paramInds)}]), 100,...
            'Normalization', 'Probability', 'FaceColor', 'k')
        set(gca,'tickdir', 'out')
        xlabel( p.Results.paramNames{setdiff([1:5], p.Results.paramInds)})
else
    numParams = length(p.Results.paramInds);
    for i = 1:numParams
        subplot(1,numParams,i); hold on;
        histogram(eval(['samples.' p.Results.paramNames{p.Results.paramInds(i)}]) , 100,...
            'Normalization', 'Probability', 'FaceColor', 'k')
        set(gca,'tickdir', 'out')
        title(p.Results.paramNames{p.Results.paramInds(i)})
    end
end
set(gcf,'Renderer', 'Painters')
    
    
    
