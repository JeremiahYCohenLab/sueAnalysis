animalName = 'ZS061';
category = 'good';
paramNames = {'lambda', 'vInit', 'omega', 'beta', 'kappa'}; % animal level
model = 'vkf_fixV_kappa';

[root, sep] = currComputer();

if contains(model,'tF')
    input = 'choice, outcome, ITI)';
else
    input = 'choice, outcome)';
end
    
%% load model fitting results
sampFile = [animalName category '_', model];
path = [root animalName sep animalName 'sorted' sep 'stan' sep 'bernoulli' sep model sep category sep];
% path = [root animalName sep 'stan' sep 'bernoulli' sep model sep];
if exist([path sampFile], 'file')
    load([path sampFile '.mat'], 'dayList');
    samples = load([path sampFile '.mat'], sampFile);
    samples = samples.(sampFile);
else
   samples = struct;
   fieldNames = [paramNames strcat('mu_', paramNames)];
   fieldNames = [fieldNames 'log_lik'];
   fieldNames = [fieldNames 'sigma'];
   fieldNames = [fieldNames 'divergent_'];
   if contains(model,'bias')
       fieldNames = [fieldNames 'bias'];   
   end
   for f = 1:length(fieldNames)
       samples.(fieldNames{f}) = [];
   end
    % read cvs files
   modelFiles = dir(path);
   outputExpression = ['^' 'output'];
   outputFiles = {modelFiles(~cellfun(@isempty, cellfun(@(x) regexp(x, outputExpression), {modelFiles.name}, 'UniformOutput', false))).name};
   for i = 1:length(outputFiles)
       currOutput = outputFiles{i};
       [nums, texts, all] = xlsread([path currOutput], currOutput(1:end-4));
        for f = 1:length(fieldNames)
            if f <= length(paramNames)
               tmpInd = find(contains(texts(39,:), [fieldNames{f}, '.']));
            else
               tmpInd = find(contains(texts(39,:), fieldNames{f}));
            end

            tmp = nums(2:end, tmpInd);
            samples.(fieldNames{f}) = [[samples.(fieldNames{f})]; tmp];
        end
   end

end

%%
pFig = figure2;
numParams = length(paramNames);
blue = [0 1 1];
purp = [0.7 0 1];
colors = cool(numParams);
for i = 1:numParams
    subplot(1,numParams,i); hold on;
    histogram(eval(['samples.mu_' paramNames{i}]) , 100,...
        'Normalization', 'Probability', 'FaceColor', colors(i,:), 'EdgeColor', 'none')
    set(gca,'tickdir', 'out') 
    title(paramNames{i})
end

%% save
saveFile = [sampFile '.mat'];
eval([sampFile,  ' = samples;']);
save([path saveFile], sampFile);
saveFigurePDF(pFig,[path model  '_posteriors.pdf'])
%%
