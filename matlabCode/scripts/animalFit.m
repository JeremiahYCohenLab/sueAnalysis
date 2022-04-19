animals = {'ZS066', 'ZS068', 'ZS069', 'ZS070', 'ZS071'};
model = '5params_2LR_k_bias';
%%
for i = 1:length(animals)-1
    stan_qLearningFit('inhibitionAll', animals{i}, 'inhibitionNrwd', 'modelName', model, 'iter', 10000);
    stan_qLearningFit('inhibitionAll', animals{i}, 'controlNrwd', 'modelName', model, 'iter', 10000);
end
%%
modelList = {'5params', '5params_k_bias', '5params_2LR_k_bias'};
for i = 1:length(animals)
    for j = 1:length(modelList)
        plotStanSessionParams(animals(i), 'beh', 'inhibitionNrwd', 'modelName', modelList{j});
    end
    for j = 1:length(modelList)
        plotStanSessionParams(animals(i), 'beh', 'controlNrwd', 'modelName', modelList{j});       
    end
end
%%