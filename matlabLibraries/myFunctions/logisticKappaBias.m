function output = logisticKappaBias(input, kappa, bias)
%% Simple logistic function
% offset is how far away from 0 and 1 the curves should saturate

output = 1 ./ (1 + exp(-1*input - kappa - bias));