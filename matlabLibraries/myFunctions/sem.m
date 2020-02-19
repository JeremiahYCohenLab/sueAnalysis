function sem_error = sem(input)
% sem   Generate standard error of the mean
%       sem_error = sem(input)
% INPUTS
%       input: vector of data
% OUTPUTS
%       sem_error: standard error of the mean for corresponding input

    sem_error = nanstd(input)./sqrt(sum(~isnan(input)));
end