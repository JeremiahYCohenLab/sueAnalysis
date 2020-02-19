function sem = sem_bernoulli(success, total)
% sem_bernoulli   Generate standard error for bernoulli distributed variables
%       sem = sem_bernoulli(success, total)
% INPUTS
%       success: number of successes
%       total: total samples
% OUTPUTS
%       sem: standard error of the mean

p = success./total;
q = 1 - p;

sem = sqrt((p.*q)./total);