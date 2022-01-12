function sem = sem_bern(events)
% sem_bernoulli   Generate standard error for bernoulli distributed variables
%       sem = sem_bern(events)
% INPUTS
%       events, 1 for success, 0 for fail
% OUTPUTS
%       sem: standard error of the mean

p = sum(events)./length(events);
q = 1 - p;

sem = sqrt((p.*q)./length(events));