function [ratioMax, x, p] = timeProjectionOpti2(ledLL, csT, startFrame, startTrial)
fun = @(x) -1*timeProjection(ledLL, csT, startFrame, startTrial, x);
startPoints = 20.5:0.1:21.5;
runs = length(startPoints);
A = [1, -1]';
b = [21.5, -20.5]';

x = zeros(size(startPoints));
p = zeros(size(startPoints));
options = optimoptions(@fmincon,'FiniteDifferenceStepSize',1e-3,'Display','off');


    parfor r = 1:runs
         [x(r), p(r)] = fmincon(fun, startPoints(r), A, b, [],[],[],[],[], options);
    end
    
    [~,B] = min(p);
    ratioMax = x(B);
    %%