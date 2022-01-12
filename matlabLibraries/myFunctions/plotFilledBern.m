function hOut = plotFilledBern(x, data, color, figHandle)
% plotFilled     Plots mean +/- SEM or std as a filled plot
%   hOut = plotFilled(x, data, color, figHandle)
%   INPUTS
%       x: x data
%       data: y data as a n x x vector (n: independent observations, x: x coordinate data)
%       color: 1x3 vector from [0 0 0] to [1 1 1]
%       figHandle: optional input
%   OUTPUTS
%       hOut: handle to plot command

if nargin < 2 || isempty(x)
    data = x;
    x = [1:size(data,2)];
end
if nargin < 3
    color = 'b';
end

if size(data,1) > 1
    y = nanmean(data);
else
    y = data;
end

tmp = sem_bernoulli(sum(data==1), sum(~isnan(data)));
yU = y + tmp;
yL = y - tmp;
yU(isnan(yU))= 0; 
yL(isnan(yL))= 0; 

if nargin > 4
    subplot(figHandle)
end
hOut = plot(x, y, 'Color', color, 'linewidth', 2); hold on;
fill([x fliplr(x)], [yU fliplr(yL)], color, 'facealpha', 0.25, 'edgecolor', 'none');