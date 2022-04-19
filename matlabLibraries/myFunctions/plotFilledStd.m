function hOut = plotFilledStd(x, data, color, figHandle)
% plotFilled     Plots mean +/- SEM as a filled plot
%   hOut = plotFilled(x, data, color, figHandle)
%   INPUTS
%       x: x data
%       data: y data as a n x x vector (n: independent observations, x: x coordinate data)
%       color: 1x3 vector from [0 0 0] to [1 1 1]
%       figHandle: optional input
%   OUTPUTS
%       hOut: handle to plot command

y = nanmean(data);
yU = y + std(data,'omitnan');
yL = y - std(data,'omitnan');

if nargin < 3
    subplot(figHandle)
end
hOut = plot(x, y, 'Color', color, 'linewidth', 2);
fill([x fliplr(x)], [yU fliplr(yL)], color, 'facealpha', 0.25, 'edgecolor', 'none')