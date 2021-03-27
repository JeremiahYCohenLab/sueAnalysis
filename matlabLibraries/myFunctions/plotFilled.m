function hOut = plotFilled(x, data, color, stdFlag, figHandle)
% plotFilled     Plots mean +/- SEM or std as a filled plot
%   hOut = plotFilled(x, data, color, figHandle)
%   INPUTS
%       x: x data
%       data: y data as a n x x vector (n: independent observations, x: x coordinate data)
%       color: 1x3 vector from [0 0 0] to [1 1 1]
%       figHandle: optional input
%   OUTPUTS
%       hOut: handle to plot command

if nargin < 4
    stdFlag = 0;
end

y = nanmean(data);

if stdFlag
    yU = y + nanstd(data);
    yL = y - nanstd(data);
else
    yU = y + sem(data);
    yL = y - sem(data);
end

if nargin > 4
    subplot(figHandle)
end
hOut = plot(x, y, 'color', color, 'linewidth', 1, 'Marker', 'none', 'LineStyle', '-'); hold on;
fill([x fliplr(x)], [yU fliplr(yL)], color, 'facealpha', 0.25, 'edgecolor', 'none')