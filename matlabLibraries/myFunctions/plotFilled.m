function hOut = plotFilled(x, data, color, lineWidth, stdFlag, figHandle)
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
    lineWidth = 1;
else
    if nargin < 5
        stdFlag = 0;
    end
end


y = nanmean(data);

if stdFlag
    yU = y + nanstd(data);
    yL = y - nanstd(data);
else
    yU = y + nanstd(data)/sqrt(size(data,1));
    yL = y - nanstd(data)/sqrt(size(data,1));
end

if nargin > 5
    subplot(figHandle)
end
hOut = plot(x, y, 'Color', color, 'linewidth', lineWidth, 'Marker', 'none', 'LineStyle', '-'); hold on;
fill([x(~isnan(yU)) fliplr(x(~isnan(yU)))], [yU(~isnan(yU)) fliplr(yL(~isnan(yU)))], color, 'facealpha', 0.25, 'edgecolor', 'none')