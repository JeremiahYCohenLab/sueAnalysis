
%% example 1
session = 'mZS061d20210424';
optoUnit = 'TT4_SS_01';
subFolder = '10ms';
sessionUnit = optoUnit;
%% example 2
session = 'mZS061d20210330';
optoUnit = 'TT4_SS_01';
sessionUnit = optoUnit;
subFolder = '20ms';
%% example 3
session = 'mZS061d20210416';
optoUnit = 'TT4_SS_01';
sessionUnit =  'TT4_SS_02';
subFolder = '10ms';

%% example increase
session = 'mZS061d20210406';
optoUnit = 'TT1_SS_01';
sessionUnit =  'TT1_SS_01';
subFolder = '20ms_1';



%%
optoIDTT(session, 'subFolder', subFolder, 'PulseWidth', str2double(subFolder(1:2)));
%%
met = getClusterMetric(session, sessionUnit, 0, 0);
%%
figure2;
plotFilled(1:101, squeeze(AD2uV*lightWaveForm(:,1,:))', 'b')
hold on;
plotFilled(1:101, squeeze(AD2uV*spontWaveForm(:,1,:))', 'k')
plot([25 41], [50 50], 'LineWidth', 2, 'Color', 'k');
plot([25 25], [50 100], 'LineWidth', 2, 'Color', 'k');
text(33, 40, '0.5 ms', 'HorizontalAlignment', 'center')
text(20, 75, '50 uV', 'HorizontalAlignment', 'center')

set(gca, 'TickDir', 'out');
set(gca, 'box', 'off')
set(gca, 'XColor', 'none')
set(gca, 'YColor', 'none')

xlim([20 90])
%%