function plotWaveforms(session, varargin)
p = inputParser;
% default parameters if none given
p.addParameter('subFolder', '')
p.addParameter('SamplingFreq', 32000);
p.parse(varargin{:});


%get session info
[root, sep] = currComputer();
pd = parseSessionString_df(session, root, sep);
sortedPath = [pd.nLynxFolder 'session' sep p.Results.subFolder sep];

%get sorted and raw ephys data
% optoFiles = dir(fullfile(dataPath,'*.rhd'));
sortedFiles = dir(fullfile(sortedPath,'*TT*.txt'));

% load(fullfile(sortedPath, 'settings.mat'), 'scaling_for_int16', 'frequency_parameters');
sampFreq = p.Results.SamplingFreq;

header = Nlx2MatCSC([sortedPath 'CSC1.ncs'], [0 0 0 0 0], 1, 1, []);

AD2uV = split(header{contains(header, '-ADBitVolts')}, 'Volts');
AD2uV = str2double(AD2uV{2})*10^6;

%%
TTprev = '';
for i = 1:length(sortedFiles)
    [cellName, ~] = strtok(sortedFiles(i).name, '.');
    spikeTimes = [load(strcat(sortedPath, sortedFiles(i).name))]';
    [TTname, unitNum] = strtok(cellName, 'SS');
    TTname = TTname(1:end-1);
    unitNum = unitNum(end);
    if strcmp(TTname, TTprev) == false % if the tetrode has changed, load a new one
        tmp_TTname = [TTname '.ntt'];
        TTdir = fullfile(sortedPath, tmp_TTname);
        [tt_ts, tt_sig] = Nlx2MatSpike(TTdir, [1 0 0 0 1], 0, 1, 1);
        TTprev = TTname;
    end
    for j = 1:4
        WaveForm{j} = AD2uV*squeeze(tt_sig(:, j, ismember(tt_ts, spikeTimes)))';
    end

    figure2();
    subplot(2,1,1); hold on; title(strcat(session, '_', cellName),'Interpreter','none')
    ylabel('Amplitude (\muV)');
    for j = 1:4
        plotFilledStd([1:32]+32*(j-1), WaveForm{j}, 'k');
    end
    line([0 128], [0 0], 'color', [0.7, 0.7, 0.7]);
    text(10, 150, sprintf('spikeNumber %d' , length(spikeTimes)));
    ylim([-200 200]);
    
    subplot(2,1,2); hold on;
    ylabel('density');
    h = histogram(diff(spikeTimes),'BinWidth',1000,'Normalization','probability');
    set(gca, 'XScale', 'log')
    line([2000 2000],[0 0.002],'color',[1 0 0]);
    formatSpec = ['RP violation = %.2f ' '%'];
    text(1000, 0.002, [sprintf(formatSpec, 100*sum(diff(spikeTimes)<2000)/(length(spikeTimes)-1)) '%']);
    title('ISI')
end
  