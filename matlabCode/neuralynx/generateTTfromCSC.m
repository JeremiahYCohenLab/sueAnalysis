
function generateTTfromCSC(session, varargin)
% written by BAB modified by ZS on 08/25/2020
[root, sep] = currComputer();
p = inputParser;
% default parameters if none given
p.addParameter('Root', root);
p.addParameter('Separator', sep)
p.addParameter('opto', false);
p.addParameter('subFolder', '')
p.addParameter('RemoveLick_Flag', false)
p.addParameter('FilterTrace_Flag', true)
p.addParameter('HighPassCutoffInHz', 300)
p.addParameter('LowPassCutoffInHz', []);
p.addParameter('SamplingFreq', 32000)
p.addParameter('ThresholdFactor', 3);
p.addParameter('RefractorySamples', 20); % Neuralynx gives 24 samples before looking for new spike
p.addParameter('AnalyzeSpecificTTs', 7);
p.addParameter('RescaleCSCs_Flag', false);
p.addParameter('changeReference',true);
p.addParameter('newReference',27)
p.addParameter('flipSign',false)
% p.addParameter('CSCscaleFactor', []);
% p.addParameter('CSCstoScale', []);

p.parse(varargin{:});

brokenChannels = [];
% brokenChannels = [11,18]; %ZS049
% brokenChannels = [2,8,12,17,32];%ZS050
% brokenChannels = [19,23];%ZS051
% brokenChannels = [9,10,11,12,29,30]; %ZS052
brokenChannels = [15]; %ZS061
% brokenChannels = [8]; %ZS059
% brokenChannels = [27]; %ZS060
pd = parseSessionString_df(session, p.Results.Root, p.Results.Separator);
if p.Results.opto
    nLynxDir = dir([pd.nLynxFolderOpto p.Results.subFolder sep]);
    pd.nLynxFolder = [pd.nLynxFolderOpto p.Results.subFolder sep];
else
    nLynxDir = dir([pd.nLynxFolderSession p.Results.subFolder sep]);
    pd.nLynxFolder = [pd.nLynxFolderSession p.Results.subFolder sep];
end

if isempty(nLynxDir)
    error('No neuralynx folder in %s', session)
end

if ~isempty(p.Results.AnalyzeSpecificTTs)
    TTnum = p.Results.AnalyzeSpecificTTs;
else
    TTmask = contains({nLynxDir.name},'TT') & contains({nLynxDir.name}, '.ntt');
    TTnum = [];
    for i = find(TTmask)
        tmpTT = nLynxDir(i).name;
        tmpTT = TTname(3:regexp(tmpTT,'.ntt') - 1);
        tmpTT = str2double(tmpTT);
        TTnum = [TTnum tmpTT];
    end
end

% test for pauses in session
if p.Results.RescaleCSCs_Flag == true
    ts = Nlx2MatCSC([pd.nLynxFolder 'CSC1.ncs'], [1 0 0 0 0], 0, 1, []);
    tsDifferences = unique(diff(ts), 'stable');
    tsDifferences = tsDifferences(2:end);
    if ~isempty(tsDifferences)
        fprintf('Session was paused %i times.\n', length(tsDifferences))
        CSCstoScale = input('Which CSCs should be scaled? Enter nothing to apply scaling uniformly. ');
        CSCscaleFactor = input('What should the scale factors be? Enter nothing to scale nothing. ');
        if ~isempty(CSCscaleFactor)
            if length(tsDifferences) ~= length(CSCscaleFactor)
                error('length(tsDifferences) (%i) and length(CSCscaleFactor) (%i) do not match.\n\n', length(tsDifferences), length(CSCscaleFactor))
            else % lengths match
                if isempty(CSCstoScale)
                    fprintf('Session paused but no specific CSCs given; assuming all will be scaled by CSCscaleFactor %s.\n\n', num2str(CSCscaleFactor))
                else
                    fprintf('Only scaling the following CSCs: %s\n\n', num2str(CSCstoScale))
                end
            end
        else
            fprintf('Session was paused but no CSCscaleFactor inputted; assuming no changes in scaling.\n\n')
        end
    else
        fprintf('No pausing in session; assuming no rescaling necessary.\n\n')
    end
end
if p.Results.changeReference
    reference = Nlx2MatCSC([pd.nLynxFolder 'CSC' num2str(p.Results.newReference) '.ncs'], [0 0 0 0 1], 0, 1, []);
    reference = reshape(reference,[],1);
end

fprintf('Analyzing %s\n', pd.sessionFolder)
for currTT = 1:length(TTnum)
    currTTnum = TTnum(currTT);
    
    chan = 4*(currTTnum-1)+1:4*currTTnum;

    fprintf('Currently on TT%0.1d: CSC %0.1d,%0.1d,%0.1d,%0.1d. ', currTTnum, chan(1), chan(2), chan(3), chan(4))
%     header = Nlx2MatSpike([pd.nLynxFolder TTname], [0 0 0 0 0], 1 , 1, []);

    [ts, samp0, header] = Nlx2MatCSC([pd.nLynxFolder 'CSC' num2str(chan(1)) '.ncs'], [1 0 0 0 1], 1, 1, []);
    [samp1] = Nlx2MatCSC([pd.nLynxFolder 'CSC' num2str(chan(2)) '.ncs'], [0 0 0 0 1], 0, 1, []);
    [samp2] = Nlx2MatCSC([pd.nLynxFolder 'CSC' num2str(chan(3)) '.ncs'], [0 0 0 0 1], 0, 1, []);
    [samp3] = Nlx2MatCSC([pd.nLynxFolder 'CSC' num2str(chan(4)) '.ncs'], [0 0 0 0 1], 0, 1, []);
    
%     AD2uV = split(header{contains(header, '-ADBitVolts')}, 'Volts');
%     AD2uV = str2num(AD2uV{2})*10^6;
    
%     samp0 = -1 * samp0;
%     samp1 = -1 * samp1;
%     samp2 = -1 * samp2;
%     samp3 = -1 * samp3;
     samp = cat(3, samp0, samp1, samp2, samp3);

    if p.Results.RescaleCSCs_Flag == true && ~isempty(tsDifferences) % any pauses in this session; indicative of change in CSC
        indMin = 1;
        for currPause = 1:length(tsDifferences)
            indMax = find(diff(ts) == tsDifferences(currPause));
            scaleFactor = CSCscaleFactor(currPause);
            if isempty(CSCstoScale) % rescale all CSCs
                samp(:, indMin:indMax,:) = scaleFactor*samp(:, inMin:Max, :);
                fprintf('\t Rescaling CSCs%i from indices %i to %i\n', chan, indMin, indMax)
            else
                if ismember(chan, CSCstoScale)
                    chanSc = find(ismember(chan, CSCstoScale));
                    samp(:,indMin:indMax,chanSc) = scaleFactor*samp(:,indMin:indMax,chanSc);
                    fprintf('\t Rescaling CSC%i from indices %i to %i\n', chanSc, indMin, indMax)
                end
            end
            indMin = indMax + 1;
        end
    end
    samp = reshape(samp, [] ,4); 
    if p.Results.changeReference
        for t = 1:4
            samp(:,t) = samp(:,t)- reference;
        end
    end
    if p.Results.FilterTrace_Flag == true
        if isempty(p.Results.LowPassCutoffInHz)
            Wn = p.Results.HighPassCutoffInHz / (p.Results.SamplingFreq/2);
            [b, a] = butter(2, Wn, 'high');
            fprintf('High pass only: %iHz. ', p.Results.HighPassCutoffInHz)
        else
            Wn = [p.Results.HighPassCutoffInHz p.Results.LowPassCutoffInHz] / (p.Results.SamplingFreq/2);
            [b, a] = butter(2, Wn, 'bandpass');
            fprintf('Bandpass: %iHz - %iHz. ', p.Results.HighPassCutoffInHz, p.Results.LowPassCutoffInHz)
        end          
        samp = filtfilt(b, a, samp);
    else
        fprintf('Not filtering the data. ');
    end

    tSamp = 1/p.Results.SamplingFreq * 1e6; % time per sample in microseconds
    
    if length(unique(diff(ts))) == 1 % no pausing
        ts_interp = ts(1):tSamp:ts(1) + tSamp*(size(samp,1) - 1);
    else % if pausing/skip due to data loss, use the proper for loop
        ts_interp = NaN(1, size(samp,1)); 
        for i = 1:length(ts)
            ts_interp(512*(i - 1) + 1:512*(i)) = ts(i):tSamp:ts(i) + tSamp*511;
        end
    end
    
    if p.Results.flipSign == true
        samp = -samp;
    end
    % threshold and median method from Rey, Pedreira, Quiroga (2015)
    thresh = p.Results.ThresholdFactor*round(median(abs(samp), 1)/0.6745);
    
    %remove broken channel
    locs = cell(size(chan));
    for channel = 1:length(chan)
        if ~ismember(chan(channel), brokenChannels) %jump broken channels
            locs{channel} = peakseek(samp(:,channel), p.Results.RefractorySamples, thresh(channel)); % look for a peak, avoid 32 samples (1ms at 32kHz)
        else
            samp(:,channel) = 0;
        end
    end

    allLocs = unique(sort([locs{1}, locs{2}, locs{3}, locs{4}]), 'stable');
    allLocs(allLocs > length(ts_interp) - p.Results.RefractorySamples) = []; % remove spikes within 1ms of the end of recording
    allLocs(allLocs < p.Results.RefractorySamples) = []; % remove spikes within 1ms of the beginning of recording
    while any(diff(allLocs) < p.Results.RefractorySamples) % while there is an overlap in peaks within 1ms, maybe means same neuron
        for i = 1:length(allLocs) - 1
            if allLocs(i + 1) < allLocs(i) + p.Results.RefractorySamples
                % save each trace and find where the best peak is
                tmp = samp(allLocs(i):min(allLocs(i) + p.Results.RefractorySamples, size(samp,1)),:);
                [~, tmpi] = max(max(tmp'));
                allLocs(i) = allLocs(i) + tmpi - 1;
                allLocs(i + 1) = allLocs(i);
            end
        end
        allLocs = unique(allLocs, 'stable');
    end

    % remove lick artifact from TT; remove the 1ms preceding every lick event; kludgy solution
    if p.Results.RemoveLick_Flag == true
        [tsTTLs, TTLs] = Nlx2MatEV([pd.nLynxFolder 'Events.nev'], [1 0 1 0 0], 0, 1);
        lickTTLs = [16 32 18 33 82 97]; % possible lick events on 295H Neuralynx rig
        tsLicks = tsTTLs(ismember(TTLs, lickTTLs));
        lickInds = find(ismember(round(tsLicks/10)*10, round(ts_interp/10)*10)); % round due to slide differences in times from .nev file and .ncs file

        for i = 1:length(lickInds)
            tmpLickArtifact = allLocs - lickInds(i);
            allLocs(tmpLickArtifact >= -p.Results.SamplingFreq*1e-3 & tmpLickArtifact < 0) = [];
        end
    end

    fprintf('\tTotal of %0.1d spikes thresholded.\n', length(allLocs));

    sampBack = round(1/3*32) - 1;
    sampFor = round(2/3*32);
    
    % in case a sample goes over the end of recording
    if allLocs
        if allLocs(end)+sampFor > size(samp,1)
            allLocs = allLocs(1:end-1);
        end

        clear sampToSave featToSave tsToSave
        sampToSave = zeros(sampFor+sampBack+1,4,length(allLocs));
%         for i = 1:length(allLocs)
%             sampToSave(:,:,i) = samp(allLocs(i) - sampBack:allLocs(i) + sampFor, :);
%         end
        
        leg = -sampBack:sampFor;
        for w = 1:length(leg)
            sampToSave(w,:,:) = samp(allLocs+leg(w), :)';
        end
            

    %     header{16} = ['-ADBitVolts ' head0{16}(13:end) ' ' head1{16}(13:end) ' ' head2{16}(13:end) ' ' head3{16}(13:end)];
    %     header{21} = ['-InputRange ' head0{21}(13:end) ' ' head1{21}(13:end) ' ' head2{21}(13:end) ' ' head3{21}(13:end)];

        featToSave(1:4,:) = squeeze(max(sampToSave,[],1));
        featToSave(5:8,:) = squeeze(min(sampToSave,[],1));

        tsToSave = ts_interp(allLocs);

        header = cell(0);

        Mat2NlxSpike([pd.nLynxFolder 'TT' num2str(currTTnum) '.ntt'], 0, 1, [], [1 1 1 1 1], tsToSave, ...
            zeros(1, length(tsToSave)), zeros(1, length(tsToSave)), featToSave, sampToSave, header);
    end
end

fprintf('Finished\n')