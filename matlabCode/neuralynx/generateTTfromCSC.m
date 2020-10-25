function generateTTfromCSC(session, varargin)
% written by BAB modified by ZS on 08/25/2020
[root, sep] = currComputer();
p = inputParser;
% default parameters if none given
p.addParameter('Root', root);
p.addParameter('Separator', sep)
p.addParameter('opto', true);
p.addParameter('subFolder', '1')
p.addParameter('RemoveLick_Flag', false)
p.addParameter('FilterTrace_Flag', true)
p.addParameter('HighPassCutoffInHz', 300)
p.addParameter('LowPassCutoffInHz', []);
p.addParameter('SamplingFreq', 32000)
p.addParameter('ThresholdFactor', 2.5);
p.addParameter('RefractorySamples', 32); % Neuralynx gives 24 samples before looking for new spike
p.addParameter('AnalyzeSpecificTTs', 1:8);
p.addParameter('RescaleCSCs_Flag', false);
% p.addParameter('CSCscaleFactor', []);
% p.addParameter('CSCstoScale', []);

p.parse(varargin{:});

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

fprintf('Analyzing %s\n', pd.sessionFolder)
for currTT = 1:length(TTnum)
    currTTnum = TTnum(currTT);
    
    chan0 = currTTnum + (3*(currTTnum - 1));
    chan1 = chan0 + 1;
    chan2 = chan0 + 2;
    chan3 = chan0 + 3;

    fprintf('Currently on TT%0.1d: CSC %0.1d,%0.1d,%0.1d,%0.1d. ', currTTnum, chan0, chan1, chan2, chan3)
%     header = Nlx2MatSpike([pd.nLynxFolder TTname], [0 0 0 0 0], 1 , 1, []);
    [ts, samp0] = Nlx2MatCSC([pd.nLynxFolder 'CSC' num2str(chan0) '.ncs'], [1 0 0 0 1], 0, 1, []);
    [samp1] = Nlx2MatCSC([pd.nLynxFolder 'CSC' num2str(chan1) '.ncs'], [0 0 0 0 1], 0, 1, []);
    [samp2] = Nlx2MatCSC([pd.nLynxFolder 'CSC' num2str(chan2) '.ncs'], [0 0 0 0 1], 0, 1, []);
    [samp3] = Nlx2MatCSC([pd.nLynxFolder 'CSC' num2str(chan3) '.ncs'], [0 0 0 0 1], 0, 1, []);
 
%     samp0 = -1 * samp0;
%     samp1 = -1 * samp1;
%     samp2 = -1 * samp2;
%     samp3 = -1 * samp3;

    if p.Results.RescaleCSCs_Flag == true && ~isempty(tsDifferences) % any pauses in this session; indicative of change in CSC
        indMin = 1;
        for currPause = 1:length(tsDifferences)
            indMax = find(diff(ts) == tsDifferences(currPause));
            scaleFactor = CSCscaleFactor(currPause);
            if isempty(CSCstoScale) % rescale all CSCs
                samp0(:, indMin:indMax) = scaleFactor*samp0(:, indMin:indMax);
                samp1(:, indMin:indMax) = scaleFactor*samp1(:, indMin:indMax);
                samp2(:, indMin:indMax) = scaleFactor*samp2(:, indMin:indMax);
                samp3(:, indMin:indMax) = scaleFactor*samp3(:, indMin:indMax);
                fprintf('\t Rescaling CSCs%i,%i,%i,%i from indices %i to %i\n', chan0, chan1, chan2, chan3, indMin, indMax)
            else
                if any(CSCstoScale == chan0)
                    samp0(:, indMin:indMax) = scaleFactor*samp0(:, indMin:indMax);
                    fprintf('\t Rescaling CSC%i from indices %i to %i\n', chan0, indMin, indMax)
                end
                if any(CSCstoScale == chan1)
                    samp1(:, indMin:indMax) = scaleFactor*samp1(:, indMin:indMax);
                    fprintf('\t Rescaling CSC%i from indices %i to %i\n', chan1, indMin, indMax)
                end
                if any(CSCstoScale == chan2)
                    samp2(:, indMin:indMax) = scaleFactor*samp2(:, indMin:indMax);
                    fprintf('\t Rescaling CSC%i from indices %i to %i\n', chan2, indMin, indMax)
                end
                if any(CSCstoScale == chan3)
                    samp3(:, indMin:indMax) = scaleFactor*samp3(:, indMin:indMax);
                    fprintf('\t Rescaling CSC%i from indices %i to %i\n', chan3, indMin, indMax)
                end
            end
            indMin = indMax + 1;
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
        samp0 = filtfilt(b, a, samp0(:));
        samp1 = filtfilt(b, a, samp1(:));
        samp2 = filtfilt(b, a, samp2(:));
        samp3 = filtfilt(b, a, samp3(:));
    else
        fprintf('Not filtering the data. ');
%         samp0 = samp0(:)*20;
%         samp1 = samp1(:)*20;
%         samp2 = samp2(:)*20;
%         samp3 = samp3(:)*20;
        samp0 = samp0(:);
        samp1 = samp1(:);
        samp2 = samp2(:);
        samp3 = samp3(:);
    end

    tSamp = 1/p.Results.SamplingFreq * 1e6; % time per sample in microseconds
    
    if length(unique(diff(ts))) == 1 % no pausing
        ts_interp = ts(1):tSamp:ts(1) + tSamp*(length(samp0) - 1);
    else % if pausing/skip due to data loss, use the proper for loop
        ts_interp = NaN(1, length(samp0)); 
        for i = 1:length(ts)
            ts_interp(512*(i - 1) + 1:512*(i)) = ts(i):tSamp:ts(i) + tSamp*511;
        end
    end

    % threshold and median method from Rey, Pedreira, Quiroga (2015)
    thresh0 = p.Results.ThresholdFactor*round(median(abs(samp0))/0.6745);
    thresh1 = p.Results.ThresholdFactor*round(median(abs(samp1))/0.6745);
    thresh2 = p.Results.ThresholdFactor*round(median(abs(samp2))/0.6745);
    thresh3 = p.Results.ThresholdFactor*round(median(abs(samp3))/0.6745);
    locs0 = peakseek(samp0, p.Results.RefractorySamples, thresh0); % look for a peak, avoid 32 samples (1ms at 32kHz)
    locs1 = peakseek(samp1, p.Results.RefractorySamples, thresh1);
    locs2 = peakseek(samp2, p.Results.RefractorySamples, thresh2);
    locs3 = peakseek(samp3, p.Results.RefractorySamples, thresh3);

    allLocs = unique(sort([locs0 locs1 locs2 locs3]), 'stable');
    allLocs(allLocs > length(ts_interp) - p.Results.RefractorySamples) = []; % remove spikes within 1ms of the end of recording
    allLocs(allLocs < p.Results.RefractorySamples) = []; % remove spikes within 1ms of the beginning of recording
    while any(diff(allLocs) < p.Results.RefractorySamples) % while there is an overlap in peaks within 1ms
        for i = 1:length(allLocs) - 1
            if allLocs(i + 1) < allLocs(i) + p.Results.RefractorySamples
                % save each trace and find where the best peak is
                tmp(1, :) = samp0(allLocs(i):allLocs(i) + p.Results.RefractorySamples);
                tmp(2, :) = samp1(allLocs(i):allLocs(i) + p.Results.RefractorySamples);
                tmp(3, :) = samp2(allLocs(i):allLocs(i) + p.Results.RefractorySamples);
                tmp(4, :) = samp3(allLocs(i):allLocs(i) + p.Results.RefractorySamples);
                [~, tmpi] = max(max(tmp));
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

    locMat = NaN(length(allLocs), p.Results.RefractorySamples);
    sampBack = round(1/3*32) - 1;
    sampFor = round(2/3*32);
    for i = 1:length(allLocs)
        locMat(i, :) = allLocs(i) - sampBack:allLocs(i) + sampFor;
    end

%     header{16} = ['-ADBitVolts ' head0{16}(13:end) ' ' head1{16}(13:end) ' ' head2{16}(13:end) ' ' head3{16}(13:end)];
%     header{21} = ['-InputRange ' head0{21}(13:end) ' ' head1{21}(13:end) ' ' head2{21}(13:end) ' ' head3{21}(13:end)];
    while any(locMat(end, :) > length(samp0)) % in case a sample goes over the end of recording
        locMat(end, :) = [];
    end

    clear sampToSave featToSave tsToSave
    sampToSave(:,1,:) = samp0(locMat)';
    sampToSave(:,2,:) = samp1(locMat)';
    sampToSave(:,3,:) = samp2(locMat)';
    sampToSave(:,4,:) = samp3(locMat)';

    featToSave(1,:) = max(samp0(locMat)');
    featToSave(2,:) = max(samp1(locMat)');
    featToSave(3,:) = max(samp2(locMat)');
    featToSave(4,:) = max(samp3(locMat)');
    featToSave(5,:) = min(samp0(locMat)');
    featToSave(6,:) = min(samp1(locMat)');
    featToSave(7,:) = min(samp2(locMat)');
    featToSave(8,:) = min(samp3(locMat)');

    tsToSave = ts_interp(allLocs);

    header = cell(0);
    Mat2NlxSpike([pd.nLynxFolder 'TT' num2str(currTTnum) '.ntt'], 0, 1, [], [1 1 1 1 1], tsToSave, ...
        zeros(1, length(tsToSave)), zeros(1, length(tsToSave)), featToSave, sampToSave, header);
end

fprintf('Finished\n')