function generateTTs(path, varargin)
p = inputParser;
    p.addParameter('HighPassCutoffInHz', 600)
    p.addParameter('LowPassCutoffInHz', 6000);
    p.addParameter('ThresholdFactor', 2.5);
    p.addParameter('RefractorySamples', 20);
    p.addParameter('medianSubtraction', 1);
    p.addParameter('flipSign', true);
    p.addParameter('AnalyzeSpecificTTs', [1:16]);
    p.addParameter('groupSize', 1); 
    p.parse(varargin{:});
    
    % mapInd = [13, 12, 11, 14, 10, 15, 9, 16, 8, 1, 7 , 2, 6, 3, 5, 4];
    mapInd = [1:16];
    session = Session(path);
    streamKey = session.recordNodes{1,1}.recordings{1,1}.continuous.keys();
    streamKey = streamKey{1};
    samples = session.recordNodes{1,1}.recordings{1,1}.continuous(streamKey).samples;
    samples = samples(mapInd, :);
    timeStamps =  session.recordNodes{1,1}.recordings{1,1}.continuous(streamKey).timestamps;
    samplingFreq = session.recordNodes{1,1}.recordings{1,1}.continuous(streamKey).metadata.sampleRate;
    numChannels = session.recordNodes{1,1}.recordings{1,1}.continuous(streamKey).metadata.numChannels;
    
    Wn = [p.Results.HighPassCutoffInHz p.Results.LowPassCutoffInHz] /(samplingFreq/2);
    [b, a] = butter(2, Wn, 'bandpass');    
    
    saveDir = [path '\' 'spikes'];
    
    if ~exist(saveDir, 'dir')
        mkdir(saveDir);
    end
    
    
    % median subtraction
    if p.Results.medianSubtraction
        referenceC = mean(samples([1:7, 9:13],:), 1);
    else
        referenceC = zeros(size(timeStamps))';
    end
    
    % loop through channels, 
    if isempty(p.Results.AnalyzeSpecificTTs)
        ttRange = 1:numChannels;
    else
        ttRange = p.Results.AnalyzeSpecificTTs;
    end
    
    for currC = ttRange
        currSamp = double(samples(currC, :)) - referenceC;

        % filtering
        currSamp = filtfilt(b, a, currSamp);
        % flip 
        if p.Results.flipSign == true
           currSamp = -currSamp;
        end     
        % amplify for spikesorting
        currSamp = 10*currSamp;
        fprintf(['Channel' num2str(currC) ' mean ' '%0.2f \n'], mean(currSamp))
        % spike detection
        % threshold and median method from Rey, Pedreira, Quiroga (2015)
        thresh = p.Results.ThresholdFactor*round(median(abs(currSamp), 2)/0.6745);
        allLocs = peakseek(currSamp, 20, thresh);
        fprintf('\tTotal of %0.1d spikes thresholded.\n', length(allLocs));
        % create ntt file
        sampBack = round(1/3*32) - 1;
        sampFor = round(2/3*32);
        
        if allLocs
            samp = [currSamp; zeros(3, length(timeStamps))]';
            % remove the last one if too long
            if allLocs(end)+sampFor > size(samp,1)
                allLocs = allLocs(1:end-1);
            end
            % remove the first one if too early
            if allLocs(1) - sampBack <= 0
                allLocs = allLocs(2:end);
            end            

            sampToSave = zeros(sampFor+sampBack+1,4,length(allLocs));


            leg = -sampBack:sampFor;
            for w = 1:length(leg)
                sampToSave(w,:,:) = samp(allLocs+leg(w), :)';
            end
            
            featToSave = [];

            featToSave(1:4,:) = squeeze(max(sampToSave,[],1));
            featToSave(5:8,:) = squeeze(min(sampToSave,[],1));

            tsToSave = round(1000*timeStamps(allLocs))';

            header = cell(0);

            Mat2NlxSpike([saveDir '\' 'TT' num2str(currC) '.ntt'], 0, 1, [], [1 1 1 1 1 1], tsToSave, ...
            zeros(1, length(tsToSave)), zeros(1, length(tsToSave)), featToSave, sampToSave, header);
            
        end        
        

    end
    
    
end