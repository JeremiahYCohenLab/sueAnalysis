function denoised = denoisingKH(rawTrace, fs, fc)
    %median filter gets rid of random/ salt and pepper noise     
    medFiltered = medfilt1(rawTrace,1); 
    %low pass filter cutoff frequency in Hz (anything lower than this passes)
    % getting rid of signal that cannot be well-described by sampling frequency
    fc = 9;
    [b,a] = butter(2, fc/(fs/2),'low'); % in KH, 9
    lowPassed = filtfilt(b,a,medFiltered); %filtfilt is zero-phase so no temporal distortion
    % photobleaching
    % [b,a] = butter(2, 0.05/(fs/2),'high');
    % bleachRemoved = filtfilt(b,a,  lowPassed, padtype='even');
    % %% exponential decay (bleaching)
    % startValue = zeros(1, 5);
    % % last 1min
    % startValue(5) = mean(lowPassed(length(lowPassed)-60*fs:end));
    % % major decay a
    % startValue(2) = (mean(lowPassed(1:60*fs)) - mean(lowPassed(length(lowPassed)-60*fs:end)))/(mean(lowPassed(60*fs:120*fs)) - mean(lowPassed(length(lowPassed)-60*fs:end)));
    % startValue(2) = log(startValue(2))/60;
    % % minor delay c
    % startValue(4) = startValue(2)/3;
    % 
    % % scaling factor
    % startValue(1) = 0.5 * (mean(lowPassed(1:60*fs)) - mean(lowPassed(length(lowPassed)-60*fs:end)));
    % startValue(3) = 0.5 * (mean(lowPassed(1:60*fs)) - mean(lowPassed(length(lowPassed)-60*fs:end)));
    % time = (1:length(lowPassed))'/fs;
    % 
    % startValue(~isreal(startValue)) = 0;
    % startValue = real(startValue);
    % startValue(startValue<0) = 0;
    % [fitresult] = doubleExpFit(time, lowPassed, startValue);
    % decay = feval(fitresult, time);
    % expFixed = lowPassed - decay;

    %% poly fit
    time = (1:length(lowPassed))'/fs;
    [p] = polyfit(time,lowPassed,4);
    polyPred = polyval(p, time);
    polyFixed = lowPassed - polyPred;
    %% high pass to get rid of fluctuations
    [b,a] = butter(2, 0.01/(fs/2), 'high');
    polyFixed = filtfilt(b, a, polyFixed);
    %% get baseline
    [b,a] = butter(2, 0.001/(fs/2), 'low');
    baseline = filtfilt(b, a, lowPassed);
    %% high-pass
    % [b,a] = butter(2, 0.005/(fs/2), 'high');
    % expFixed = filtfilt(b, a, expFixed);
    %%
    dFF = polyFixed;
    %%
    b_percentile = 0.7;
    sorted = sort(dFF);
    b_median = median(sorted(1:round(length(sorted) * b_percentile)));
    denoised = dFF - b_median;

    denoised = dFF;
end