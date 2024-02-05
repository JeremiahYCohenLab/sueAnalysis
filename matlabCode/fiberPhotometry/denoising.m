function denoised = denoising(rawTrace, fs, fc)
    dn1_R1_sig = medfilt1(rawTrace,1); %median filter gets rid of random/ salt and pepper noise 
     %cutoff frequency in Hz (anything lower than this passes)
%     fs = 40; %sampling rate in Hz
    [b,a] = butter(2,fc/(fs/2),'low');
%     [b,a] = butter(2,0.99999999,'low'); % getting rid of signal that cannot be well-described by sampling frequency
    dn2_R1_sig = filtfilt(b,a,dn1_R1_sig); %filtfilt is zero-phase so no temporal distortion
    %% remove exp
    % exponential decay (bleaching)
    % 'a*exp(-b*x)+c*exp(-d*x)+e'
    bleachRemoved = dn2_R1_sig;
    startValue = zeros(1, 5);
    % last 1min
    startValue(5) = mean(bleachRemoved(length(bleachRemoved)-60*fs:end));
    % major decay a
    startValue(2) = (mean(bleachRemoved(1:60*fs)) - mean(bleachRemoved(length(bleachRemoved)-60*fs:end)))/(mean(bleachRemoved(60*fs:120*fs)) - mean(bleachRemoved(length(bleachRemoved)-60*fs:end)));
    startValue(2) = log(startValue(2))/60;
    % minor delay c
    startValue(4) = startValue(2)/3;
    
    % scaling factor
    startValue(1) = 0.5 * (mean(bleachRemoved(1:60*fs)) - mean(bleachRemoved(length(bleachRemoved)-60*fs:end)));
    startValue(3) = 0.5 * (mean(bleachRemoved(1:60*fs)) - mean(bleachRemoved(length(bleachRemoved)-60*fs:end)));
    time = (1:length(bleachRemoved))'/fs;
    
    startValue(~isreal(startValue)) = 0;
    startValue = real(startValue);
    startValue(startValue<0) = 0;
    [fitresult] = doubleExpFit(time, bleachRemoved, startValue);
    decay = feval(fitresult, time);
    dn2_R1_sig = bleachRemoved - decay;
    %% High-pass filter method of accounting for slow decay 
    fc2 = 0.01; %cutoff frequency in Hz (anything higher than this f passes)
    [b2,a2] = butter(2,fc2/(fs/2),'high'); %this accounts for slow time course changes on the order of ~16 minutes
    denoised = filtfilt(b2,a2,dn2_R1_sig) + 5; %filtfilt is zero-phase %add a constant so the trace doesn't go down to zero on y-axis
    % SR = 30;
    denoised = dn2_R1_sig;
end