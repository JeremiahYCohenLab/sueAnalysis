% 40 Hz sine
Fs = 20000;
sinHz = 40;
rampTail = .1; % ramp down at last 100-ms

% 1100-ms
duration = 1100;

dur = duration/1000;
t = 1/Fs:1/Fs:dur;
y = (sin(2*pi*sinHz*t-pi/2)+1)/2;

rampFun = ones(1,numel(t));
rampOnset = (dur - rampTail)*Fs;
rampFun(round((1:Fs*rampTail)+rampOnset)) = 1-(1:Fs*rampTail)/(Fs*rampTail);
y = y.*rampFun;
audiowrite('40HzSine_1100ms.wav',y,Fs)

% 1325-ms
duration = 1325;

dur = duration/1000;
t = 1/Fs:1/Fs:dur;
y = (sin(2*pi*sinHz*t-pi/2)+1)/2;

rampFun = ones(1,numel(t));
rampOnset = (dur - rampTail)*Fs;
rampFun(round((1:Fs*rampTail)+rampOnset)) = 1-(1:Fs*rampTail)/(Fs*rampTail);
y = y.*rampFun;
audiowrite('40HzSine_1325ms.wav',y,Fs)

% 800-ms
duration = 800;

dur = duration/1000;
t = 1/Fs:1/Fs:dur;
y = (sin(2*pi*sinHz*t-pi/2)+1)/2;

rampFun = ones(1,numel(t));
rampOnset = (dur - rampTail)*Fs;
rampFun(round((1:Fs*rampTail)+rampOnset)) = 1-(1:Fs*rampTail)/(Fs*rampTail);
y = y.*rampFun;
audiowrite('40HzSine_800ms.wav',y,Fs)

% 400-ms
duration = 400;

dur = duration/1000;
t = 1/Fs:1/Fs:dur;
y = (sin(2*pi*sinHz*t-pi/2)+1)/2;

rampFun = ones(1,numel(t));
rampOnset = (dur - rampTail)*Fs;
rampFun(round((1:Fs*rampTail)+rampOnset)) = 1-(1:Fs*rampTail)/(Fs*rampTail);
y = y.*rampFun;
audiowrite('40HzSine_400ms.wav',y,Fs)

% 500-ms
duration = 500;

dur = duration/1000;
t = 1/Fs:1/Fs:dur;
y = (sin(2*pi*sinHz*t-pi/2)+1)/2;

rampFun = ones(1,numel(t));
rampOnset = (dur - rampTail)*Fs;
rampFun(round((1:Fs*rampTail)+rampOnset)) = 1-(1:Fs*rampTail)/(Fs*rampTail);
y = y.*rampFun;
audiowrite('40HzSine_500ms.wav',y,Fs)

% 175-ms
duration = 175;

dur = duration/1000;
t = 1/Fs:1/Fs:dur;
y = (sin(2*pi*sinHz*t-pi/2)+1)/2;

rampFun = ones(1,numel(t));
rampOnset = (dur - rampTail)*Fs;
rampFun(round((1:Fs*rampTail)+rampOnset)) = 1-(1:Fs*rampTail)/(Fs*rampTail);
y = y.*rampFun;
audiowrite('40HzSine_175ms.wav',y,Fs)

% 300-ms
duration = 300;

dur = duration/1000;
t = 1/Fs:1/Fs:dur;
y = (sin(2*pi*sinHz*t-pi/2)+1)/2;

rampFun = ones(1,numel(t));
rampOnset = (dur - rampTail)*Fs;
rampFun(round((1:Fs*rampTail)+rampOnset)) = 1-(1:Fs*rampTail)/(Fs*rampTail);
y = y.*rampFun;
audiowrite('40HzSine_300ms.wav',y,Fs)

% 1300-ms
duration = 1300;

dur = duration/1000;
t = 1/Fs:1/Fs:dur;
y = (sin(2*pi*sinHz*t-pi/2)+1)/2;

rampFun = ones(1,numel(t));
rampOnset = (dur - rampTail)*Fs;
rampFun(round((1:Fs*rampTail)+rampOnset)) = 1-(1:Fs*rampTail)/(Fs*rampTail);
y = y.*rampFun;
audiowrite('40HzSine_1300ms.wav',y,Fs)

% 125-ms
duration = 125;

dur = duration/1000;
t = 1/Fs:1/Fs:dur;
y = (sin(2*pi*sinHz*t-pi/2)+1)/2;

rampFun = ones(1,numel(t));
rampOnset = (dur - rampTail)*Fs;
rampFun(round((1:Fs*rampTail)+rampOnset)) = 1-(1:Fs*rampTail)/(Fs*rampTail);
y = y.*rampFun;
audiowrite('40HzSine_125ms.wav',y,Fs)

% 1200-ms
duration = 1200;

dur = duration/1000;
t = 1/Fs:1/Fs:dur;
y = (sin(2*pi*sinHz*t-pi/2)+1)/2;

rampFun = ones(1,numel(t));
rampOnset = (dur - rampTail)*Fs;
rampFun(round((1:Fs*rampTail)+rampOnset)) = 1-(1:Fs*rampTail)/(Fs*rampTail);
y = y.*rampFun;
audiowrite('40HzSine_1200ms.wav',y,Fs)

%% 50-ms ramp down
Fs = 20000;
sinHz = 40;
rampTail = .05; % ramp down at last 100-ms

% 225-ms
duration = 225;

dur = duration/1000;
t = 1/Fs:1/Fs:dur;
y = (sin(2*pi*sinHz*t-pi/2)+1)/2;

rampFun = ones(1,numel(t));
rampOnset = (dur - rampTail)*Fs;
rampFun(round((1:Fs*rampTail)+rampOnset)) = 1-(1:Fs*rampTail)/(Fs*rampTail);
y = y.*rampFun;
audiowrite('40HzSine_225ms_down50ms.wav',y,Fs)

% 125-ms
duration = 125;

dur = duration/1000;
t = 1/Fs:1/Fs:dur;
y = (sin(2*pi*sinHz*t-pi/2)+1)/2;

rampFun = ones(1,numel(t));
rampOnset = (dur - rampTail)*Fs;
rampFun(round((1:Fs*rampTail)+rampOnset)) = 1-(1:Fs*rampTail)/(Fs*rampTail);
y = y.*rampFun;
audiowrite('40HzSine_125ms_down50ms.wav',y,Fs)

%%  square pulse with ramp down

Fs = 20000;
flatT = 3;  %time in s
rampT = 0.5;

% 400ms flat with 100ms ramp down
flat = repmat(1, flatT*Fs, 1);
ramp = linspace(0, 1, rampT*Fs);
y = [flat; fliplr(ramp)'];
figure; plot(y)

txt = ['flat' num2str(flatT*1000) 'msRamp' num2str(rampT*1000) 'ms.wav']
audiowrite(txt,y,Fs)




