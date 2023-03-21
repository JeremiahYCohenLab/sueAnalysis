figure2;
subplot(4,1,1);
scatter(currTrial_lick, ones(size(currTrial_lick)));
%
x = zeros(1, max(currTrial_lick)+10);
x(currTrial_lick + 3) = 1;
subplot(4,1,2);
plot(x);
title('raw')
%
fs = 1000;
n = length(x); 
ySig = fft(x);% number of samples
f = (0:n-1)*(fs/n);     % frequency range
powerSig = abs(ySig).^2/n;    % power of the DFT
subplot(4,1,3);
plot(f,powerSig)
%
fc = [0.1 10];
[b, a] = butter(2,fc/(fs/2),'bandpass');
sigFilt = filtfilt(b, a, x);
subplot(4,1,4);
plot(sigFilt);
%%
fo = 8;
q = 2;
bw = (fo/(fs/2))/q;
[b,a] = iircomb(fs/fo,bw,'peak');

sigCom = filter(b, a, x);
subplot(4,1,4);
plot(sigCom);
%%