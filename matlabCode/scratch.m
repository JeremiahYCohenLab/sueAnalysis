signal = GSig(:,1);
time = 1:(length(signal(:,1))/samplingFreq);
[p, f] = pspectrum(signal, samplingFreq);
signal = IsoSig(:,1);
time = 1:(length(signal(:,1))/samplingFreq);
[pIso, f] = pspectrum(signal, samplingFreq);
figure2;
hold on;
plot(f,p/p(1))
plot(f,pIso/pIso(1))
legend({'G', 'Iso'})
set(gca,'xscale','log')
%%
 stan_qLearningFit('grabNE', '676353', 'goodBeh', 'modelName', '5params', 'iter', 10000);
 stan_qLearningFit('grabNE', '673595', 'goodBeh', 'modelName', '5params', 'iter', 10000);
 %%