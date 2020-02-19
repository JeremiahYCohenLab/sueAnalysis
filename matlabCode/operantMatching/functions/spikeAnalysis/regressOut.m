fr = maxFRcs(1,s.responseInds)';
beta_mod = fitlm(1:numel(beta), beta);
fr_mod = fitlm(1:numel(fr), fr);

beta_regressed = beta - beta_mod.predict([1:numel(beta)]');
fr_regressed = fr - fr_mod.predict([1:numel(fr)]');

myMod = fitlm(beta_regressed, fr_regressed)

beta_bins = linspace(min(beta_regressed), max(beta_regressed), 10);
fr_binned = arrayfun(@(i,j) fr_regressed(beta_regressed >= i & beta_regressed < j), ...
    beta_bins(1:end - 1), beta_bins(2:end), 'uniformoutput', false);


figure; 
subplot(211); hold on
scatter(beta_regressed, fr_regressed)

subplot(212); hold on
for i = 1:numel(fr_binned)
    errorbar(beta_bins(i), mean(fr_binned{i}), std(fr_binned{i})/sqrt(numel(fr_binned{i})),'k', 'linewidth', 2)
end