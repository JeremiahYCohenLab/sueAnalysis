function plotROC(xF, xT)
% xF data from real false
% xT data from real true
if size(xF, 2) ~= 1
    xF = xF';
    xT = xT';
end
allX = [xF; xT];
threshs = flip(linspace(min(allX), max(allX), 50));

TP = zeros(size(threshs));
FP = zeros(size(threshs));

for i = 1:length(threshs)
    TP(i) = sum(xT>=threshs(i))/length(xT);
    FP(i) = sum(xF>=threshs(i))/length(xF);
end

auc = trapz(FP, TP);

figure2Wide;
edges = linspace(min(allX) - 0.0001, max(allX) + 0.0001, 20);
subplot(1, 2, 1)
hold on;
histogram(xF, edges, 'Normalization', 'probability');
histogram(xT, edges, 'Normalization', 'probability');
legend({'False', 'True'})
subplot(1, 2, 2)
plot(FP, TP)
xlim([0 1])
ylim([0 1])
xlabel('FP')
ylabel('recall (TP)')

title(sprintf(['AUC = %0.2f'], auc))

end