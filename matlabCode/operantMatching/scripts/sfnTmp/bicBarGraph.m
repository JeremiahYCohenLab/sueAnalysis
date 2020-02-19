%for use with compareLH

figure;  hold on;
for i = 1:4
    barAvg(i) = mean(bic(:,i));
    barSem(i) = sem(bic(:,i));
    bar(i, barAvg(i), 'FaceColor', colors(i,:))
    errorbar(i, barAvg(i), barSem(i), '.', 'Color', 'k')
end

xticks([])
xticklabels('')
set(gca, 'tickdir', 'out')
set(gcf, 'renderer', 'painters')