function [preFractions, postFractions] = compareChoiceRwdFractions_opMD(xlFile, sheet, pre, post)

preFractions = choiceRwdFractions_opMD(xlFile, sheet, pre);
postFractions = choiceRwdFractions_opMD(xlFile, sheet, post);

v1 = [0 0 0]; 
v2 = [1 1 0];
for i = 1:length(preFractions)
    pt = [preFractions(i,:) 0]; 
    a = v1 - v2; 
    b = pt - v2;
    preDist(i) = norm(cross(a,b)) / norm(a);
end

for i = 1:length(postFractions)
    pt = [postFractions(i,:) 0]; 
    a = v1 - v2; 
    b = pt - v2;
    postDist(i) = norm(cross(a,b)) / norm(a);
end

xVals = linspace(0, 1, 10);
preCoeffs = polyfit(preFractions(:,1), preFractions(:,2), 1);
preY = polyval(preCoeffs, xVals);
postCoeffs = polyfit(postFractions(:,1), postFractions(:,2), 1);
postY = polyval(postCoeffs, xVals);

figure; hold on; 
scatter(preFractions(:,1), preFractions(:,2), 'MarkerEdgeColor', 'c');
scatter(postFractions(:,1), postFractions(:,2), 'MarkerEdgeColor', 'm');
plot([0 1], [0 1], '-k', 'linewidth', 2);
plot(xVals, preY, '--c')
plot(xVals, postY, '--m')
preLeg = ['pre, mean dist: ' num2str(nanmean(preDist))];
postLeg = ['post, mean dist: ' num2str(nanmean(postDist))];
legend(preLeg, postLeg)
set(gca, 'tickdir', 'out')
xlabel('Blockwise Reward Fraction')
ylabel('Blockwise Choice Fraction')
title(sheet)

end