function myROC = auROCZS(x,y)
% calculate auROC of x and y
% default: au>0.5, y>x
%% TPs and FPs
thresholds = linspace(min([x;y])-0.01, max([x;y])+0.01, 50);
thresholds = thresholds(1:end-1);
thresholds = unique([x;y]);
TPs = zeros(size(thresholds));
FPs = zeros(size(thresholds));
for i = 1:length(thresholds)
   TPs(i) = sum(y >= thresholds(end+1-i))/length(y);
   FPs(i) = sum(x >= thresholds(end+1-i))/length(x);
end
FPs = [0; FPs];
TPs = [0; TPs];
myROC = cumtrapz(FPs, TPs);
myROC = myROC(end);
end
%%