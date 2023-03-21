% data prep
len = 5000;
x = rand(1, len);
y = rand(1, len);
diffXY = reshape(y'-x, [], 1);
%% epsilon greedy
epsilon = 0:0.1:0.5;
colors = repmat(linspace(0, 0.7, length(epsilon)), 3, 1)';
binNum = 20;
edges = linspace(-1, 1, binNum+1);
figure2;
hold on;
for j = 1:length(epsilon)
    currE = epsilon(j);
    randInd = randperm(length(diffXY), round(currE*length(diffXY)));
    choices = diffXY > 0;
    choices(randInd) = rand(round(currE*length(diffXY)), 1)>0.5;
    meanDiff = NaN(1,binNum);
    meanChoices = NaN(1,binNum);
    for i = 1:binNum
        meanDiff(i) = mean(diffXY(diffXY>=edges(i) & diffXY<edges(i+1)));
        meanChoices(i) = mean(choices(diffXY>=edges(i) & diffXY<edges(i+1)));
    end
    hold on; 
    plot(meanDiff, meanChoices, 'LineWidth', 1.5, 'Color', colors(j,:))
end
 legend(cellfun(@num2str, num2cell(epsilon), 'UniformOutput', false))
%% softmax
beta = 8:-1:3;
colors = repmat(linspace(0, 0.7, length(beta)), 3, 1)';
binNum = 20;
edges = linspace(-1, 1, binNum+1);
figure2;
hold on;
for j = 1:length(beta)
    currB = beta(j);
    probs = logistic(currB*diffXY);
    meanDiff = NaN(1,binNum);
    meanChoices = NaN(1,binNum);
    for i = 1:binNum
        meanDiff(i) = mean(diffXY(diffXY>=edges(i) & diffXY<edges(i+1)));
        meanChoices(i) = mean(probs(diffXY>=edges(i) & diffXY<edges(i+1)));
    end
    hold on; 
    plot(meanDiff, meanChoices, 'LineWidth', 1.5, 'Color', colors(j,:))
end
 legend(cellfun(@num2str, num2cell(beta), 'UniformOutput', false))
%% thompson
uncertainty = 0.1:0.05:0.3;
colors = repmat(linspace(0, 0.7, length(uncertainty)), 3, 1)';
binNum = 20;
edges = linspace(-1, 1, binNum+1);
figure2;
hold on;
for j = 1:length(uncertainty)
    xNoise = x + normrnd(0, uncertainty(j), size(x));
    yNoise = y + normrnd(0, uncertainty(j), size(x));
    diffXYnoise = reshape(yNoise'-xNoise, [], 1);
    choices = diffXYnoise > 0;
    meanDiff = NaN(1,binNum);
    meanChoices = NaN(1,binNum);
    for i = 1:binNum
        meanDiff(i) = mean(diffXY(diffXY>=edges(i) & diffXY<edges(i+1)));
        meanChoices(i) = mean(choices(diffXY>=edges(i) & diffXY<edges(i+1)));
    end
    hold on; 
    plot(meanDiff, meanChoices, 'LineWidth', 1.5, 'Color', colors(j,:))    
end
legend(cellfun(@num2str, num2cell(uncertainty), 'UniformOutput', false))
%%



