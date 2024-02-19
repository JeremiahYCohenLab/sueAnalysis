dayList = getDayList('npOpto', 'optoTagging', 'withUnits');
waveformRange = {[1:7], [1:60]}; % 2:6, 5:15
%% combine session by session data
root = 'F:\npOptoRecordings\withOpto\';
allWFs = [];
loc = [];
for i = 1:length(dayList)
    session = dayList{i};
    wfPath = [root session '\sorted\wfLoc.mat'];
    load(wfPath, 'maxChannelOpto','waveforms2DAligned', 'waveformsAligned', 'waveformsOpto');
    currWFLin = waveforms2DAligned(:, waveformRange{1}, waveformRange{2});
    currWFLin = permute(currWFLin, [3,2,1]);
    baseline = mean(currWFLin(1:3,:,:), 1);
    currWFLin = currWFLin - baseline;
    currWFLin = reshape(currWFLin, [], size(waveforms2DAligned, 1))';
    currWFLin = currWFLin./max(squeeze(currWFLin(:, 181:240)), [], 2);
    allWFs = [allWFs; currWFLin];
    loc = [loc; [maxChannelOpto - mean(maxChannelOpto)]/(max(maxChannelOpto) - min(maxChannelOpto))];
end
%%

meanWF = mean(allWFs, 1, 'omitmissing');
stdWF = std(allWFs, 0, 1, 'omitmissing');

allWFzs = (allWFs - meanWF)./stdWF;

[coeff,score,latent, ~, explained, mu] = pca(allWFs);


%% clustering
numCat = 3;
indAll = {};
dis = {};
for a = 1:10
    [indAll{a}, ~, dis{a}] = kmeans([score(:, 1:5)], numCat);
end
[~,optiInds] = min(cellfun(@mean, dis));
ind = indAll{optiInds};
%% plotting
Nstep = 1000;
Pstep = 1000;
myMap = [[linspace(0, 1, Nstep)' linspace(0, 1, Nstep)' ones(Nstep, 1)];...
[ones(Pstep, 1) linspace(1, 0, Pstep)' linspace(1, 0, Pstep)']];
for i = 1:5
    temp = coeff(:,i);
    figure2;
    imagesc(reshape(temp, [], length(waveformRange{1}))');
    colormap(myMap)
    clim([-max(abs(temp)) max(abs(temp))])
    colorbar
end

figure2;
hold on
for i = 1:numCat
scatter3(score(ind==i,4), score(ind==i,2), score(ind==i,3))
end
xlabel('1')
ylabel('2')
zlabel('3')

figure2;
hold on
for i = 1:numCat
    plot(mean(allWFs(ind==i, :)))
end

figure;
hold on
edges = linspace(min(loc)-0.01, max(loc)+0.01, 15);
for i = 1:numCat
    histogram(loc(ind==i), edges, 'Normalization', 'probability')
end
   
%% 1D clustering
root = 'F:\npOptoRecordings\withOpto\';
allWFs = [];
loc = [];
for i = 1:length(dayList)
    session = dayList{i};
    wfPath = [root session '\sorted\wfLoc.mat'];
    load(wfPath, 'maxChannelOpto','waveforms2DAligned', 'waveformsAligned', 'waveformsOpto');
    currWFLin = squeeze(mean(waveforms2DAligned(:, waveformRange{1}, waveformRange{2}), 2));
    % normalizing
    currWFLin = currWFLin - mean(currWFLin(:, 1:10),2);
    currWFLin = currWFLin./(max(currWFLin, [], 2));
    allWFs = [allWFs; currWFLin];
    loc = [loc; [maxChannelOpto - mean(maxChannelOpto)]/(max(maxChannelOpto) - min(maxChannelOpto))];
end
%%

meanWF = mean(allWFs, 1, 'omitmissing');
stdWF = std(allWFs, 0, 1, 'omitmissing');

allWFzs = (allWFs - meanWF)./stdWF;

[coeff,score,latent, ~, explained, mu] = pca(allWFzs);

%% clustering
numCat = 4;
indAll = {};
dis = {};
for a = 1:10
    [indAll{a}, ~, dis{a}] = kmeans([score(:, :)], numCat);
end
[~,optiInds] = min(cellfun(@mean, dis));
ind = indAll{optiInds};
%% plotting

figure2;
hold on
for i = 1:numCat
scatter3(score(ind==i,4), score(ind==i,5), score(ind==i,3), 'filled')
end
figure2;
hold on
time = ([1:420] - 22)/30000 * 1000;
for i = 1:numCat
    plot(time, mean(allWFs(ind==i, :)))
end

figure;
hold on
edges = linspace(min(loc)-0.01, max(loc)+0.01, 15);
for i = 1:numCat
    histogram(loc(ind==i), edges, 'Normalization', 'probability')
end

legend
    

Nstep = 1000;
Pstep = 1000;
myMap = [[linspace(0, 1, Nstep)' linspace(0, 1, Nstep)' ones(Nstep, 1)];...
[ones(Pstep, 1) linspace(1, 0, Pstep)' linspace(1, 0, Pstep)']];

for i = 1:numCat
    figure2; 
    temp = squeeze(mean(allWFs2D(ind==i, :, :)));
    imagesc(temp);
    title([num2str(i)]);
    colormap(myMap)
    clim([-max(abs(temp), [], 'all') max(abs(temp), [], 'all')])
    colorbar
end

%% plot data
%% combine session by session data
root = 'F:\npOptoRecordings\withOpto\';
allWFs2D = [];
loc = [];
allSessions = [];
for i = 1:length(dayList)
    session = dayList{i};
    wfPath = [root session '\sorted\wfLoc.mat'];
    load(wfPath, 'maxChannelOpto','waveforms2DAligned', 'waveformsAligned', 'waveformsOpto');
    allWFs2D = cat(1, allWFs2D, waveforms2DAligned(:, waveformRange{1}, waveformRange{2}));
    sessionRep = repmat(i, size(waveforms2DAligned, 1), 1);
    allSessions = [allSessions; sessionRep];
    loc = [loc; maxChannelOpto-median(maxChannelOpto)];
end
%%

 pTov = allWFs2D(:, 4, 16) - min(allWFs2D(:, :, 4:16), [], [2,3]) - mean(allWFs2D(:,:,1), 2, 'omitmissing');

%%
Nstep = 1000;
Pstep = 1000;
myMap = [[linspace(0, 1, Nstep)' linspace(0, 1, Nstep)' ones(Nstep, 1)];...
[ones(Pstep, 1) linspace(1, 0, Pstep)' linspace(1, 0, Pstep)']];
for i = 98:115
    figure2;
    imagesc(squeeze(allWFs2D(i, :, :)));
    title([num2str(i) '--' num2str(allSessions(i)) 'C' num2str(ind(i)) 'D' num2str(loc(i)) 'PV' num2str(pTov(i))]);
    colormap(myMap)
    temp = reshape(squeeze(allWFs2D(i, :, :)), [], 1);
    clim([-max(abs(temp)) max(abs(temp))])
    colorbar
end
%% ica
% find independend features

Mdl = rica(allWFs, 5);
score = transform(Mdl, allWFs);

%%
Nstep = 1000;
Pstep = 1000;
myMap = [[linspace(0, 1, Nstep)' linspace(0, 1, Nstep)' ones(Nstep, 1)];...
[ones(Pstep, 1) linspace(1, 0, Pstep)' linspace(1, 0, Pstep)']];
for i = 98:115
    temp = Mdl.TransformWeights(:,i);
    figure2;
    imagesc(reshape(temp, [], length(waveformRange{1}))');
    colormap(myMap)
    clim([-max(abs(temp)) max(abs(temp))])
    colorbar
    title(num2str(i))
end

%% realign loc
locNew = loc;
for i = 1:length(dayList)

    % fprintf(num2str(i))
    currLoc = loc(allSessions == i);
    currInd = ind(allSessions == i);

    tempCat = currInd;
    tempCat(currInd == 1 | currInd == 4) = 1;
    tempCat(currInd == 2 | currInd == 3) = 0;

    glm = fitglm(currLoc, tempCat, 'Distribution', 'binomial');
    threshold = (0.5 - glm.Coefficients.Estimate(1))/glm.Coefficients.Estimate(2);
    
    locNew(allSessions == i) = loc(allSessions == i) - threshold;

end

%% perm
locPerm = loc;
for i = 1:length(dayList)

    % fprintf(num2str(i))
    currLoc = loc(allSessions == i);
    currInd = ind(allSessions == i);
    currLoc = currLoc(randperm(length(currLoc)));

    tempCat = currInd;
    tempCat(currInd == 1 | currInd == 4) = 1;
    tempCat(currInd == 2 | currInd == 3) = 0;

    glm = fitglm(currLoc, tempCat, 'Distribution', 'binomial');
    threshold = (0.5 - glm.Coefficients.Estimate(1))/glm.Coefficients.Estimate(2);
    
    locPerm(allSessions == i) = loc(allSessions == i) - threshold;

end

%% ratio
edges = linspace(min(locNew)-0.1, max(locNew)+0.1, 6);
ratio = NaN(1, length(edges)-1);
typeI = histcounts(locNew(ind==2 | ind==3), edges);
typeII = histcounts(locNew(ind==1), edges);

ratioI = typeI./histcounts(locNew, edges);
ratioII = typeII./histcounts(locNew, edges);

figure2;
hold on;
plot(ratioI)
plot(ratioII)
legend('I', 'II')

%% ratio real

locNewReal = loc(~isnan(ind)) * 10;
indNewReal = ind(~isnan(ind));
pI = sum(indNewReal==2 | indNewReal==1)/length(indNewReal);
pII = sum(indNewReal==3 | indNewReal==4)/length(indNewReal);
edges = linspace(min(locNewReal)-0.1, max(locNewReal)+0.1, 6);
disReal = 0.5 *(edges(1:end-1) + edges(2:end));
ratio = NaN(1, length(edges)-1);
typeI = histcounts(locNewReal(indNewReal==1 | indNewReal==3), edges);
typeII =     histcounts(locNewReal(indNewReal==2 | indNewReal==4), edges);
all = histcounts(locNewReal, edges);
ratioI = typeI./all;
ratioII = typeII./all;

figure2;
hold on;
plot(disReal, ratioI, 'Color', 'k', 'LineWidth', 2)
plot(disReal, ratioII, 'Color', [0.7 0.7 0.7], 'LineWidth', 2)
legend('I', 'II')
xlabel('Centered distance (um)')
ylabel('Proprotion')
%% ratio 4
locNewReal = 10*locNew;
edges = linspace(min(locNewReal)-0.1, max(locNewReal)+0.1, 6);
disReal = 0.5*(edges(1:end-1) + edges(2:end));
ratio = NaN(1, length(edges)-1);
typeRatio = cell(numCat, 1);
for i = 1:numCat
    typeRatio{i} = histcounts(locNewReal(ind==i), edges)./histcounts(locNewReal, edges);

end

figure2;
hold on;
for i = 1:numCat
    plot(disReal, typeRatio{i});
end
legend
%%
figure2;
edges = linspace(min(locNewReal)-0.001, max(locNewReal)+0.001, 12);
histogram(locNewReal(indNewReal==1|indNewReal==4), edges, 'Normalization', 'probability', 'FaceColor', 'k', 'EdgeColor', 'k')
hold on;
histogram(locNewReal(indNewReal==3|indNewReal==2), edges, 'Normalization', 'probability', 'FaceColor', [0.6 0.6 0.6], 'EdgeColor', [0.6 0.6 0.6])
%% GMM clustering
% run the EM-algorithm
scoreTemp = score(~isnan(sum(score, 2)),:);
% initialization

[expectations, theta] = fitGaussMixture(scoreTemp(:, 3:5),4,'kmeans');

confThresh = 0.4;
[confG, indG] = max(expectations.z, [], 2);
indG(confG<confThresh) = NaN;
indNew = NaN(1,size(score, 1));
indNew(~isnan(sum(score, 2))) = indG;

y = tsne(score);
figure2;
gscatter(y(~isnan(indG),1), y(~isnan(indG),2), indNew(~isnan(indNew)));
%%