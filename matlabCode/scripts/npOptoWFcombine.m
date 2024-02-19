dayList = getDayList('npOpto', 'optoTagging', 'withUnits');
waveformRange = {[1:7], [1:30]}; % 2:6, 5:15
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
    allWFs = [allWFs; currWFLin];
    loc = [loc; [maxChannelOpto - mean(maxChannelOpto)]/(max(maxChannelOpto) - min(maxChannelOpto))];
end
%%

meanWF = mean(allWFs, 1, 'omitmissing');
stdWF = std(allWFs, 0, 1, 'omitmissing');

allWFzs = (allWFs - meanWF)./stdWF;

[coeff,score,latent, ~, explained, mu] = pca(allWFs);


%% clustering
numCat = 4;
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
scatter3(score(ind==i,1), score(ind==i,2), score(ind==i,4))
end
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
    loc = [loc; maxChannelOpto];
end
%%

 pTov = allWFs2D(:, 4, 16) - min(allWFs2D(:, :, 4:16), [], [2,3]) - mean(allWFs2D(:,:,1), 2, 'omitmissing');

%%
Nstep = 1000;
Pstep = 1000;
myMap = [[linspace(0, 1, Nstep)' linspace(0, 1, Nstep)' ones(Nstep, 1)];...
[ones(Pstep, 1) linspace(1, 0, Pstep)' linspace(1, 0, Pstep)']];
for i = 1:length(ind)
    figure2;
    imagesc(squeeze(allWFs2D(i, :, :)));
    title([num2str(i) '--' num2str(allSessions(i)) 'C' num2str(label(i)) 'D' num2str(loc(i)) 'PV' num2str(pTov(i))]);
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
for i = 1:5
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
typeII = histcounts(locNew(ind==1 | ind==4), edges);

ratioI = typeI./histcounts(locNew, edges);
ratioII = typeII./histcounts(locNew, edges);

figure2;
hold on;
plot(ratioI)
plot(ratioII)
legend('I', 'II')
%% ratio 4
edges = linspace(min(locNew)-0.1, max(locNew)+0.1, 6);
ratio = NaN(1, length(edges)-1);
typeRatio = cell(numCat, 1);
for i = 1:numCat
    typeRatio{i} = histcounts(locNew(ind==i), edges)./histcounts(locNew, edges);

end

figure2;
hold on;
for i = 1:numCat
    plot(typeRatio{i});
end
legend
%%