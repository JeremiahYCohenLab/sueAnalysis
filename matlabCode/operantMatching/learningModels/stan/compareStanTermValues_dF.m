function compareStanTermValues_dF(xlFile, sheet, category, beh, varargin)

p = inputParser;
% default parameters if none given
p.addParameter('modelName', 'sixParam_absPePeAN_bi')
p.addParameter('bernFlag', 1)
p.addParameter('varNames', {'peBar' 'aN' 'pe'})
p.parse(varargin{:});

[root, sep] = currComputer();

[~, dayList, ~] = xlsread(xlFile, sheet);
[~,col] = find(~cellfun(@isempty,strfind(dayList, category)) == 1);
dayList = dayList(2:end, col);
endInd = find(cellfun(@isempty,dayList),1);
if ~isempty(endInd)
    dayList = dayList(1:endInd-1,:);
end

varStruct = [];
for currSesh = 1:length(dayList)
    sessionName = dayList{currSesh};
    [animal, ~] = strtok(sessionName, 'd'); 
    animal = animal(2:end);

    if p.Results.bernFlag
        modelPath = [root animal sep animal 'sorted' sep 'stan' sep 'bernoulli' sep p.Results.modelName sep animal...
                beh '_' p.Results.modelName '.mat'];
    else
        modelPath = [root animal sep animal 'sorted' sep 'stan' sep p.Results.modelName sep animal...
                    beh '_' p.Results.modelName '.mat'];
    end
    [t, o] = generateStanModelTerms_opMD(p.Results.modelName, modelPath, sessionName, 0);
    
    if sum(cellfun(@any, regexp(p.Results.varNames, 'pe'))) > 0
        tmp = eval(['[t.' p.Results.varNames{1} ' t.' p.Results.varNames{2} ' abs(t.' p.Results.varNames{3} ')]']);
    else
        tmp = eval(['[t.' p.Results.varNames{1} ' t.' p.Results.varNames{2} ' t.' p.Results.varNames{3} ']']);
    end
    tmp = tmp(~abs(o.allRewards),:);
    varStruct = [varStruct; tmp];
    
end

varStruct(:,2) = log(varStruct(:,2)); %use the ln of the learning rate
varStruct(varStruct(:,2)==-Inf, :) = [];  %replace -Inf learning rates w/ 0
varStruct(varStruct(:,1)==0,:) = [];  %ignore times when peBar == 0

[probMap, edges, ~, loc] = histcn(varStruct);


for i =1:length(varStruct)
    colorInds(i) = probMap(loc(i,1), loc(i,2), loc(i,3));
end

colorMap = jet(max(colorInds)+1);
colors = colorMap(colorInds,:);
figure; scatter3(varStruct(:,1), varStruct(:,2), varStruct(:,3), [], colors, 'filled')
titles = generateParamTitles(p.Results.varNames);
xlabel(titles{1}, 'interpreter', 'latex')
ylabel(['ln(' titles{2} ')'])
zlabel('|RPE|')
set(gca, 'FontSize', 15)

x = [-20:-1:-380];
y = ones(length(x), 1) * 10;
angles = [x' y];
CaptureFigVid(angles,'C:\Users\cooper\Desktop\meow')


%plot 2d heat map of two variables
[probMap, xEdges, yEdges] = histcounts2(varStruct(:,1), varStruct(:,2), 50, 'Normalization', 'Probability');
probMap = rot90(probMap);

figure; 
imagesc(probMap);
colormap(gray(256));
ax = gca;
cb = colorbar('Peer', ax);
cb.Ticks = [];
caxis(ax, [0 max(max(probMap))]);
ylabel(cb, 'probability')
set(ax, 'XTick', [], 'YTick', [], 'FontSize', 15)
xlabel(ax, titles{1}, 'interpreter', 'latex')
ylabel(ax, ['ln(' titles{2} ')'])
set(gcf, 'renderer', 'painters')


