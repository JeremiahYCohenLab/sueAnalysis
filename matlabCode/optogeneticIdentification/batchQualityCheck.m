animalName = 'ZS061';
[root, sep] = currComputer();
xlFile = [animalName '.xlsx'];
savePath = [root animalName sep animalName 'sorted' sep 'optoFiles' sep];
if ~exist(savePath, 'dir')
    mkdir(savePath)
end
% get corresponding unit files
[nums, unitsInfo,~] = xlsread([root xlFile], 'neurons');
sessionList = unitsInfo(2:end,1); 
unitList = unitsInfo(2:end,2);
quality = nums(:,1);
sessionList = sessionList(quality<=0.05);
unitList = unitList(quality<=0.05);

for i = 1:length(sessionList)
    rasters = qualityCheck(sessionList{i},1,'unit',unitList{i});
    saveFigurePDF(rasters, [savePath sep sessionList{i} '_' unitList{i} '.pdf']);
end  