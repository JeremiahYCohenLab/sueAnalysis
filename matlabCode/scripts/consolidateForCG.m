desinationRoot = 'F:\forCooper\';
dayList = getDayList('PSpupil', 'PS01pupil', 'all');
bigT = 8000000000;
for i = 1:length(dayList)
    session = dayList{i};
    pd = parseSessionString_df(session, desinationRoot, '\');
    allFiles = dir(pd.videopath);
    allNames = {allFiles.name};
    allSizes = [allFiles.bytes]; 
    videoInd = contains(allNames, '.avi');
    allNames = allNames(videoInd);
    for j = 1:length(allNames)
        delete([pd.videopath '\' allNames{j}]);
    end

end
%%