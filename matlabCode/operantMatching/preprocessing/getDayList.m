function dayList = getDayList(xlFile, sheet, col)
    [root, sep] = currComputer();
    [~, dayList, ~] = xlsread([root xlFile '.xlsx'], sheet);
    cate = cell(1,size(dayList,2));
    cate(:) = {col};
    cate = cellfun(@strcmp, dayList(1,:), cate)>0;
    % [~,col] = find(~cellfun(@isempty,strfind(dayList, category)) == 1);
    dayList = dayList(2:end,cate);
    endInd = find(cellfun(@isempty,dayList),1);
    if ~isempty(endInd)
        dayList = dayList(1:endInd-1,:);
    end
end