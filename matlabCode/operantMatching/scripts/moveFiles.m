animals = {'ZS066', 'ZS068', 'ZS069', 'ZS070', 'ZS071'};
[root, sep] = currComputer();
newDir = 'H:\';
states = struct;
for ani = 1:length(animals)
    aniDir = [root animals{ani} sep];
    list = dir(aniDir);
    expression = ['^' 'm' animals{ani} 'd'];
    allFolders = list(~cellfun(@isempty, cellfun(@(x) regexp(x, expression), {list.name}, 'UniformOutput', false)));
    allFolders = {allFolders.name};
    currStates = zeros(length(allFolders),1);
    for ses = 1:length(allFolders)
        source = [aniDir allFolders{ses} sep 'sorted'];
        desti = [newDir animals{ani} sep allFolders{ses} sep 'sorted'];
        if exist(source, 'dir')
           currStates(ses) = copyfile(source, desti);
        end
    end
    states.(animals{ani}).sessions = allFolders';
    states.(animals{ani}).states = currStates;
end
%%