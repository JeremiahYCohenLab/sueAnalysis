animals = {'689515'};
[root, sep] = currComputer();
root = 'I:\sueData\';
newDir = 'F:\';
folders = {'behavior', 'sorted', 'pupil', 'figures'};
states = struct;
for ani = 1:length(animals)
    aniDir = [root animals{ani} sep];
    list = dir(aniDir);
    expression = ['^' 'm' animals{ani} 'd'];
    allFolders = list(~cellfun(@isempty, cellfun(@(x) regexp(x, expression), {list.name}, 'UniformOutput', false)));
    allFolders = {allFolders.name};
    currStates = zeros(length(allFolders), length(folders));
    for ses = 1:length(allFolders)
        for f = 1:length(folders)
            source = [aniDir allFolders{ses} sep folders{f}];
            desti = [newDir animals{ani} sep allFolders{ses} sep folders{f}];
            if exist(source, 'dir')
                if ~exist(desti, 'dir')
                    mkdir(desti)
                end
                currStates(ses, f) = copyfile(source, desti);
            end
        end
    end
    states.(['m' animals{ani}]).sessions = allFolders';
    states.(['m' animals{ani}]).states = currStates;
end
%%