% Get all text files in the current folder
path = 'C:\Users\zhixi\Documents\dlc\tongueTrackingStraight-ZS-2022-12-01\labeled-data';
allDir = dir(path);
ind = ~contains({allDir.name}, 'labeled');
allFolders = {allDir.name};
allFolders = allFolders(ind);
allFolders = allFolders(3:end);
%%
for i = 1:length(allFolders)
    pathTmp = ['C:\Users\zhixi\Documents\dlc\tongueTrackingStraight-ZS-2022-12-01\labeled-data\' allFolders{i}];
    files = dir(pathTmp);
    files = {files.name};
    ind = contains(files, 'Collected') & contains(files, 'inhibition');
    files = files(ind);
    if sum(ind)>0
        % Loop through each file 
        for id = 1:length(files)
            % Get the file name 
            [~, f,ext] = fileparts(files{id});
            rename = ['CollectedData_ZS' ext]; 
            movefile([pathTmp '\' files{id}], [pathTmp '\' rename]); 
        end
    end
end
%%