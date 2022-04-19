function pathData = parseSessionString_df(fileOrFolder, root, sep)
% parseSessionString_df    Parses input string to generate corresponding pathData outputs
%   INPUTS
%       fileOrFolder: sessionName or name of .asc file
%           e.g.: 'mBB041d20161006' or 'mBB041d20161006.asc'
%       root: root folder
%           e.g.: 'G:\'
%       sep: separator
%           e.g.: '\' or '/'
%   OUTPUTS
%       pathData
%           Structure with sessionFolder, sortedFolder, etc...

if contains(fileOrFolder,'.asc') % input is .asc file
    filename = fileOrFolder;
    [animalName, date] = strtok(filename, 'd'); 
    animalName = animalName(2:end);
    date = date(1:9);
    sessionFolder = ['m' animalName date];
    behavioralDataPath = [root animalName sep sessionFolder sep 'behavior' sep filename];
    suptitleName = filename(1:strfind(filename,'.asc')-1);
    saveFigName = suptitleName;
    videopath = [root animalName sep sessionFolder sep 'pupil'];
else % input is the folder
    [animalName, date] = strtok(fileOrFolder, 'd'); 
    animalName = animalName(2:end); 
    date = date(2:9);
    sessionFolder = ['m' animalName 'd' date];
    filepath = [root animalName sep sessionFolder sep 'behavior' sep];
    allFiles = dir(filepath);
    fileInd = contains({allFiles.name},[fileOrFolder '.asc']) | contains({allFiles.name},[fileOrFolder '.asc']);
    behavioralDataPath = [filepath allFiles(fileInd).name];
    if any(fileInd)
        suptitleName = allFiles(fileInd).name(1:end-4);
    else % if looking at a folder w/o behavioral data (e.g. optoID)
        suptitleName = [];
    end
    saveFigName = sessionFolder(~(sessionFolder==sep));
    videopath = [root animalName sep sessionFolder sep 'pupil'];
end
sortedFolderLocation = [root animalName sep sessionFolder sep 'sorted' sep];

% append path information
pathData.suptitleName = suptitleName;
pathData.sessionFolder = sessionFolder;
pathData.sortedFolder = sortedFolderLocation;
pathData.animalName = animalName;
pathData.saveFigName = saveFigName;
pathData.saveFigFolder = [root animalName sep sessionFolder sep 'figures' sep];
pathData.baseFolder = [root animalName sep sessionFolder sep];
pathData.behavioralDataPath = behavioralDataPath;
pathData.date = date;
pathData.videopath = videopath;

if isdir([pathData.baseFolder 'neuralynx'])
    pathData.nLynxFolder = [pathData.baseFolder 'neuralynx' sep];
    pathData.nLynxFolderOpto = [pathData.baseFolder 'neuralynx' sep 'opto' sep];
    pathData.nLynxFolderSession = [pathData.baseFolder 'neuralynx' sep 'session' sep];
end