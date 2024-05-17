function pathData = parseSessionString_allen(session, root, sep)
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

% parse behavior dir
tempString = split(session, '_');
aniID = tempString{1};
date = tempString{2};
time = tempString{3};

sessionFolder = [root aniID sep session sep];
sortedFolder = [sessionFolder 'sorted' sep];


% append path information
pathData.aniName = aniID;
pathData.sessionFolder = sessionFolder;
pathData.sortedFolder = sortedFolder;
pathData.saveFigFolder = [sessionFolder sep 'figures' sep];
pathData.date = date;
