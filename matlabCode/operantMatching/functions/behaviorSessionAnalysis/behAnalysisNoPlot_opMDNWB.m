function s = behAnalysisNoPlot_opMDNWB(session, varargin)

p = inputParser;
% default parameters if none given
p.addParameter('simpleFlag', 0)
p.parse(varargin{:});
[root, sep] = currComputerAcute();
% parse behavior dir
tempString = split(session, '_');
aniID = tempString{1};
date = tempString{2};
time = tempString{3};

sessionFolder = [root aniID sep session sep];
sortedFolder = [sessionFolder 'sorted' sep];

% load nwb
allFiles = {dir(sortedFolder).name};
nwbInd = contains(allFiles, session) & contains(allFiles, '.nwb');  
if sum(nwbInd)==0
    fprintf([session ' no nwb file. \n'])
    s = [];
    return
end
nwbName = allFiles{nwbInd};
nwb = nwbRead([sortedFolder nwbName], 'ignorecache');
sessionData = nwb.intervals_trials.toTable;
% output choices and outcomes
s.responseInds = find(sessionData.animal_response == 1 | sessionData.animal_response == 0);
s.allChoices = 2 * sessionData.animal_response(s.responseInds) - 1;
leftRewards = sessionData.rewarded_historyL(s.responseInds);
rightRewards = sessionData.rewarded_historyR(s.responseInds);  
s.allRewards = -leftRewards + rightRewards;
s.sessionData = sessionData;