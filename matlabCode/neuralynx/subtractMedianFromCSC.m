function subtractMedianFromCSC(sessionString, varargin)

tic
[root, sep] = currComputer_operantMatching();
p = inputParser;
% default parameters if none given
p.addParameter('Root', root);
p.addParameter('Separator', sep)
p.addParameter('NlynxFolder', []);
p.parse(varargin{:});

if isempty(p.Results.NlynxFolder)
    pd = parseSessionString_oM(sessionString, p.Results.Root, p.Results.Separator);
    nLynxDir = dir(pd.nLynxFolder);
    if isempty(nLynxDir)
        error('No neuralynx folder in %s', sessionString)
    end
    nLynxFolder = pd.nLynxFolder;
else
    nLynxDir = dir(p.Results.NlynxFolder);
    nLynxFolder = p.Results.NlynxFolder;
end
numCSC = sum(contains({nLynxDir.name}, '.ncs'));
if ~ismember(numCSC, [32 64])
    error('CSC number does not equal 32 or 64')
end

% median for median filter subtraction
[timestamps, channelNumbers, sampleFreq, numValidSamp, tempCSC, header{1}] = Nlx2MatCSC([nLynxFolder 'CSC1.ncs'], [1 1 1 1 1], 1, 1, []);
CSCarray = NaN(numCSC, length(tempCSC(:)));
CSCarray(1, :) = tempCSC(:);
fprintf('Loaded CSC1.ncs\n')
parfor i = 2:numCSC
    [tempCSC header{i}] = Nlx2MatCSC([nLynxFolder 'CSC' num2str(i) '.ncs'], [0 0 0 0 1], 1, 1, []);
    CSCarray(i, :) = tempCSC(:);
    fprintf(['Loaded CSC' num2str(i) '.ncs\n']);
end
CSCmedian = median(CSCarray);
CSCarray = CSCarray - CSCmedian;
fprintf('Finished median subtraction\n');

lengthCSC = size(CSCarray, 2)/512;
parfor i = 1:numCSC
    Mat2NlxCSC([nLynxFolder 'CSC' num2str(i) '.ncs'], 0, 1, [], [1 1 1 1 1 1], ...
        timestamps, channelNumbers, sampleFreq, numValidSamp, reshape(CSCarray(i, :), 512, lengthCSC), header{i});
    fprintf(['Saved CSC' num2str(i) '.ncs\n']);
end
save([nLynxFolder 'CSCmedian.mat'], 'CSCmedian')
fprintf('Finished\n')
toc