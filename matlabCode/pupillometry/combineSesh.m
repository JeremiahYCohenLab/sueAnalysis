function [cmbSesh, startTrials, allChoices] = combineSesh(xlFile, animal, category, varargin)
%task and model parameters
p = inputParser;
% default parameters if none given
p.addParameter('numBins', 10)
p.addParameter('plotFlag', 1)
p.addParameter('maxTrials', 1000)

p.parse(varargin{:})

%%
[root, sep] = currComputer();

[~ , dayList, ~] = xlsread([root xlFile], animal);
[~,col] = find(~cellfun(@isempty,strfind(dayList, category)) == 1);
dayList = dayList(2:end,col);
endInd = find(cellfun(@isempty,dayList),1);
if ~isempty(endInd)
    dayList = dayList(1:endInd-1,:);
end

startTrials = [];
%% start wtih first day
for i = 1:length(dayList)
    sessionName = dayList{i};
    [animalName, date] = strtok(sessionName, 'd'); 
    animalName = animalName(2:end);
    date = date(1:9);
    sessionFolder = ['m' animalName date];

    if isstrprop(sessionName(end), 'alpha')
        sessionDataPath = [root animalName sep sessionFolder sep 'sorted' sep 'session ' sessionName(end) sep sessionName '_sessionData_behav.mat'];
    else
        sessionDataPath = [root animalName sep sessionFolder sep 'sorted' sep 'session' sep sessionName '_sessionData_behav.mat'];
    end
    
    if exist(sessionDataPath,'file')
        load(sessionDataPath);
    else                                %otherwise generate the struct
        [behSessionData, ~, ~, ~] = generateSessionData_operantMatchingDecoupledRwdDelay(sessionName);
    end
    
    behSessionData = behSessionData(1:min(p.Results.maxTrials, length(behSessionData)));
    
    if i == 1
        fields = fieldnames(behSessionData);
        fields = fields';
        fields{2,8} = NaN;
        fields{2,9} = NaN;
        fields{2,10} = NaN;
        gap = struct(fields{:});
        gap(2:10) =gap(1);
        cmbSesh = behSessionData;
        startTrials = 1;
    else
        startTrials = [startTrials length(cmbSesh)+11];
        cmbSesh = [cmbSesh, gap, behSessionData];
    end
    
    
end
