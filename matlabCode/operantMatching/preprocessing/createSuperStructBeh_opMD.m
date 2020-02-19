function createSuperStructBeh_opMD(xlFile, exp, expName, varargin)

p = inputParser;
% default parameters if none given
p.addParameter('category', 'all');
p.addParameter('revForFlag',0);
p.addParameter('intanFlag',0);
p.parse(varargin{:});

[root, sep] = currComputer();
s = struct;

if iscell(exp)
    animals = exp;
else
    [~, animals, ~] = xlsread('expLists.xlsx', exp);
end
    


for aI = 1:length(animals)
    aN = animals{aI};
    s.(aN) = struct;
    [~, dayList, ~] = xlsread(xlFile, animals{aI});
    [~,col] = find(~cellfun(@isempty,strfind(dayList, p.Results.category)) == 1);
    dayList = dayList(2:end,col);
    endInd = find(cellfun(@isempty,dayList),1);
    if ~isempty(endInd)
        dayList = dayList(1:endInd-1,:);
    end
    
    for sI = 1:length(dayList)
        sN = dayList{sI}
        s.(aN).(sN) = struct;
        
        sessionName = char(dayList{sI});
        [animalName, date] = strtok(sessionName, 'd'); 
        animalName = animalName(2:end);
        date = date(1:9);
        sessionFolder = ['m' animalName date];

        if isstrprop(sessionName(end), 'alpha')
            sessionDataPath = [root animalName sep sessionFolder sep 'sorted' sep 'session ' sessionName(end) sep sessionName '_sessionData_behav.mat'];
        else
            sessionDataPath = [root animalName sep sessionFolder sep 'sorted' sep 'session' sep sessionName '_sessionData_behav.mat'];
        end

        if p.Results.revForFlag
            if exist(sessionDataPath,'file')
                load(sessionDataPath)
                behSessionData = sessionData;
            else
                [behSessionData, blockSwitch, blockProbs] = generateSessionData_behav_operantMatching(sessionName);
            end
        else
            if exist(sessionDataPath,'file')
                load(sessionDataPath)
            else
                [behSessionData, blockSwitch, blockSwitchL, blockSwitchR] = generateSessionData_operantMatchingDecoupled(sessionName);
            end
        end
        
        s.(aN).(sN).s = behSessionData;
        
    end
end

savePath = ['C:\Users\cooper_PC\Desktop\data\' expName sep];
if ~exist(savePath)
    mkdir(savePath)
end
save([savePath expName '.mat'], 's')


end