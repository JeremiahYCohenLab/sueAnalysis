function createSuperStruct_opMD(xlFile, animals, expName, varargin)

p = inputParser;
% default parameters if none given
p.addParameter('category', 'all');
p.addParameter('revForFlag',0);
p.addParameter('intanFlag',0);
p.parse(varargin{:});

[root, sep] = currComputer();
s = struct;


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
        sN = dayList{sI};
        s.(aN).(sN) = struct;
        
        [bD] = behDataForStruct_opMD(sN, 'revForFlag', p.Results.revForFlag);
        s.(aN).(sN).s = bD.behSessionData;
        s.(aN).(sN).b = rmfield(bD, 'behSessionData');
        
        if isstrprop(sN(end), 'alpha')
            sPath = [root aN sep sN(1:end-1) sep 'sorted' sep 'session ' sN(end) sep];
            s.(aN).(sN).p.sessionFolder = sN(1:end-1);
            s.(aN).(sN).p.saveFigFolder = [root aN sep sN(1:end-1) sep 'figures' sep];
        else
            sPath = [root aN sep sN sep 'sorted' sep 'session' sep];
            s.(aN).(sN).p.sessionFolder = sN;
            s.(aN).(sN).p.saveFigFolder = [root aN sep sN sep 'figures' sep];
        end
        sortedFolder = dir(sPath);
        
        s.(aN).(sN).p.animalName = aN;
        s.(aN).(sN).p.sessionName = sN;
        s.(aN).(sN).p.sortedFolder = sPath;
        
            
        if any(~cellfun(@isempty,strfind({sortedFolder.name},'_intan.mat'))) == 1
            sessionDataInd = ~cellfun(@isempty,strfind({sortedFolder.name},'_intan.mat'));
            load([sPath sortedFolder(sessionDataInd).name])
            s.(aN).(sN).e = sessionData;
        elseif any(~cellfun(@isempty,strfind({sortedFolder.name},'_nL.mat'))) == 1
            sessionDataInd = ~cellfun(@isempty,strfind({sortedFolder.name},'_nL.mat'));
            load([sPath sortedFolder(sessionDataInd).name])
            s.(aN).(sN).e = sessionData;
        end
    end
end

savePath = ['C:\Users\cooper_PC\Desktop\data\' expName sep];
if ~exist(savePath)
    mkdir(savePath)
end
save([savePath 's.mat'], 's')


end