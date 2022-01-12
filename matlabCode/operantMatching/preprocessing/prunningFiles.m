function files = prunningFiles(xlFile, animal,category)
% clean up neuralynx raw recording files; including session recordings

[root, sep] = currComputer();

[~, dayList, ~] = xlsread([root xlFile], animal);
[~,col] = find(contains(dayList, category));
dayList = dayList(2:end,col);
endInd = find(cellfun(@isempty,dayList),1);
if ~isempty(endInd)
    dayList = dayList(1:endInd-1,:);
end
files = struct;
for i = 1:length(dayList)
    sessionName = dayList{i};
    name = sessionName;
    files(i).name = name;
    files(i).rmvSession = 0;
    [animalName, date] = strtok(sessionName, 'd'); 
    animalName = animalName(2:end);
    date = date(1:9);
    sessionFolder = ['m' animalName date];
    
    if isstrprop(sessionName(end), 'alpha')
        dataPath = [root animalName sep sessionFolder sep'];
    else
        dataPath = [root animalName sep sessionFolder sep];
    end
    % initializing as zeros
    % behavior
    % pupil
    % recordings
    recFiles = [dataPath 'neuralynx'];
    if exist([recFiles sep 'session'], 'dir')
        sessionRec = dir([recFiles sep 'session']);
        rawFileIDs = contains({sessionRec.name},'.nrd');
        if sum(rawFileIDs)>0
            sizes = [sessionRec(rawFileIDs).bytes];
            if max(sizes) > 10000000000
                filename = sessionRec(rawFileIDs).name;
                delete([recFiles sep 'session' sep filename]);  
                files(i).rmvSession = 1;
            end
        else
            %in old version raw data were in subfolders
            subFolder = {sessionRec([sessionRec.isdir]).name};
            subFolder = subFolder(~contains(subFolder, '.'));
            for j = 1:length(subFolder)
                subFolderTemp = dir([recFiles sep 'session' sep subFolder{j}]);
                rawFileIDsSub = contains({subFolderTemp.name},'.nrd');
                if sum(rawFileIDsSub)>0
                   sizes = subFolderTemp(rawFileIDsSub).bytes;
                    if max(sizes) > 10000000000
                        filename = subFolderTemp(rawFileIDsSub).name;
                        delete([recFiles sep 'session' sep subFolder{j} sep filename]);
                        files(i).rmvSession = 1;
                        break
                    end
                end
            end
        end
    end
    
%     if exist([recFiles sep 'opto'], 'dir')
%         optoRec = dir([recFiles sep 'opto']);
%         for j = 1:length(optoRec)
%             if contains(optoRec(j).name,'ms')
%                 currOpto = dir([optoRec(j).folder sep optoRec(j).name]);
%                 rawFileIDs = contains({currOpto.name},'.nrd');
%                 if sum(rawFileIDs)>0
%                     sizes = currOpto(rawFileIDs).bytes;
%                     if max(sizes) > 100000000
%                         files(i).opto = 1;
%                         break
%                     end   
%                 end
%             end
%         end
%     end
                
        
    
end

