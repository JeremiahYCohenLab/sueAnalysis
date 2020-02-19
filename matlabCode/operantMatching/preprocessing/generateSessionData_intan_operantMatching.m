function [sessionData] = generateSessionData_intan_operantMatching(sessionName)

% Determine if computer is PC or Mac and set roots and separators appropriately
[root,sep] = currComputer();

%% Assign variables
% Path
[animalName] = strtok(sessionName, 'd');
animalName = animalName(2:end);
if isstrprop(sessionName(end), 'alpha')
    ephysPath = [root animalName sep sessionName(1:end-1) sep 'ephys' sep 'spikeSort_session ' sessionName(end) sep];
    sortedPath = [root animalName sep sessionName(1:end-1) sep 'sorted' sep 'session ' sessionName(end) sep];
else
    ephysPath = [root animalName sep sessionName sep 'ephys' sep 'spikeSort_session' sep];
    sortedPath = [root animalName sep sessionName sep 'sorted' sep 'session' sep];
end

% DI numbers
waterR = 0;
waterL = 1;
CSminus = 2;
CSplus = 3;
lickR = 4;
lickL = 5;
diArray = [waterR waterL CSminus CSplus lickR lickL];


%% Import event TTLs, TTL timestamps, and cluster data (if it exists)

%extract and organize digital input instances
load([ephysPath 'events.mat']);
allEvents = [];
for diInd = 1:length(diArray)
    inputTmp = eval(['events.DIN0', num2str(diArray(diInd))]);
    inputTmp = [inputTmp; ones(1,length(inputTmp))*diArray(diInd)];
    allEvents = [allEvents inputTmp];
end
[~,sortInds] = sort(allEvents,2);
allEvents = allEvents(:, sortInds(1,:));
allEvents(1,:) = round(allEvents(1,:));
    
% Clustered data
sortedFiles = dir(fullfile(ephysPath, '*.txt'));
for i = 1:length(sortedFiles)
    load([ephysPath sortedFiles(i).name]);
end

%% Generate structure

% Initialize structure
sessionData.trialType = [];
sessionData.trialEnd = [];
sessionData.CSon = [];
sessionData.licksL = [];
sessionData.licksR = [];
sessionData.rewardL = [];
sessionData.rewardR = [];
sessionData.rewardTime = [];
sessionData.allSpikes = [];

responseWindow = 2000;

% Find tetrode data (files with SS and TT in the name)
allVar = who;
allClust = find(~cellfun(@isempty,strfind(allVar,'C_')) | ~cellfun(@isempty,strfind(allVar,'TT')))';
for m = 1:length(allClust)
    sessionData(m).allSpikes = transpose(round(eval(allVar{allClust(m)})));
end


trial = 1;
for i = 1:length(allEvents)
    if allEvents(2,i) == CSplus || allEvents(2,i) == CSminus % if this is a trial begin and continue to end of trial
        if allEvents(2,i) == CSplus
            sessionData(trial).trialType = 'CSplus';
        elseif allEvents(2,i) == CSminus
            sessionData(trial).trialType = 'CSminus';
        end
        sessionData(trial).CSon = allEvents(1,i);
        sessionData(trial).rewardR = NaN; 
        sessionData(trial).rewardL = NaN;
        sessionData(trial).rewardTime = NaN;
        
        currInd = i+1; % begin at the next index
        if currInd <= length(allEvents)
            tEndFlag = false; % end when this is true
            while tEndFlag == false
                if allEvents(2,currInd) == lickR % right lick
                    sessionData(trial).licksR = [sessionData(trial).licksR allEvents(1,currInd)];
                    currInd = currInd + 1;
                elseif allEvents(2,currInd) == lickL % left lick
                    sessionData(trial).licksL = [sessionData(trial).licksL allEvents(1,currInd)];
                    currInd = currInd + 1;
                elseif allEvents(2,currInd) == waterR & allEvents(2,currInd-1) == lickR % R reward; waterR by itself could be evap code
                    if allEvents(1,currInd) - allEvents(1,currInd-1) <= responseWindow
                        sessionData(trial).rewardR = 1;
                        sessionData(trial).rewardTime = allEvents(1,currInd);
                        currInd = currInd + 1;
                    end                  
                elseif allEvents(2,currInd) == waterL & allEvents(2,currInd-1) == lickL % L reward; waterL by itself could be evap code
                    if allEvents(1,currInd) - allEvents(1,currInd-1) <= responseWindow
                        sessionData(trial).rewardL = 1;
                        sessionData(trial).rewardTime = allEvents(1,currInd);
                        currInd = currInd + 1;
                    end
                elseif allEvents(2,currInd) == waterL || allEvents(2,currInd) == waterR
                    currInd = currInd + 1;
                end

                if currInd <= length(allEvents) % if not at the end of the file, stop when index is on the next trial
                    if allEvents(2,currInd) == CSplus || allEvents(2,currInd) == CSminus
                        sessionData(trial).trialEnd = allEvents(1,currInd);
                        tEndFlag = true;
                    end
                else % if at the end of the file, escape
                    sessionData(trial).trialEnd = NaN;
                    tEndFlag = true;
                end
            end
        else
            sessionData(trial).trialEnd = NaN;
        end

        % whole bunch of code to figure out when to indicate no reward given
        if isnan(sessionData(trial).rewardR) && isnan(sessionData(trial).rewardL) % if fields are default
            if ~isempty(sessionData(trial).licksR) && ~isempty(sessionData(trial).licksL) % L and R licks
                lickRlat = sessionData(trial).licksR(1) - sessionData(trial).CSon;
                lickLlat = sessionData(trial).licksL(1) - sessionData(trial).CSon;
                if lickRlat <= responseWindow && lickLlat <= responseWindow % L and R within lick window; shouldn't ever happen
                    if lickRlat < lickLlat
                        sessionData(trial).rewardR = 0;
                        sessionData(trial).rewardTime = sessionData(trial).licksR(1);
                    else
                        sessionData(trial).rewardL = 0;
                        sessionData(trial).rewardTime = sessionData(trial).licksL(1);
                    end
                elseif lickRlat <= responseWindow % R lick within lick window
                    sessionData(trial).rewardR = 0;
                    sessionData(trial).rewardTime = sessionData(trial).licksR(1);
                elseif lickLlat <= responseWindow % L lick within lick window
                    sessionData(trial).rewardL = 0;
                    sessionData(trial).rewardTime = sessionData(trial).licksL(1);
                end
            elseif ~isempty(sessionData(trial).licksR) && isempty(sessionData(trial).licksL) % R lick window
                lickRlat = sessionData(trial).licksR(1) - sessionData(trial).CSon;
                if lickRlat <= responseWindow % R lick within lick window
                    sessionData(trial).rewardR = 0;
                    sessionData(trial).rewardTime = sessionData(trial).licksR(1);
                end
            elseif isempty(sessionData(trial).licksR) && ~isempty(sessionData(trial).licksL) % L lick window
                lickLlat = sessionData(trial).licksL(1) - sessionData(trial).CSon;
                if lickLlat <= responseWindow % L lick within lick window
                    sessionData(trial).rewardL = 0;
                    sessionData(trial).rewardTime = sessionData(trial).licksL(1);
                end
            end
        end

        
        % remove bounced lick signals (those that occur within 55ms of one another)
        if ~isempty(sessionData(trial).licksR)
            sessionData(trial).licksR = removeBadLicks(sessionData(trial).licksR);
        end
        if ~isempty(sessionData(trial).licksL)
            sessionData(trial).licksL = removeBadLicks(sessionData(trial).licksL);
        end
        
        % add spike data
        L = sessionData(trial).CSon;
        U = sessionData(trial).trialEnd;
        if isnan(U)
            U = allEvents(1,end);
        end
        for l = allClust
            eval(sprintf('sessionData(trial).%s = %s(%s >= L & %s < U)'';', ...
            allVar{l}, allVar{l}, allVar{l}, allVar{l}))
        end
        
        trial = trial+1;
    end
end

% save the data
if isempty(dir(sortedPath))
    mkdir(sortedPath)
end

% sessionFolder(~(sessionFolder==sep)) removes the separator (/ or \) in the filename before saving
save([sortedPath sep sessionName(~(sessionName==sep)) '_sessionData_intan.mat'], 'sessionData')