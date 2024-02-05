% behavior realign
function behPhotometryPreprocessingAllen(session)
%%
    [root, sep] = currComputer();
    pd = parseSessionString_df(session, root, sep);
     
    p = parseSessionString_df(session, root, sep);
    fpDir = [p.baseFolder sep 'photometry' sep];
    
    allFiles = dir(fpDir);
    allFiles = {allFiles([allFiles.isdir]==0).name}';
    % load ttl file
    ttlInd = contains(allFiles, 'TTL') & contains(allFiles, 'csv');
    ttlFile = allFiles{ttlInd};
    ttlSigInd = contains(allFiles, 'TTL') & ~contains(allFiles, 'csv');
    ttlSignalFile = allFiles{ttlSigInd};
    ttlTime = readmatrix([fpDir ttlFile]); % time in ms, 1 Hz, updates each second
    fileID = fopen([fpDir ttlSignalFile]);
    ttlSignal = fread(fileID, 'float64'); % sampling rate 1kHz
    fclose(fileID);
    ttlSig1 = ttlSignal(1:3:end);
    ttlSig2 = ttlSignal(2:3:end);
    ttlSig3 = ttlSignal(3:3:end);
    timeTS = linspace(ttlTime(1), ttlTime(end)+1000, length(ttlSig1)+1) - 1000; % correction for signal lag

    

    %%  extract square waves
    ttlSig = ttlSig1;
    ttlSig = ttlSig >= 3;
    ttlDiff = [0; ttlSig(2:end) - ttlSig(1:end-1)];
    ttl_p = [];
    ttl_l = [];
    
    for ii = 1:length(ttlSig)
        if ttlDiff(ii) == 1
            for jj = 1:120
                if ii+jj> length(ttlSig)-1
                    break
                end
                if ttlDiff(ii+jj) == -1
                    ttl_p = [ttl_p ii];
                    ttl_l = [ttl_l jj];
                    break
                end
            end              
        end
    end
    
    %%
    behSessionData.trialType = [];
    behSessionData.trialEnd = [];
    behSessionData.CSon = [];
    behSessionData.licksL = [];
    behSessionData.licksR = [];
    behSessionData.rewardL = [];
    behSessionData.rewardR = [];
    behSessionData.respondTime = [];
    behSessionData.rewardTime = [];
    behSessionData.rewardProbL = [];
    behSessionData.rewardProbR = [];
    behSessionData.allLicks = [];
    behSessionData.laser = [];
    
    % find trials starts
    trialOnsInds = find(ttl_l==1);
   
    % check if last trial was finished
    trialEnds = find(ttl_l==40);
    if isempty(find(trialEnds > trialOnsInds(end), 1))
        trialOnsInds = trialOnsInds(1:end-1);
    end

    for i = 1:length(trialOnsInds)-1
        behSessionData(i).CSon = NaN;
        behSessionData(i).trialEnd = NaN;
        behSessionData(i).rewardL = NaN;
        behSessionData(i).rewardR = NaN;
        behSessionData(i).respondTime = NaN;
        behSessionData(i).rewardTime = NaN;
        behSessionData(i).laser = NaN;
    
        tmpChoice = NaN;
        tmpRwd = NaN;
        behSessionData(i).trialType = 'CSplus';
        behSessionData(i).CSon = timeTS(ttl_p(trialOnsInds(i)));

        if ttl_l(trialOnsInds(i) + 1) == 2
            tmpChoice = 0;
            if ttl_l(trialOnsInds(i) + 2) <= 30 && ttl_l(trialOnsInds(i) + 2) > 28
                tmpRwd = 1;
            else
                tmpRwd = 0;
            end
        else 
            if ttl_l(trialOnsInds(i) + 1) == 3
                tmpChoice = 1;
                if ttl_l(trialOnsInds(i) + 2) <= 30 && ttl_l(trialOnsInds(i) + 2) > 28
                    tmpRwd = 1;
                else
                    tmpRwd = 0;
                end          
            end 
        end
            
        if tmpChoice == 0
            behSessionData(i).rewardL = tmpRwd;
            behSessionData(i).respondTime = timeTS(ttl_p(trialOnsInds(i)+1));
            behSessionData(i).rewardTime = timeTS(ttl_p(trialOnsInds(i)+1));
        else 
            if tmpChoice == 1
                behSessionData(i).rewardR = tmpRwd;
                behSessionData(i).respondTime = timeTS(ttl_p(trialOnsInds(i)+1));
                behSessionData(i).rewardTime = timeTS(ttl_p(trialOnsInds(i)+1));
            end
        end
        
        if i<length(trialOnsInds)
            behSessionData(i).trialEnd = timeTS(ttl_p(trialOnsInds(i+1)));
        end
    end
    if ~exist(pd.sortedFolder, 'dir')
        mkdir([pd.sortedFolder])
    end
    save([pd.sortedFolder session '_sessionData_behav.mat'], 'behSessionData');
%%

