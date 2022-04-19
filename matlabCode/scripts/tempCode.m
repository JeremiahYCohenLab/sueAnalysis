dayList =  getDayList('allDBh-cre', 'all-DBh', 'good');
[root, sep] = currComputer();
for i = 1:length(dayList)
    session = dayList{i};
    pd = parseSessionString_df(session, root, sep);
    if isstrprop(session(end), 'alpha')
        behSessionDataPath = [pd.sortedFolder 'session ' session(end) sep session '_sessionData_behav.mat'];
        pupilAlignPath = [pd.sortedFolder 'session ' session(end) sep session '_pupil.mat'];
    else
        behSessionDataPath = [pd.sortedFolder 'session' sep session '_sessionData_behav.mat'];
        pupilAlignPath = [pd.sortedFolder 'session' sep session '_pupil.mat'];
    end
    if exist(pupilAlignPath, 'file')
        fprintf([session ' aligned \n'])
        load(pupilAlignPath);
        load(behSessionDataPath);
        % append pupil diameter of all timestamps
        list = dir(pd.videopath);
        expressionSkeleton = ['^' session 'DLC' '\w*' '_skeleton.csv' '$']; % any skeleton is good
        skeleton = list(~cellfun(@isempty, cellfun(@(x) regexp(x, expressionSkeleton), {list.name}, 'UniformOutput', false)));
        skeleton = skeleton(1).name;
        diaRaw = csvread([pd.videopath sep skeleton], 2, 0);
        qualF = diaRaw(:, 4)>0.9;
        qualF(end) = 1; % make sure not having a lot of NaN
        dia = interp1(find(qualF>0),diaRaw(qualF,2),1:length(diaRaw(:,2)));
        % put realigned pupil frame back to linear time
        diaRealign = NaN(1,cueFT(1)-1 + ceil(FR/1000 * (behSessionData(end).CSon + 10000 - behSessionData(1).CSon)));
        diaRealign(1:cueFT(1)-1) = dia(1:cueFT(1)-1); % quality before cue is usually good
        for j = 1:length(cueFT)
            if qualInd(j)
                startF = cueFT(1) + round(FR/1000*(behSessionData(j).CSon - behSessionData(1).CSon));
                if j~=length(cueFT)
                    endF = cueFT(1) + round(FR/1000*(behSessionData(j+1).CSon - behSessionData(1).CSon));
                else
                    endF = min(cueFT(1) + round(FR/1000*(behSessionData(j).CSon + 10000 - behSessionData(1).CSon)), length(diaRealign));
                end
                
                if cueFT(j)+endF-startF > length(dia) % in case pupil ended early
                    endF = length(dia)+startF-cueFT(j);
                end
                diaRealign(startF:endF) = dia(cueFT(j):cueFT(j)+endF-startF); 
            end      
        end
        
        save(pupilAlignPath,'dia', 'diaRealign','-append');
        fprintf([session ' pupil time reorganized \n'])
    else
        fprintf([session ' not aligned \n'])
    end

end
%%
dayList = getDayList('dynaVol.xlsx', 'ZS076', 'gap');
for i = 1:length(dayList)
    generateSessionData_behav_operantMatching_RwdDelay_dynaVol(dayList{i});
    analyzeBehavioralData_operantMatching_RwdDelay(dayList{i});
end
%%

%%




