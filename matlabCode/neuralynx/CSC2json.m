% loop through all sessions to get recording stard and end time
aniName = 'ZS062';
dayList = getDayList(aniName, aniName, 'all');
[root, sep] = currComputer();
for ind = 1:length(dayList)
    session = dayList{ind};
    pathData = parseSessionString_df(session, root, sep);
    if isfield(pathData, 'nLynxFolderSession')
        file = [pathData.nLynxFolderSession 'CSC1.ncs'];
        file_save = [pathData.nLynxFolderSession 'rec_info.json'];
        if exist(file, 'file')
            [Timestamps, Header] = Nlx2MatCSC(file, [1 0 0 0 0], 1, 1, []);
            % Initialize containers.Map for dictionary behavior
            headerStruct = struct();
        
            for i = 1:length(Header)
                line = strtrim(Header{i});
                if startsWith(line, '-')
                    parts = regexp(line, '^-([\w]+)\s+(.*)', 'tokens', 'once');
                    if ~isempty(parts)
                        key = matlab.lang.makeValidName(parts{1});
                        valStr = strtrim(parts{2});
                        valNum = str2double(valStr);
                        if ~isnan(valNum)
                            value = valNum;
                        else
                            value = valStr;
                        end
                        headerStruct.(key) = value;
                    end
                end
            end
            headerStruct.start = Timestamps(1);
            headerStruct.end = Timestamps(end);
            jsonStr = jsonencode(headerStruct);
            fid = fopen(file_save, 'w');
            fwrite(fid, jsonStr, 'char');
            fclose(fid);
            fprintf([session ' finished \n'])
        else
            fprintf([session ' no csc1.ncs \n'])
        end
    else
        fprintf([session ' no session recording \n'])
    end
end

    