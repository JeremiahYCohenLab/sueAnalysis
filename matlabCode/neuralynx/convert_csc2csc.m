function convert_csc2csc(path_nlyx, path_save)
    file = [path_nlyx 'CSC1.ncs'];
    if path_save == 0
        path_save = path_nlyx;
    end
    file_save = [path_save '\csc_data.hdf5'];
    % Load everything
    % [Timestamps, ChannelNumbers, SampleFrequencies,NumberOfValidSamples, Samples, Header] = Nlx2MatCSC('test.ncs', [1 1 1 1 1], 1, 1, [] );
    [Timestamps, Header] = Nlx2MatCSC(file, [1 0 0 0 0], 1, 1, [] );
    % [Timestamps, Header] = Nlx2MatNRD(file, 0, [1 0], 1, 1, []);
    
    %% save metadata
    headerStruct = struct();
    
    % Parse into structure and auto-convert numeric fields
    for i = 1:length(Header)
        line = strtrim(Header{i});
        if startsWith(line, '-')
            splitIdx = regexp(line, '\s', 'once');
            if ~isempty(splitIdx)
                key = strrep(line(1:splitIdx-1), '-', '');
                key = matlab.lang.makeValidName(key); % ensure valid field name
                valueStr = strtrim(line(splitIdx+1:end));
    
                % Try converting to number
                valueNum = str2double(valueStr);
                if ~isnan(valueNum) && isfinite(valueNum)
                    value = valueNum;
                else
                    value = valueStr;
                end
    
                headerStruct.(key) = value;
            end
        end
    end
    
    %% save timestamps
    % === 3. Write timestamps ===
    timestamps_sec = double(Timestamps) * 1e-6;  % microseconds → seconds
    h5create(file_save, '/timestamps', size(timestamps_sec), 'Datatype', 'double');
    h5write(file_save, '/timestamps', timestamps_sec);
    
    %% 
    
    % === 5. Write header as attributes ===
    headerFields = fieldnames(headerStruct);
    for i = 1:numel(headerFields)
        key = headerFields{i};
        value = headerStruct.(key);
        h5writeatt(file_save, '/', key, value);
    end
    
    
    
    %% 
    all_samples = cell(32);
    parfor channel = 1:32
        curr_path = [path_nlyx sprintf('CSC%d.ncs', channel)];
        samples = Nlx2MatCSC(curr_path, [0 0 0 0 1], 0, 1, []);
        samples = double(samples);  % ensure HDF5-compatible
        samples = reshape(samples, [], 1);
        all_samples{channel} = samples;
    end
    all_matrix = cell2mat(all_samples');
    all_matrix = double(all_matrix);

    h5create(file_save, sprintf('/samples'), size(all_matrix), 'Datatype', 'double');
    h5write(file_save, sprintf('/samples'), all_matrix)
%%
end






