function convert_nrd2csc(path_nlyx, path_save)
    % Get list of all .nrd files
    nrd_files = dir(fullfile(path_nlyx, '*.nrd'));
    
    % Filter files larger than 0.5 GB (0.5 * 1024^3 bytes)
    large_nrd_files = nrd_files([nrd_files.bytes] > 0.5 * 1024^3);
    
    % If any files matched
    if ~isempty(large_nrd_files)
        % Find the index of the largest file
        [~, idx_largest] = max([large_nrd_files.bytes]);
    
        % Get the full path of the largest file
        file = fullfile(path_nlyx, large_nrd_files(idx_largest).name);
    else
        fprintf('No .nrd files larger than 0.5 GB found in %s.\n', path_nlyx);
        return 
    end

    if path_save == 0
        path_save = path_nlyx;
    end
    file_save = [path_save '\raw_data.hdf5'];
    % remove if file_save exists

    if isfile(file_save)
        delete(file_save);
        fprintf('Deleted existing file: %s\n', file_save);
    end

    % Load everything
    [Timestamps, Header] = Nlx2MatNRD(file, 0, [1 0], 1, 1, []);
    
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
    Timestamps = double(Timestamps) * 1e-6;  % microseconds → seconds
    h5create(file_save, '/timestamps', size(Timestamps), 'Datatype', 'double');
    h5write(file_save, '/timestamps', Timestamps);
    % clear Timestamps
    
    %% 
    
    % === 5. Write header as attributes ===
    headerFields = fieldnames(headerStruct);
    for i = 1:numel(headerFields)
        key = headerFields{i};
        value = headerStruct.(key);
        h5writeatt(file_save, '/', key, value);
    end
    
    
    
    %%
    % Create dataset with known size beforehand
    % sample_info = Nlx2MatNRD(file, 0, [0 1], 0, 1, []);
    % Get number of samples
    n_samples = length(Timestamps);
    clear Timestamps
    
    % Create HDF5 dataset: rows = samples, cols = channels
    h5create(file_save, '/samples', [n_samples, 32], 'Datatype', 'int32');
    
    % Read each channel and write as one column
    % parfor channel = 1:32
    %     samples = Nlx2MatNRD(file, channel-1, [0 1], 0, 1, []);
    %     samples = int16(reshape(samples, [], 1));  % Make sure it's a column vector [n_samples x 1]
    %     h5write(file_save, '/samples', samples, [1, channel], size(samples));
    %     % fprintf('Channel %d written, %d samples\n', channel, length(samples));
    % end
    all_samples = int32(zeros(n_samples, 32));
    parfor channel = 1:32
        samples = int32(Nlx2MatNRD(file, channel-1, [0 1], 0, 1, []));
        samples = reshape(samples, [], 1);  % Make sure it's a column vector [n_samples x 1]
        all_samples(:, channel) = samples;
        % fprintf('Channel %d written, %d samples\n', channel, length(samples));
    end
    % h5write(file_save, '/samples', samples, [1, channel], size(samples));
    h5write(file_save, '/samples', all_samples);

    %%
end
