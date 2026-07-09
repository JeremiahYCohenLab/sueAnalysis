function [all_min, all_max] = read_range(path_nlyx)
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
        all_min = [];
        all_max = [];
        return;
    end
    
    % Load timestamps and header
    [Timestamps, Header] = Nlx2MatNRD(file, 0, [1 0], 1, 1, []);
    
    % Preallocate output arrays
    all_min = zeros(32, 1);
    all_max = zeros(32, 1);


    start = Timestamps(1);

    % Loop over 32 channels
    parfor channel = 1:32
        % Correct MATLAB indexing for timestamps
        samples = Nlx2MatNRD(file, channel-1, [0 1], 0, 4, [start, start+1000000]);
        
        all_min(channel) = min(samples(:));
        all_max(channel) = max(samples(:));
    end
end
