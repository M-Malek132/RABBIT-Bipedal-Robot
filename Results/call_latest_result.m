function call_latest_result()
    % CALL_LATEST_RESULT Automatically finds the most recent simulation
    % results in the /result folder and calls illustrate_results()
    
    result_dir = fullfile(pwd, 'Results');
    
    % Get all .mat files in the directory
    files = dir(fullfile(result_dir, '*.mat'));
    
    if isempty(files)
        fprintf('No result files found in /result directory.\n');
        return;
    end
    
    % Sort files by date modified
    [~, idx] = sort([files.datenum], 'descend');
    latest_file = fullfile(files(idx(1)).folder, files(idx(1)).name);
    
    fprintf('Found %d results. Loading most recent: %s\n', length(files), files(idx(1)).name);
    
    % Call the illustration function
    illustrate_results(latest_file);
end
