clc;
clear;

% Get current script directory
current_dir = fileparts(mfilename('fullpath'));

fprintf('Folder Hierarchy:\n\n');

% Start recursive display
displayFolderTree(current_dir, 1);

%% =========================
% Recursive Function
%% =========================
function displayFolderTree(folder, level)

    % Indentation based on hierarchy level
    indent = repmat('    ', 1, level);

    % Get folder name
    [~, folder_name, ~] = fileparts(folder);

    % Print current folder
    fprintf('%s|-- %s\n', indent, folder_name);

    % Get contents
    contents = dir(folder);

    % Remove "." and ".."
    contents = contents(~ismember({contents.name}, {'.', '..'}));

    % Separate folders and files
    dirs = contents([contents.isdir]);
    files = contents(~[contents.isdir]);

    % Print files
    for i = 1:length(files)
        fprintf('%s    |-- %s\n', indent, files(i).name);
    end

    % Recursive call for subfolders
    for i = 1:length(dirs)
        subfolder = fullfile(folder, dirs(i).name);
        displayFolderTree(subfolder, level + 1);
    end
end
