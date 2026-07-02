%% startup.m
% Initialize RABBIT Robot Project

fprintf('Initializing RABBIT Robot Project...\n');

%% Reset MATLAB path to avoid conflicts
restoredefaultpath;
rehash toolboxcache;

%% Get project root directory
project_root = fileparts(mfilename('fullpath'));

%% Add project folders
addpath(genpath(project_root));

%% Create Results folder
results_dir = fullfile(project_root,'Results');

if ~exist(results_dir,'dir')
    mkdir(results_dir);
    fprintf('Created Results folder.\n');
end

%% Graphics defaults
set(groot,'DefaultFigureColor','w');
set(groot,'DefaultAxesFontSize',12);
set(groot,'DefaultLineLineWidth',1.5);

%% Print project info
fprintf('Project root:\n%s\n',project_root);
fprintf('All subfolders added to MATLAB path.\n');

%% Verify important functions
required_functions = {
    'rabbit_dynamics'
    'simulate_one_step'
    'animate_rabbit'
    'rabbit_controller'
};

for i = 1:length(required_functions)

    f = which(required_functions{i});

    if ~isempty(f)
        fprintf('[OK] %s found\n',required_functions{i});
    else
        fprintf('[MISSING] %s not found\n',required_functions{i});
    end

end

%% Initialize configuration
config('init');

fprintf('Startup complete.\n');
