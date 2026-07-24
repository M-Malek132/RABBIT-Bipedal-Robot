%% startup.m
% Initialize RABBIT Robot Project

fprintf('Initializing RABBIT Robot Project...\n');

%% 1. Reset and Add Paths
restoredefaultpath;
project_root = fileparts(mfilename('fullpath'));

% Add everything recursively, but then we handle potential 
% missing folders or conflicts explicitly if needed.
addpath(genpath(project_root));

%% 2. Setup Results Folder
results_dir = fullfile(project_root, 'Results');
if ~exist(results_dir, 'dir')
    mkdir(results_dir);
    fprintf('Created Results folder.\n');
end

%% 3. Graphics Defaults
set(0, 'DefaultFigureColor', 'w', ...
       'DefaultAxesFontSize', 12, ...
       'DefaultLineLineWidth', 1.5);

%% 4. Verify Critical Functions
% These are the ACTUAL HZD-pipeline entry points. The old legacy demo helpers
% (rabbit_dynamics / rabbit_controller / main_demo) are intentionally NOT checked
% here: they still call functions that no longer exist, yet `which` finds the
% files and would report a misleading [OK]. See README "Broken legacy entry
% points".
required = {'simulate_hzd_gait', 'hzd_constraints', 'hzd_cost', ...
            'rabbit_impact_map', 'rabbit_reset_map', 'theta_of_q', ...
            'animate_hzd_result'};

fprintf('\nVerifying Dependencies:\n');
for i = 1:length(required)
    if isempty(which(required{i}))
        fprintf('  [MISSING] %s\n', required{i});
    else
        fprintf('  [OK]      %s\n', required{i});
    end
end

%% 5. Initialization
% If this failed before, ensure config.m is in the root 
% or a folder added by genpath.
if exist('config', 'file')
    config('init');
else
    warning('config.m not found. Skipping initialization.');
end

fprintf('Startup complete.\n');
