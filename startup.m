%% startup.m
% Initialize The-RABBIT-Robot project

clc;
clear;
close all;

fprintf('Initializing RABBIT Robot Project...\n');

%% Get project root directory
project_root = fileparts(mfilename('fullpath'));

%% Add all folders and subfolders to MATLAB path
addpath(genpath(project_root));

%% Create Results folder if missing
results_dir = fullfile(project_root, 'Results');

if ~exist(results_dir, 'dir')
    mkdir(results_dir);
    fprintf('Created Results folder.\n');
end

%% Optional graphics settings
set(0, 'DefaultFigureColor', 'w');
set(0, 'DefaultAxesFontSize', 12);
set(0, 'DefaultLineLineWidth', 1.5);

%% Display project structure loaded
fprintf('Project root:\n%s\n', project_root);
fprintf('All subfolders added to MATLAB path.\n');

%% Test basic function accessibility
required_functions = {
    'rabbit_dynamics'
    'simulate_one_step'
    'animate_rabbit'
    'rabbit_controller'
};

for i = 1:length(required_functions)
    if exist(required_functions{i}, 'file')
        fprintf('[OK] %s found\n', required_functions{i});
    else
        fprintf('[MISSING] %s not found\n', required_functions{i});
    end
end

fprintf('Startup complete.\n');
