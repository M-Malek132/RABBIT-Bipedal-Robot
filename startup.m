%% startup.m
% Initialize RABBIT Robot Project

fprintf('Initializing RABBIT Robot Project...\n');

%% 1. Reset and Add Paths
restoredefaultpath;
project_root = fileparts(mfilename('fullpath'));

% Add everything recursively -- EXCEPT hidden directories.
%
% genpath walks folders whose names begin with a dot, and .claude/worktrees
% holds COMPLETE CHECKOUTS of this repo (that is what a git worktree is). Those
% copies contain a ch3_*.m for nearly every function here, and '.claude' sorts
% BEFORE 'Chapter3', so an unfiltered genpath puts them first on the path and
% every Chapter 3 call silently resolves to whatever commit the worktree is
% parked at. Measured: which('ch3_continuation') returned the worktree copy at
% f363689, four commits behind main, while the edited file sat third in the
% list -- so edits appeared to have no effect, and a function deleted on main
% still ran because a worktree still had it.
%
% Dropping any entry with a dot-prefixed component fixes this for every future
% worktree too, without disturbing the worktrees themselves.
parts = strsplit(genpath(project_root), pathsep);
parts = parts(~cellfun(@isempty, parts));
parts = parts(~contains(parts, [filesep '.']));
addpath(strjoin(parts, pathsep));

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
