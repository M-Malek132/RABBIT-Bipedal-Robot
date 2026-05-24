function main_demo_stepping_stone()
%MAIN_DEMO_STEPPING_STONE  Main entry point for the RABBIT stepping-stone demo.
%
% This script:
%   1) loads model parameters,
%   2) initializes the state,
%   3) runs a multi-step hybrid simulation,
%   4) saves and visualizes results.
%
% Expected compatible helpers:
%   - parameters()
%   - make_initial_state()
%   - simulate_n_steps()
%   - rabbit_clf_cbf_controller()
%   - save_results()           % optional, if you have it
%   - illustrate_results()     % optional, if you have it

clc;
close all;

fprintf('========================================================\n');
fprintf(' RABBIT Biped Robot: Stepping-Stone Simulation\n');
fprintf('========================================================\n');

% Add project folders to path
addpath(genpath(pwd));

% Load robot parameters
params = parameters();

% Initialize state
[x0, ~] = make_initial_state();

% Number of walking steps equals number of stones
nSteps = size(params.stones, 1);
fprintf('Simulation setup: %d steps\n', nSteps);

% Global step index used by the controller wrapper
global CURRENT_STEP;
CURRENT_STEP = 1;

% Controller wrapper:
% simulate_n_steps expects a controller handle of the form controller(t, x)
controller = @(t, x, params) execution_wrapper(t, x, params);

% Run hybrid multi-step simulation
fprintf('Starting multi-step simulation...\n');
[t_all, x_all, impact_log] = simulate_n_steps(x0, params, nSteps, controller);

% Save results if your project provides a saver
if exist('save_results', 'file') == 2
    save_results(t_all, x_all, params, impact_log);
end

% Visualize if your project provides a visualization routine
if exist('illustrate_results', 'file') == 2
    call_latest_result();
end

fprintf('Simulation complete.\n');
end

% -------------------------------------------------------------------------
% Local controller wrapper
% -------------------------------------------------------------------------
function u = execution_wrapper(t, x, params)
%EXECUTION_WRAPPER  Bridges simulate_n_steps() and the step-indexed controller.
%
% This wrapper reads CURRENT_STEP and passes it to rabbit_clf_cbf_controller.

global CURRENT_STEP;

% Clamp to valid range just in case
step_idx = min(max(CURRENT_STEP, 1), size(params.stones, 1));

% Call your step-aware controller.
% Adjust the output signature here if your controller returns more values.
[u, ~, ~] = rabbit_clf_cbf_controller(t, x, params, step_idx);
end
