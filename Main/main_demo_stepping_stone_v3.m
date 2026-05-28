function main_demo_stepping_stone_v3()
%MAIN_DEMO_STEPPING_STONE  Main entry point for the RABBIT stepping-stone demo.
clc; clear; close all;

fprintf('========================================================\n');
fprintf(' RABBIT Biped Robot: Stepping-Stone Simulation\n');
fprintf('========================================================\n');

% 1. Setup Environment
addpath(genpath(pwd));
params = parameters();

% 2. Initialize State
[x0, ~] = make_initial_state();

% 3. Multi-step Configuration
nSteps = size(params.stones, 1);
fprintf('Simulation setup: %d steps\n', nSteps);

global CURRENT_STEP;
CURRENT_STEP = 1;

% 4. Run Simulation
% Pass the execution wrapper which handles the controller calls
controller_handle = @(t, x, p) execution_wrapper(t, x, p);

fprintf('Starting multi-step simulation...\n');
% simulate_n_steps should handle the impact dynamics and step incrementing
[t_all, x_all, impact_log] = simulate_n_steps(x0, params, 2, controller_handle);

% 5. Visualization & Results
if exist('save_results', 'file') == 2
    save_results(t_all, x_all, params, impact_log);
end

if exist('illustrate_results', 'file') == 2
    % Call your specific visualization function
    % ensure illustrate_results(t_all, x_all, params) is compatible
    call_latest_result()
end

fprintf('Simulation complete.\n');
end

% -------------------------------------------------------------------------
% Execution Wrapper
% -------------------------------------------------------------------------
function u = execution_wrapper(t, x, params)
    global CURRENT_STEP;
    [u, ~] = rabbit_clf_controller_v3(t, x, params);
    
    % Optional: Stop simulation if it detects a fall
    if x(2) < 0.3 % If torso height drops too low
        error('Robot fell at t=%.3f during step %d', t, CURRENT_STEP);
    end
end
