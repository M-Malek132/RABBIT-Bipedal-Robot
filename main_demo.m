%% MAIN_DEMO.M
% =========================================================
% RABBIT Robot Main Demonstration Script
%
% This script:
%   1. Initializes the project
%   2. Loads robot parameters
%   3. Selects controller mode
%   4. Runs hybrid walking simulation
%   5. Plots simulation results
%   6. Animates the robot
%
% =========================================================

clc;
clear;
close all;

fprintf('\n=========================================\n');
fprintf('      RABBIT ROBOT MAIN DEMO\n');
fprintf('=========================================\n');

%% --------------------------------------------------------
% 1. Initialize Project
% ---------------------------------------------------------

fprintf('\n[1] Initializing project...\n');

startup;

%% --------------------------------------------------------
% 2. Load Robot Parameters
% ---------------------------------------------------------

fprintf('[2] Loading robot parameters...\n');

params = parameters();

%% --------------------------------------------------------
% 3. Select Controller Mode
% ---------------------------------------------------------

fprintf('[3] Selecting controller mode...\n');

% Available modes:
%   'passive'
%   'pd'
%   'hzd'

controller_mode = 'pd';

fprintf('Controller Mode: %s\n', upper(controller_mode));

%% --------------------------------------------------------
% 4. Define Initial Conditions
% ---------------------------------------------------------

fprintf('[4] Defining initial conditions...\n');

% ---------------------------------------------------------
% State Vector Format
%
% x = [q ; dq]
%
% q  = joint positions
% dq = joint velocities
%
% Assuming 5 DOF robot:
%   q  -> 5x1
%   dq -> 5x1
%
% Total state dimension = 10
% ---------------------------------------------------------

nq = 5;

q0  = zeros(nq,1);
dq0 = zeros(nq,1);

% Small forward lean for walking initialization
q0(3) = 0.2;

% Initial state
x0 = [q0; dq0];

%% --------------------------------------------------------
% 5. Simulation Settings
% ---------------------------------------------------------

fprintf('[5] Configuring simulation...\n');

sim.nSteps = 5;
sim.tMax   = 10;
sim.dt     = 0.001;

fprintf('Number of steps : %d\n', sim.nSteps);
fprintf('Maximum time    : %.2f sec\n', sim.tMax);

%% --------------------------------------------------------
% 6. Run Simulation
% ---------------------------------------------------------

fprintf('[6] Running simulation...\n');

try

    switch lower(controller_mode)

        case 'passive'

            [t, x, impacts] = simulate_n_steps( ...
                x0, ...
                params, ...
                sim.nSteps);

        case 'pd'

            controller = @(t,x) pd_controller(t, x, params);

            [t, x, impacts] = simulate_n_steps( ...
                x0, ...
                params, ...
                sim.nSteps, ...
                controller);

        case 'hzd'

            controller = @(t,x) hzd_controller(t, x, params);

            [t, x, impacts] = simulate_n_steps( ...
                x0, ...
                params, ...
                sim.nSteps, ...
                controller);

        otherwise

            error('Unknown controller mode.');

    end

    fprintf('Simulation completed successfully.\n');

catch ME

    fprintf('\nSimulation failed.\n');
    fprintf('Error Message:\n%s\n', ME.message);

    rethrow(ME);

end

%% --------------------------------------------------------
% 7. Basic Validation
% ---------------------------------------------------------

fprintf('[7] Checking simulation results...\n');

if isempty(t)
    error('Simulation returned empty time vector.');
end

if any(isnan(x(:)))
    error('NaN detected in state trajectory.');
end

fprintf('Trajectory validation passed.\n');

%% --------------------------------------------------------
% 8. Plot Results
% ---------------------------------------------------------

fprintf('[8] Plotting results...\n');

figure('Name','Joint Positions');

plot(t, x(:,1:nq), 'LineWidth', 1.5);

xlabel('Time [s]');
ylabel('Joint Angles [rad]');
title('Joint Positions');
grid on;

legend( ...
    'q_1', ...
    'q_2', ...
    'q_3', ...
    'q_4', ...
    'q_5');

figure('Name','Joint Velocities');

plot(t, x(:,nq+1:end), 'LineWidth', 1.5);

xlabel('Time [s]');
ylabel('Joint Velocities [rad/s]');
title('Joint Velocities');
grid on;

legend( ...
    'dq_1', ...
    'dq_2', ...
    'dq_3', ...
    'dq_4', ...
    'dq_5');

%% --------------------------------------------------------
% 9. Energy Analysis
% ---------------------------------------------------------

fprintf('[9] Computing energy...\n');

try

    energy = zeros(length(t),1);

    for k = 1:length(t)

        energy(k) = rabbit_energy_model( ...
            x(k,:)', ...
            params);

    end

    figure('Name','Robot Energy');

    plot(t, energy, 'LineWidth', 2);

    xlabel('Time [s]');
    ylabel('Energy [J]');
    title('Total Mechanical Energy');
    grid on;

catch

    fprintf('Energy model unavailable or failed.\n');

end

%% --------------------------------------------------------
% 10. Animate Robot
% ---------------------------------------------------------

fprintf('[10] Launching animation...\n');

try

    animate_rabbit(t, x, params);

catch ME

    fprintf('Animation failed:\n%s\n', ME.message);

end

%% --------------------------------------------------------
% 11. Save Results
% ---------------------------------------------------------

fprintf('[11] Saving results...\n');

results.time      = t;
results.state     = x;
results.params    = params;
results.controller= controller_mode;

save('Results/demo_results.mat', 'results');

fprintf('Results saved:\n');
fprintf('Results/demo_results.mat\n');

%% --------------------------------------------------------
% 12. Simulation Summary
% ---------------------------------------------------------

fprintf('\n=========================================\n');
fprintf('         SIMULATION COMPLETE\n');
fprintf('=========================================\n');

fprintf('Controller Mode : %s\n', upper(controller_mode));
fprintf('Simulation Time : %.2f sec\n', t(end));
fprintf('Total Samples   : %d\n', length(t));
fprintf('Hybrid Steps    : %d\n', sim.nSteps);

fprintf('=========================================\n\n');
