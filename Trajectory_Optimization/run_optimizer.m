% run_optimizer.m
clear; clc; close all;
if exist('startup.m','file'), startup; end

% First generate initial guess (only needed once)
if ~exist(fullfile('Results', 'initial_guess.mat'), 'file')
    make_initial_guess;
end

% Setup
params = parameters();
ctrl = RabbitController.default();

% Create optimizer with good initial guess
optimizer = RabbitGaitOptimizer(params, ctrl, 0.5, ...
    fullfile('Results', 'initial_guess.mat'));

[CP_opt, fval, exitflag] = optimizer.optimize();

% Save and test if successful
if exitflag > 0
    save(fullfile('Results', 'optimized_gait.mat'), 'CP_opt');
    
    ctrl.trajectory.CP = CP_opt;
    [x0, ~, ~] = make_initial_state();
    handle = ctrl.to_function_handle();
    [t_all, x_all, ~] = simulate_n_steps(x0, params, 5, handle);
    
    if size(x_all,1)==14, x_traj = x_all'; else, x_traj = x_all; end
    valid = x_traj(:,2) > 0.3;
    
    figure;
    subplot(2,2,1); plot(t_all(valid), x_traj(valid,4:7)); title('Joints'); grid on;
    subplot(2,2,2); plot(x_traj(valid,1), x_traj(valid,2)); title('Hip'); grid on;
    subplot(2,2,3); plot(t_all(valid), x_traj(valid,8)); title('Speed'); grid on;
    subplot(2,2,4); ctrl.trajectory.plot(); title('Optimized');
    
    if sum(valid) > 10
        animate_rabbit_stepping_stones(x_traj(valid,:)', params);
    end
end