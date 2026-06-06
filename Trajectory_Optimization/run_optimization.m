% run_optimization.m
% Find periodic walking gait via optimization

clear; clc; close all;
if exist('startup.m','file'), startup; end

%% Setup
[x0, params, ~] = make_initial_state();
ctrl = RabbitController.default();

% Adjust phase bounds
ctrl.trajectory.theta0 = x0(3);
ctrl.trajectory.thetaf = x0(3) + 0.2;

%% Create optimizer
target_speed = 0.5;  % m/s
optimizer = RabbitGaitOptimizer(params, ctrl, target_speed);

%% Run optimization
fprintf('============================================\n');
fprintf('  OPTIMIZING PERIODIC GAIT\n');
fprintf('  Target speed: %.2f m/s\n', target_speed);
fprintf('  Variables: %d\n', ctrl.trajectory.num_variables());
fprintf('============================================\n\n');

[CP_opt, fval, exitflag] = optimizer.optimize();

%% Save result
save('optimized_gait.mat', 'CP_opt', 'fval', 'exitflag', 'target_speed');
fprintf('\nSaved to optimized_gait.mat\n');

%% Test optimized gait
if exitflag > 0
    fprintf('\n============================================\n');
    fprintf('  TESTING OPTIMIZED GAIT\n');
    fprintf('============================================\n');
    
    ctrl.trajectory.CP = CP_opt;
    controller_handle = ctrl.to_function_handle();
    
    [t_all, x_all, ~] = simulate_n_steps(x0, params, 5, controller_handle);
    
    % Convert if needed
    if size(x_all, 1) == 14
        x_traj = x_all';
    else
        x_traj = x_all;
    end
    
    % Plot
    figure('Name', 'Optimized Gait');
    subplot(2,2,1);
    plot(t_all, x_traj(:,4:7)); xlabel('Time'); ylabel('Angle');
    title('Joint Angles'); grid on; legend('q1','q2','q3','q4');
    
    subplot(2,2,2);
    plot(x_traj(:,1), x_traj(:,2)); xlabel('x'); ylabel('z');
    title('Hip Trajectory'); grid on;
    
    subplot(2,2,3);
    plot(t_all, x_traj(:,8)); xlabel('Time'); ylabel('v_x');
    title('Forward Speed'); grid on;
    
    subplot(2,2,4);
    ctrl.trajectory.plot();
    title('Optimized B-Spline Trajectories');
    
    % Animate if valid
    valid_idx = x_traj(:,2) > 0.3;
    if sum(valid_idx) > 10
        x_valid = x_traj(valid_idx, :)';
        params.speed = 0.5;
        animate_rabbit_stepping_stones(x_valid, params);
    end
    
    % Summary
    fprintf('\nFinal speed: %.3f m/s\n', (x_traj(end,1)-x_traj(1,1))/t_all(end));
    fprintf('Distance: %.3f m\n', x_traj(end,1)-x_traj(1,1));
end