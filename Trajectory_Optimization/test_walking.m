% test_walking.m
% Test B-spline walking controller with animation

clear; clc; close all;
if exist('startup.m','file'), startup; end

%% Create controller
ctrl = RabbitController.default();

%% Get robot parameters and initial state
[x0, params, ~] = make_initial_state();

% Adjust phase bounds
ctrl.trajectory.theta0 = x0(3);
ctrl.trajectory.thetaf = x0(3) + 0.2;

%% Get function handle
controller_handle = ctrl.to_function_handle();

%% Simulate
fprintf('Simulating with default B-spline trajectory...\n');
n_steps = 5;
[t_all, x_all, impact_log] = simulate_n_steps(x0, params, n_steps, controller_handle);

%% Fix orientation: simulator returns 14xN, we need Nx14
if size(x_all, 1) == 14
    x_traj = x_all';  % Convert to Nx14
else
    x_traj = x_all;
end

%% ============================================================
%  PLOTS
%  ============================================================

% Figure 1: State trajectories
figure('Name', 'RABBIT B-Spline Walking', 'Position', [50, 50, 1400, 800]);

% Joint angles
subplot(2,3,1);
plot(t_all, x_traj(:, 4:7), 'LineWidth', 1.5);
xlabel('Time (s)'); ylabel('Angle (rad)');
legend('q1 (St Hip)', 'q2 (St Knee)', 'q3 (Sw Hip)', 'q4 (Sw Knee)', 'Location', 'best');
title('Actuated Joint Angles'); grid on;

% Torso (phase variable)
subplot(2,3,2);
plot(t_all, x_traj(:, 3), 'b-', 'LineWidth', 1.5);
hold on;
yline(ctrl.trajectory.theta0, 'r--', '\theta_0');
yline(ctrl.trajectory.thetaf, 'g--', '\theta_f');
xlabel('Time (s)'); ylabel('Torso (rad)');
title('Phase Variable'); grid on;

% Hip trajectory
subplot(2,3,3);
plot(x_traj(:, 1), x_traj(:, 2), 'k-', 'LineWidth', 1.5);
xlabel('x (m)'); ylabel('z (m)');
title('Hip Trajectory'); grid on; axis equal;

% Forward velocity
subplot(2,3,4);
plot(t_all, x_traj(:, 8), 'b-', 'LineWidth', 1.5);
xlabel('Time (s)'); ylabel('v_x (m/s)');
title('Forward Velocity'); grid on;

% Virtual constraint errors
subplot(2,3,5);
y_errors = zeros(length(t_all), 4);
for i = 1:length(t_all)
    [y, ~] = ctrl.trajectory.virtual_constraint(x_traj(i,1:7)', x_traj(i,8:14)');
    y_errors(i,:) = y';
end
plot(t_all, y_errors, 'LineWidth', 1.5);
xlabel('Time (s)'); ylabel('Error (rad)');
legend('y1','y2','y3','y4', 'Location', 'best');
title('Virtual Constraint Errors'); grid on;

% Normalized phase
subplot(2,3,6);
s_traj = zeros(length(t_all), 1);
for i = 1:length(t_all)
    s_traj(i) = ctrl.trajectory.phase(x_traj(i,1:7)');
end
plot(t_all, s_traj, 'k-', 'LineWidth', 1.5);
xlabel('Time (s)'); ylabel('Phase s');
title('Normalized Phase'); grid on;

sgtitle('RABBIT Walking with B-Spline Virtual Constraints');

% Figure 2: Virtual constraints in phase domain
figure('Name', 'Virtual Constraints Tracking', 'Position', [100, 100, 1000, 600]);
joint_names = {'Stance Hip (q1)', 'Stance Knee (q2)', 'Swing Hip (q3)', 'Swing Knee (q4)'};

for j = 1:4
    subplot(2, 2, j);
    
    % Desired
    s_plot = linspace(0, 1, 200);
    hd_plot = zeros(1, 200);
    for k = 1:200
        hd = ctrl.trajectory.evaluate(s_plot(k));
        hd_plot(k) = hd(j);
    end
    plot(s_plot, hd_plot, 'b-', 'LineWidth', 2); hold on;
    
    % Control points
    cp_s = linspace(0, 1, ctrl.trajectory.n+1);
    plot(cp_s, ctrl.trajectory.CP(:, j), 'bo', 'MarkerSize', 10, 'MarkerFaceColor', 'b');
    
    % Actual
    for i = 1:length(t_all)
        s = ctrl.trajectory.phase(x_traj(i,1:7)');
        if s >= 0 && s <= 1
            plot(s, x_traj(i, 3+j), 'r.', 'MarkerSize', 3);
        end
    end
    
    xlabel('Phase s'); ylabel('Angle (rad)');
    title(joint_names{j}); grid on;
    legend('Desired', 'CPs', 'Actual', 'Location', 'best');
end

sgtitle('Virtual Constraints: Desired vs Actual');

%% ============================================================
%  ANIMATION
%  ============================================================

% Filter out bad states (hip below ground, or joint angles > 2*pi)
valid_idx = x_traj(:,2) > 0.3 & max(abs(x_traj(:,4:7)), [], 2) < 2*pi;

if sum(valid_idx) > 10
    fprintf('\nAnimating valid portion of trajectory...\n');
    x_valid = x_traj(valid_idx, :)';
    
    % Slow down animation
    if ~isfield(params, 'speed')
        params.speed = 0.5;
    end
    
    animate_rabbit_stepping_stones(x_valid', params);
else
    fprintf('\nNot enough valid states for animation.\n');
end

%% Summary
fprintf('\n============================================\n');
fprintf('  SIMULATION SUMMARY\n');
fprintf('============================================\n');
fprintf('Steps attempted: %d\n', n_steps);
fprintf('Total time: %.3f s\n', t_all(end));
fprintf('Distance: %.3f m\n', x_traj(end,1) - x_traj(1,1));

% Find where robot fell (hip height < 0.3)
fell_idx = find(x_traj(:,2) < 0.3, 1);
if ~isempty(fell_idx)
    fprintf('Robot fell at t = %.3f s\n', t_all(fell_idx));
end

fprintf('============================================\n');