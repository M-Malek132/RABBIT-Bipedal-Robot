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

try
    [t_all, x_all, impact_log] = simulate_n_steps(x0, params, n_steps, controller_handle);
catch ME
    fprintf('Simulation crashed: %s\n', ME.message);
    
    % Try to use whatever data we got
    if exist('t_all', 'var') && ~isempty(t_all)
        fprintf('Using partial trajectory data...\n');
    else
        fprintf('No trajectory data available.\n');
        return;
    end
end

%% Fix orientation: simulator returns 14xN or Nx14
if size(x_all, 1) == 14 && size(x_all, 2) > 14
    x_traj = x_all';  % Convert to Nx14
else
    x_traj = x_all;
end

% Filter valid states (hip above ground, reasonable angles)
valid = x_traj(:,2) > 0.3 & max(abs(x_traj(:,4:7)), [], 2) < 5;
x_valid = x_traj(valid, :);
t_valid = t_all(valid);

if isempty(x_valid)
    fprintf('No valid states for plotting.\n');
    return;
end

%% ============================================================
%  PLOTS
%  ============================================================

figure('Name', 'RABBIT B-Spline Walking', 'Position', [50, 50, 1400, 800]);

% Joint angles
subplot(2,3,1);
plot(t_valid, x_valid(:, 4:7), 'LineWidth', 1.5);
xlabel('Time (s)'); ylabel('Angle (rad)');
legend('q1', 'q2', 'q3', 'q4', 'Location', 'best');
title('Joint Angles'); grid on;

% Torso (phase variable)
subplot(2,3,2);
plot(t_valid, x_valid(:, 3), 'b-', 'LineWidth', 1.5);
hold on;
yline(ctrl.trajectory.theta0, 'r--');
yline(ctrl.trajectory.thetaf, 'g--');
xlabel('Time (s)'); ylabel('Torso (rad)');
title('Phase Variable'); grid on;

% Hip trajectory
subplot(2,3,3);
plot(x_valid(:, 1), x_valid(:, 2), 'k-', 'LineWidth', 1.5);
xlabel('x (m)'); ylabel('z (m)');
title('Hip Trajectory'); grid on; axis equal;

% Forward velocity
subplot(2,3,4);
plot(t_valid, x_valid(:, 8), 'b-', 'LineWidth', 1.5);
xlabel('Time (s)'); ylabel('v_x (m/s)');
title('Forward Speed'); grid on;

% Virtual constraint errors
subplot(2,3,5);
y_err = zeros(length(t_valid), 4);
for i = 1:length(t_valid)
    [y, ~] = ctrl.trajectory.virtual_constraint(x_valid(i,1:7)', x_valid(i,8:14)');
    y_err(i,:) = y';
end
plot(t_valid, y_err, 'LineWidth', 1.5);
xlabel('Time (s)'); ylabel('Error (rad)');
title('VC Errors'); grid on;

% Normalized phase
subplot(2,3,6);
s_traj = zeros(length(t_valid), 1);
for i = 1:length(t_valid)
    s_traj(i) = ctrl.trajectory.phase(x_valid(i,1:7)');
end
plot(t_valid, s_traj, 'k-', 'LineWidth', 1.5);
xlabel('Time (s)'); ylabel('Phase s');
title('Normalized Phase'); grid on;

sgtitle('RABBIT with B-Spline Virtual Constraints');

%% ============================================================
%  ANIMATION
%  ============================================================

fprintf('\nAnimating valid portion of trajectory...\n');

% Use only valid states
x_anim = x_valid';

% Slow down animation
if ~isfield(params, 'speed')
    params.speed = 0.3;
end

try
    animate_rabbit_stepping_stones(x_anim', params);
catch ME
    fprintf('Animation error: %s\n', ME.message);
end

%% Summary
fprintf('\n============================================\n');
fprintf('  SIMULATION SUMMARY\n');
fprintf('============================================\n');
fprintf('Total time: %.3f s\n', t_all(end));
fprintf('Distance: %.3f m\n', x_traj(end,1) - x_traj(1,1));

if size(x_valid,1) < size(x_traj,1)
    fprintf('Filtered %d/%d invalid states\n', ...
        size(x_traj,1)-size(x_valid,1), size(x_traj,1));
end
fprintf('============================================\n');