% test_bspline_controller_full.m
% Full integration test with RABBIT simulator

clear; clc; close all;

%% Initialize B-spline controller parameters
params_ctrl = init_bspline_params();

% Adjust phase bounds based on expected step
params_ctrl.theta0 = -0.15;  
params_ctrl.thetaf =  0.25;

% Adjust PD gains (may need tuning)
params_ctrl.Kp = diag([300, 300, 200, 200]);
params_ctrl.Kd = diag([40,  40,  25,  25 ]);

fprintf('============================================\n');
fprintf('  B-SPLINE VIRTUAL CONSTRAINT CONTROLLER\n');
fprintf('============================================\n');
fprintf('B-spline: n=%d, p=%d (%d control points/joint)\n', ...
    params_ctrl.n, params_ctrl.p, params_ctrl.n+1);
fprintf('Phase range: [%.3f, %.3f]\n', params_ctrl.theta0, params_ctrl.thetaf);
fprintf('Gains: Kp=[%d,%d,%d,%d], Kd=[%d,%d,%d,%d]\n\n', ...
    diag(params_ctrl.Kp), diag(params_ctrl.Kd));

%% Load robot parameters (use your existing params)
% Assuming you have a function or script that loads robot params
% params = load_rabbit_params();  % <-- replace with your actual call
% For now, we need to merge params_ctrl into your params structure

% TEMPORARY: Create a minimal params structure
% REPLACE THIS with your actual robot parameters
params = struct();
params.g = 9.81;
% Add all your required fields here...
% We'll use packParameters(params) so it must have all robot fields

% Merge controller params into params structure
params.ctrl = params_ctrl;

%% Define controller handle
% This matches your simulator's calling convention:
% controller(t, x, param)
controller_handle = @(t, x, param) rabbit_controller_bspline(t, x, param.ctrl);

%% Initial state
% Start of step configuration
q0 = [
    0.0;                % x: hip horizontal position
    0.85;               % z: hip vertical position  
    params_ctrl.theta0; % qt: torso angle = theta0 (start of step)
    -0.3;               % q1: stance knee (bent)
    0.4;                % q2: stance hip
    -0.2;               % q3: swing knee (bent behind)
    0.5                 % q4: swing hip
];

dq0 = [
    0.8;    % dx: forward velocity
    0.0;    % dz    0.15;   % dqt: torso angular velocity
    -0.5;   % dq1
    1.0;    % dq2
    -0.5;   % dq3
    1.0     % dq4
];

x0 = [q0; dq0];

fprintf('Initial state:\n');
fprintf('  q  = [%.3f, %.3f, %.3f, %.3f, %.3f, %.3f, %.3f]\n', q0);
fprintf('  dq = [%.3f, %.3f, %.3f, %.3f, %.3f, %.3f, %.3f]\n\n', dq0);

%% Simulate
n_steps = 3;
fprintf('Starting %d-step simulation...\n\n', n_steps);

try
    [t_all, x_all, impact_log] = simulate_n_steps(x0, params, n_steps, controller_handle);
    
    fprintf('\n============================================\n');
    fprintf('  SIMULATION COMPLETE\n');
    fprintf('============================================\n');
    
catch ME
    fprintf('\n============================================\n');
    fprintf('  SIMULATION FAILED\n');
    fprintf('============================================\n');
    fprintf('Error: %s\n', ME.message);
    fprintf('In: %s (line %d)\n', ME.stack(1).name, ME.stack(1).line);
    return;
end

%% ============================================================
%  ANALYSIS AND PLOTTING
%  ============================================================

if isempty(t_all)
    fprintf('No data to plot.\n');
    return;
end

%% Figure 1: State trajectories
figure('Name', 'RABBIT B-Spline Controller - States', 'Position', [50, 50, 1200, 800]);

% Joint angles
subplot(3,3,1);
plot(t_all, x_all(:, 4:7), 'LineWidth', 1.5);
xlabel('Time (s)'); ylabel('Angle (rad)');
legend('q1 (st knee)', 'q2 (st hip)', 'q3 (sw knee)', 'q4 (sw hip)', 'Location', 'best');
title('Actuated Joint Angles');
grid on;

% Torso angle (phase variable)
subplot(3,3,2);
plot(t_all, x_all(:, 3), 'b-', 'LineWidth', 1.5);
hold on;
yline(params_ctrl.theta0, 'r--', '\theta_0');
yline(params_ctrl.thetaf, 'g--', '\theta_f');
xlabel('Time (s)'); ylabel('Torso Angle (rad)');
title('Phase Variable (Torso)');
grid on;

% Hip trajectory
subplot(3,3,3);
plot(x_all(:, 1), x_all(:, 2), 'k-', 'LineWidth', 1.5);
xlabel('x (m)'); ylabel('z (m)');
title('Hip Trajectory');
grid on; axis equal;

% Joint velocities
subplot(3,3,4);
plot(t_all, x_all(:, 11:14), 'LineWidth', 1.5);
xlabel('Time (s)'); ylabel('Velocity (rad/s)');
legend('dq1', 'dq2', 'dq3', 'dq4', 'Location', 'best');
title('Actuated Joint Velocities');
grid on;

% Forward velocity
subplot(3,3,5);
v_forward = x_all(:, 8);
plot(t_all, v_forward, 'b-', 'LineWidth', 1.5);
xlabel('Time (s)'); ylabel('Velocity (m/s)');
title('Forward Velocity');
grid on;

% Torso angular velocity
subplot(3,3,6);
plot(t_all, x_all(:, 10), 'r-', 'LineWidth', 1.5);
xlabel('Time (s)'); ylabel('Angular Velocity (rad/s)');
title('Torso Angular Velocity');
grid on;

% Phase variable over time
subplot(3,3,7);
s_traj = (x_all(:, 3) - params_ctrl.theta0) / (params_ctrl.thetaf - params_ctrl.theta0);
s_traj = min(max(s_traj, 0), 1);
plot(t_all, s_traj, 'k-', 'LineWidth', 1.5);
xlabel('Time (s)'); ylabel('Phase s');
title('Normalized Phase Variable');
grid on;

% Virtual constraint errors
subplot(3,3,8);
y_errors = zeros(size(t_all,1), 4);
for i = 1:length(t_all)
    x_i = x_all(i, :)';
    s = (x_i(3) - params_ctrl.theta0) / (params_ctrl.thetaf - params_ctrl.theta0);
    s = min(max(s, 0), 1);
    hd = desired_gait_bspline(s, params_ctrl.ControlPoints, params_ctrl.n, params_ctrl.p);
    y_errors(i,:) = x_i(4:7)' - hd';
end
plot(t_all, y_errors, 'LineWidth', 1.5);
xlabel('Time (s)'); ylabel('Error (rad)');
legend('y1', 'y2', 'y3', 'y4', 'Location', 'best');
title('Virtual Constraint Errors');
grid on;

% Control torques
subplot(3,3,9);
u_history = zeros(length(t_all), 4);
for i = 1:length(t_all)
    u_history(i,:) = controller_handle(t_all(i), x_all(i,:)', params)';
end
plot(t_all, u_history, 'LineWidth', 1.5);
xlabel('Time (s)'); ylabel('Torque (Nm)');
legend('u1', 'u2', 'u3', 'u4', 'Location', 'best');
title('Control Torques');
grid on;

sgtitle('RABBIT Walking with B-Spline Virtual Constraints');

%% Figure 2: Virtual constraints in phase domain
figure('Name', 'Virtual Constraints - Phase Domain', 'Position', [100, 100, 1000, 600]);

joint_names = {'Stance Knee (q1)', 'Stance Hip (q2)', 'Swing Knee (q3)', 'Swing Hip (q4)'};
ylabels = {'q1 (rad)', 'q2 (rad)', 'q3 (rad)', 'q4 (rad)'};

for j = 1:4
    subplot(2, 2, j);
    
    % Plot desired trajectory
    s_plot = linspace(0, 1, 100);
    hd_plot = zeros(1, 100);
    for k = 1:100
        hd_plot(k) = desired_gait_bspline(s_plot(k), params_ctrl.ControlPoints, params_ctrl.n, params_ctrl.p);
        hd_plot(k) = hd_plot(j);  % Extract joint j
    end
    plot(s_plot, hd_plot, 'b-', 'LineWidth', 2);
    hold on;
    
    % Plot control points
    cp_s = linspace(0, 1, params_ctrl.n+1);
    plot(cp_s, params_ctrl.ControlPoints(:, j), 'bo', 'MarkerSize', 8, 'MarkerFaceColor', 'b');
    
    % Plot actual trajectory
    s_actual = (x_all(:, 3) - params_ctrl.theta0) / (params_ctrl.thetaf - params_ctrl.theta0);
    s_actual = min(max(s_actual, 0), 1);
    q_actual = x_all(:, 3+j);
    plot(s_actual, q_actual, 'r.', 'MarkerSize', 3);
    
    xlabel('Phase s'); ylabel(ylabels{j});
    title(joint_names{j});
    legend('Desired', 'Control Points', 'Actual', 'Location', 'best');
    grid on;
end

sgtitle('Virtual Constraints: Desired vs Actual');

%% Figure 3: Animation-like stick figure at key frames
figure('Name', 'Key Poses', 'Position', [150, 150, 800, 400]);

% Sample a few key poses
n_poses = 5;
indices = round(linspace(1, length(t_all), n_poses));

for i = 1:n_poses
    subplot(1, n_poses, i);
    % Call your drawing function if available
    % draw_rabbit(x_all(indices(i), 1:7)', params);
    % For now, just show the joint angles as text
    q_i = x_all(indices(i), 1:7)';
    bar(q_i(4:7));
    title(sprintf('t=%.2f s', t_all(indices(i))));
    xlabel('Joint'); ylabel('Angle (rad)');
    set(gca, 'XTickLabel', {'q1','q2','q3','q4'});
    ylim([-1, 1]);
end

sgtitle('Key Poses During Walking');

%% Print summary
fprintf('\n============================================\n');
fprintf('  SIMULATION SUMMARY\n');
fprintf('============================================\n');
fprintf('Steps completed: %d\n', length(impact_log));
fprintf('Total time: %.3f s\n', t_all(end));
fprintf('Distance traveled: %.3f m\n', x_all(end,1) - x_all(1,1));
fprintf('Avg forward speed: %.3f m/s\n', (x_all(end,1) - x_all(1,1)) / t_all(end));

if ~isempty(impact_log)
    fprintf('\nImpact details:\n');
    for i = 1:min(length(impact_log), 5)
        if ~isempty(impact_log(i).impact_time)
            fprintf('  Step %d: t=%.3f s\n', i, impact_log(i).impact_time);
        end
    end
end

% Check periodicity
fprintf('\nPeriodicity check (compare first and last state):\n');
if length(impact_log) >= 2
    q_first = x_all(1, 1:7)';
    q_last = x_all(end, 1:7)';
    fprintf('  First q: [%.3f, %.3f, %.3f, %.3f, %.3f, %.3f, %.3f]\n', q_first);
    fprintf('  Last  q: [%.3f, %.3f, %.3f, %.3f, %.3f, %.3f, %.3f]\n', q_last);
    
    % After leg swap, compare appropriate joints
    % q1_first should match q3_last (swing becomes stance)
    % q2_first should match q4_last
    % q3_first should match q1_last
    % q4_first should match q2_last
    fprintf('  |q1_first - q3_last| = %.4f\n', abs(q_first(4) - q_last(6)));
    fprintf('  |q2_first - q4_last| = %.4f\n', abs(q_first(5) - q_last(7)));
    fprintf('  |q3_first - q1_last| = %.4f\n', abs(q_first(6) - q_last(4)));
    fprintf('  |q4_first - q2_last| = %.4f\n', abs(q_first(7) - q_last(5)));
end

fprintf('============================================\n');