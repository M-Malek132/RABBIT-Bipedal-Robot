% test_bspline_walking.m
% Full integration test: B-spline virtual constraints on RABBIT

clear; clc; close all;

%% Get initial state and robot parameters
[x0, params, ~] = make_initial_state();

%% Add controller parameters
ctrl = init_bspline_params_v2(x0, params);
params.ctrl = ctrl;

fprintf('\n============================================\n');
fprintf('  B-SPLINE VIRTUAL CONSTRAINT CONTROLLER\n');
fprintf('============================================\n');
fprintf('Robot: mT=%.1f, m1=%.1f, m2=%.1f kg\n', params.mT, params.m1, params.m2);
fprintf('       l1=%.2f, l2=%.2f, lt=%.2f m\n', params.l1, params.l2, params.lt);
fprintf('B-spline: n=%d, p=%d (%d CPs/joint)\n', ctrl.n, ctrl.p, ctrl.n+1);
fprintf('Phase range: [%.3f, %.3f]\n', ctrl.theta0, ctrl.thetaf);
fprintf('PD Gains: Kp=[%d,%d,%d,%d], Kd=[%d,%d,%d,%d]\n', ...
    diag(ctrl.Kp), diag(ctrl.Kd));

%% Controller handle
controller_handle = @(t, x, param) rabbit_controller_bspline(t, x, param.ctrl);

%% Display initial state
fprintf('\nInitial state:\n');
fprintf('  q  = [%.4f, %.4f, %.4f, %.4f, %.4f, %.4f, %.4f]\n', x0(1:7));
fprintf('  dq = [%.4f, %.4f, %.4f, %.4f, %.4f, %.4f, %.4f]\n', x0(8:14));

%% Check initial phase
qt0 = x0(3);
s0 = (qt0 - ctrl.theta0) / (ctrl.thetaf - ctrl.theta0);
fprintf('  Initial phase s = %.4f (qt=%.4f)\n', s0, qt0);

if s0 < 0 || s0 > 1
    fprintf('  WARNING: Initial torso outside phase range!\n');
    fprintf('  Adjusting theta0 to %.4f\n', qt0 - 0.05);
    params.ctrl.theta0 = qt0 - 0.05;
end

%% Simulate
n_steps = 3;
fprintf('\nStarting %d-step simulation...\n', n_steps);

try
    [t_all, x_all, impact_log] = simulate_n_steps(x0, params, n_steps, controller_handle);
    
    if isempty(t_all)
        error('No trajectory data returned.');
    end
    
    fprintf('\n============================================\n');
    fprintf('  SIMULATION COMPLETE\n');
    fprintf('============================================\n');
    
catch ME
    fprintf('\n============================================\n');
    fprintf('  SIMULATION FAILED\n');
    fprintf('============================================\n');
    fprintf('Error ID: %s\n', ME.identifier);
    fprintf('Error msg: %s\n', ME.message);
    
    % Show stack trace
    for i = 1:length(ME.stack)
        fprintf('  In %s (line %d)\n', ME.stack(i).name, ME.stack(i).line);
    end
    
    % Try to show what state caused the error
    fprintf('\nCheck your controller and dynamics functions.\n');
    return;
end

%% ============================================================
%  ANALYSIS
%  ============================================================

% Compute virtual constraint errors
y_errors = zeros(length(t_all), 4);
for i = 1:length(t_all)
    x_i = x_all(:, i);
    s = (x_i(3) - params.ctrl.theta0) / ...
        (params.ctrl.thetaf - params.ctrl.theta0);
    s = min(max(s, 0), 1);
    hd = desired_gait_bspline(s, params.ctrl.ControlPoints, ...
                              params.ctrl.n, params.ctrl.p);
    y_errors(i,:) = x_i(4:7)' - hd';
end

% Compute control torques
u_history = zeros(length(t_all), 4);
for i = 1:length(t_all)
    u_history(i,:) = controller_handle(t_all(i), x_all(:, i), params)';
end

%% ============================================================
%  PLOTS
%  ============================================================

figure('Name', 'RABBIT Walking', 'Position', [50, 50, 1400, 900]);

joint_labels = {'q1 (St Hip)', 'q2 (St Knee)', 'q3 (Sw Hip)', 'q4 (Sw Knee)'};

% Row 1: Joint angles, Torso, Hip XY, Forward velocity
subplot(3,4,1);
plot(t_all, x_all(4:7, :), 'LineWidth', 1.5);
xlabel('Time (s)'); ylabel('Angle (rad)');
legend(joint_labels, 'Location', 'best');
title('Joint Angles'); grid on;

subplot(3,4,2);
plot(t_all, x_all(3, :), 'b-', 'LineWidth', 1.5);
hold on;
yline(params.ctrl.theta0, 'r--');
yline(params.ctrl.thetaf, 'g--');
xlabel('Time (s)'); ylabel('Torso (rad)');
title('Phase Variable (Torso)'); grid on;

subplot(3,4,3);
plot(x_all(1, :), x_all(2, :), 'k-', 'LineWidth', 1.5);
xlabel('x (m)'); ylabel('z (m)');
title('Hip Trajectory'); grid on; axis equal;

subplot(3,4,4);
plot(t_all, x_all(8, :), 'b-', 'LineWidth', 1.5);
xlabel('Time (s)'); ylabel('v_x (m/s)');
title('Forward Velocity'); grid on;

% Row 2: Joint velocities, Phase s, Torso ang vel, Hip height
subplot(3,4,5);
plot(t_all, x_all(11:14, :), 'LineWidth', 1.5);
xlabel('Time (s)'); ylabel('Vel (rad/s)');
title('Joint Velocities'); grid on;

subplot(3,4,6);
s_traj = (x_all(3, :) - params.ctrl.theta0) / ...
         (params.ctrl.thetaf - params.ctrl.theta0);
s_traj = min(max(s_traj, 0), 1);
plot(t_all, s_traj, 'k-', 'LineWidth', 1.5);
xlabel('Time (s)'); ylabel('Phase s');
title('Normalized Phase'); grid on;

subplot(3,4,7);
plot(t_all, x_all(10, :), 'r-', 'LineWidth', 1.5);
xlabel('Time (s)'); ylabel('dqt/dt (rad/s)');
title('Torso Ang Vel'); grid on;

subplot(3,4,8);
plot(t_all, x_all(2, :), 'm-', 'LineWidth', 1.5);
xlabel('Time (s)'); ylabel('z (m)');
title('Hip Height'); grid on;

% Row 3: Virtual constraint errors
for j = 1:4
    subplot(3,4,8+j);
    plot(t_all, y_errors(:,j), 'LineWidth', 1.5);
    xlabel('Time (s)'); ylabel('Error (rad)');
    title(['VC Error: ' joint_labels{j}]); grid on;
end

sgtitle('RABBIT Walking with B-Spline Virtual Constraints');

%% Figure 2: Virtual constraints tracking
figure('Name', 'Virtual Constraints', 'Position', [100, 100, 1000, 600]);

for j = 1:4
    subplot(2, 2, j);
    
    % Desired trajectory
    s_plot = linspace(0, 1, 200);
    hd_plot = zeros(1, 200);
    for k = 1:200
        hd_full = desired_gait_bspline(s_plot(k), params.ctrl.ControlPoints, ...
                                       params.ctrl.n, params.ctrl.p);
        hd_plot(k) = hd_full(j);
    end
    plot(s_plot, hd_plot, 'b-', 'LineWidth', 2); hold on;
    
    % Control points
    cp_s = linspace(0, 1, params.ctrl.n+1);
    plot(cp_s, params.ctrl.ControlPoints(:, j), 'bo', ...
         'MarkerSize', 10, 'MarkerFaceColor', 'b');
    
    % Actual trajectory
    s_actual = (x_all(:, 3) - params.ctrl.theta0) / ...
               (params.ctrl.thetaf - params.ctrl.theta0);
    s_actual = min(max(s_actual, 0), 1);
    q_actual = x_all(:, 3+j);
    plot(s_actual, q_actual, 'r.', 'MarkerSize', 4);
    
    xlabel('Phase s'); ylabel('Angle (rad)');
    title(joint_labels{j});
    legend('Desired', 'CPs', 'Actual', 'Location', 'best');
    grid on;
end

sgtitle('Virtual Constraints: Desired vs Actual');

%% Figure 3: Torques
figure('Name', 'Control Torques', 'Position', [150, 150, 800, 400]);

for j = 1:4
    subplot(2, 2, j);
    plot(t_all, u_history(:,j), 'LineWidth', 1.5);
    xlabel('Time (s)'); ylabel('Torque (Nm)');
    title(['u_' num2str(j) ': ' joint_labels{j}]);
    grid on;
end

sgtitle('Control Torques');

%% Summary
fprintf('\n============================================\n');
fprintf('  SIMULATION SUMMARY\n');
fprintf('============================================\n');

if ~isempty(impact_log)
    fprintf('Steps completed:   %d\n', length(impact_log));
    
    % Filter valid impacts (non-empty)
    valid_impacts = 0;
    for i = 1:length(impact_log)
        if ~isempty(impact_log(i).impact_time)
            valid_impacts = valid_impacts + 1;
        end
    end
    fprintf('Valid impacts:     %d\n', valid_impacts);
end

fprintf('Total time:        %.3f s\n', t_all(end));
fprintf('Distance traveled: %.3f m\n', x_all(end,1) - x_all(1,1));
fprintf('Avg speed:         %.3f m/s\n', ...
    (x_all(end,1) - x_all(1,1)) / t_all(end));

fprintf('\nVirtual Constraint RMS Errors:\n');
for j = 1:4
    fprintf('  %s: %.4f rad\n', joint_labels{j}, rms(y_errors(:,j)));
end

fprintf('\nTorque Stats (max |u|):\n');
for j = 1:4
    fprintf('  u_%d: %.1f Nm\n', j, max(abs(u_history(:,j))));
end

fprintf('============================================\n');