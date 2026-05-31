% test_gait_library.m
% Load gait library from Results/, visualize, interpolate, and test

clear; clc; close all;
if exist('startup.m','file'), startup; end

results_dir = fullfile(pwd, 'Results');

%% Load library
library_file = fullfile(results_dir, 'gait_library.mat');

if exist(library_file, 'file')
    load(library_file);
    fprintf('Loaded gait library from Results/ with %d gaits.\n\n', length(gait_library));
    
    % Print library summary
    fprintf('%-8s %-10s %-10s %-10s %-10s\n', 'Speed', 'T_step', 'Length', 'Actual', 'Fit Err');
    fprintf('%-8s %-10s %-10s %-10s %-10s\n', '------', '------', '------', '------', '------');
    for i = 1:length(gait_library)
        fprintf('%-8.2f %-10.3f %-10.3f %-10.3f %-10.4f\n', ...
            gait_library(i).speed, ...
            gait_library(i).T_step, ...
            gait_library(i).step_length, ...
            gait_library(i).actual_speed, ...
            gait_library(i).fit_error);
    end
    
else
    fprintf('No gait library found. Generating...\n');
    gait_library = generate_gait_library();
    if isempty(gait_library)
        error('Failed to generate gait library.');
    end
end

if length(gait_library) < 2
    error('Need at least 2 gaits. Only have %d.', length(gait_library));
end

%% ============================================================
%  PLOT GAIT LIBRARY
%  ============================================================
figure('Name', 'Gait Library', 'Position', [50, 50, 1200, 800]);
joint_names = {'Stance Hip (q1)', 'Stance Knee (q2)', ...
               'Swing Hip (q3)', 'Swing Knee (q4)'};
colors = lines(length(gait_library));

for j = 1:4
    subplot(2, 2, j); hold on;
    for i = 1:length(gait_library)
        CP = gait_library(i).CP;
        n = gait_library(i).n;
        p = gait_library(i).p;
        
        s_plot = linspace(0, 1, 100);
        hd_plot = zeros(1, 100);
        for k = 1:100
            N = BSpline(n, p, s_plot(k));
            hd_plot(k) = N * CP(:, j);
        end
        
        plot(s_plot, hd_plot, 'Color', colors(i,:), 'LineWidth', 2, ...
             'DisplayName', sprintf('%.2f m/s', gait_library(i).actual_speed));
    end
    xlabel('Phase s'); ylabel('Angle (rad)');
    title(joint_names{j}); grid on;
    if j == 1, legend('Location', 'best'); end
end
sgtitle('Gait Library: B-Spline Virtual Constraints');

% Save figure
saveas(gcf, fullfile(results_dir, 'gait_library.png'));

%% ============================================================
%  INTERPOLATION TEST
%  ============================================================
fprintf('\n============================================\n');
fprintf('  GAIT INTERPOLATION\n');
fprintf('============================================\n');

speeds = [gait_library.actual_speed];
v_min = min(speeds);
v_max = max(speeds);

% Test interpolation at mid-speed
v_desired = (v_min + v_max) / 2;

% Find bracketing gaits
[~, idx_low] = max(speeds(speeds <= v_desired));
[~, idx_high] = min(speeds(speeds >= v_desired));

if isempty(idx_low), idx_low = idx_high; end
if isempty(idx_high), idx_high = idx_low; end

alpha = (v_desired - speeds(idx_low)) / (speeds(idx_high) - speeds(idx_low));
alpha = max(0, min(1, alpha));  % clamp

CP_interp = (1-alpha)*gait_library(idx_low).CP + alpha*gait_library(idx_high).CP;

fprintf('Interpolating between:\n');
fprintf('  Low:  %.2f m/s (idx %d)\n', speeds(idx_low), idx_low);
fprintf('  High: %.2f m/s (idx %d)\n', speeds(idx_high), idx_high);
fprintf('  Target: %.2f m/s (alpha = %.2f)\n', v_desired, alpha);

% Plot interpolation
figure('Name', 'Gait Interpolation', 'Position', [100, 100, 1000, 500]);

for j = 1:4
    subplot(2, 2, j); hold on;
    
    s_plot = linspace(0, 1, 100);
    n = gait_library(1).n;
    p = gait_library(1).p;
    
    % Low speed
    hd_low = zeros(1, 100);
    for k = 1:100
        N = BSpline(n, p, s_plot(k));
        hd_low(k) = N * gait_library(idx_low).CP(:, j);
    end
    
    % High speed
    hd_high = zeros(1, 100);
    for k = 1:100
        N = BSpline(n, p, s_plot(k));
        hd_high(k) = N * gait_library(idx_high).CP(:, j);
    end
    
    % Interpolated
    hd_interp = zeros(1, 100);
    for k = 1:100
        N = BSpline(n, p, s_plot(k));
        hd_interp(k) = N * CP_interp(:, j);
    end
    
    plot(s_plot, hd_low, 'b-', 'LineWidth', 2);
    plot(s_plot, hd_high, 'r-', 'LineWidth', 2);
    plot(s_plot, hd_interp, 'k--', 'LineWidth', 2);
    
    xlabel('Phase s'); ylabel('Angle (rad)');
    title(joint_names{j}); grid on;
    
    if j == 1
        legend(sprintf('%.2f m/s', speeds(idx_low)), ...
               sprintf('%.2f m/s', speeds(idx_high)), ...
               sprintf('%.2f m/s (interp)', v_desired), ...
               'Location', 'best');
    end
end
sgtitle(sprintf('Gait Interpolation at %.2f m/s', v_desired));

% Save figure
saveas(gcf, fullfile(results_dir, 'gait_interpolation.png'));

%% ============================================================
%  WALKING TEST WITH INTERPOLATED GAIT
%  ============================================================
fprintf('\n============================================\n');
fprintf('  WALKING TEST\n');
fprintf('============================================\n');

% Create trajectory
traj = BSplineTrajectory(gait_library(1).n, gait_library(1).p, ...
                          gait_library(1).theta0, gait_library(1).thetaf);
traj.CP = CP_interp;

% Create controller
Kp = diag([200, 200, 150, 150]);
Kd = diag([30,  30,  20,  20 ]);
ctrl = RabbitController(traj, Kp, Kd);

% Get initial state
[x0, params, ~] = make_initial_state();
ctrl.trajectory.theta0 = x0(3);
ctrl.trajectory.thetaf = x0(3) + 0.2;

controller_handle = ctrl.to_function_handle();

% Simulate
[t_all, x_all, impact_log] = simulate_n_steps(x0, params, 5, controller_handle);

% Process data
if size(x_all, 1) == 14
    x_traj = x_all';
else
    x_traj = x_all;
end

valid = x_traj(:,2) > 0.3 & max(abs(x_traj(:,4:7)), [], 2) < 5;
t_valid = t_all(valid);
x_valid = x_traj(valid, :);

if isempty(x_valid)
    fprintf('No valid states for plotting.\n');
    return;
end

% Plot walking results
figure('Name', 'Walking Test', 'Position', [150, 150, 1200, 800]);

subplot(2,3,1);
plot(t_valid, x_valid(:,4:7), 'LineWidth', 1.5);
xlabel('Time (s)'); ylabel('Angle (rad)');
legend('q1','q2','q3','q4', 'Location', 'best');
title('Joint Angles'); grid on;

subplot(2,3,2);
plot(t_valid, x_valid(:,3), 'LineWidth', 1.5);
xlabel('Time (s)'); ylabel('Torso (rad)');
title('Phase Variable'); grid on;

subplot(2,3,3);
plot(x_valid(:,1), x_valid(:,2), 'LineWidth', 1.5);
xlabel('x (m)'); ylabel('z (m)');
title('Hip Trajectory'); grid on; axis equal;

subplot(2,3,4);
plot(t_valid, x_valid(:,8), 'LineWidth', 1.5);
xlabel('Time (s)'); ylabel('v_x (m/s)');
title('Forward Speed'); grid on;

subplot(2,3,5);
y_err = zeros(length(t_valid), 4);
for i = 1:length(t_valid)
    [y, ~] = ctrl.trajectory.virtual_constraint(x_valid(i,1:7)', x_valid(i,8:14)');
    y_err(i,:) = y';
end
plot(t_valid, y_err, 'LineWidth', 1.5);
xlabel('Time (s)'); ylabel('Error (rad)');
title('VC Tracking Errors'); grid on;

subplot(2,3,6);
traj.plot();
title('B-Spline Virtual Constraints');

sgtitle(sprintf('Walking Test at %.2f m/s (Interpolated)', v_desired));

% Save figure
saveas(gcf, fullfile(results_dir, 'walking_test.png'));

% Save walking data
timestamp = datestr(now, 'yyyy-mm-dd_HH-MM-SS');
walking_data.t_all = t_all;
walking_data.x_all = x_all;
walking_data.v_desired = v_desired;
walking_data.CP_interp = CP_interp;
save(fullfile(results_dir, ['walking_test_', timestamp, '.mat']), 'walking_data');

%% Summary
fprintf('\n============================================\n');
fprintf('  RESULTS SUMMARY\n');
fprintf('============================================\n');
fprintf('Target speed: %.2f m/s\n', v_desired);
fprintf('Steps completed: %d\n', length(impact_log));
fprintf('Total time: %.3f s\n', t_all(end));
fprintf('Distance: %.3f m\n', x_traj(end,1) - x_traj(1,1));
fprintf('Avg speed: %.3f m/s\n', (x_traj(end,1) - x_traj(1,1))/t_all(end));
fprintf('\nAll results saved to Results/\n');
fprintf('============================================\n');

%% Animate
if size(x_valid, 1) > 10
    fprintf('\nStarting animation...\n');
    params.speed = 0.3;
    animate_rabbit_stepping_stones(x_valid, params);
end