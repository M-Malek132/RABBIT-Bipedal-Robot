% test_gait_library.m
% Load gait library, visualize, interpolate, and test walking
% All outputs saved to Results/

clear; clc; close all;
if exist('startup.m','file'), startup; end

%% Setup
results_dir = fullfile(pwd, 'Results');
if ~exist(results_dir, 'dir')
    mkdir(results_dir);
end

%% Load library
library_file = fullfile(results_dir, 'gait_library.mat');

if exist(library_file, 'file')
    load(library_file);
    fprintf('Loaded gait library with %d gaits.\n\n', length(gait_library));
    
    % Print summary (handle missing fields)
    fprintf('%-6s %-8s %-8s %-8s', 'Speed', 'T_step', 'Length', 'Actual');
    if isfield(gait_library, 'fit_error')
        fprintf(' %-8s', 'Fit Err');
    end
    fprintf('\n');
    fprintf('%-6s %-8s %-8s %-8s', '------', '------', '------', '------');
    if isfield(gait_library, 'fit_error')
        fprintf(' %-8s', '------');
    end
    fprintf('\n');
    
    for i = 1:length(gait_library)
        fprintf('%-6.2f %-8.3f %-8.3f %-8.3f', ...
            gait_library(i).speed, ...
            gait_library(i).T_step, ...
            gait_library(i).step_length, ...
            gait_library(i).actual_speed);
        if isfield(gait_library, 'fit_error')
            fprintf(' %-8.4f', gait_library(i).fit_error);
        end
        fprintf('\n');
    end
else
    fprintf('No gait library found. Run generate_gait_library first.\n');
    return;
end

%% Plot 1: Gait Library
fig1 = figure('Name', 'Gait Library', 'Position', [50, 50, 1200, 800]);
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

saveas(fig1, fullfile(results_dir, 'gait_library.png'));
fprintf('Saved: gait_library.png\n');

%% Interpolation
speeds = [gait_library.actual_speed];
v_desired = mean([min(speeds), max(speeds)]);

[~, idx_low] = max(speeds(speeds <= v_desired));
[~, idx_high] = min(speeds(speeds >= v_desired));

if isempty(idx_low) || isempty(idx_high) || idx_low == idx_high
    idx_low = 1;
    idx_high = min(2, length(gait_library));
end

alpha = (v_desired - speeds(idx_low)) / (speeds(idx_high) - speeds(idx_low));
alpha = max(0, min(1, alpha));
CP_interp = (1-alpha)*gait_library(idx_low).CP + alpha*gait_library(idx_high).CP;

fprintf('\nInterpolation:\n');
fprintf('  %.2f (idx %d) + %.2f (idx %d) -> %.2f m/s (alpha=%.2f)\n', ...
    speeds(idx_low), idx_low, speeds(idx_high), idx_high, v_desired, alpha);

%% Plot 2: Interpolation
fig2 = figure('Name', 'Interpolation', 'Position', [100, 100, 1000, 500]);

for j = 1:4
    subplot(2, 2, j); hold on;
    s_plot = linspace(0, 1, 100);
    n = gait_library(1).n;
    p = gait_library(1).p;
    
    hd_low = zeros(1, 100); hd_high = zeros(1, 100); hd_interp = zeros(1, 100);
    for k = 1:100
        N = BSpline(n, p, s_plot(k));
        hd_low(k) = N * gait_library(idx_low).CP(:, j);
        hd_high(k) = N * gait_library(idx_high).CP(:, j);
        hd_interp(k) = N * CP_interp(:, j);
    end
    
    plot(s_plot, hd_low, 'b-', 'LineWidth', 2);
    plot(s_plot, hd_high, 'r-', 'LineWidth', 2);
    plot(s_plot, hd_interp, 'k--', 'LineWidth', 2);
    xlabel('Phase s'); ylabel('Angle (rad)');
    title(joint_names{j}); grid on;
end
sgtitle(sprintf('Interpolation at %.2f m/s', v_desired));
legend('Low', 'High', 'Interp', 'Location', 'best');

saveas(fig2, fullfile(results_dir, 'gait_interpolation.png'));
fprintf('Saved: gait_interpolation.png\n');

%% Walking Test
fprintf('\n============================================\n');
fprintf('  WALKING TEST at %.2f m/s\n', v_desired);
fprintf('============================================\n');

traj = BSplineTrajectory(gait_library(1).n, gait_library(1).p, ...
                          gait_library(1).theta0, gait_library(1).thetaf);
traj.CP = CP_interp;

Kp = diag([200, 200, 150, 150]);
Kd = diag([30,  30,  20,  20 ]);
ctrl = RabbitController(traj, Kp, Kd);

[x0, params, ~] = make_initial_state();
ctrl.trajectory.theta0 = x0(3);
ctrl.trajectory.thetaf = x0(3) + 0.2;

controller_handle = ctrl.to_function_handle();

[t_all, x_all, impact_log] = simulate_n_steps(x0, params, 5, controller_handle);

if size(x_all, 1) == 14
    x_traj = x_all';
else
    x_traj = x_all;
end

valid = x_traj(:,2) > 0.3;
t_valid = t_all(valid);
x_valid = x_traj(valid, :);

%% Plot 3: Walking Results
fig3 = figure('Name', 'Walking Test', 'Position', [150, 150, 1200, 800]);

subplot(2,3,1);
plot(t_valid, x_valid(:,4:7), 'LineWidth', 1.5);
xlabel('Time (s)'); ylabel('Angle (rad)');
legend('q1','q2','q3','q4'); title('Joint Angles'); grid on;

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
title('VC Errors'); grid on;

subplot(2,3,6);
traj.plot();
title('B-Spline Constraints');

sgtitle(sprintf('Walking Test at %.2f m/s', v_desired));

saveas(fig3, fullfile(results_dir, 'walking_test.png'));
fprintf('Saved: walking_test.png\n');

%% Save walking data
timestamp = datestr(now, 'yyyy-mm-dd_HH-MM-SS');
walking_data.t_all = t_all;
walking_data.x_all = x_all;
walking_data.v_desired = v_desired;
walking_data.CP_interp = CP_interp;
save(fullfile(results_dir, ['walking_test_', timestamp, '.mat']), 'walking_data');
fprintf('Saved: walking_test_%s.mat\n', timestamp);

%% Summary
fprintf('\n============================================\n');
fprintf('  RESULTS SAVED TO Results/\n');
fprintf('============================================\n');
fprintf('  gait_library.png\n');
fprintf('  gait_interpolation.png\n');
fprintf('  walking_test.png\n');
fprintf('  walking_test_%s.mat\n', timestamp);
fprintf('============================================\n');

%% Animate
if size(x_valid, 1) > 10
    fprintf('\nAnimating...\n');
    animate_rabbit_stepping_stones(x_valid, params);
end