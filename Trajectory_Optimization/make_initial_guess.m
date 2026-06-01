% make_initial_guess.m - Using x_hip as phase variable
clear; clc;
if exist('startup.m','file'), startup; end

params = parameters();

%% Setup
nq = 7; x0 = zeros(2*nq, 1);
qt=0.1; q1=-0.3; q2=0.6; q3=-1.0; q4=0.6;
l1=params.l1; l2=params.l2;
px=l1*sin(qt+q1)+l2*sin(qt+q1+q2);
pz=l1*cos(qt+q1)+l2*cos(qt+q1+q2);
x0(1:nq)=[px;pz;qt;q1;q2;q3;q4];
dq0=zeros(nq,1); dq0(3)=0.3;
J=J_stance(x0(1:nq),packParameters(params));
dq0_c=(eye(7)-pinv(J)*J)*dq0;
x0(nq+1:end)=dq0_c;

%% Simulate to steady state
controller = @rabbit_controller;
[t_all, x_all, impact_log] = simulate_n_steps(x0, params, 5, controller);
if size(x_all,1)==14, x_all=x_all'; end

%% Extract step 3
impact_t = [];
for i = 1:length(impact_log)
    if ~isempty(impact_log(i).impact_time)
        impact_t = [impact_t; impact_log(i).impact_time];
    end
end

t_start = impact_t(2); t_end = impact_t(3);
step_idx = t_all >= t_start & t_all <= t_end;
t_step = t_all(step_idx); t_step = t_step - t_step(1);
x_step = x_all(step_idx, :);

fprintf('Step 3: %.3f s, %d samples\n', t_step(end), sum(step_idx));

%% Extract and fit using x_hip as phase
q_act = x_step(:, 4:7);
theta_traj = x_step(:, 1);  % X_HIP!

theta0_fit = theta_traj(1);
thetaf_fit = theta_traj(end);
s_traj = (theta_traj - theta0_fit) / (thetaf_fit - theta0_fit);

fprintf('Phase (x_hip): [%.4f, %.4f] m, Δ=%.4f m\n', ...
    theta0_fit, thetaf_fit, thetaf_fit-theta0_fit);
fprintf('s range: [%.4f, %.4f]\n', min(s_traj), max(s_traj));

%% Fit B-splines
n_cp = 8; p_deg = 3; n = n_cp - 1;
CP = zeros(n_cp, 4);
s_target = linspace(0, 1, n_cp)';

for j = 1:4
    [s_unique, idx_unique] = unique(s_traj);
    q_unique = q_act(idx_unique, j);
    if length(s_unique) >= 2
        CP(:, j) = interp1(s_unique, q_unique, s_target, 'pchip', 'extrap');
    else
        CP(:, j) = linspace(q_act(1,j), q_act(end,j), n_cp)';
    end
end

%% Verify
s_check = linspace(0, 1, 100)';
hd_check = zeros(100, 4);
for k = 1:100
    N = BSpline(n, p_deg, s_check(k));
    hd_check(k, :) = (N * CP)';
end
q_check = zeros(100, 4);
for j = 1:4
    [s_unique, idx_unique] = unique(s_traj);
    q_unique = q_act(idx_unique, j);
    if length(s_unique) >= 2
        q_check(:, j) = interp1(s_unique, q_unique, s_check, 'pchip', 'extrap');
    end
end
fit_err = rms(q_check - hd_check, 'all');
fprintf('Fit error: %.4f rad\n', fit_err);

%% Save
step_length = x_step(end,1) - x_step(1,1);
save(fullfile('Results', 'initial_guess.mat'), ...
    'CP', 'theta0_fit', 'thetaf_fit', 'fit_err', ...
    'step_length', 'n', 'p_deg', 't_step');
fprintf('Saved. Step: %.3f m, speed: %.3f m/s\n', ...
    step_length, step_length/t_step(end));

%% Plot
figure('Name', 'B-Spline Fit (Phase = x_{hip})', 'Position', [50, 50, 1200, 800]);
joint_names = {'q1', 'q2', 'q3', 'q4'};
for j = 1:4
    subplot(2,4,j); hold on;
    plot(s_traj, q_act(:,j), 'b.', 'MarkerSize', 4);
    plot(s_check, hd_check(:,j), 'r-', 'LineWidth', 2);
    plot(s_target, CP(:,j), 'ko', 'MarkerSize', 10, 'MarkerFaceColor', 'k');
    xlabel('s'); ylabel('rad'); title(joint_names{j}); grid on;
    
    subplot(2,4,4+j); hold on;
    plot(t_step, q_act(:,j), 'b-', 'LineWidth', 1.5);
    xlabel('Time (s)'); ylabel('rad'); title([joint_names{j} ' (time)']); grid on;
end
sgtitle(sprintf('B-Spline Fit, Phase = x_{hip} (err=%.4f rad)', fit_err));

%% Animate
fprintf('\nStarting animation...\n');
valid = x_all(:,2) > 0.5 & max(abs(x_all(:,4:7)), [], 2) < 5;
if sum(valid) > 10
    params.speed = 0.3;
    animate_rabbit_stepping_stones(x_all(valid,:), params);
end