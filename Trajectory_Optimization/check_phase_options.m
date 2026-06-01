% check_phase_options.m
% Check which phase variable is monotonic

clear; clc;
if exist('startup.m','file'), startup; end

params = parameters();

%% Setup and simulate
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

[t_all, x_all, impact_log] = simulate_n_steps(x0, params, 5, @rabbit_controller);
if size(x_all,1)==14, x_all=x_all'; end

%% Find step 3 using impact times
impact_t = [];
for i = 1:length(impact_log)
    if ~isempty(impact_log(i).impact_time)
        impact_t = [impact_t; impact_log(i).impact_time];
    end
end

if length(impact_t) >= 3
    t_start = impact_t(2);
    t_end = impact_t(3);
else
    % Fallback: use approximate times from previous runs
    t_start = 1.07;
    t_end = 1.82;
end

step_idx = t_all >= t_start & t_all <= t_end;
t_s = t_all(step_idx); t_s = t_s - t_s(1);
x_s = x_all(step_idx, :);

%% Compute phase variable candidates
q = x_s(:, 1:7);

theta_old  = q(:,4) + 0.5*q(:,6);     % original: q1 + 0.5*q3
theta_x    = q(:,1);                    % hip x position
theta_qt   = q(:,3);                    % torso angle
theta_q1   = q(:,4);                    % stance hip only
theta_sub  = q(:,4) - 0.5*q(:,6);      % q1 - 0.5*q3

%% Check monotonicity
fprintf('\n=== Phase Variable Candidates (Step 3) ===\n');
fprintf('%-25s  %-8s  %-8s  %-8s  %s\n', 'Variable', 'Min', 'Max', 'Δ', 'Monotonic?');
fprintf('%-25s  %-8s  %-8s  %-8s  %s\n', '---------', '---', '---', '---', '----------');

vars = {'theta = q1+0.5*q3', 'theta = x_hip', 'theta = qt', ...
        'theta = q1', 'theta = q1-0.5*q3'};
datas = {theta_old, theta_x, theta_qt, theta_q1, theta_sub};

for i = 1:length(vars)
    d = datas{i};
    is_mono = all(diff(d) >= -1e-6) || all(diff(d) <= 1e-6);
    fprintf('%-25s  %-8.4f  %-8.4f  %-8.4f  %s\n', ...
        vars{i}, min(d), max(d), max(d)-min(d), check(is_mono));
end

%% Plot all candidates
figure('Name', 'Phase Variable Candidates', 'Position', [50, 50, 1200, 800]);

for i = 1:5
    subplot(2, 3, i);
    plot(t_s, datas{i}, 'b-', 'LineWidth', 2);
    xlabel('Time (s)'); ylabel('\theta');
    title(sprintf('%s  [Δ=%.4f]', vars{i}, max(datas{i})-min(datas{i})));
    grid on;
end

subplot(2, 3, 6);
plot(t_s, q(:,4:7), 'LineWidth', 1.5);
xlabel('Time (s)'); ylabel('Angle (rad)');
legend('q1','q2','q3','q4', 'Location', 'best');
title('Joint Angles (context)'); grid on;

sgtitle('Phase Variable Candidates for B-Spline Virtual Constraints');

%% Also show full 5 steps for the best candidate
figure('Name', 'Best Candidate - All Steps', 'Position', [100, 100, 1000, 400]);

theta_best = x_all(:, 1);  % x_hip (probably the best)

subplot(1,2,1);
plot(t_all, theta_best, 'b-', 'LineWidth', 1.5);
xlabel('Time (s)'); ylabel('x_{hip} (m)');
title('Hip X Position (All 5 Steps)'); grid on;

subplot(1,2,2);
plot(t_all, x_all(:,4:7), 'LineWidth', 1.5);
xlabel('Time (s)'); ylabel('Angle (rad)');
legend('q1','q2','q3','q4', 'Location', 'best');
title('Joint Angles'); grid on;

function result = check(condition)
    if condition, result = '✅'; else, result = '❌'; end
end