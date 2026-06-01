% make_initial_guess.m
% Generate a good initial guess from the working controller

clear; clc;
if exist('startup.m','file'), startup; end

params = parameters();

% Use the EXACT same setup as main_demo
nq = 7;
x0 = zeros(2*nq, 1);

qt = 0.1; q1 = -0.3; q2 = 0.6; q3 = -1.0; q4 = 0.6;

l1 = params.l1; l2 = params.l2;
px = l1*sin(qt+q1) + l2*sin(qt+q1+q2);
pz = l1*cos(qt+q1) + l2*cos(qt+q1+q2);

q0 = [px; pz; qt; q1; q2; q3; q4];
x0(1:nq) = q0;

dq0 = zeros(nq, 1);
dq0(3) = 0.3;

J = J_stance(q0, packParameters(params));
dq0_corrected = (eye(7) - pinv(J)*J) * dq0;
x0(nq+1:end) = dq0_corrected;

% Simulate ONE step with the working controller
controller = @rabbit_controller;
[t_step, x_step, impact_info] = simulate_one_step(x0, params, controller);

fprintf('Step time: %.3f s, samples: %d\n', t_step(end), length(t_step));

% Extract actuated joint trajectories
q_traj = x_step(:, 1:7);
q_act = q_traj(:, 4:7);  % [q1, q2, q3, q4]

% Phase variable (same as rabbit_controller)
theta_traj = q_traj(:, 4) + 0.5 * q_traj(:, 6);
theta0 = theta_traj(1);
thetaf = theta_traj(end);
s_traj = (theta_traj - theta0) / (thetaf - theta0);

fprintf('Phase range: [%.4f, %.4f]\n', theta0, thetaf);

% Fit B-splines
n_cp = 8; p_deg = 3;
CP = zeros(n_cp, 4);
s_target = linspace(0, 1, n_cp)';

for j = 1:4
    CP(:, j) = interp1(s_traj, q_act(:, j), s_target, 'linear', 'extrap');
end

% Verify fit
n = n_cp - 1;
s_check = linspace(0, 1, 100)';
hd_check = zeros(100, 4);
for k = 1:100
    N = BSpline(n, p_deg, s_check(k));
    hd_check(k, :) = (N * CP)';
end
q_check = interp1(s_traj, q_act, s_check);
fit_err = rms(q_check - hd_check);
fprintf('B-spline fit RMS error: %.4f rad\n', fit_err);

% Save initial guess
save(fullfile('Results', 'initial_guess.mat'), 'CP', 'theta0', 'thetaf', 'fit_err');
fprintf('Saved initial guess to Results/initial_guess.mat\n');

% Plot
figure;
for j = 1:4
    subplot(2,2,j); hold on;
    plot(s_traj, q_act(:,j), 'b.', 'MarkerSize', 3);
    plot(s_check, hd_check(:,j), 'r-', 'LineWidth', 2);
    plot(s_target, CP(:,j), 'ko', 'MarkerSize', 8);
    xlabel('Phase s'); ylabel('Angle (rad)');
    title(sprintf('q%d', j)); grid on;
    legend('Recorded', 'B-spline', 'CPs');
end
sgtitle('Initial Guess from Working Controller');