% test_initial_guess.m
% Verify the B-spline initial guess actually walks before optimizing

clear; clc; close all;
if exist('startup.m','file'), startup; end

%% Load the initial guess
data = load(fullfile('Results', 'initial_guess.mat'));
fprintf('Loaded initial guess:\n');
fprintf('  CP size: %dx%d\n', size(data.CP,1), size(data.CP,2));
fprintf('  Phase range: [%.4f, %.4f]\n', data.theta0, data.thetaf);
fprintf('  Fit error: %.4f\n', data.fit_err);

%% Create controller with these B-splines
traj = BSplineTrajectory(7, 3, data.theta0, data.thetaf);
traj.CP = data.CP;

Kp = diag([200, 200, 150, 150]);
Kd = diag([30,  30,  20,  20 ]);
ctrl = RabbitController(traj, Kp, Kd);

%% Test 1: Use the EXACT same initial state as main_demo
fprintf('\n--- Test 1: Same state as main_demo ---\n');

params = parameters();
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

% Check phase
theta0_actual = x0(4) + 0.5*x0(6);
fprintf('  Initial state theta = %.4f\n', theta0_actual);
fprintf('  B-spline theta0    = %.4f\n', data.theta0);
fprintf('  Match: %s\n', check(abs(theta0_actual - data.theta0) < 0.01));

% Try one step
controller_handle = ctrl.to_function_handle();
[t_step, x_step, impact_info] = simulate_one_step(x0, params, controller_handle);

fprintf('  Step time: %.3f s\n', t_step(end));
fprintf('  Samples: %d\n', length(t_step));
fprintf('  Impact detected: %d\n', impact_info.detected);

% Check final phase
theta_final = x_step(end,4) + 0.5*x_step(end,6);
fprintf('  Final theta = %.4f (B-spline thetaf = %.4f)\n', theta_final, data.thetaf);

%% Test 2: Try multi-step
fprintf('\n--- Test 2: Multi-step ---\n');
[t_all, x_all, ~] = simulate_n_steps(x0, params, 5, controller_handle);

if size(x_all,1) == 14, x_traj = x_all'; else, x_traj = x_all; end
valid = x_traj(:,2) > 0.3;
fprintf('  Valid states: %d / %d\n', sum(valid), length(t_all));

%% Plot if we have data
if sum(valid) > 10
    figure('Name', 'Initial Guess Test');
    subplot(2,2,1);
    plot(t_all(valid), x_traj(valid,4:7)); xlabel('Time'); ylabel('Angle');
    title('Joint Angles'); grid on; legend('q1','q2','q3','q4');
    
    subplot(2,2,2);
    plot(x_traj(valid,1), x_traj(valid,2)); xlabel('x'); ylabel('z');
    title('Hip Trajectory'); grid on; axis equal;
    
    subplot(2,2,3);
    plot(t_all(valid), x_traj(valid,8)); xlabel('Time'); ylabel('v_x');
    title('Forward Speed'); grid on;
    
    subplot(2,2,4);
    traj.plot();
    title('B-Spline Virtual Constraints');
    
    % Animate
    if sum(valid) > 20
        params.speed = 0.3;
        animate_rabbit_stepping_stones(x_traj(valid,:), params);
    end
end

function result = check(condition)
    if condition, result = '✅'; else, result = '❌ MISMATCH'; end
end