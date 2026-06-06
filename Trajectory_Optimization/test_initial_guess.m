% test_initial_guess.m - FIXED
% Set theta0 to current x_hip at each step start

clear; clc; close all;
if exist('startup.m','file'), startup; end

%% Load B-spline fit
data = load(fullfile('Results', 'initial_guess.mat'));
fprintf('Loaded initial guess:\n');
fprintf('  CP: %dx%d\n', size(data.CP,1), size(data.CP,2));
fprintf('  Phase (x_hip): [%.4f, %.4f] (Δ=%.3f m)\n', ...
    data.theta0_fit, data.thetaf_fit, data.thetaf_fit - data.theta0_fit);

% Use RELATIVE phase: always 0 to step_length
step_length = data.step_length;
CP = data.CP;
n = data.n;
p_deg = data.p_deg;

%% Setup initial state (same as main_demo)
params = parameters();
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

%% Create controller with phase relative to current hip position
% theta0 = current x_hip, thetaf = current x_hip + step_length
x_hip_start = x0(1);
traj = BSplineTrajectory(n, p_deg, x_hip_start, x_hip_start + step_length);
traj.CP = CP;

Kp = diag([200, 200, 150, 150]);
Kd = diag([30,  30,  20,  20 ]);
ctrl = RabbitController(traj, Kp, Kd);

fprintf('Initial x_hip = %.4f, phase range = [%.4f, %.4f]\n', ...
    x_hip_start, x_hip_start, x_hip_start + step_length);

%% Simulate
controller_handle = ctrl.to_function_handle();
fprintf('Simulating 5 steps...\n');
[t_all, x_all, impact_log] = simulate_n_steps(x0, params, 5, controller_handle);

if size(x_all,1)==14, x_traj = x_all'; else, x_traj = x_all; end
valid = x_traj(:,2) > 0.3 & max(abs(x_traj(:,4:7)),[],2) < 5;
fprintf('Valid states: %d/%d\n', sum(valid), length(t_all));

%% Plot
if sum(valid) > 10
    t_valid = t_all(valid);
    x_valid = x_traj(valid, :);
    
    figure('Name', 'B-Spline Controller', 'Position', [50, 50, 1200, 800]);
    
    subplot(2,3,1);
    plot(t_valid, x_valid(:,4:7), 'LineWidth', 1.5);
    xlabel('Time'); ylabel('Angle'); legend('q1','q2','q3','q4');
    title('Joint Angles'); grid on;
    
    subplot(2,3,2);
    s_traj = zeros(size(t_valid));
    for i = 1:length(t_valid)
        s_traj(i) = ctrl.trajectory.phase(x_valid(i,1:7)');
    end
    plot(t_valid, s_traj, 'LineWidth', 1.5);
    xlabel('Time'); ylabel('s');
    title('Normalized Phase'); grid on;
    
    subplot(2,3,3);
    plot(x_valid(:,1), x_valid(:,2), 'LineWidth', 1.5);
    xlabel('x'); ylabel('z'); title('Hip'); grid on; axis equal;
    
    subplot(2,3,4);
    plot(t_valid, x_valid(:,8), 'LineWidth', 1.5);
    xlabel('Time'); ylabel('v_x'); title('Speed'); grid on;
    
    subplot(2,3,5);
    y_err = zeros(length(t_valid), 4);
    for i = 1:length(t_valid)
        [y, ~] = ctrl.trajectory.virtual_constraint(x_valid(i,1:7)', x_valid(i,8:14)');
        y_err(i,:) = y';
    end
    plot(t_valid, y_err, 'LineWidth', 1.5);
    xlabel('Time'); ylabel('Error'); title('VC Errors'); grid on;
    
    subplot(2,3,6);
    traj.plot();
    title('B-Spline Trajectories');
    
    sgtitle('B-Spline Controller (Phase = relative x_{hip})');
    
    % Animate
    if sum(valid) > 20
        animate_rabbit_stepping_stones(x_valid, params);
    end
end