% test_initial_guess_v2.m
clear; clc; close all;
if exist('startup.m','file'), startup; end

data = load(fullfile('Results', 'initial_guess.mat'));
fprintf('Loaded. Phase: [%.1f, %.1f]\n', data.theta0, data.thetaf);

%% Create controller (SAME phase as rabbit_controller)
traj = BSplineTrajectory(data.n, data.p_deg, data.theta0, data.thetaf);
traj.CP = data.CP;

Kp = diag([200, 200, 150, 150]);
Kd = diag([30,  30,  20,  20 ]);
ctrl = RabbitController(traj, Kp, Kd);

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

%% Simulate
controller_handle = ctrl.to_function_handle();
fprintf('Simulating...\n');
[t_all, x_all, ~] = simulate_n_steps(x0, params, 5, controller_handle);

if size(x_all,1)==14, x_traj = x_all'; else, x_traj = x_all; end
valid = x_traj(:,2) > 0.5;
fprintf('Valid: %d/%d\n', sum(valid), length(t_all));

if sum(valid) > 10
    t_v = t_all(valid); x_v = x_traj(valid,:);
    figure('Name', 'B-Spline = desired\_gait');
    subplot(2,2,1); plot(t_v, x_v(:,4:7)); title('Joints'); grid on;
    subplot(2,2,2); plot(x_v(:,1), x_v(:,2)); title('Hip'); grid on;
    subplot(2,2,3); plot(t_v, x_v(:,8)); title('Speed'); grid on;
    subplot(2,2,4); traj.plot(); title('B-Splines');
    
    if sum(valid) > 20
        animate_rabbit_stepping_stones(x_v, params);
    end
end