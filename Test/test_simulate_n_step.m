clc;
clear;
close all;

fprintf("Testing simulate_one_step...\n");

%% Initialize project
startup;

%% Load robot parameters
params = parameters();

%% Number of coordinates
nq = 7;

%% Initial state
x0 = zeros(2*nq,1);

%% Choose initial configuration

qt = 0.1;     % torso

q1 = -0.3;    % stance hip
q2 = 0.6;     % stance knee

q3 = -1.0;    % swing hip
q4 = 0.6;     % swing knee

%% Link lengths
l1 = params.l1;
l2 = params.l2;

%% Solve base position so stance foot is on ground

x = l1*sin(qt+q1) + l2*sin(qt+q1+q2);
z = l1*cos(qt+q1) + l2*cos(qt+q1+q2);

q0 = [x; z; qt; q1; q2; q3; q4];

x0(1:nq) = q0;
x0(8) = x0(8) * 0; % Increase forward velocity by 20%

%% Initial velocities
x0(nq+1:end) = zeros(nq,1);

% small torso velocity to break symmetry
x0(nq+3) = 0.3;

J = J_stance(x0(1:7),packParameters(params));
initial_foot_vel = J * x0(8:14);
dq0_corrected = (eye(7) - pinv(J) * J) * x0(8:14);
initial_foot_vel_check = J * dq0_corrected;
disp('Corrected foot velocity:');
disp(initial_foot_vel_check);

x0(nq+1:end) = dq0_corrected;

%% Check initial foot positions

[p_st,p_sw,~,~,~,~] = rabbit_kinematics(q0,packParameters(params));

fprintf("Stance foot position: [%f  %f]\n",p_st(1),p_st(2));
fprintf("Swing foot position:  [%f  %f]\n",p_sw(1),p_sw(2));

% [t,x,impact] = simulate_one_step(x0, params, controller);

nSteps = 100;
controller = @rabbit_controller;

[t_all, x_all, impact_log] = simulate_n_steps( ...
    x0, ...
    params, ...
    nSteps, ...
    controller);

animate_rabbit(x_all,params)
