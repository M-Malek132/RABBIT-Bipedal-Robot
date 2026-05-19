function main_demo()
clc;
clear;
close all;

fprintf("=====================================\n");
fprintf(" RABBIT Biped Robot Demo\n");
fprintf(" Hybrid Walking Simulation\n");
fprintf("=====================================\n");

%% Initialize project
startup;

%% Load robot parameters
params = parameters();
p = packParameters(params);

%% Number of coordinates
nq = 7;

%% Initial state vector
x0 = zeros(2*nq,1);

%% Initial configuration (angles)

qt = 0.1;      % torso angle

q1 = -0.3;     % stance hip
q2 = 0.6;      % stance knee

q3 = -1.0;     % swing hip
q4 = 0.6;      % swing knee

%% Link lengths
l1 = params.l1;
l2 = params.l2;

%% Solve base position so stance foot touches ground

px = l1*sin(qt+q1) + l2*sin(qt+q1+q2);
pz = l1*cos(qt+q1) + l2*cos(qt+q1+q2);

q0 = [px; pz; qt; q1; q2; q3; q4];

x0(1:nq) = q0;

%% Initial velocities
dq0 = zeros(nq,1);

% small torso velocity to break symmetry
dq0(3) = 0.3;

%% Enforce stance foot contact constraint

J = J_stance(q0,p);

dq0_corrected = (eye(7) - pinv(J)*J) * dq0;

x0(nq+1:end) = dq0_corrected;

%% Check initial foot positions

[p_st,p_sw,~,~,~,~] = rabbit_kinematics(q0,p);

fprintf("Initial stance foot: [%f , %f]\n",p_st(1),p_st(2));
fprintf("Initial swing  foot: [%f , %f]\n",p_sw(1),p_sw(2));

%% Simulation settings

nSteps = 20;
controller = @rabbit_controller;

fprintf("Running walking simulation...\n");

[t_all, x_all, impact_log] = simulate_n_steps( ...
    x0, ...
    params, ...
    nSteps, ...
    controller);

fprintf("Simulation finished.\n");

%% Animate robot

animate_rabbit_stepping_stones(x_all,params);

fprintf("Animation complete.\n");

end
