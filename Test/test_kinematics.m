clc;
clear;
close all;

fprintf("Testing simulate_one_step...\n");

%% Initialize project
startup;

%% Load robot parameters
p = parameters();

%% Initial state
q0 = [0; 1; 0.1; -0.3; 0.6; 0; 0];

% ----- check Jacobian finite difference -----

eps = 1e-6;
dq_test = randn(7,1);

q_test = q0 + eps*dq_test;

[p1,~,~,~,~,~] = rabbit_kinematics(q0,packParameters(p));
[p2,~,~,~,~,~] = rabbit_kinematics(q_test,packParameters(p));

fd = (p2 - p1)/eps;      % finite difference velocity

J = J_stance(q0,packParameters(p));
Jdq = J*dq_test;

disp('finite difference vs Jacobian')
disp([fd Jdq fd-Jdq])


% ----- check constraint acceleration equation -----

% sample state
q = q0;
dq = randn(7,1)*0.1;

D = D_matrix(q,packParameters(p));
C = C_vector(q,dq,packParameters(p));
G = G_vector(q,packParameters(p));

J = J_stance(q,packParameters(p));
Jdotdq = Jdotdq_stance(q,dq,packParameters(p));

B = input_matrix();
u = zeros(4,1);

A = [D -J';
     J zeros(2,2)];

rhs = [B*u - C - G;
       -Jdotdq];

sol = A\rhs;

qdd = sol(1:7);

disp('constraint acceleration (should be ~0)')
disp(J*qdd + Jdotdq)


% ----- check stance foot acceleration directly -----

dt = 1e-6;

q_next  = q + dt*dq + 0.5*dt^2*qdd;
dq_next = dq + dt*qdd;

[p1,~,~,~,~,~] = rabbit_kinematics(q,packParameters(p));
[p2,~,~,~,~,~] = rabbit_kinematics(q_next,packParameters(p));

acc_fd = (p2 - p1 - dt*(J*dq)) / (0.5*dt^2);

disp('foot acceleration finite difference')
disp(acc_fd)
