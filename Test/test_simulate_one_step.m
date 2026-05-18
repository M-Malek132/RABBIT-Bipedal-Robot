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

qt = 0.2;     % torso

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

%% Run simulation

controller = @rabbit_controller;

[t,x,impact] = simulate_one_step(x0,params,controller);

%% Plot joint angles

figure
plot(t,x(:,1:nq),'LineWidth',1.5)
title("Joint Angles")
xlabel("time (s)")
ylabel("rad")
grid on

%% Plot joint velocities

figure
plot(t,x(:,nq+1:end),'LineWidth',1.5)
title("Joint Velocities")
xlabel("time (s)")
ylabel("rad/s")
grid on

%% Display impact information

disp("Impact information:")
disp(impact)

if isempty(impact.time)
    disp("No impact detected during simulation.")
else
    fprintf("Impact time: %f seconds\n",impact.time);
end

%% Compute controller effort (torques)

nu = 4;   % RABBIT has 4 actuators
u = zeros(length(t),nu);

for k = 1:length(t)
    u(k,:) = controller(t(k),x(k,:)',params)';
end

%% Show maximum torque for debugging

fprintf("Maximum controller torque: %f Nm\n",max(abs(u(:))));

%% Plot controller torques

figure
plot(t,u,'LineWidth',1.5)
title("Controller Effort (Torques)")
xlabel("time (s)")
ylabel("Torque (Nm)")
legend("Stance Hip","Stance Knee","Swing Hip","Swing Knee")
grid on

%% Animate robot motion
figure
animate_rabbit(x',params)

%% Final stance foot verification

[p_st,~,~,~,~,~] = rabbit_kinematics(q0,packParameters(params));

disp("Initial stance foot position:")
disp(p_st)

[stance,swing,~,~,~,~] = rabbit_kinematics(x(end,1:7),packParameters(params))
