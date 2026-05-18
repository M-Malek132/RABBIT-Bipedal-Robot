clear; clc; close all;
% Get the directory of the current script
current_dir = fileparts(mfilename('fullpath'));

% Get the parent directory
parent_dir = fileparts(current_dir);

% Add parent directory and all its subfolders to MATLAB path
addpath(genpath(parent_dir));

disp('RABBIT 5-Link Passive Walker');

params.num_steps = input("How many steps: ");
params.anim = input("Show animation? (1/0): ");

params.alpha = deg2rad(3);   % slope

nq = 7;

% initial state
q0 = [0; 1; -0.2; 0.4; 0.1; -0.3; 0.6];
qd0 = zeros(nq,1);

x0 = [q0; qd0];

t_traj = [];
x_traj = [];

for i = 1:params.num_steps

    if i > 1
        x_impact = rabbit_impact_map(x_end,params);
        x0 = rabbit_reset_map(x_impact,params);
    end

    tspan = 0:0.02:10;

    options = odeset('Events',@(t,x) rabbit_impact_event(t,x,params));

    [t,x] = ode45(@(t,x) rabbit_ode(t,x,params),tspan,x0,options);

    t = t';
    x = x';

    t_traj = [t_traj t];
    x_traj = [x_traj x];

    x_end = x(:,end);

end

if params.anim
    animate_rabbit(x_traj,params)
end
