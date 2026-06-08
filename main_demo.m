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

[x0, params, p] = make_initial_state();

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

animate_rabbit_stepping_stones(x_all',params);

fprintf("Animation complete.\n");

end
