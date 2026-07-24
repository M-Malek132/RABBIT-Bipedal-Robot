function main_demo()
%MAIN_DEMO  Legacy stepping-stone walking demo  (NON-RUNNABLE REFERENCE).
%
% ==> NON-RUNNABLE. This predates the HZD pipeline and calls interfaces that no
%     longer exist: make_initial_state with 3 outputs (it returns 1),
%     simulate_n_steps with 4 args (it takes 3), and
%     animate_rabbit_stepping_stones (does not exist). Kept for reference only.
%
%     Working entry point:
%       Trajectory_Optimization/rabbit_hzd_trajectory_optimization.m
%     See README "Broken legacy entry points".

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
