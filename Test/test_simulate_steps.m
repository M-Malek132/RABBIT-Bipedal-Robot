%% SETUP & INITIALIZATION
clc; clear; close all;
startup; 
config('init');

% No params needed here
x0 = make_initial_state();
animate_rabbit(x0)
% Controller is now decoupled from the workspace state (pass [] = passive/no
% torque, so this is a smoke test of the hybrid integration + impact/reset).
[t, x] = simulate_n_steps(x0, 10, []);

%% VISUALIZATION
% No params passed to animation
animate_rabbit(x)