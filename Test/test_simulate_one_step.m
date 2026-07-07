%% SETUP & INITIALIZATION
clc; clear; close all;
startup; 
config('init');

% No p
% arams needed here
x0 = make_initial_state(); 
% animate_rabbit(x0);

% Controller is now decoupled from the workspace state
[t, x, impact] = simulate_n_steps(x0, 5, @rabbit_controller);

%% VISUALIZATION
% No params passed to animation
animate_rabbit(x)