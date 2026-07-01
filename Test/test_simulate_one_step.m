%% SETUP & INITIALIZATION
clc; clear; close all;
startup; 
config('init');

% No p
% arams needed here
x0 = make_initial_state(); 

% Controller is now decoupled from the workspace state
[t, x, impact] = simulate_one_step(x0, @rabbit_controller);

%% VISUALIZATION
% No params passed to animation
animate_rabbit_stepping_stones(x'); 
