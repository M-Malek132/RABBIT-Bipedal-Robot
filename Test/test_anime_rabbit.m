%% Main Test Script for RABBIT Animation
clear; clc; close all;
% Get the directory of the current script
current_dir = fileparts(mfilename('fullpath'));

% Get the parent directory
parent_dir = fileparts(current_dir);

% Add parent directory and all its subfolders to MATLAB path
addpath(genpath(parent_dir));

% Number of frames
N = 100;
t = linspace(0, 1, N);

% 1. Hip Position (Moving forward and slightly bobbing)
px = linspace(0, 0.8, N);
pz = 0.75 + 0.05 * sin(pi * t); % slight vertical oscillation

% 2. Stance Leg (q1: thigh, q2: knee)
% Stance leg starts forward and moves backward relative to hip
q1 = linspace(0.4, -0.4, N); 
q2 = 0.2 * ones(1, N); % slightly bent knee

% 3. Torso (q3)
q3 = 0.05 * sin(2 * pi * t); % slight swaying

% 4. Swing Leg (q4: thigh, q5: knee)
% Swing leg starts back and swings forward
q4 = linspace(-0.5, 0.5, N);
% Knee bends in the middle of the swing to clear the ground
q5 = 0.8 * sin(pi * t); 

% Combine into state trajectory matrix [7 x N]
% q = [px; pz; q1; q2; q3; q4; q5]
q_traj = [px; pz; q1; q2; q3; q4; q5];

% Add dummy velocities (zeros) to match the 14-state format [q; qdot]
x_traj = [q_traj; zeros(7, N)];

% Run Animation
disp('Starting Animation Test...');
animate_rabbit_stepping_stones(x_traj',parameters());
disp('Animation Finished.');
