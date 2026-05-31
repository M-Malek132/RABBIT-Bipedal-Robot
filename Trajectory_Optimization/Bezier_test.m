% test_virtual_constraints.m
% Test the updated virtual constraints with torso phase variable

clear; clc; close all;

%% Parameters
M = 5;              % 5th order Bézier (6 coefficients per output)
theta0 = -0.1;      % Torso angle at start of step
thetaf = 0.3;       % Torso angle at end of step

%% Create virtual constraints object
vc = bezier_virtual_constraints(M, theta0, thetaf);

%% Define boundary conditions for actuated joints
% Start of step (s=0)
qa_start = [-0.4;    % q1: stance knee (slightly bent)
             0.5;    % q2: stance hip 
            -0.4;    % q3: swing knee (slightly bent behind)
             0.5];   % q4: swing hip

% End of step (s=1, just before impact)
qa_end   = [ 0.3;    % q1: stance knee extended
            -0.2;    % q2: stance hip behind
             0.3;    % q3: swing knee extended forward
            -0.2];   % q4: swing hip forward

%% Initialize coefficients with linear interpolation
alpha = vc.initialize_coefficients(qa_start, qa_end);

%% Plot trajectories
vc.plot_control_points(alpha);

%% Test evaluation at different states
fprintf('Testing virtual constraint evaluation:\n\n');

% Test at beginning of step
q_start = [0; 0.85; theta0; qa_start];
dq_start = [0.5; 0; 0.2; 0.5; -1.0; 0.5; -1.0];  % example velocities

[s, ds] = vc.compute_phase(q_start, dq_start);
fprintf('At s=0:\n');
fprintf('  s = %.3f (should be ~0)\n', s);
fprintf('  ds/dt = %.3f\n\n', ds);

[y, dy] = vc.virtual_constraint(q_start, dq_start, alpha);
fprintf('  y = [%.3f, %.3f, %.3f, %.3f]\n', y);
fprintf('  (should be near zero at s=0)\n\n');

% Test at end of step
q_end = [0.5; 0.85; thetaf; qa_end];
dq_end = [1.0; 0; 0.3; -1.0; 1.0; -1.0; 1.0];

[s, ds] = vc.compute_phase(q_end, dq_end);
fprintf('At s=1:\n');
fprintf('  s = %.3f (should be ~1)\n', s);
fprintf('  ds/dt = %.3f\n\n', ds);

[y, dy] = vc.virtual_constraint(q_end, dq_end, alpha);
fprintf('  y = [%.3f, %.3f, %.3f, %.3f]\n', y);
fprintf('  (should be near zero at s=1)\n');

%% Print coefficient structure
fprintf('\nCoefficient matrix structure:\n');
fprintf('  Size: %d x %d\n', size(alpha));
fprintf('  Each column j: coefficients for output j\n');
fprintf('  Row k (k=1...%d): coefficient for B_{k-1}(s)\n', M+1);
fprintf('  Total optimization variables: %d\n', vc.num_coefficients());