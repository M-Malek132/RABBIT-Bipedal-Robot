clear; clc; close all;

% Get the directory of the current script
current_dir = fileparts(mfilename('fullpath'));

% Get the parent directory
parent_dir = fileparts(current_dir);

% Add parent directory and all its subfolders to MATLAB path
addpath(genpath(parent_dir));

%% ============================================================
%  HZD Gait Optimization Template
%  ------------------------------------------------------------
%  This script optimizes virtual constraints for a walking robot
%  using the Hybrid Zero Dynamics framework.
%% ============================================================

%% Robot and gait parameters

model.nq = 5;          % number of generalized coordinates
model.nu = 4;          % number of actuators
model.g  = 9.81;

% Desired average walking speed
opt.v_des = 0.8;       % m/s

% Bezier polynomial degree
opt.M = 6;

% Number of outputs
opt.ny = model.nu;

% Number of Bezier coefficients per output
opt.nCoeff = opt.M + 1;

% Total number of alpha parameters
opt.nAlpha = opt.ny * opt.nCoeff;

% Friction coefficient
opt.mu = 0.6;

% Minimum hip height
opt.hipHeightMin = 0.55;

% Step time bounds
opt.Tmin = 0.25;
opt.Tmax = 1.0;

% Torque bounds
opt.uMax = 80;         % Nm
opt.uMin = -80;        % Nm

% Joint limits, example values
opt.qMin = [-pi; -pi; -pi; -pi; -pi];
opt.qMax = [ pi;  pi;  pi;  pi;  pi];

%% Initial guess

% Initial Bezier coefficients
alpha0 = zeros(opt.ny, opt.nCoeff);

% Example initial guess
% You should replace this with a reasonable pose trajectory
alpha0(1,:) = linspace( 0.2, -0.2, opt.nCoeff);
alpha0(2,:) = linspace(-0.3,  0.3, opt.nCoeff);
alpha0(3,:) = linspace(-0.5, -0.2, opt.nCoeff);
alpha0(4,:) = linspace(-0.2, -0.6, opt.nCoeff);

alpha0_vec = alpha0(:);

% Initial condition after impact
q0 = [0.1; -0.2; -0.3; -0.4; 0.2];
dq0 = [0.5; -0.2; 0.1; -0.1; 0.3];

% Initial step time
T0 = 0.6;

% Decision variable:
% z = [alpha(:); q0; dq0; T]
z0 = [alpha0_vec; q0; dq0; T0];

%% Bounds

lb_alpha = -5 * ones(opt.nAlpha,1);
ub_alpha =  5 * ones(opt.nAlpha,1);

lb_q  = opt.qMin;
ub_q  = opt.qMax;

lb_dq = -10 * ones(model.nq,1);
ub_dq =  10 * ones(model.nq,1);

lb_T = opt.Tmin;
ub_T = opt.Tmax;

lb = [lb_alpha; lb_q; lb_dq; lb_T];
ub = [ub_alpha; ub_q; ub_dq; ub_T];

%% Optimization options

options = optimoptions('fmincon', ...
    'Algorithm', 'sqp', ...
    'Display', 'iter', ...
    'MaxIterations', 300, ...
    'MaxFunctionEvaluations', 2e5, ...
    'ConstraintTolerance', 1e-5, ...
    'OptimalityTolerance', 1e-5, ...
    'StepTolerance', 1e-8);

%% Run optimization

problem.objective = @(z) objectiveHZD(z, model, opt);
problem.nonlcon   = @(z) constraintsHZD(z, model, opt);
problem.x0        = z0;
problem.lb        = lb;
problem.ub        = ub;
problem.solver    = 'fmincon';
problem.options   = options;

[zStar, JStar, exitflag, output] = fmincon(problem);

disp('Optimization complete.');
disp(['Cost = ', num2str(JStar)]);
disp(['Exit flag = ', num2str(exitflag)]);

%% Extract solution

[alphaStar, q0Star, dq0Star, TStar] = unpackDecisionVariables(zStar, model, opt);

disp('Optimized step time:');
disp(TStar);

disp('Optimized Bezier coefficients:');
disp(alphaStar);

%% Simulate optimized gait

x0Star = [q0Star; dq0Star];

simOpt.alpha = alphaStar;
simOpt.T = TStar;
simOpt.Kp = 100;
simOpt.Kd = 20;

nSteps = 5;

[tAll, xAll, uAll] = simulateHybridWalking(x0Star, model, opt, simOpt, nSteps);

%% Plot results

plotHZDResults(tAll, xAll, uAll, model, opt);
