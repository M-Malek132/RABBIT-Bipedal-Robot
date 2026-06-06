%% ============================================================
%  main_HZD_optimization.m
%  HZD Gait Optimization for RABBIT (7-DOF floating-base model)
%
%  Coordinates:  q = [px, pz, qt, q1, q2, q3, q4]'
%    px, pz  : base (torso) position
%    qt      : torso absolute angle
%    q1, q2  : stance hip / knee angles
%    q3, q4  : swing  hip / knee angles
%
%  Actuation:    u = [u1, u2, u3, u4]'  (q1..q4 are actuated)
%
%  Dependencies (all already in your repo):
%    Model/parameters.m
%    Dynamics/D_matrix.m, C_vector.m, G_vector.m, input_matrix.m
%    Model/rabbit_kinematics.m
%    Contact/rabbit_impact_map.m, rabbit_impact_event.m
%    Reset_Map/rabbit_reset_map.m
%    Controller/rabbit_virtual_constraints.m  (if you want to reuse it)
%
%  NEW files this script depends on (provided alongside this file):
%    hzd_closedLoopDynamics.m
%    hzd_simulateOneStep.m
%    hzd_virtualConstraints.m
%    hzd_h0_outputs.m
%    hzd_phaseVariable.m
%    hzd_objectiveHZD.m
%    hzd_constraintsHZD.m
%    hzd_unpackDecisionVars.m
%    hzd_bezier.m
%    hzd_plotResults.m
%% ============================================================

clear; clc; close all;

% ---- make sure repo is on path --------------------------------
% (startup.m already does this; call it if needed)
if ~exist('parameters', 'file')
    run(fullfile(fileparts(which('main_HZD_optimization')), '..', 'startup.m'));
end

%% ---- Robot & model dimensions --------------------------------
params      = parameters();      % load RABBIT physical parameters
model.nq    = 7;                 % q = [px pz qt q1 q2 q3 q4]
model.nu    = 4;                 % actuated: q1 q2 q3 q4
model.g     = params.g;
model.params = params;           % pass through to all sub-functions

%% ---- Optimization settings -----------------------------------
opt.v_des          = 0.5;        % desired average walking speed [m/s]
opt.M              = 5;          % Bezier polynomial degree
opt.ny             = model.nu;   % number of virtual constraints (= nu)
opt.nCoeff         = opt.M + 1;
opt.nAlpha         = opt.ny * opt.nCoeff;

opt.mu             = 0.6;        % ground friction coefficient
opt.hipHeightMin   = 0.55;       % minimum hip height [m]

opt.Tmin           = 0.25;       % step duration bounds [s]
opt.Tmax           = 1.0;

opt.uMax           =  80;        % joint torque limits [Nm]
opt.uMin           = -80;

% Phase variable limits (theta at start / end of step)
% theta = c'*q;  see hzd_phaseVariable.m for definition.
% These should bracket one full step; tune if optimiser fails.
opt.thetaStart     = -0.20;
opt.thetaEnd       =  0.20;

% Joint angle limits (for all 7 DOF; px/pz limits are large)
opt.qMin = [-5; -5;  -pi; -pi/2; -pi; -pi/2; -pi];
opt.qMax = [ 5;  5;   pi;  pi/2;   0;  pi/2;    0];
% NOTE: knee angles (q2, q4) are negative on RABBIT when bent.
%       Adjust the above to match your sign convention.

% Controller gains (for the inner feedback linearisation loop)
opt.Kp = 200;
opt.Kd =  30;

%% ---- Initial guess -------------------------------------------
% Reasonable flat-footed standing pose for RABBIT
q0_init  = [0.00;   % px  – doesn't matter for periodicity
            0.90;   % pz  – torso ~0.9 m above ground
           -0.05;   % qt  – torso slightly backward
           -0.20;   % q1  – stance hip flexed
           -0.40;   % q2  – stance knee bent
            0.20;   % q3  – swing  hip extended
           -0.30];  % q4  – swing  knee slightly bent

dq0_init = [0.50;   % px_dot  – forward walking speed
            0.00;
            0.00;
            0.00;
            0.00;
            0.00;
            0.00];

% Initial Bezier coefficients: shape a smooth ramp for each output
alpha0 = zeros(opt.ny, opt.nCoeff);
alpha0(1,:) = linspace(-0.20,  0.20, opt.nCoeff);  % q1 trajectory
alpha0(2,:) = linspace(-0.40, -0.25, opt.nCoeff);  % q2 trajectory
alpha0(3,:) = linspace( 0.20, -0.20, opt.nCoeff);  % q3 trajectory
alpha0(4,:) = linspace(-0.30, -0.25, opt.nCoeff);  % q4 trajectory

T0 = 0.50;   % initial step duration guess

%% ---- Decision variable vector --------------------------------
% z = [alpha(:); q0(1..nq); dq0(1..nq); T]
alpha0_vec = alpha0(:);
z0 = [alpha0_vec; q0_init; dq0_init; T0];

%% ---- Bounds --------------------------------------------------
lb_alpha = -8 * ones(opt.nAlpha, 1);
ub_alpha =  8 * ones(opt.nAlpha, 1);

lb_q  = opt.qMin;
ub_q  = opt.qMax;

lb_dq = -15 * ones(model.nq, 1);
ub_dq =  15 * ones(model.nq, 1);

lb = [lb_alpha; lb_q; lb_dq; opt.Tmin];
ub = [ub_alpha; ub_q; ub_dq; opt.Tmax];

%% ---- fmincon options -----------------------------------------
options = optimoptions('fmincon', ...
    'Algorithm',              'sqp', ...
    'Display',                'iter', ...
    'MaxIterations',          500, ...
    'MaxFunctionEvaluations', 5e5, ...
    'ConstraintTolerance',    1e-4, ...
    'OptimalityTolerance',    1e-5, ...
    'StepTolerance',          1e-9, ...
    'FiniteDifferenceStepSize', 1e-6);

%% ---- Run optimisation ----------------------------------------
fprintf('\n=== Starting HZD gait optimisation ===\n');
fprintf('Decision vars: %d  (alpha:%d  q0:%d  dq0:%d  T:1)\n', ...
    length(z0), opt.nAlpha, model.nq, model.nq);

problem.objective = @(z) hzd_objectiveHZD(z, model, opt);
problem.nonlcon   = @(z) hzd_constraintsHZD(z, model, opt);
problem.x0        = z0;
problem.lb        = lb;
problem.ub        = ub;
problem.solver    = 'fmincon';
problem.options   = options;

[zStar, JStar, exitflag, output] = fmincon(problem);

fprintf('\n=== Optimisation complete ===\n');
fprintf('Cost      = %.6f\n', JStar);
fprintf('Exit flag = %d\n',   exitflag);
fprintf('Message   : %s\n',   output.message);

%% ---- Extract & display solution ------------------------------
[alphaStar, q0Star, dq0Star, TStar] = hzd_unpackDecisionVars(zStar, model, opt);

fprintf('\nOptimised step duration: %.4f s\n', TStar);
fprintf('Optimised Bezier coefficients:\n');
disp(alphaStar);

%% ---- Simulate optimised gait (5 steps) -----------------------
simOpt.alpha  = alphaStar;
simOpt.T      = TStar;
simOpt.Kp     = opt.Kp;
simOpt.Kd     = opt.Kd;

nSteps = 5;
x0Star = [q0Star; dq0Star];

[tAll, xAll, uAll] = hzd_simulateNSteps(x0Star, model, opt, simOpt, nSteps);

%% ---- Plot ----------------------------------------------------
hzd_plotResults(tAll, xAll, uAll, model, opt);

%% ---- Save results --------------------------------------------
timestamp = datestr(now, 'yyyy-mm-dd_HH-MM-SS');
save(fullfile('..','Results', ['hzd_result_', timestamp, '.mat']), ...
    'alphaStar', 'q0Star', 'dq0Star', 'TStar', 'tAll', 'xAll', 'uAll', ...
    'model', 'opt', 'JStar');
fprintf('Results saved.\n');