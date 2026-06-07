%% ============================================================
%  main_HZD_optimization.m
%  HZD Gait Optimization for RABBIT  –  B-spline virtual constraints
%
%  Coordinates:  q = [px, pz, qt, q1, q2, q3, q4]'   (nq = 7)
%  Actuated:     u = [u1, u2, u3, u4]'               (nu = 4)
%
%  Decision variables:
%    z = [CP_vec; q0; dq0; T]
%    CP_vec  : (n+1)*4  – B-spline control-point matrix CP flattened
%              CP is (n+1)×4,  reshape(CP_vec, n+1, 4) recovers it
%    q0      : 7×1  initial state (configuration)
%    dq0     : 7×1  initial state (velocity)
%    T       : 1    step duration
%% ============================================================

clear; clc; close all;

if ~exist('parameters', 'file')
    run(fullfile(fileparts(which('main_HZD_optimization')), '..', 'startup.m'));
end

%% ---- Robot model ---------------------------------------------
params      = parameters();
model.nq    = 7;           % [px pz qt q1 q2 q3 q4]
model.nu    = 4;           % actuated: q1 q2 q3 q4
model.g     = params.g;
model.params = params;

%% ---- B-spline settings --------------------------------------
opt.n_bs   = 5;            % number of data points  (control points = n+1 = 6)
opt.p_bs   = 3;            % B-spline degree (cubic)
opt.ny     = model.nu;     % 4 outputs

opt.nCP_vars = (opt.n_bs + 1) * opt.ny;   % 6*4 = 24  (same count as Bezier M=5)

% Phase variable:  theta = q(1) = px  (increases monotonically as robot walks)
% Set these to the expected px range over one step.
% For v_des = 0.5 m/s, T ~ 0.5 s  =>  step_length ~ 0.25 m
% Start at 0, end ~ step_length.
opt.thetaStart =  0.00;    % px at step start [m]  — adjust after first run
opt.thetaEnd   =  0.30;    % px at step end   [m]  — adjust after first run

%% ---- Walking targets ----------------------------------------
opt.v_des        = 0.5;    % desired average walking speed [m/s]
opt.Tmin         = 0.25;
opt.Tmax         = 1.0;

%% ---- Physical constraints -----------------------------------
opt.mu           = 0.6;
opt.hipHeightMin = 0.55;
opt.uMax         =  150;
opt.uMin         = -150;

opt.qMin = [-5; -5;  -pi;  -pi/2; -pi;  -pi/2; -pi];
opt.qMax = [ 5;  5;   pi;   pi/2;   0;   pi/2;    0];

%% ---- Controller gains ---------------------------------------
opt.Kp = 200;
opt.Kd =  30;

%% ---- Initial guess -------------------------------------------
% Control points: (n+1) rows × 4 columns  [q1 q2 q3 q4]
% Initialise with a linear ramp from a plausible start to end pose.
q_act_start = [-0.20; -0.40;  0.20; -0.30];   % [q1 q2 q3 q4] at step start
q_act_end   = [ 0.20; -0.25; -0.20; -0.25];   % [q1 q2 q3 q4] at step end

n_pts = opt.n_bs + 1;
CP0 = zeros(n_pts, opt.ny);
for j = 1:opt.ny
    CP0(:, j) = linspace(q_act_start(j), q_act_end(j), n_pts);
end
CP0_vec = CP0(:);   % column-major, length = (n+1)*4

% Initial robot state
q0_init  = [0.00; 0.90; -0.05; -0.20; -0.40;  0.20; -0.30];
dq0_init = [0.50; 0.00;  0.00;  0.00;  0.00;  0.00;  0.00];
T0       = 0.50;

z0 = [CP0_vec; q0_init; dq0_init; T0];

%% ---- Variable bounds ----------------------------------------
lb_CP = -pi * ones(opt.nCP_vars, 1);
ub_CP =  pi * ones(opt.nCP_vars, 1);

lb = [lb_CP;  opt.qMin; -15*ones(model.nq,1); opt.Tmin];
ub = [ub_CP;  opt.qMax;  15*ones(model.nq,1); opt.Tmax];

%% ---- fmincon ------------------------------------------------
options = optimoptions('fmincon', ...
    'Algorithm',              'sqp', ...
    'Display',                'iter', ...
    'MaxIterations',          500, ...
    'MaxFunctionEvaluations', 5e5, ...
    'ConstraintTolerance',    1e-4, ...
    'OptimalityTolerance',    1e-5, ...
    'StepTolerance',          1e-9, ...
    'FiniteDifferenceStepSize', 1e-6);

fprintf('\n=== HZD B-spline gait optimisation ===\n');
fprintf('Decision vars: %d  (CP:%d  q0:%d  dq0:%d  T:1)\n', ...
    length(z0), opt.nCP_vars, model.nq, model.nq);

problem.objective = @(z) hzd_objectiveHZD(z, model, opt);
problem.nonlcon   = @(z) hzd_constraintsHZD(z, model, opt);
problem.x0        = z0;
problem.lb        = lb;
problem.ub        = ub;
problem.solver    = 'fmincon';
problem.options   = options;

[zStar, JStar, exitflag, output] = fmincon(problem);

fprintf('\n=== Done ===\n');
fprintf('Cost = %.6f,  exit = %d\n', JStar, exitflag);
fprintf('%s\n', output.message);

%% ---- Extract solution ---------------------------------------
[CPstar, q0Star, dq0Star, TStar] = hzd_unpackDecisionVars(zStar, model, opt);

fprintf('\nStep duration: %.4f s\n', TStar);
fprintf('Control points (rows=knot, cols=joints):\n');
disp(CPstar);

%% ---- Multi-step simulation ----------------------------------
simOpt.CP = CPstar;
simOpt.T  = TStar;
simOpt.Kp = opt.Kp;
simOpt.Kd = opt.Kd;

[tAll, xAll, uAll] = hzd_simulateNSteps([q0Star; dq0Star], model, opt, simOpt, 5);

%% ---- Plot ---------------------------------------------------
hzd_plotResults(tAll, xAll, uAll, model, opt);

%% ---- Save ---------------------------------------------------
ts = datestr(now, 'yyyy-mm-dd_HH-MM-SS');
save(fullfile('..','Results',['hzd_bspline_',ts,'.mat']), ...
    'CPstar','q0Star','dq0Star','TStar','tAll','xAll','uAll','model','opt','JStar');
fprintf('Saved.\n');
