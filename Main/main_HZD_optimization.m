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

% Phase variable:
%   theta = -q(3) - q(4) - 0.5*q(5)
%
% These values are set later after q0_init and q_act_end are defined.
% They are not step lengths.
opt.thetaStart = [];
opt.thetaEnd   = [];

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
q_act_start = [-0.20; -0.40;  0.20; -0.30];
q_act_end   = [ 0.20; -0.25; -0.20; -0.25];

n_pts = opt.n_bs + 1;
CP0 = zeros(n_pts, opt.ny);

for j = 1:opt.ny
    CP0(:, j) = linspace(q_act_start(j), q_act_end(j), n_pts);
end

CP0_vec = CP0(:);

% Initial robot configuration.
q0_init = [0.00; 0.90; -0.05; q_act_start];

% Set thetaStart from the actual initial configuration.
theta0_init = -q0_init(3) - q0_init(4) - 0.5*q0_init(5);
opt.thetaStart = theta0_init;

% Estimate thetaEnd from approximate final shape.
qt_end_guess = -0.05;
qEnd_shape_guess = [0.30; 0.90; qt_end_guess; q_act_end];

thetaEnd_guess = -qEnd_shape_guess(3) ...
                 -qEnd_shape_guess(4) ...
                 -0.5*qEnd_shape_guess(5);

opt.thetaEnd = thetaEnd_guess;

% Step duration guess.
T0 = 0.50;

% Initial velocity guess consistent with dy = 0.
[~, dhd_dtheta0] = hzd_evalBSpline(opt.thetaStart, CP0, opt);

% Choose phase speed.
% If thetaEnd < thetaStart, use negative dtheta.
if opt.thetaEnd < opt.thetaStart
    dtheta0 = -1.0;
else
    dtheta0 =  1.0;
end

dq_act0 = dhd_dtheta0 * dtheta0;

dpx0 = opt.v_des;
dpz0 = 0.0;

% Because theta_dot = -dq(3) - dq(4) - 0.5*dq(5)
dqt0 = -(1 + dhd_dtheta0(1) + 0.5*dhd_dtheta0(2)) * dtheta0;

dq0_init = [dpx0; dpz0; dqt0; dq_act0];

z0 = [CP0_vec; q0_init; dq0_init; T0];

%% ---- Variable bounds ----------------------------------------

% Control-point bounds based on actuated joint limits.
lb_CP_mat = repmat(opt.qMin(4:7).', n_pts, 1);
ub_CP_mat = repmat(opt.qMax(4:7).', n_pts, 1);

lb_CP = lb_CP_mat(:);
ub_CP = ub_CP_mat(:);

% Initial configuration bounds.
q0_lb = opt.qMin;
q0_ub = opt.qMax;

% Fix horizontal translation gauge.
q0_lb(1) = 0.0;
q0_ub(1) = 0.0;

% Floating-base height.
q0_lb(2) = opt.hipHeightMin;
q0_ub(2) = 1.15;

% Torso angle.
q0_lb(3) = -0.50;
q0_ub(3) =  0.50;

% Initial velocity bounds.
dq0_lb = -10 * ones(model.nq, 1);
dq0_ub =  10 * ones(model.nq, 1);

% Forward velocity.
dq0_lb(1) = 0.05;
dq0_ub(1) = 1.50;

% Vertical velocity.
dq0_lb(2) = -1.00;
dq0_ub(2) =  1.00;

% Torso angular velocity.
dq0_lb(3) = -6.00;
dq0_ub(3) =  6.00;

lb = [lb_CP; q0_lb; dq0_lb; opt.Tmin];
ub = [ub_CP; q0_ub; dq0_ub; opt.Tmax];

fprintf('\nInitial phase diagnostics:\n');
fprintf('thetaStart = %.6f\n', opt.thetaStart);
fprintf('thetaEnd   = %.6f\n', opt.thetaEnd);
fprintf('Delta theta = %.6f\n', opt.thetaEnd - opt.thetaStart);
fprintf('Initial dtheta = %.6f\n', dtheta0);

x0_init = [q0_init; dq0_init];
[y0, dy0, ~] = hzd_virtualConstraints(x0_init, CP0, model, opt);

fprintf('Initial ||y||  = %.3e\n', norm(y0));
fprintf('Initial ||dy|| = %.3e\n', norm(dy0));

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
