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
opt.qMax = [5; 5; pi; 1.5; 1.5; 1.5; 1.5];

%% ---- Controller gains ---------------------------------------
opt.Kp = 200;
opt.Kd =  30;

%% ---- Initial guess -------------------------------------------
% Control points: (n+1) rows × 4 columns  [q1 q2 q3 q4]
q_act_start =    [-0.316969066015193
    0.740545951688271
   -0.737791239686444
    0.516687166027459];
q_act_end   = [ 0.20; -0.25; -0.20; -0.25];

n_pts = opt.n_bs + 1;
CP0 = zeros(n_pts, opt.ny);

for j = 1:opt.ny
    CP0(:, j) = linspace(q_act_start(j), q_act_end(j), n_pts);
end

CP0_vec = CP0(:);

%% ---- Initial guess from provided state -----------------------

x0_guess = [
    0.184960894261270
    0.913697468187791
    0.146428704138018
   -0.316969066015193
    0.740545951688271
   -0.737791239686444
    0.516687166027459
    0.676639780570254
   -0.209372372209695
    0.879620500001476
   -0.319771891690171
    0.392229976631628
    0.963605877665479
   -0.0270056036752210
];

q0_init  = x0_guess(1:model.nq);
dq0_init = x0_guess(model.nq+1:2*model.nq);

% Update actuated start to match the provided initial condition.
q_act_start = q0_init(4:7);

% Keep your chosen end shape, or make it approximately relabelled/opposite.
q_act_end = [ ...
    0.3000
   -1.0000
    0.6000
   -0.3000
];

% Rebuild B-spline control points using the provided initial q_act_start.
n_pts = opt.n_bs + 1;
CP0 = zeros(n_pts, opt.ny);

for j = 1:opt.ny
    CP0(:, j) = linspace(q_act_start(j), q_act_end(j), n_pts);
end

CP0_vec = CP0(:);

% Phase values from actual start and guessed end.
opt.thetaStart = -q0_init(3) - q0_init(4) - 0.5*q0_init(5);

qEnd_shape_guess = [
    q0_init(1) + opt.v_des * 0.50
    q0_init(2)
    q0_init(3)
    q_act_end
];

opt.thetaEnd = -qEnd_shape_guess(3) ...
               -qEnd_shape_guess(4) ...
               -0.5*qEnd_shape_guess(5);

% Step duration guess.
T0 = 0.50;

z0 = [CP0_vec; q0_init; dq0_init; T0];

fprintf('\nInitial guess check:\n');
disp('q0_init =');  disp(q0_init);
disp('dq0_init ='); disp(dq0_init);

fprintf('thetaStart = %.6f\n', opt.thetaStart);
fprintf('thetaEnd   = %.6f\n', opt.thetaEnd);
fprintf('Delta theta = %.6f\n', opt.thetaEnd - opt.thetaStart);

[y0, dy0, ~] = hzd_virtualConstraints(q0_init, dq0_init, CP0, model, opt);
fprintf('Initial ||y||  = %.3e\n', norm(y0));
fprintf('Initial ||dy|| = %.3e\n', norm(dy0));

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
q0_lb(2) = max(opt.hipHeightMin, q0_init(2) - 0.15);
q0_ub(2) = q0_init(2) + 0.15;

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
% fprintf('Initial dtheta = %.6f\n', dtheta0);

x0_init = [q0_init; dq0_init];
[y0, dy0, ~] = hzd_virtualConstraints(q0_init, dq0_init, CP0, model, opt);

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

function theta = hzd_theta(q)
%HZD_THETA Phase variable used in this optimization.
%
% theta = -q(3) - q(4) - 0.5*q(5)

    theta = -q(3) - q(4) - 0.5*q(5);
end


function pStance = getStanceFootPosition(q, model)
%GETSTANCEFOOTPOSITION Return stance-foot Cartesian position.
%
% This wrapper handles both possible rabbit_kinematics styles:
%
%   1. Struct output:
%        kin = rabbit_kinematics(q, params)
%        kin.stanceFoot
%
%   2. Multiple outputs:
%        [stanceFoot, swingFoot, ...] = rabbit_kinematics(q, p)
%
% If your rabbit_kinematics has a different output order, adjust this
% function only.

    params = model.params;

    try
        kin = rabbit_kinematics(q, params);

        if isstruct(kin)
            if isfield(kin, 'stanceFoot')
                pStance = kin.stanceFoot(:);
                return;
            elseif isfield(kin, 'stance_foot')
                pStance = kin.stance_foot(:);
                return;
            elseif isfield(kin, 'stance')
                pStance = kin.stance(:);
                return;
            else
                error('Struct output has no recognizable stance-foot field.');
            end
        end

    catch
        % Fall through to multiple-output version.
    end

    p = packParameters(params);

    [stanceFoot, ~, ~, ~, ~, ~] = rabbit_kinematics(q, p);

    pStance = stanceFoot(:);
end


function J = numericalFootJacobian(footFcn, q)
%NUMERICALFOOTJACOBIAN Numerically compute foot-position Jacobian.
%
% footFcn : function handle returning 2x1 foot position
% q       : nq x 1 configuration
%
% J       : 2 x nq Jacobian

    q = q(:);

    nq = numel(q);
    f0 = footFcn(q);
    nf = numel(f0);

    J = zeros(nf, nq);

    epsFD = 1e-6;

    for i = 1:nq
        dq = zeros(nq, 1);
        dq(i) = epsFD;

        fp = footFcn(q + dq);
        fm = footFcn(q - dq);

        J(:, i) = (fp - fm) / (2 * epsFD);
    end
end
