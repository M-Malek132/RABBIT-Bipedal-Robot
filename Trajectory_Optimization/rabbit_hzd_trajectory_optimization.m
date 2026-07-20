%% =========================================================
% HZD PARAMETERS
%% =========================================================
p.nu            = 4;
p.n_coeffs      = 6;
p.nq            = 7;
p.bs_degree     = 3;
p.Kp = 400;  p.Kd = 40;
p.T_max         = 3;        % reduced for faster debugging
p.theta_plus    = 0.3;
p.min_theta_gap = 0.15;     % reject x_start too close to theta_plus (denominator blowup)
p.ground_tol    = 1e-3;     % ground penetration / foot-height tolerance

% Westervelt NIC3 / thesis Table 3.1 "mid-step swing foot clearance": the
% swing foot must stay this far above the ground through the middle of the
% step, so the impact surface is reached only at the END of the step.
% (Thesis uses 0.1 m on ATRIAS. Started at 0.05 here and the gait could only
% reach ~0.014, so that was infeasible from the initial guess; 0.02 is a
% reachable first target. Raise it once the gait closes.)
p.swing_clearance_min = 0.02;   % [m]

% Cost is torque^2 PER STEP LENGTH (Westervelt eq. 6.43). Below this length
% the divisor is clamped and a smooth quadratic penalty takes over -- see
% hzd_cost.m for why this must NOT be a hard cutoff.
p.min_step_length   = 0.05;     % [m]
p.short_step_weight = 1e5;      % weight on (min_step_length - L_step)^2

% DESIRED STEP LENGTH lower bound (Westervelt NEC1 / thesis "desired step
% length"), enforced as a hard inequality in hzd_constraints. Without it the
% optimizer collapses to a ~0.015 m scuffing "step in place". Earlier runs
% reached L_step ~ 0.49 m easily, so 0.30 is comfortably achievable.
p.step_length_min = 0.30;       % [m]

% DESIRED WALKING SPEED (Westervelt NEC1 average walking rate, eq. 6.50):
% enforced as an EQUALITY  L_step / T_step = v_des  in hzd_constraints, which
% PINS the walking speed instead of leaving it a byproduct of the dynamics.
% NOTE: set this near the speed the UNCONSTRAINED gait already walks at, or
% the equality will fight the rest of the problem. inspect_solution now prints
% the achieved speed -- run it on the current z_opt first and match v_des to
% that (thesis uses 0.5 m/s for ATRIAS).
p.v_des = 0.5;                  % [m/s]

%% =========================================================
% INITIAL STATE — natural post-impact start-of-step pose
% (px pz qt q1 q2 q3 q4  dpx dpz dqt dq1 dq2 dq3 dq4)
% Constructed so that: stance foot on ground; SWING FOOT BEHIND the
% stance foot (-0.28 m) with 6 cm clearance; theta = -0.15 rad (hip
% behind stance foot, so phase sweeps -0.15 -> theta_plus=0.30 forward);
% both knees bent naturally; stance-foot velocity = 0 at t=0.
%% =========================================================
x0 = [ +0.0000; -0.9310; +0.2999; -0.7934; +0.6869; -0.6318; +0.9810; ...
       +0.3952; -0.0419; +0.1847; +0.1847; +0.1045; +0.0000; +0.0000];

r  = check_ground_validity(x0(1:7), x0(8:14));
if ~r.is_valid
    disp(r.violations);
end

%% =========================================================
% INITIAL SPLINE COEFFICIENTS  (see make_coeffs.m for rationale)
%% =========================================================
coeffs0 = make_coeffs(x0, p);

%% =========================================================
% SETTLE THE INITIAL STATE (Poincare / fixed-point iteration)
% DISABLED (n_settle = 0) -- tried and it makes things WORSE.
%
% Iterating the step map only converges if the gait is already near a
% STABLE periodic orbit. coeffs0 is a crude linear interpolation, far
% from periodic, so the iteration diverges instead of settling: theta
% blew up to +/-2..3 rad (robot flipping over) within 2 steps and the
% resulting x0 violated the bounds and ground constraints.
%
% Worth revisiting ONLY once fmincon can produce a near-periodic gait --
% then settling could polish it. Set n_settle > 0 to re-enable.
%% =========================================================
n_settle = 0;
if n_settle > 0
    [x0, settle_hist] = settle_initial_state(coeffs0, x0, p, n_settle);
    coeffs0 = make_coeffs(x0, p);   % rebuild so spline(s=0) == settled joints
end

z0 = [coeffs0(:); x0];

%% =========================================================
% WARM START (optional but strongly recommended once a gait exists)
% Set warm_start_file to a previously saved result to start fmincon from
% that gait instead of the cold x0. Cold starts from x0 are a gamble (even
% the first success wandered for 100+ iterations); warm starting from a
% known periodic gait is what makes speed continuation / a gait library
% work: solve one gait, nudge p.v_des, re-solve from the previous solution.
%
% The saved p may predate newer fields (e.g. v_des); we keep the CURRENT p
% so all constraints use current settings, and only borrow z_opt as z0.
%% =========================================================
warm_start_file = 'hzd_result_2026-07-20_12-06-36.mat';   % '' = cold start from x0
if ~isempty(warm_start_file)
    % Resolve relative to the repo's Results/ folder, independent of cwd
    % (same convention as the SAVE block below).
    ws_path = fullfile(fileparts(mfilename('fullpath')), '..', 'Results', warm_start_file);
    S  = load(ws_path);
    z0 = S.z_opt;
    fprintf('Warm start from %s\n', ws_path);
    fprintf('--- starting gait evaluated under CURRENT p (check its speed!) ---\n');
    inspect_solution(z0, p);
end

%% =========================================================
% VARIABLE BOUNDS — with no bounds fmincon is free to drive the torso
% and legs into self-colliding or otherwise nonphysical poses. (These
% also used to guard against theta_of_q's atan2 branch cut; that is no
% longer a concern now that the phase variable is linear, but keeping
% the robot in a sane configuration range is still worth doing.)
%% =========================================================
n_cp = p.nu * p.n_coeffs;
lb_cp = -2.0 * ones(n_cp, 1);
ub_cp =  2.0 * ones(n_cp, 1);

% x_start = [px pz qt q1 q2 q3 q4  dpx dpz dqt dq1 dq2 dq3 dq4]
lb_x0 = [   0; -1.4; -0.6; -2.0; -2.0; -2.0; -2.0;  -6; -6; -6; -6; -6; -6; -6];
ub_x0 = [   0; -0.5;  0.6;  2.0;  2.0;  2.0;  2.0;   6;  6;  6;  6;  6;  6;  6];

lb = [lb_cp; lb_x0];
ub = [ub_cp; ub_x0];

%% =========================================================
% OPTIMIZATION
%% =========================================================
% The cost/constraints run through an event-triggered ode45 (impact
% detection), so a small perturbation of z shifts the impact time and the
% integrator step pattern -- the DEFAULT forward finite differences see
% that as noise and return bad gradients, which makes SQP thrash (feasibility
% goes non-monotonic, optimality explodes). Central differences with a step
% large enough to average over that noise give much cleaner gradients.
options = optimoptions('fmincon', ...
    'Display','iter', ...
    'Algorithm','sqp', ...
    'MaxFunctionEvaluations', 10000, ...
    'MaxIterations', 2000, ...
    'OptimalityTolerance', 1e-3, ...
    'ConstraintTolerance', 1e-2, ...
    'StepTolerance', 1e-6, ...
    'ScaleProblem', false, ...
    'FiniteDifferenceType', 'central', ...
    'FiniteDifferenceStepSize', 1e-3);
% ScaleProblem turned OFF deliberately: with it on, fmincon derived its
% constraint scaling from the failure sentinel in hzd_constraints and then
% reported "constraints satisfied" while the TRUE max violation was ~1.57
% (157x ConstraintTolerance). All decision variables here are O(1), so
% scaling buys little, and off means the reported Feasibility is directly
% comparable to ConstraintTolerance and to inspect_solution's numbers.

[z_opt, fval, exitflag] = fmincon(...
    @(z) hzd_cost(z, p), ...
    z0, [], [], [], [], lb, ub, ...
    @(z) hzd_constraints(z, p), ...
    options);

report = inspect_solution(z_opt, p);

%% =========================================================
% SAVE RESULT
%% =========================================================
results_dir = fullfile(fileparts(mfilename('fullpath')), '..', 'Results');
if ~exist(results_dir, 'dir')
    mkdir(results_dir);
end
ts = datestr(now, 'yyyy-mm-dd_HH-MM-SS');
result_file = fullfile(results_dir, ['hzd_result_', ts, '.mat']);
save(result_file, 'z_opt', 'p', 'fval', 'exitflag', 'report');
fprintf('Result saved to: %s\n', result_file);

animate_hzd_result(z_opt, p, 10);
