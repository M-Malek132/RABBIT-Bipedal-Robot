function [p, lb, ub, options] = hzd_problem_setup()
%HZD_PROBLEM_SETUP Single source of truth for the HZD gait-optimization
% problem: parameters p, decision-variable bounds lb/ub, and fmincon options.
% Used by BOTH rabbit_hzd_trajectory_optimization.m (single solve) and
% build_gait_library.m (speed continuation) so they cannot drift apart.

    %% ---- HZD PARAMETERS ----
    p.nu            = 4;
    p.n_coeffs      = 6;
    p.nq            = 7;
    p.bs_degree     = 3;
    p.Kp = 400;  p.Kd = 40;
    p.T_max         = 3;
    p.theta_plus    = 0.3;
    p.min_theta_gap = 0.15;     % reject x_start too close to theta_plus (denominator blowup)
    p.ground_tol    = 1e-3;     % ground penetration / foot-height tolerance

    % Baumgarte stabilization of the stance-foot contact constraint. The
    % acceleration-only constraint is a double integrator with no restoring
    % term, so the foot drifts ~0.1-0.3 m per step; these gains make position
    % and velocity violations decay as c_ddot + 2*alpha*c_dot + beta^2*c = 0
    % (critically damped, time constant 1/alpha = 0.05 s -- fast vs a ~1 s
    % step, but resolvable at the ode45 MaxStep of 0.01 s).
    p.baumgarte_alpha = 0;   % Baumgarte DISABLED (0 => original dynamics).
    p.baumgarte_beta  = 0;   % It fought the fixed virtual-constraint PD law.

    % Westervelt NIC3 / thesis "mid-step swing foot clearance": swing foot must
    % stay this far above ground through the middle of the step (thesis uses
    % 0.1 m on ATRIAS; 0.02 is a reachable first target here -- raise later).
    p.swing_clearance_min = 0.02;   % [m]

    % Cost is torque^2 PER STEP LENGTH (Westervelt eq. 6.43). Below this length
    % the divisor is clamped + a smooth quadratic penalty takes over (hzd_cost).
    p.min_step_length   = 0.05;     % [m]
    p.short_step_weight = 1e5;      % weight on (min_step_length - L_step)^2

    % DESIRED STEP LENGTH lower bound (hard inequality) -- without it the
    % optimizer collapses to a ~0.015 m scuffing "step in place".
    p.step_length_min = 0.30;       % [m]

    % DESIRED WALKING SPEED (Westervelt NEC1, eq. 6.50): equality
    % L_step / T_step = v_des. Default is the measured natural speed of the
    % seed gait; callers override p.v_des for continuation.
    p.v_des = 0.355;                % [m/s]

    %% ---- VARIABLE BOUNDS ----
    % Keep the robot in a sane configuration range (no self-collision / absurd
    % poses). CP = B-spline control points; x0 = initial state.
    n_cp  = p.nu * p.n_coeffs;
    lb_cp = -2.0 * ones(n_cp, 1);
    ub_cp =  2.0 * ones(n_cp, 1);

    % x_start = [px pz qt q1 q2 q3 q4  dpx dpz dqt dq1 dq2 dq3 dq4]
    lb_x0 = [ 0; -1.4; -0.6; -2.0; -2.0; -2.0; -2.0;  -6; -6; -6; -6; -6; -6; -6];
    ub_x0 = [ 0; -0.5;  0.6;  2.0;  2.0;  2.0;  2.0;   6;  6;  6;  6;  6;  6;  6];

    lb = [lb_cp; lb_x0];
    ub = [ub_cp; ub_x0];

    %% ---- FMINCON OPTIONS ----
    % Central finite differences (event-ODE makes forward FD noisy); scaling
    % OFF so reported Feasibility is comparable to ConstraintTolerance.
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
end
