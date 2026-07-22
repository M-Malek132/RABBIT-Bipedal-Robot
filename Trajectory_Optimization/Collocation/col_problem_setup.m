function [p, lb, ub, options] = col_problem_setup(N, p_over)
%COL_PROBLEM_SETUP  Parameters, bounds and fmincon options for the HZD
% direct-collocation solver. Inherits every field of hzd_problem_setup (so the
% two solvers share cost weights, speed target, physical limits, spline
% degree, ...) and adds the node count p.N plus the collocation-vector bounds.
%
%   [p, lb, ub, options] = col_problem_setup()          % default N = 20 nodes
%   [p, lb, ub, options] = col_problem_setup(N)
%   [p, lb, ub, options] = col_problem_setup(N, p_over) % override p fields
%
% p_over is an optional struct whose fields overwrite the inherited p BEFORE
% the bounds are built -- so toggling e.g. enforce_torque also updates the
% torque box on U. Decision-vector order matches col_pack: [X;U;Lam;coeffs;T].

    if nargin < 1 || isempty(N), N = 20; end

    p   = hzd_problem_setup();     % single source of truth for shared params

    if nargin >= 2 && ~isempty(p_over)
        fn = fieldnames(p_over);
        for k = 1:numel(fn), p.(fn{k}) = p_over.(fn{k}); end
    end
    p.N = N;                       % N always wins over any p_over.N

    nq = p.nq;  nu = p.nu;  nc = p.n_coeffs;  nx = 2*nq;

    % ---- per-node state bounds [q; dq] ----------------------------------
    %  q = [hip_x hip_z torso q1 q2 q3 q4]; dq likewise. hip_x is free to
    %  advance during the step (unlike the shooting x0 which fixed it at 0).
    lb_x = [-1.0; -1.4; -0.6; -2.0; -2.0; -2.0; -2.0;  -12*ones(nq,1)];
    ub_x = [ 1.0; -0.5;  0.6;  2.0;  2.0;  2.0;  2.0;   12*ones(nq,1)];

    % ---- per-node torque bounds (this is the Table-3.1 torque box) -------
    if isfield(p,'enforce_torque') && p.enforce_torque
        ulim = p.u_max;
    else
        ulim = 1000;                 % effectively unbounded
    end
    lb_u = -ulim*ones(nu,1);   ub_u = ulim*ones(nu,1);

    % ---- per-node contact-force bounds [Fx; Fz] -------------------------
    lb_l = [-3000; -1000];     ub_l = [3000; 4000];

    % ---- spline coeffs + step time --------------------------------------
    lb_c = -2.0*ones(nu*nc,1); ub_c = 2.0*ones(nu*nc,1);
    lb_T = 0.10;               ub_T = p.T_max;

    lb = [repmat(lb_x,N,1); repmat(lb_u,N,1); repmat(lb_l,N,1); lb_c; lb_T];
    ub = [repmat(ub_x,N,1); repmat(ub_u,N,1); repmat(ub_l,N,1); ub_c; ub_T];

    % ---- fmincon options -------------------------------------------------
    % Central FD (the transcription is smooth, so this is accurate) and
    % ScaleProblem off so reported feasibility matches ConstraintTolerance.
    % No ODE in the loop, so each evaluation is cheap despite the large vector.
    options = optimoptions('fmincon', ...
        'Display','iter', ...
        'Algorithm','sqp', ...
        'MaxFunctionEvaluations', 3e5, ...
        'MaxIterations', 1500, ...
        'OptimalityTolerance', 1e-3, ...
        'ConstraintTolerance', 1e-4, ...
        'StepTolerance', 1e-8, ...
        'ScaleProblem', false, ...
        'FiniteDifferenceType', 'central', ...
        'FiniteDifferenceStepSize', 1e-5);
end
