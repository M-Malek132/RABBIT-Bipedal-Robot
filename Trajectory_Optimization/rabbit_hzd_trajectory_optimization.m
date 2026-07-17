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

%% =========================================================
% INITIAL STATE — random, projected + rejection-sampled until valid
%% =========================================================
x0 = make_random_initial_state(p);
r  = check_ground_validity(x0(1:7), x0(8:14));
if ~r.is_valid
    disp(r.violations);
end

%% =========================================================
% INITIAL SPLINE COEFFICIENTS
% Smooth interpolation from x0's actuated joints to a placeholder
% touchdown pose. Also guarantees spline(s=0) == x0's joints exactly,
% since a clamped B-spline's first control point IS its s=0 value.
%% =========================================================
q_desired = zeros(p.nu, p.n_coeffs);
for i = 1:p.nu
    q_start_i = x0(i+3);          % act_idx = 4:7, offset by 3
    q_end_i   = q_start_i * 0.5;  % placeholder target — tune per joint
    q_desired(i,:) = linspace(q_start_i, q_end_i, p.n_coeffs);
end
coeffs0 = q_desired;
z0 = [coeffs0(:); x0];

%% =========================================================
% OPTIMIZATION
%% =========================================================
options = optimoptions('fmincon', ...
    'Display','iter', ...
    'Algorithm','sqp', ...
    'MaxFunctionEvaluations', 10000, ...
    'MaxIterations', 2000, ...
    'OptimalityTolerance', 1e-3, ...
    'ConstraintTolerance', 1e-2, ...
    'StepTolerance', 1e-6, ...
    'ScaleProblem', true);

[z_opt, fval, exitflag] = fmincon(...
    @(z) hzd_cost(z, p), ...
    z0, [], [], [], [], [], [], ...
    @(z) hzd_constraints(z, p), ...
    options);

animate_hzd_result(z_opt, p, 10);
