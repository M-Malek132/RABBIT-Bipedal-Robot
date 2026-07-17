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
% INITIAL STATE — hand-picked natural-looking standing pose
% (px pz qt q1 q2 q3 q4  dpx dpz dqt dq1 dq2 dq3 dq4)
%% =========================================================
x0 = [ 0; -0.9163; 0.2355; -0.5144; 0.7922; -0.9752; 0.8873; ...
      -0.1879; -0.0102; -0.7345; 0.5003; 0.0611; 0.4137; -0.3175];

r  = check_ground_validity(x0(1:7), x0(8:14));
if ~r.is_valid
    disp(r.violations);
end

%% =========================================================
% INITIAL SPLINE COEFFICIENTS
% Smooth interpolation from x0's actuated joints to a periodicity-
% motivated touchdown target: for a periodic step, the pose right
% before impact (pre-relabel) should look like the current pose with
% stance/swing roles swapped -- exactly what relabel_state computes.
% Also guarantees spline(s=0) == x0's joints exactly, since a clamped
% B-spline's first control point IS its s=0 value.
%% =========================================================
x0_relabeled = relabel_state(x0);
q_end_guess  = x0_relabeled(4:7);

q_desired = zeros(p.nu, p.n_coeffs);
for i = 1:p.nu
    q_start_i = x0(i+3);          % act_idx = 4:7, offset by 3
    q_end_i   = q_end_guess(i);
    q_desired(i,:) = linspace(q_start_i, q_end_i, p.n_coeffs);
end
coeffs0 = q_desired;
z0 = [coeffs0(:); x0];

%% =========================================================
% VARIABLE BOUNDS — theta_of_q is atan2-based (wraps at +/-pi), and
% with no bounds fmincon is free to drive the torso/legs into
% self-colliding or wrapped-around poses, which is what produced the
% huge theta_end constraint violation and exploding gradients seen
% without bounds. These are generous (they comfortably contain what
% make_random_initial_state actually samples) but keep the search out
% of physically nonsensical territory.
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
