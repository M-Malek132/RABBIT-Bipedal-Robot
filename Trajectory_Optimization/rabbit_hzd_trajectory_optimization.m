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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [coeffs, x_start] = unpack_z(z, p)
    n1 = p.nu * p.n_coeffs;
    coeffs  = reshape(z(1:n1), [p.nu, p.n_coeffs]);
    x_start = z(n1+1:end);
end

function [value,isterminal,direction] = impact_event_wrapper(t,xi)
    nq = 7;
    [value,isterminal,direction] = rabbit_impact_event(t, xi(1:2*nq));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% COST — penalizes only genuine simulation failure. Equality constraints
% are handled by fmincon's nonlcon, not duplicated here.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function J = hzd_cost(z, p)
    [coeffs, x_start] = unpack_z(z, p);
    [~, total_torque_sq, ~, status] = simulate_hzd_gait(coeffs, x_start, p);

    if status < 0
        J = 1e6;   % failed / incomplete step — steer optimizer away
    else
        J = total_torque_sq;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% COST WITH PENALTY FOR CONSTRAINT VIOLATIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function J = hzd_cost_with_penalty(z, p)
    [coeffs, x_start] = unpack_z(z, p);
    [~, total_torque_sq, status] = simulate_hzd_gait(coeffs, x_start, p);
    
    % Compute constraint violation penalty
    [c, ceq] = hzd_constraints(z, p);
    penalty = 1000 * sum(ceq.^2);  % Quadratic penalty for constraint violations
    
    % Add penalty for failed simulations
    if status < 0
        J = 1e6 + penalty;
    else
        J = total_torque_sq + penalty;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [coeffs, x_start] = unpack_z(z, p)
    n1 = p.nu * p.n_coeffs;
    coeffs  = reshape(z(1:n1), [p.nu, p.n_coeffs]);
    x_start = z(n1+1:end);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PHASE VARIABLE: relative stance-foot -> hip vector angle
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function th = theta_of_q(q)
    hip  = q(1:2);
    foot = P_st(q);
    rel  = hip - foot;
    th   = atan2(rel(1), rel(2));
end

function g = dtheta_dq_of(q)
    hip  = q(1:2);
    foot = P_st(q);
    rel  = hip - foot;
    r2   = rel(1)^2 + rel(2)^2;
    
    % Avoid division by zero
    if r2 < 1e-10
        r2 = 1e-10;
    end

    dhip_dq = [1 0 zeros(1,5); 0 1 zeros(1,5)];
    Jst     = J_st(q);              % 2x7, real repo function
    drel_dq = dhip_dq - Jst;

    g = ( rel(2)*drel_dq(1,:) - rel(1)*drel_dq(2,:) ) / r2;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% B-SPLINE EVAL
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [b, db] = bspline_eval(c, s, p)
    n  = p.n_coeffs - 1;
    % Clamp s to [0,1] to avoid numerical issues
    s = max(0, min(1, s));
    
    try
        N  = BSpline(n, p.bs_degree, s);
        dN = BSpline_derivative(n, p.bs_degree, s);
    catch
        % Fallback for numerical issues
        N = zeros(n+1, 1);
        dN = zeros(n+1, 1);
    end
    
    b  = c(:).' * N(:);
    db = c(:).' * dN(:);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GAIT SIMULATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [x_end, total_torque_sq, max_penetration] = simulate_hzd_gait(coeffs, x_start, p)
    nq = p.nq;
    theta_minus = theta_of_q(x_start(1:nq));
    theta_plus  = p.theta_plus;

    xi0 = [x_start; 0];
    opts = odeset('RelTol',1e-6,'AbsTol',1e-6,'MaxStep',0.01, ...
                  'Events', @impact_event_wrapper);

    sol = ode45(@(t,xi) hzd_closed_loop_ode(t,xi,coeffs,theta_minus,theta_plus,p), ...
                [0 p.T_max], xi0, opts);

    xi_end          = sol.y(:,end);
    x_end           = xi_end(1:2*nq);
    total_torque_sq = xi_end(end);

    % Dense-sample the WHOLE step to find worst ground clipping, not just endpoints
    t_samples  = linspace(sol.x(1), sol.x(end), 60);
    xi_samples = deval(sol, t_samples);

    max_penetration = -inf;
    for k = 1:length(t_samples)
        q_k  = xi_samples(1:nq, k);
        pts  = get_body_points(q_k);
        max_penetration = max([max_penetration, pts.torso_top(2), ...
                                pts.stance_knee(2), pts.swing_knee(2), pts.swing_foot(2)]);
    end
end

function [value,isterminal,direction] = impact_event_wrapper(t,xi)
    nq = 7;
    [value,isterminal,direction] = rabbit_impact_event(t, xi(1:2*nq));
end

function dxi = hzd_closed_loop_ode(t,xi,coeffs,theta_minus,theta_plus,p)
    nq = p.nq;
    x  = xi(1:2*nq);
    q  = x(1:nq);
    dq = x(nq+1:end);

    theta     = theta_of_q(q);
    s         = (theta - theta_minus) / (theta_plus - theta_minus);
    s         = max(0, min(1, s));
    ds_dtheta = 1/(theta_plus - theta_minus);
    dtheta_dt = dtheta_dq_of(q) * dq;
    
    % Clamp dtheta_dt
    if abs(dtheta_dt) > 100
        dtheta_dt = sign(dtheta_dt) * 100;
    end

    act_idx = 4:7;
    y  = q(act_idx);
    dy = dq(act_idx);

    yd  = zeros(p.nu,1);
    dyd = zeros(p.nu,1);
    for i = 1:p.nu
        [b, db] = bspline_eval(coeffs(i,:), s, p);
        yd(i)  = b;
        dyd(i) = db * ds_dtheta * dtheta_dt;
    end

    e  = y - yd;
    de = dy - dyd;
    
    % Clamp errors
    e = max(-10, min(10, e));
    de = max(-100, min(100, de));
    
    tau = -p.Kp.*e - p.Kd.*de;
    tau = max(-1000, min(1000, tau));

    ddq = rabbit_constrained_dynamics(q, dq, tau);
    ddq = max(-1000, min(1000, ddq));

    dxi = [dq; ddq; sum(tau.^2)];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLAYBACK VERSION — same dynamics as simulate_hzd_gait, but returns
% the full trajectory instead of just the final state
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [t_out, x_out] = simulate_hzd_gait_full(coeffs, x_start, p)
    nq = p.nq;
    theta_minus = theta_of_q(x_start(1:nq));
    theta_plus  = p.theta_plus;

    xi0 = [x_start; 0];
    opts = odeset('RelTol',1e-6,'AbsTol',1e-6,'MaxStep',0.01, ...
                  'Events', @impact_event_wrapper);

    [t_out, XI] = ode45(@(t,xi) hzd_closed_loop_ode(t,xi,coeffs,theta_minus,theta_plus,p), ...
                         [0 p.T_max], xi0, opts);

    x_out = XI(:, 1:2*nq);   % drop the appended torque-cost column
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ANIMATE THE OPTIMIZED GAIT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function animate_hzd_result(z_opt, p, nSteps)
    if nargin < 3, nSteps = 10; end

    [coeffs, x_current] = unpack_z(z_opt, p);

    t_all = [];
    x_all = [];
    time_offset = 0;

    for step = 1:nSteps
        [t_step, x_step] = simulate_hzd_gait_full(coeffs, x_current, p);

        if isempty(t_all)
            t_all = t_step;
            x_all = x_step;
        else
            t_all = [t_all; t_step(2:end) + time_offset];
            x_all = [x_all; x_step(2:end, :)];
        end

        x_minus        = x_step(end, :)';
        x_after_impact = rabbit_impact_map(x_minus);
        x_current      = rabbit_reset_map(x_after_impact);

        time_offset = t_all(end);
    end

    animate_rabbit(x_all);   % your existing function — expects 14 x N or N x 14, handles both
end

function pts = get_body_points(q)
    % Returns struct of all key body points' [x; z] world positions,
    % using the same transform chain as animate_rabbit.m's rabbit_points
    pos_hip   = Tt(q) * [0 0 0 1]';
    P_torso   = Tt(q) * [0 -0.75 0 1]';
    P_knee_st = T2(q) * [0 0 0 1]';
    P_knee_sw = T4(q) * [0 0 0 1]';

    pts.hip         = [pos_hip(1,1);   pos_hip(3,1)];
    pts.torso_top   = [P_torso(1,1);   P_torso(3,1)];
    pts.stance_knee = [P_knee_st(1,1); P_knee_st(3,1)];
    pts.swing_knee  = [P_knee_sw(1,1); P_knee_sw(3,1)];
    pts.stance_foot = P_st(q);
    pts.swing_foot  = P_sw(q);
end

function report = check_ground_validity(q, dq, tol)
    if nargin < 3, tol = 1e-3; end
    pts = get_body_points(q);

    % Recall: downward-positive pz axis. Ground = 0. Above ground = pz <= 0.
    % A point with pz > tol is PENETRATING the ground.
    fields = {'torso_top','stance_knee','swing_knee','swing_foot'};
    report.max_penetration = -inf;
    report.violations = {};
    for i = 1:numel(fields)
        h = pts.(fields{i})(2);
        if h > tol
            report.violations{end+1} = sprintf('%s penetrating ground: z=%.4f', fields{i}, h);
        end
        report.max_penetration = max(report.max_penetration, h);
    end

    report.stance_foot_height = pts.stance_foot(2);
    if abs(report.stance_foot_height) > tol
        report.violations{end+1} = sprintf('stance foot NOT on ground: z=%.4f', report.stance_foot_height);
    end

    if nargin >= 2 && ~isempty(dq)
        J = J_st(q);
        v = J*dq;
        report.stance_foot_velocity = norm(v);
        if norm(v) > 1e-2
            report.violations{end+1} = sprintf('stance foot slipping: |v|=%.4f', norm(v));
        end
    end

    report.is_valid = isempty(report.violations);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function J = hzd_cost(z, p)
    [coeffs, x_start] = unpack_z(z, p);
    [~, total_torque_sq, ~] = simulate_hzd_gait(coeffs, x_start, p);
    J = total_torque_sq;
end

function [c, ceq] = hzd_constraints(z, p)
    [coeffs, x_start] = unpack_z(z, p);
    nq = p.nq;

    [x_end, ~, max_penetration] = simulate_hzd_gait(coeffs, x_start, p);

    x_next = rabbit_reset_map(rabbit_impact_map(x_end));
    idx_periodic = 2:14;
    ceq_periodicity = x_next(idx_periodic) - x_start(idx_periodic);

    q_start  = x_start(1:nq);
    dq_start = x_start(nq+1:end);
    foot0 = P_st(q_start);
    ceq_foot_height = foot0(2);
    Jst = J_st(q_start);
    ceq_foot_vel = Jst * dq_start;

    q_end = x_end(1:nq);
    ceq_theta_end = theta_of_q(q_end) - p.theta_plus;

    ceq = [ceq_periodicity; ceq_foot_height; ceq_foot_vel; ceq_theta_end];

    % NEW: inequality — nothing may penetrate the ground during the step
    tol = 1e-3;
    c = max_penetration - tol;   % fmincon requires c <= 0
end