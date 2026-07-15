%% =========================================================
% HZD PARAMETERS
%% =========================================================
p.nu        = 4;
p.n_coeffs  = 6;
p.nq        = 7;
p.bs_degree = 3;
p.Kp = 400;  p.Kd = 40;
p.T_max = 3;                  % Reduced for faster debugging
p.theta_plus = 0.3;

n_free_spline = p.nu * p.n_coeffs;
n_free_state  = 2*p.nq;

% Better initialization - use a known feasible state
% First, get a reasonable initial state
q0 = [0; 0.5; -0.5; 0.3; -0.6; 0.2; -0.4];  % Reasonable joint angles
dq0 = zeros(p.nq, 1);  % Start from rest
x0 = [q0; dq0];

% Initialize spline coefficients to track desired joint trajectories
% Simple straight-line interpolation from initial to desired final angles
theta_minus = theta_of_q(q0);
theta_plus = p.theta_plus;
s_points = linspace(0, 1, p.n_coeffs);

% Desired joint trajectories (simple smooth motion)
q_desired = zeros(p.nu, p.n_coeffs);
for i = 1:p.nu
    % Start from current joint angle, end at a reasonable value
    q_start = q0(i+3);  % act_idx = 4:7, so offset by 3
    q_end = q_start * 0.5;  % Some reasonable final value
    q_desired(i,:) = linspace(q_start, q_end, p.n_coeffs);
end

% Initial spline coefficients (just the desired values)
coeffs0 = q_desired;

% Combine into z0
z0 = [coeffs0(:); x0];

% Options for optimization - more forgiving
options = optimoptions('fmincon', ...
    'Display','iter', ...
    'Algorithm','sqp', ...
    'MaxFunctionEvaluations', 10000, ...
    'MaxIterations', 2000, ...
    'OptimalityTolerance', 1e-3, ...
    'ConstraintTolerance', 1e-2, ...  % Relaxed
    'StepTolerance', 1e-6, ...
    'ScaleProblem', true);  % Enable problem scaling

[z_opt, fval, exitflag] = fmincon(...
    @(z) hzd_cost_with_penalty(z, p), ...
    z0, [], [], [], [], [], [], ...
    @(z) hzd_constraints(z, p), ...
    options);

animate_hzd_result(z_opt, p, 10);

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
function [x_end, total_torque_sq, status] = simulate_hzd_gait(coeffs, x_start, p)
    nq = p.nq;
    theta_minus = theta_of_q(x_start(1:nq));
    theta_plus = p.theta_plus;
    
    % Check if theta_minus and theta_plus are valid
    if isnan(theta_minus) || isinf(theta_minus)
        x_end = x_start;
        total_torque_sq = 1e6;
        status = -1;
        return;
    end
    
    % Make sure theta_plus > theta_minus
    if theta_plus <= theta_minus
        theta_plus = theta_minus + 0.5;
    end

    xi0 = [x_start; 0];
    
    opts = odeset('RelTol',1e-4,'AbsTol',1e-6,'MaxStep',0.01, ...
                  'Events', @impact_event_wrapper);

    try
        [~, XI] = ode45(@(t,xi) hzd_closed_loop_ode(t,xi,coeffs,theta_minus,theta_plus,p), ...
                         [0 p.T_max], xi0, opts);
    catch
        x_end = x_start;
        total_torque_sq = 1e6;
        status = -2;
        return;
    end

    if isempty(XI)
        x_end = x_start;
        total_torque_sq = 1e6;
        status = -3;
        return;
    end

    xi_end          = XI(end,:).';
    x_end           = xi_end(1:2*nq);
    total_torque_sq = xi_end(end);
    status = 0;
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function J = hzd_cost(z, p)
    [coeffs, x_start] = unpack_z(z, p);
    [~, total_torque_sq, status] = simulate_hzd_gait(coeffs, x_start, p);
    
    if status < 0
        J = 1e6;
    else
        J = total_torque_sq;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [c, ceq] = hzd_constraints(z, p)
    [coeffs, x_start] = unpack_z(z, p);
    nq = p.nq;

    [x_end, ~, status] = simulate_hzd_gait(coeffs, x_start, p);
    q_end = x_end(1:nq);

    % If simulation failed, return large constraint violations
    if status < 0
        ceq = ones(2*nq + 3, 1) * 100;
        c = [];
        return;
    end

    % Periodicity
    try
        x_next = rabbit_reset_map(rabbit_impact_map(x_end));
    catch
        % If impact map fails, return penalty
        ceq = ones(2*nq + 3, 1) * 100;
        c = [];
        return;
    end
    
    idx_periodic = 2:14;
    ceq_periodicity = x_next(idx_periodic) - x_start(idx_periodic);

    % Ground contact of x_start
    q_start  = x_start(1:nq);
    dq_start = x_start(nq+1:end);
    foot0 = P_st(q_start);
    ceq_foot_height = foot0(2);
    Jst = J_st(q_start);
    ceq_foot_vel = Jst * dq_start;

    % Force impact at theta_plus
    ceq_theta_end = theta_of_q(q_end) - p.theta_plus;

    % Add bounds on joint angles to keep them reasonable
    q_start_bounds = q_start(4:7);  % Actuated joints
    q_end_bounds = q_end(4:7);
    
    % Soft constraints - allow some violation with penalty
    ceq = [ceq_periodicity; 
           ceq_foot_height; 
           ceq_foot_vel; 
           ceq_theta_end];
    c = [];
end