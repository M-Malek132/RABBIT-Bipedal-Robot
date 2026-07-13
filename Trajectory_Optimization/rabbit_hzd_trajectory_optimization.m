%% =========================================================
% HZD PARAMETERS  (these are YOUR optimization hyperparameters,
% separate from the robot's own compiled physics functions)
%% =========================================================
p.nu        = 4;
p.n_coeffs  = 6;
p.nq        = 7;
p.bs_degree = 3;
p.Kp = 400;  p.Kd = 40;
p.T_max = 10;                 % match simulate_one_step's own span
p.theta_plus = 0.3;

n_free_spline = p.nu * p.n_coeffs;
n_free_state  = 2*p.nq;
z0 = [randn(n_free_spline,1); make_random_initial_state(p)];

options = optimoptions('fmincon', 'Display','iter', 'Algorithm','sqp');
[z_opt, fval] = fmincon(...
    @(z) hzd_cost(z, p), ...
    z0, [], [], [], [], [], [], ...
    @(z) hzd_constraints(z, p), ...
    options);

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

    dhip_dq = [1 0 zeros(1,5); 0 1 zeros(1,5)];
    Jst     = J_st(q);              % 2x7, real repo function
    drel_dq = dhip_dq - Jst;

    g = ( rel(2)*drel_dq(1,:) - rel(1)*drel_dq(2,:) ) / r2;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% B-SPLINE EVAL (your repo's own BSpline / BSpline_derivative)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [b, db] = bspline_eval(c, s, p)
    n  = p.n_coeffs - 1;
    N  = BSpline(n, p.bs_degree, s);
    dN = BSpline_derivative(n, p.bs_degree, s);
    b  = c(:).' * N(:);
    db = c(:).' * dN(:);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GAIT SIMULATION — mirrors simulate_one_step, but with HZD control inline
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [x_end, total_torque_sq] = simulate_hzd_gait(coeffs, x_start, p)
    nq = p.nq;
    theta_minus = theta_of_q(x_start(1:nq));   % derived fresh from THIS x_start

    % Need theta_plus to define the phase window: since it's not physically
    % meaningful to fix in advance when x_start is free too, derive it
    % from the spline's own natural endpoint — see note below.
    theta_plus = p.theta_plus;   % TODO — see note under the code block

    xi0 = [x_start; 0];
    opts = odeset('RelTol',1e-6,'AbsTol',1e-6,'MaxStep',0.01, ...
                  'Events', @impact_event_wrapper);

    [~, XI] = ode45(@(t,xi) hzd_closed_loop_ode(t,xi,coeffs,theta_minus,theta_plus,p), ...
                     [0 p.T_max], xi0, opts);

    xi_end          = XI(end,:).';
    x_end           = xi_end(1:2*nq);
    total_torque_sq = xi_end(end);
end

function [value,isterminal,direction] = impact_event_wrapper(t,xi)
    nq = 7;
    [value,isterminal,direction] = rabbit_impact_event(t, xi(1:2*nq));  % real signature, no p
end

function dxi = hzd_closed_loop_ode(t,xi,coeffs,theta_minus,theta_plus,p)
    nq = p.nq;
    x  = xi(1:2*nq);
    q  = x(1:nq);
    dq = x(nq+1:end);

    theta     = theta_of_q(q);
    s         = (theta - theta_minus) / (theta_plus - theta_minus);
    s         = min(max(s,0),1);
    ds_dtheta = 1/(theta_plus - theta_minus);
    dtheta_dt = dtheta_dq_of(q) * dq;

    act_idx = 4:7;     % confirmed by relabel_state: stance_hip, stance_knee, swing_hip, swing_knee
    y  = q(act_idx);  dy = dq(act_idx);

    yd  = zeros(p.nu,1);
    dyd = zeros(p.nu,1);
    for i = 1:p.nu
        [b, db] = bspline_eval(coeffs(i,:), s, p);
        yd(i)  = b;
        dyd(i) = db * ds_dtheta * dtheta_dt;
    end

    e  = y - yd;  de = dy - dyd;
    tau = -p.Kp.*e - p.Kd.*de;

    ddq = rabbit_constrained_dynamics(q, dq, tau);   % real signature: returns [ddq, lambda], no p

    dxi = [dq; ddq; sum(tau.^2)];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function J = hzd_cost(z, p)
    [coeffs, x_start] = unpack_z(z, p);
    [~, total_torque_sq] = simulate_hzd_gait(coeffs, x_start, p);
    J = total_torque_sq;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [c, ceq] = hzd_constraints(z, p)
    [coeffs, x_start] = unpack_z(z, p);
    nq = p.nq;

    [x_end, ~] = simulate_hzd_gait(coeffs, x_start, p);
    q_end = x_end(1:nq);

    % Periodicity (unchanged, excludes world px)
    x_next = rabbit_reset_map(rabbit_impact_map(x_end));
    idx_periodic = 2:14;
    ceq_periodicity = x_next(idx_periodic) - x_start(idx_periodic);

    % Ground contact of x_start
    q_start  = x_start(1:nq);
    dq_start = x_start(nq+1:end);
    foot0 = P_st(q_start);
    ceq_foot_height = foot0(2);
    Jst = J_st(q_start);
    ceq_foot_vel = Jst * dq_start;

    % NEW: force impact to occur exactly at theta_plus
    % (i.e., spline's s=1 endpoint coincides with the physical touchdown)
    ceq_theta_end = theta_of_q(q_end) - p.theta_plus;

    ceq = [ceq_periodicity; ceq_foot_height; ceq_foot_vel; ceq_theta_end];
    c = [];
end