%% =========================================================
% HZD PARAMETERS
%% =========================================================
p.nu       = 4;      % actuated joints: q1,q2,q3,q4
p.n_coeffs = 6;       % Bezier control points per joint (degree 5)
p.nq       = 7;
p.theta_plus = ...;   % FIXED design target: desired theta(q) at impact (sets nominal step size — tune this)
p.T_max      = 2.0;   % max integration time per step (safety cap for the event detector)
p.Kp = 400;  p.Kd = 40;  % PD gains tracking the virtual constraints (tune)
% ... plus whatever fields your dynamics/kinematics functions need (link lengths, masses, etc.)

n_free_bezier = p.nu * p.n_coeffs;
n_free_state  = 2*p.nq;                 % x_start is now part of z
z0 = [randn(n_free_bezier,1); x0_guess]; % x0_guess: reuse your main_demo.m style initial guess

%% =========================================================
% OPTIMIZATION
%% =========================================================
options = optimoptions('fmincon', 'Display','iter', 'Algorithm','sqp');
[z_opt, fval] = fmincon(...
    @(z) hzd_cost(z, p), ...
    z0, [], [], [], [], [], [], ...
    @(z) hzd_constraints(z, p), ...
    options);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% UNPACK z
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [coeffs, x_start] = unpack_z(z, p)
    n1 = p.nu * p.n_coeffs;
    coeffs  = reshape(z(1:n1), [p.nu, p.n_coeffs]);
    x_start = z(n1+1:end);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PHASE VARIABLE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function th = theta_of_q(q)
    % Absolute angle of the line from the (pinned, origin) stance foot to the hip
    th = atan2(q(1), q(2));   % q(1)=px, q(2)=pz
end

function g = dtheta_dq_of(q)
    % Row vector d(theta)/dq, nonzero only in px,pz
    px = q(1); pz = q(2); r2 = px^2 + pz^2;
    g = zeros(1, numel(q));
    g(1) =  pz / r2;
    g(2) = -px / r2;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% BEZIER (degree n_coeffs-1)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [b, db] = bezier_eval(c, s)
    c = c(:).';
    M = numel(c) - 1;
    k = 0:M;
    binom = arrayfun(@(kk) nchoosek(M,kk), k);
    b = sum(c .* binom .* s.^k .* (1-s).^(M-k));

    if M >= 1
        dc = diff(c);
        Mm1 = M - 1;
        km1 = 0:Mm1;
        binom2 = arrayfun(@(kk) nchoosek(Mm1,kk), km1);
        db = M * sum(dc .* binom2 .* s.^km1 .* (1-s).^(Mm1-km1));
    else
        db = 0;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% GAIT SIMULATION (event-based, one full swing phase)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [x_start_out, x_end, total_torque_sq] = simulate_hzd_gait(coeffs, x_start, p)
    nq = p.nq;
    x_start_out = x_start;  % passed through unchanged, kept for signature symmetry with hzd_cost/hzd_constraints

    theta_minus = theta_of_q(x_start(1:nq));
    theta_plus  = p.theta_plus;

    xi0 = [x_start; 0];  % append running torque-cost integral

    opts = odeset('Events', @(t,xi) impact_event_wrapper(t,xi,p), ...
                  'RelTol',1e-8, 'AbsTol',1e-9);

    [~, XI] = ode45(@(t,xi) hzd_closed_loop_ode(t,xi,coeffs,theta_minus,theta_plus,p), ...
                     [0 p.T_max], xi0, opts);

    xi_end          = XI(end,:).';
    x_end           = xi_end(1:2*nq);
    total_torque_sq = xi_end(end);
end

function [value,isterminal,direction] = impact_event_wrapper(t,xi,p)
    nq = p.nq;
    x = xi(1:2*nq);
    [value,isterminal,direction] = rabbit_impact_event(t,x,p);  % reuse your existing swing-foot-height event fn
end

function dxi = hzd_closed_loop_ode(t,xi,coeffs,theta_minus,theta_plus,p)
    nq = p.nq;
    x  = xi(1:2*nq);
    q  = x(1:nq);
    dq = x(nq+1:end);

    theta      = theta_of_q(q);
    s          = (theta - theta_minus) / (theta_plus - theta_minus);
    s          = min(max(s,0),1);           % clamp — guards overshoot right at impact
    ds_dtheta  = 1/(theta_plus - theta_minus);
    dtheta_dt  = dtheta_dq_of(q) * dq;       % scalar

    act_idx = 4:7;                            % q1,q2,q3,q4
    y  = q(act_idx);
    dy = dq(act_idx);

    yd  = zeros(p.nu,1);
    dyd = zeros(p.nu,1);
    for i = 1:p.nu
        [b, db] = bezier_eval(coeffs(i,:), s);
        yd(i)  = b;
        dyd(i) = db * ds_dtheta * dtheta_dt;   % chain rule
    end

    e  = y  - yd;
    de = dy - dyd;
    tau = -p.Kp.*e - p.Kd.*de;

    % TODO: match this to your actual constrained-dynamics function signature
    ddq = rabbit_constrained_dynamics(q, dq, tau, p);

    dxi = [dq; ddq; sum(tau.^2)];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% HZD COST FUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function J = hzd_cost(z, p)
    [coeffs, x_start] = unpack_z(z, p);
    [~, ~, total_torque_sq] = simulate_hzd_gait(coeffs, x_start, p);
    J = total_torque_sq;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% HZD CONSTRAINTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [c, ceq] = hzd_constraints(z, p)
    [coeffs, x_start] = unpack_z(z, p);
    nq = p.nq;

    [~, x_end, ~] = simulate_hzd_gait(coeffs, x_start, p);

    % --- Periodicity (Poincare map): post-impact relabeled state == x_start
    x_next = rabbit_reset_map(rabbit_impact_map(x_end, p), p);
    ceq_periodicity = x_next - x_start;

    % --- Ground-contact validity of x_start (needed now that x_start is free)
    q_start  = x_start(1:nq);
    dq_start = x_start(nq+1:end);

    foot_pos = foot_positions(q_start, p);        % TODO: match actual name/signature
    ceq_foot_pos = foot_pos(1:2);                  % stance foot must sit at [0;0]

    J_stance = stance_foot_jacobian(q_start, p);   % TODO: match actual name/signature
    ceq_foot_vel = J_stance * dq_start;            % stance foot velocity must be zero

    ceq = [ceq_periodicity; ceq_foot_pos; ceq_foot_vel];

    % --- Stability (optional, left for later)
    c = [];
end