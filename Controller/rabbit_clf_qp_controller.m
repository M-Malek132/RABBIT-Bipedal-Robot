function [u, info] = rabbit_clf_qp_controller(q, dq, coeffs, theta_minus, theta_plus, p)
%RABBIT_CLF_QP_CONTROLLER  Control-Lyapunov-function quadratic-program controller.
%
%   [u, info] = rabbit_clf_qp_controller(q, dq, coeffs, theta_minus, theta_plus, p)
%
% Implements the input-output-linearizing CLF-QP controller of Chapter 3
% (Hybrid Zero Dynamics + control-Lyapunov-function QP). It is a drop-in
% alternative to the fixed-gain PD law in hzd_control_and_dynamics.m: same
% virtual-constraint outputs, same constrained single-support dynamics, but the
% torque is chosen by a quadratic program so it can respect actuator and
% friction limits that a fixed PD law cannot.
%
% ---------------------------------------------------------------------------
% OUTPUTS AND RELATIVE DEGREE
% ---------------------------------------------------------------------------
% The controlled outputs are the HZD virtual constraints (thesis Eq. 3.2),
%
%       y(q) = q(4:7) - yd(s),
%
% with q(4:7) the four actuated leg joints and yd(s) the B-spline desired
% evolution in the phase variable s = (theta - theta_minus)/(theta_plus-...).
% h0(q)=q(4:7) is linear in q, so y has relative degree 2 and
%
%       ydot = Jy * dq,           Jy = H - kappa*(dyd/ds)*c,
%       yddot = L_f^2 y + (L_g L_f y) u,
%
% where H = [0 I] selects the actuated joints, c = dtheta/dq (constant row,
% from dtheta_dq_of), kappa = ds/dtheta = 1/(theta_plus-theta_minus), and,
% because both H and c are constant,
%
%       L_f^2 y  = Jy*ddq0 - kappa^2 * (d^2 yd/ds^2) * (c*dq)^2,
%       L_g L_f y = Jy * (d ddq / du)            (the decoupling matrix, m x m).
%
% ddq0 and d(ddq)/du come from the acceleration-level constrained dynamics
% (stance foot pinned): the KKT solve makes ddq affine in u, ddq = ddq0 + Mu*u,
% which mirrors rabbit_constrained_dynamics.m with the Baumgarte term off.
%
% ---------------------------------------------------------------------------
% RES-CLF (thesis Eqs. 3.12-3.18)
% ---------------------------------------------------------------------------
% In the transverse coordinates eta = [y; ydot], the I/O-linearized system is
% etadot = F*eta + G*mu with mu := yddot. F,G are in controllable canonical
% form, so the CARE F'P + PF - PGG'P + I = 0 (Q=R=I) has the closed-form
% block solution P = [sqrt(3) I, I; I, sqrt(3) I]. The rapidly-exponentially-
% stabilizing CLF uses the epsilon-scaling I_eps = diag((1/eps)I, I):
%
%       V(eta) = eta' * P_eps * eta,     P_eps = I_eps * P * I_eps,
%
% and the QP enforces Vdot <= -(gamma/eps) V, gamma = lambda_min(Q)/lambda_max(P).
% The eps factor sets the convergence rate (needed to dominate the impact-map
% expansion; see thesis Sec. 3.3).
%
% ---------------------------------------------------------------------------
% QUADRATIC PROGRAM
% ---------------------------------------------------------------------------
%   decision z = [u (m x1); delta (relaxation, scalar)]
%   min   1/2 u'Ru + p_relax*delta^2
%   s.t.  (LgV*Adec) u - delta <= -(gamma/eps)V - LfV - LgV*L_f^2y   (CLF)
%         -u_max <= u <= u_max                                        (torque)
%         [optional friction cone |Fx| <= mu*Fz on the contact force]
% The relaxation delta keeps the QP feasible when the torque box would
% otherwise make the CLF inequality infeasible (standard Ames-CLF practice).
%
% Tunables (all read from p with the defaults shown):
%   p.clf_eps    (0.10)  CLF convergence rate factor (smaller = faster)
%   p.clf_R      (1)     weight on ||u||^2 in the QP
%   p.clf_relax  (1e3)   penalty on the CLF-constraint relaxation delta
%   p.u_max      (120)   |torque| box bound [Nm] (shared with the optimizer)
%   p.clf_enforce_torque   (true)  apply the torque box in the QP
%   p.clf_enforce_friction (false) add the friction-cone inequality
%   p.mu_max     (0.4)   friction coefficient for the cone
%   p.clf_Kp,p.clf_Kd    PD gains for the feedback-linearizing fallback
%
% info fields: ddq, tau(=u), lambda, V, Vdot, delta, s, y, dy, Adec, Lf2y,
%              exitflag (quadprog flag, or -99 if the fallback law was used).

    nq = p.nq;
    m  = p.nu;
    act_idx = 4:nq;                      % actuated joints q(4:7)

    % ---- defaults -------------------------------------------------------
    eps_c   = getdef(p, 'clf_eps', 0.10);
    Rwt     = getdef(p, 'clf_R',   1);
    p_relax = getdef(p, 'clf_relax', 1e3);
    u_max   = getdef(p, 'u_max',   120);
    enf_tau = getdef(p, 'clf_enforce_torque',   true);
    enf_fri = getdef(p, 'clf_enforce_friction', false);
    mu_max  = getdef(p, 'mu_max',  0.4);
    Kp_fb   = getdef(p, 'clf_Kp',  getdef(p, 'Kp', 400));
    Kd_fb   = getdef(p, 'clf_Kd',  getdef(p, 'Kd', 40));

    % ---- phase variable -------------------------------------------------
    theta   = theta_of_q(q);
    kappa   = 1 / (theta_plus - theta_minus);        % ds/dtheta
    s       = max(0, min(1, (theta - theta_minus) * kappa));
    c_row   = dtheta_dq_of(q);                        % 1 x nq  (constant)
    thd     = c_row * dq;                             % theta_dot = c*dq

    % ---- desired output and its s-derivatives ---------------------------
    yd = zeros(m,1); dyd_ds = zeros(m,1); d2yd_ds2 = zeros(m,1);
    for i = 1:m
        [b, db, ddb] = bspline_eval(coeffs(i,:), s, p);
        yd(i)       = b;
        dyd_ds(i)   = db;
        d2yd_ds2(i) = ddb;
    end

    % ---- outputs y, ydot ------------------------------------------------
    Hsel = zeros(m, nq);
    for i = 1:m, Hsel(i, act_idx(i)) = 1; end        % selects q(4:7)
    Jy = Hsel - kappa * (dyd_ds * c_row);            % dy/dq   (m x nq)
    y  = q(act_idx) - yd;
    dy = Jy * dq;

    % ---- constrained dynamics, affine in u:  ddq = ddq0 + Mu*u ----------
    [ddq0, Mu, lam0, Mlam] = constrained_dyn_affine(q, dq);

    % ---- Lie derivatives (relative degree 2) ----------------------------
    Lf2y = Jy * ddq0 - kappa^2 * d2yd_ds2 * (thd^2);  % L_f^2 y   (m x 1)
    Adec = Jy * Mu;                                   % L_g L_f y (m x m)

    % ---- RES-CLF ---------------------------------------------------------
    Im  = eye(m);
    a   = sqrt(3);                                    % CARE P block value
    P   = [a*Im, Im; Im, a*Im];                       % P (2m x 2m)
    Ie  = blkdiag((1/eps_c)*Im, Im);
    Pe  = Ie * P * Ie;                                % P_eps
    Fcl = [zeros(m), Im; zeros(m), zeros(m)];
    Gcl = [zeros(m); Im];
    eta = [y; dy];

    Vval  = eta' * Pe * eta;
    LfV   = eta' * (Fcl'*Pe + Pe*Fcl) * eta;          % scalar
    LgV   = 2 * eta' * Pe * Gcl;                       % 1 x m
    gamma = 1 / (a + 1);                               % lam_min(Q)/lam_max(P)

    % ---- assemble & solve the QP ----------------------------------------
    % CLF (soft):  (LgV*Adec) u - delta <= b_clf
    Aclf = [LgV * Adec, -1];
    bclf = -(gamma/eps_c)*Vval - LfV - LgV*Lf2y;

    Aineq = Aclf;
    bineq = bclf;

    if enf_fri
        % Contact force lambda = [Fx; Fz] = lam0 + Mlam*u, affine in u.
        % Cone |Fx| <= mu*Fz  ->  Fx - mu Fz <= 0 and -Fx - mu Fz <= 0.
        r1 = Mlam(1,:) - mu_max*Mlam(2,:);
        r2 = -Mlam(1,:) - mu_max*Mlam(2,:);
        Aineq = [Aineq; r1, 0; r2, 0];
        bineq = [bineq; -(lam0(1) - mu_max*lam0(2)); -(-lam0(1) - mu_max*lam0(2))];
    end

    Hqp = blkdiag(Rwt*Im, 2*p_relax);                 % 1/2 z'Hz -> p_relax*delta^2
    fqp = zeros(m+1, 1);

    if enf_tau
        lb = [-u_max*ones(m,1); 0];
        ub = [ u_max*ones(m,1); inf];
    else
        lb = [-inf(m,1); 0];
        ub = [ inf(m,1); inf];
    end

    persistent qpopts
    if isempty(qpopts)
        qpopts = optimoptions('quadprog', 'Display', 'off');
    end

    % try/catch so the controller never throws inside the optimizer's step
    % simulation (fmincon probes odd states); on any failure fall through to
    % the feedback-linearizing PD law below. warning() is silenced so a soft
    % quadprog "infeasible" message does not spam thousands of ODE evaluations.
    ws = warning('off', 'all');
    try
        [z, ~, exitflag] = quadprog(Hqp, fqp, Aineq, bineq, [], [], lb, ub, [], qpopts);
    catch
        z = []; exitflag = -10;
    end
    warning(ws);

    if isempty(z) || exitflag <= 0
        % Fallback: feedback-linearize + PD virtual control (thesis Eq. 3.16),
        % then saturate. Keeps the sim alive if the QP fails a step.
        mu_pd = -Kp_fb*y - Kd_fb*dy;
        u = Adec \ (mu_pd - Lf2y);
        u = max(-u_max, min(u_max, u));
        delta = NaN;
        exitflag = -99;
    else
        u     = z(1:m);
        delta = z(end);
    end

    % ---- outputs for the caller -----------------------------------------
    ddq    = ddq0 + Mu   * u;
    lambda = lam0 + Mlam * u;
    if nargout > 1
        info.ddq      = ddq;
        info.tau      = u;
        info.lambda   = lambda;
        info.V        = Vval;
        info.Vdot     = LfV + LgV*(Lf2y + Adec*u);
        info.delta    = delta;
        info.s        = s;
        info.y        = y;
        info.dy       = dy;
        info.Adec     = Adec;
        info.Lf2y     = Lf2y;
        info.exitflag = exitflag;
    end
end

% =========================================================================
function [ddq0, Mu, lam0, Mlam] = constrained_dyn_affine(q, dq)
% Acceleration-level single-support dynamics (stance foot pinned), written as
% an affine map in u. Mirrors rabbit_constrained_dynamics.m (Baumgarte OFF):
%   [M -J'; J 0] [ddq; -lambda] = [B u - C - G; -Jdot dq].
% Solving once against b0 (u=0) and against the input columns Bu gives
%   ddq    = ddq0 + Mu  * u,
%   lambda = lam0 + Mlam* u.
    Mm  = M(q);
    Cv  = V([q; dq]);
    Gv  = G(q);
    B   = input_matrix();
    J   = J_st(q);
    Jdq = Jdotdq_st(q, dq);
    nc  = size(J, 1);
    nq  = numel(q);

    A    = [Mm, -J'; J, zeros(nc)];
    Ainv = A \ eye(size(A));
    b0   = [-Cv - Gv; -Jdq];
    Bu   = [B; zeros(nc, size(B, 2))];

    sol0 = Ainv * b0;
    solU = Ainv * Bu;

    ddq0 = sol0(1:nq);        Mu   = solU(1:nq, :);
    lam0 = sol0(nq+1:end);    Mlam = solU(nq+1:end, :);
end

% =========================================================================
function v = getdef(s, f, d)
% Field f of struct s, or default d when the field is absent/empty.
    if isfield(s, f) && ~isempty(s.(f)), v = s.(f); else, v = d; end
end
