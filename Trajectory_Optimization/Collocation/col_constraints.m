function [c, ceq] = col_constraints(z, p)
%COL_CONSTRAINTS  Hermite-Simpson HZD direct-collocation constraints.
%
%   [c, ceq] = col_constraints(z, p)
%
% EQUALITIES (ceq = 0):
%   * Hermite-Simpson dynamics defects between consecutive nodes (no ODE solve).
%   * Stance-foot pinned: P_st(q_k) = P_st(q_1) at every node (kills drift), the
%     foot on the ground at node 1 (P_st(q_1)_z = 0), and zero foot velocity at
%     node 1 (J_st(q_1) v_1 = 0). Position-at-all-nodes + velocity-at-node-1 is
%     the index-reduced holonomic form.
%   * VIRTUAL CONSTRAINTS embedded: y(q_k) = q_k(4:7) - yd(s_k, coeffs) = 0 at
%     EVERY node, and ydot = 0 at node 1. This forces the trajectory onto the
%     hybrid zero-dynamics manifold, so `coeffs` (the B-spline / Bezier params)
%     parametrizes the gait -- exactly the thesis formulation.
%   * Terminal impact surface (swing-foot height 0), theta(q_N) = theta_plus.
%   * Periodicity through impact + relabel: reset(impact(x_N)) matches x_1 on
%     all coordinates except the horizontal hip position (which advances).
%   * Average walking speed L_step / T = v_des.
%
% INEQUALITIES (c <= 0):
%   * No ground penetration + mid-step swing-foot clearance (NIC3).
%   * Step-length floor.
%   * Table-3.1 physical limits, each gated: friction cone and min vertical GRF
%     on the explicit contact force Lam, impact impulse at touchdown. (The
%     torque box is enforced as simple bounds on U in col_problem_setup, not
%     here.)
%
% s_k is NOT clamped here (unlike the controller): the terminal theta_end
% equality pins s_N = 1 and monotone theta keeps s in [0,1], so leaving it
% smooth gives fmincon clean derivatives.

    [X, U, Lam, coeffs, T] = col_unpack(z, p);
    N  = p.N;   nq = p.nq;   nu = p.nu;   nx = 2*nq;
    h  = T/(N-1);
    theta_plus = p.theta_plus;

    q1 = X(1:nq,1);   v1 = X(nq+1:nx,1);
    theta_minus = theta_of_q(q1);
    denom = theta_plus - theta_minus;
    kappa = 1/denom;
    c_row = dtheta_dq_of(q1);              % constant 1 x nq
    Hsel  = [zeros(nu,3), eye(nu)];        % selects actuated joints q(4:7)

    % ---- f at every node -------------------------------------------------
    F = zeros(nx, N);
    for k = 1:N
        F(:,k) = col_dynamics(X(:,k), U(:,k), Lam(:,k), p);
    end

    % ---- Hermite-Simpson defects ----------------------------------------
    Def = zeros(nx, N-1);
    for k = 1:N-1
        xmid = 0.5*(X(:,k) + X(:,k+1)) + (h/8)*(F(:,k) - F(:,k+1));
        umid = 0.5*(U(:,k) + U(:,k+1));
        lmid = 0.5*(Lam(:,k) + Lam(:,k+1));
        fmid = col_dynamics(xmid, umid, lmid, p);
        Def(:,k) = X(:,k+1) - X(:,k) - (h/6)*(F(:,k) + 4*fmid + F(:,k+1));
    end

    % ---- foot no-drift + virtual constraints at every node ---------------
    foot1 = P_st(q1);
    ceq_footdrift = zeros(2, N-1);
    ceq_vc        = zeros(nu, N);
    for k = 1:N
        qk  = X(1:nq,k);
        sk  = (theta_of_q(qk) - theta_minus) / denom;
        yd  = zeros(nu,1);
        for i = 1:nu, yd(i) = bspline_eval(coeffs(i,:), sk, p); end
        ceq_vc(:,k) = qk(4:7) - yd;
        if k >= 2
            ceq_footdrift(:,k-1) = P_st(qk) - foot1;
        end
    end

    % ---- node-1 anchors (velocity level) --------------------------------
    ceq_footZ1   = foot1(2);               % stance foot on ground
    ceq_footvel1 = J_st(q1) * v1;          % 2
    dyd_ds1 = zeros(nu,1);
    for i = 1:nu, [~, db] = bspline_eval(coeffs(i,:), 0, p); dyd_ds1(i) = db; end
    Jy1 = Hsel - kappa*(dyd_ds1 * c_row);
    ceq_ydot1 = Jy1 * v1;                  % nu

    % ---- terminal + periodicity + speed ---------------------------------
    qN = X(1:nq,N);
    swN = P_sw(qN);   stN = P_st(qN);
    ceq_impactsurf = swN(2);               % swing foot strikes ground
    ceq_thetaend   = theta_of_q(qN) - theta_plus;

    xplus = rabbit_reset_map(rabbit_impact_map(X(:,N)));
    ceq_period = xplus(2:14) - X(2:14,1);  % 13 (hip x, index 1, advances)

    L_step   = swN(1) - stN(1);
    ceq_speed = L_step/T - p.v_des;

    ceq = [ Def(:); ceq_footdrift(:); ceq_vc(:);
            ceq_footZ1; ceq_footvel1; ceq_ydot1;
            ceq_impactsurf; ceq_thetaend; ceq_period; ceq_speed ];

    % ======================= INEQUALITIES ================================
    % The NUMBER of inequalities must NOT depend on the state (fmincon requires
    % a fixed count). So the mid-step swing-foot clearance uses a fixed
    % NODE-INDEX band [k_lo,k_hi] (phase is monotone in the node index) rather
    % than a state-dependent phase test.
    c = [];
    k_lo = round(0.2*(N-1)) + 1;
    k_hi = round(0.8*(N-1)) + 1;
    for k = 2:N-1
        sw = P_sw(X(1:nq,k));
        c(end+1,1) = -sw(2);                              %#ok<AGROW> height >= 0
        if k >= k_lo && k <= k_hi
            c(end+1,1) = p.swing_clearance_min - sw(2);   %#ok<AGROW>
        end
    end

    % Step-length floor
    c(end+1,1) = p.step_length_min - L_step;

    % Table 3.1 physical limits (gated). NOTE the sign convention of Lam(2)
    % (vertical GRF) is whatever J_st' produces; check col_verify_seed output
    % and flip if the seed shows compression as negative before enforcing.
    gated = @(f) isfield(p,f) && p.(f);
    if gated('enforce_friction')
        for k = 1:N
            c(end+1,1) =  Lam(1,k) - p.mu_max*Lam(2,k);   %#ok<AGROW>
            c(end+1,1) = -Lam(1,k) - p.mu_max*Lam(2,k);   %#ok<AGROW>
        end
    end
    if gated('enforce_grf')
        for k = 1:N
            c(end+1,1) = p.grf_min - Lam(2,k);            %#ok<AGROW>
        end
    end
    if gated('enforce_impulse')
        [~, imp] = rabbit_impact_map(X(:,N));
        c(end+1,1) = norm(imp) - p.impulse_max;           %#ok<AGROW>
    end
end
