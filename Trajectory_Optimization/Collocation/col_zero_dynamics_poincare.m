function [rho_zd, thdot_fix, R_at_fix, info] = col_zero_dynamics_poincare(col_file)
%COL_ZERO_DYNAMICS_POINCARE  Intrinsic (zero-dynamics) orbital stability.
%
%   [rho_zd, thdot_fix, R_at_fix, info] = col_zero_dynamics_poincare(col_file)
%
% Computes the eigenvalue of the HYBRID-ZERO-DYNAMICS restricted Poincare map
% for a collocation gait -- the orbital stability the gait would have under a
% controller that renders the virtual-constraint manifold invariant (the CLF-QP
% ideal). Unlike poincare_stability (full-order, closed-loop, under a specific
% controller, ~26 shooting sims), this integrates only the 1-DOF underactuated
% ("zero") dynamics, so it is cheap and CONTROLLER-INDEPENDENT.
%
% Why it matters: the full-order rho under PD tells you the closed loop is
% unstable, but not WHY. If rho_zd < 1 the ORBIT is intrinsically stabilizable
% (PD is just too weak; CLF-QP should do better); if rho_zd > 1 the orbit
% itself must be reshaped -- i.e. stability must enter the gait OPTIMIZATION.
%
% Method. On the manifold the whole state is an explicit function of the phase
% theta and its rate thetadot (reconstruct_state below). The reduced Poincare
% section is theta = theta_minus (step start, post-impact); the scalar return
% map is
%     R : thetadot_minus  --(flow of the zero dynamics to theta_plus)-->
%         thetadot at impact  --(impact + relabel)-->  thetadot_minus(next).
% rho_zd = |dR/dthetadot| at the gait's own thetadot (central difference).
%
%   rho_zd     : |dR/dthetadot|  (< 1 => intrinsically orbitally stable)
%   thdot_fix  : the gait's start-of-step thetadot (nominal fixed point)
%   R_at_fix   : R(thdot_fix); |R_at_fix - thdot_fix| small => consistent orbit
%   info       : struct (theta_minus, theta_plus, coeffs, p)

    here = fileparts(mfilename('fullpath'));
    root = fileparts(fileparts(here));
    addpath(genpath(root));

    if nargin < 1 || isempty(col_file)
        d = dir(fullfile(root,'Results','col_result_*.mat'));
        [~,i] = max([d.datenum]);
        col_file = fullfile(d(i).folder, d(i).name);
    elseif ~isfile(col_file)
        col_file = fullfile(root,'Results', col_file);
    end
    S = load(col_file);
    p = S.p;
    [X, ~, ~, coeffs, ~] = col_unpack(S.z_opt, p);

    q0  = X(1:p.nq, 1);
    dq0 = X(p.nq+1:end, 1);
    theta_minus = theta_of_q(q0);
    theta_plus  = p.theta_plus;
    c_row = dtheta_dq_of(q0);                 % constant row
    thdot_fix = c_row * dq0;                   % nominal start-of-step thetadot

    ctx.p = p; ctx.coeffs = coeffs; ctx.c_row = c_row;
    ctx.theta_minus = theta_minus; ctx.theta_plus = theta_plus;
    ctx.x_foot = subsref(P_st(q0), struct('type','()','subs',{{1}}));  % stance foot x

    % central difference of the scalar return map
    h = 1e-4 * max(1, abs(thdot_fix));
    Rp = return_map(thdot_fix + h, ctx);
    Rm = return_map(thdot_fix - h, ctx);
    R_at_fix = return_map(thdot_fix, ctx);
    rho_zd = abs((Rp - Rm) / (2*h));

    info = struct('theta_minus',theta_minus,'theta_plus',theta_plus, ...
                  'coeffs',coeffs,'p',p);

    fprintf('%s\n', col_file);
    fprintf('  zero-dynamics restricted Poincare:\n');
    fprintf('    thetadot* (fixed point)   : %.4f rad/s\n', thdot_fix);
    fprintf('    R(thetadot*)              : %.4f rad/s  (consistency |.-.|=%.2e)\n', ...
            R_at_fix, abs(R_at_fix - thdot_fix));
    fprintf('    rho_zd = |dR/dthetadot|   : %.4f   %s\n', rho_zd, ...
            tern(rho_zd < 1, '[intrinsically STABLE]', '[intrinsically UNSTABLE]'));
end

% =========================================================================
function thd_next = return_map(thd_start, ctx)
% One step of the reduced return map: flow the zero dynamics from theta_minus
% to theta_plus starting at thetadot=thd_start, then apply impact + relabel.
    [q0, dq0] = reconstruct_state(ctx.theta_minus, thd_start, ctx);

    odef = @(t, z) [ z(2); theta_ddot(z(1), z(2), ctx) ];
    opts = odeset('RelTol',1e-8,'AbsTol',1e-10, ...
                  'Events', @(t,z) reach_theta_plus(t, z, ctx.theta_plus), ...
                  'MaxStep', 0.02);

    z0 = [ctx.theta_minus; thd_start];
    sol = ode45(odef, [0 5], z0, opts);
    if isempty(sol.ie)
        thd_next = NaN; return;      % never reached theta_plus (fell / stalled)
    end
    theta_end = sol.ye(1);   thdot_end = sol.ye(2);

    [qm, dqm] = reconstruct_state(theta_end, thdot_end, ctx);
    xplus = rabbit_reset_map(rabbit_impact_map([qm; dqm]));
    thd_next = ctx.c_row * xplus(8:14);
end

% =========================================================================
function thdd = theta_ddot(theta, thdot, ctx)
% Zero-dynamics acceleration: reconstruct the manifold state, apply the
% feedback-linearizing torque that holds the outputs (ydot,yddot = 0), and read
% off thetaddot = c_row * qddot. Reuses rabbit_constrained_dynamics (affine in
% u) exactly as the CLF-QP controller does.
    p = ctx.p; nu = p.nu; nq = p.nq;
    [q, dq] = reconstruct_state(theta, thdot, ctx);

    kappa = 1/(ctx.theta_plus - ctx.theta_minus);
    s = (theta - ctx.theta_minus) * kappa;

    dyd_ds = zeros(nu,1);  d2yd = zeros(nu,1);
    for i = 1:nu
        [~, db, ddb] = bspline_eval(ctx.coeffs(i,:), s, p);
        dyd_ds(i) = db;  d2yd(i) = ddb;
    end
    Hsel = [zeros(nu,3), eye(nu)];
    Jy   = Hsel - kappa*(dyd_ds * ctx.c_row);

    % ddq = ddq0 + Mu*u  (rabbit_constrained_dynamics is affine in u)
    ddq0 = rabbit_constrained_dynamics(q, dq, zeros(nu,1));
    Mu = zeros(nq, nu);
    for i = 1:nu
        e = zeros(nu,1); e(i) = 1;
        Mu(:,i) = rabbit_constrained_dynamics(q, dq, e) - ddq0;
    end
    Adec = Jy * Mu;
    Lf2y = Jy * ddq0 - kappa^2 * d2yd * (thdot^2);
    u_ff = -Adec \ Lf2y;                       % hold yddot = 0

    ddq = rabbit_constrained_dynamics(q, dq, u_ff);
    thdd = ctx.c_row * ddq;
end

% =========================================================================
function [q, dq] = reconstruct_state(theta, thdot, ctx)
% Explicit manifold reconstruction (q,dq) from (theta, thetadot):
%   actuated joints slaved to the spline, torso from theta's definition, base
%   (hip) from the pinned stance foot. Both legs length 1/2, so
%   hip - P_st = cos(q5/2) * [sin(theta); cos(theta)].
    p = ctx.p; nu = p.nu;
    kappa = 1/(ctx.theta_plus - ctx.theta_minus);
    s = (theta - ctx.theta_minus) * kappa;

    yd = zeros(nu,1);  dyd_ds = zeros(nu,1);
    for i = 1:nu
        [b, db] = bspline_eval(ctx.coeffs(i,:), s, p);
        yd(i) = b;  dyd_ds(i) = db;
    end

    q = zeros(7,1);
    q(4:7) = yd;                                % actuated joints
    q(3)   = theta - q(4) - 0.5*q(5);           % torso from theta = q3+q4+0.5 q5
    ee = cos(q(5)/2) * [sin(theta); cos(theta)];
    q(1:2) = [ctx.x_foot; 0] + ee;              % hip from pinned foot P_st=[x_foot;0]

    dq = zeros(7,1);
    dq(4:7) = dyd_ds * kappa * thdot;           % ydot = 0
    dq(3)   = thdot - dq(4) - 0.5*dq(5);        % thetadot = dq3+dq4+0.5 dq5
    % d/dt [ cos(q5/2)*(sin theta; cos theta) ]
    dee = -sin(q(5)/2)*(dq(5)/2) * [sin(theta); cos(theta)] + ...
           cos(q(5)/2) * [cos(theta); -sin(theta)] * thdot;
    dq(1:2) = dee;
end

function [value, isterminal, direction] = reach_theta_plus(~, z, theta_plus)
% Stop the zero-dynamics flow when theta reaches theta_plus (rising).
    value = z(1) - theta_plus;
    isterminal = 1;
    direction = 1;
end

function out = tern(c,a,b), if c, out=a; else, out=b; end, end
