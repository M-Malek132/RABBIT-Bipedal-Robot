function out = ch6_step(x0, alpha, p)
%CH6_STEP  One walking step under the Chapter-6 controller, with diagnostics.
%
%   out = ch6_step(x0, alpha, p)
%
% Integrates the closed loop from x0 to swing-foot strike, applies Delta, and
% records what the chapter's figures need: the barrier values along the step
% (Figs. 6.13-6.15, 6.25b), the contact forces and torques (Fig. 6.26), and the
% step length at impact (Figs. 6.5, 6.10, 6.25a).
%
% ============================================ EVERY CONTROLLER IS SAMPLED, ALWAYS
% p.control_dt must be > 0 here, and this is a deliberate restriction rather
% than an omission.
%
% The performance reason is Chapter 3's: a QP control law is only piecewise
% smooth in x, the active set changes as rows engage, and an adaptive explicit
% solver treats every kink as a failed error test. ch3_params records the
% measurement -- 51908 RHS evaluations to advance 0.5% of a step. Chapter 6
% makes it worse, because barrier rows switch more often than a torque box does.
%
% The FAIRNESS reason matters more, and it is why the restriction applies to the
% PD and plain-CLF baselines too. Table 6.1 compares three controllers. A
% continuously evaluated controller enjoys an advantage no digital
% implementation of it has, so letting the smooth baselines run continuously
% while the QP is held at 1 kHz would make the comparison meaningless in the
% baselines' favour. ch3_compare_controllers samples all three for the same
% reason.
%
% ------------------------------------------------------------ the t = 0 guard
% At the start of a step the swing foot is the foot that just landed, so its
% height is exactly zero and the guard is already satisfied. The event function
% holds a positive constant until p.guard_min_time, by which point the foot has
% lifted. Same construction as ch3_step; see there for why no spurious crossing
% is created.
%
% Inputs
%   x0    : 14x1 start-of-step state
%   alpha : ny x n_ctrl coefficients for this step
%   p     : parameter struct, with p.stone already resolved for this step
%
% Outputs
%   out : struct
%     .t .x                trajectory (control samples plus solver points)
%     .x_end .x_next       pre- and post-impact states
%     .T                   step duration [s]
%     .L_step              stance-foot advance [m]
%     .l_s                 STEP LENGTH AT IMPACT, l_f at strike            (6.13)
%     .h_f_end             swing-foot height at termination [m]
%     .in_window           l_min(T) <= l_s <= l_max(T)
%     .stone_end           the foothold window at the impact instant
%     .impulse             2x1 impact impulse [Ns]
%     .ok                  the guard fired
%     .log                 per-control-sample diagnostics (see below)
%     .adm                 ch6_admissible at x0
%
%   out.log fields, each 1 x n_samples unless noted:
%     .t .h (nb x n) .h_cbf (nb x n) .margin (nb x n) .cbf_active (nb x n)
%     .u (nu x n) .Fx .Fz .mu_fric .delta .V .feasible .n_dead .l_f .h_f
%     .l_min .l_max
%
% See also CH6_CONTROL, CH6_SIMULATE, CH6_ADMISSIBLE, CH3_IMPACT.

if ~isfield(p, 'control_dt') || p.control_dt <= 0
    error('ch6_step:control_dt', ...
          ['Chapter 6 requires p.control_dt > 0 (sampled-data control). ' ...
           'See the header of this file: a QP law stalls an adaptive solver, ' ...
           'and sampling the baselines at the same rate is what makes ' ...
           'Table 6.1 a comparison rather than a handicap.']);
end

x0 = x0(:);

out = struct();
out.adm = ch6_admissible(x0, alpha, p, 0);

opts = odeset('RelTol',  p.ode_reltol, ...
              'AbsTol',  p.ode_abstol, ...
              'MaxStep', p.ode_maxstep, ...
              'Events',  @(t,x) guard_event(t, x, p));

T_cap = p.T_max * 2;
dt    = p.control_dt;

t_all = 0;
X_all = x0;
tk    = 0;
xk    = x0;
fired = false;

L = init_log(p);

while tk < T_cap - eps(T_cap)

    [u, ci] = ch6_control(tk, xk, alpha, p);
    L = push_log(L, tk, u, ci, p);

    te = min(tk + dt, T_cap);
    s  = ode45(@(~, xx) zoh_rhs(xx, u, p), [tk te], xk, opts);

    t_all = [t_all, s.x(2:end)];        %#ok<AGROW>
    X_all = [X_all, s.y(:, 2:end)];     %#ok<AGROW>

    tk = s.x(end);
    xk = s.y(:, end);

    if ~isempty(s.ie)
        fired = true;
        break;
    end
    if ~all(isfinite(xk))
        break;
    end
end

out.t     = t_all;
out.x     = X_all;
out.x_end = X_all(:, end);
out.T     = t_all(end);
out.ok    = fired;
out.log   = finish_log(L);

%% ------------------------------------------------------- what landed where
q_end = out.x_end(1:p.nq);
r_end = P_sw(q_end) - P_st(q_end);

out.l_s     = r_end(1);                  % step length at impact          (6.13)
out.h_f_end = r_end(2);

foot_st_0  = P_st(x0(1:p.nq));
foot_sw_e  = P_sw(q_end);
out.L_step = foot_sw_e(1) - foot_st_0(1);

% The window is evaluated AT THE IMPACT INSTANT, which for a moving stone is
% not the window the step started with -- that is the whole content of Section
% 6.2, so reading it at t = 0 would report success against a stone that is no
% longer there.
% TOLERANCE: 1 micron. A barrier that is doing its job lands the foot ON the
% edge of a short stone, so this test is decided at the boundary and a 1e-9
% tolerance would be asking floating point to settle it. 1e-6 m is far below
% anything physical and far above the integrator's own noise.
lv = ch6_stone_level(p.stone, out.T);
out.stone_end = lv;
out.in_window = (out.l_s >= lv.l_min - 1e-6) && (out.l_s <= lv.l_max + 1e-6);

[out.x_next, out.impulse] = ch3_impact(out.x_end, p);

end

% ---------------------------------------------------------------------------
function xdot = zoh_rhs(x, u, p)
[f, g] = ch3_control_affine(x, p);
xdot   = f + g*u;
end

% ---------------------------------------------------------------------------
function [value, isterminal, direction] = guard_event(t, x, p)
if t < p.guard_min_time
    value = 1;
else
    value = ch3_guard(x, p);
end
isterminal = 1;
direction  = -1;
end

% ---------------------------------------------------------------------------
function L = init_log(p)
%INIT_LOG  Preallocate generously; trimmed in finish_log.
%
% The sample count is known within one: T_cap/dt. Growing the arrays instead
% would reallocate a few thousand times per step inside the integration loop.
n = ceil(2*p.T_max / p.control_dt) + 2;

L = struct();
L.n  = 0;
L.t  = nan(1, n);
L.u  = nan(p.nu, n);
L.Fx = nan(1, n);   L.Fz = nan(1, n);   L.mu_fric = nan(1, n);
L.delta = nan(1, n);  L.V = nan(1, n);
L.feasible = true(1, n);   L.n_dead = zeros(1, n);
L.l_f = nan(1, n);  L.h_f = nan(1, n);
L.l_min = nan(1, n); L.l_max = nan(1, n);
L.h = []; L.h_cbf = []; L.margin = []; L.cbf_active = [];
L.nmax = n;
end

% ---------------------------------------------------------------------------
function L = push_log(L, t, u, ci, p)
i = L.n + 1;
if i > L.nmax, return; end       % T_cap reached; the loop is about to stop
L.n = i;

L.t(i)    = t;
L.u(:, i) = u;

% GRF from the SAME KKT solve the controller used -- lambda = lam_drift +
% lam_in u is exact, so this is the force the simulated robot actually applied,
% not a reconstruction of it.
lam = ci.aux.lam_drift + ci.aux.lam_in * u;
L.Fx(i) = lam(1);
L.Fz(i) = lam(2);
L.mu_fric(i) = abs(lam(1)) / max(abs(lam(2)), eps);

L.delta(i)    = ci.delta;
L.V(i)        = ci.V;
L.feasible(i) = ci.qp.feasible;
L.n_dead(i)   = ci.qp.n_dead;

nb = numel(ci.B);
if isempty(L.h)
    L.h          = nan(nb, L.nmax);
    L.h_cbf      = nan(nb, L.nmax);
    L.margin     = nan(nb, L.nmax);
    L.cbf_active = false(nb, L.nmax);
end
for j = 1:nb
    L.h(j, i)          = ci.B(j).g;
    L.h_cbf(j, i)      = p.cbf.gamma_b * ci.B(j).g + ci.B(j).gdot;
    L.margin(j, i)     = ci.qp.margin(j);
    L.cbf_active(j, i) = ci.qp.cbf_active(j);
end

if strcmpi(p.cbf.problem, 'stones')
    L.l_f(i)   = ci.geom.l_f;
    L.h_f(i)   = ci.geom.h_f;
    L.l_min(i) = ci.geom.l_min;
    L.l_max(i) = ci.geom.l_max;
elseif strcmpi(p.cbf.problem, 'obstacle')
    L.l_f(i) = ci.geom.l_H;
    L.h_f(i) = ci.geom.h_H;
end
end

% ---------------------------------------------------------------------------
function L = finish_log(L)
n = L.n;
f = {'t','Fx','Fz','mu_fric','delta','V','feasible','n_dead', ...
     'l_f','h_f','l_min','l_max'};
for i = 1:numel(f), L.(f{i}) = L.(f{i})(1:n); end
L.u = L.u(:, 1:n);
for f2 = {'h','h_cbf','margin','cbf_active'}
    if ~isempty(L.(f2{1})), L.(f2{1}) = L.(f2{1})(:, 1:n); end
end
L = rmfield(L, {'n','nmax'});
end
