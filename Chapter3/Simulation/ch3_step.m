function out = ch3_step(x0, alpha, p)
%CH3_STEP  Integrate ONE step of the hybrid system, up to the guard.
%
%   out = ch3_step(x0, alpha, p)
%
% Integrates xdot = f + g u(x) from x0 until the swing foot strikes the
% ground (the switching surface S), then applies Delta.
%
% THE t = 0 GUARD PROBLEM.  At the start of a step the "swing" foot is the
% foot that just landed, so its height is EXACTLY zero and the guard is
% already satisfied at t = 0. Terminating there would give a zero-length step.
% The event function therefore reports a positive constant until
% p.guard_min_time, by which point the foot has lifted; from then on it
% reports the true height. Because the held value is positive and the true
% height is also positive at the switch-over, no spurious crossing is created.
%
% Inputs
%   x0    : 14x1 start-of-step state
%   alpha : ny x n_ctrl coefficients
%   p     : parameter struct
%
% Outputs
%   out : struct with
%           .t       1 x nt time grid
%           .x       14 x nt state trajectory
%           .x_end   14x1 pre-impact state
%           .x_next  14x1 post-Delta state (start of the next step)
%           .T       step duration [s]
%           .impulse 2x1 impact impulse [Ns]
%           .L_step  step length [m] (advance of the stance foot)
%           .ok      false if the guard never fired
%           .sol     the raw ode45 solution struct, so callers can resample
%                    with deval at the solver's own accuracy rather than
%                    re-interpolating the output grid
%
% See also CH3_ODE_RHS, CH3_IMPACT, CH3_SIMULATE.

opts = odeset('RelTol',   p.ode_reltol, ...
              'AbsTol',   p.ode_abstol, ...
              'MaxStep',  p.ode_maxstep, ...
              'Events',   @(t,x) guard_event(t, x, p));

T_cap = p.T_max * 2;      % hard cap so a non-striking gait cannot run forever

if isfield(p, 'control_dt') && p.control_dt > 0
    [t, X, fired] = integrate_zoh(x0, alpha, p, T_cap, opts);
    sol = [];             % see note in integrate_zoh on dense output
else
    sol   = ode45(@(t,x) ch3_ode_rhs(t, x, alpha, p), [0 T_cap], x0(:), opts);
    t     = sol.x;
    X     = sol.y;
    fired = ~isempty(sol.ie);
end

out = struct();
out.t     = t;
out.x     = X;
out.x_end = X(:, end);
out.T     = t(end);
out.ok    = fired;
out.sol   = sol;

% stance-foot position at the start, swing-foot position at strike: their
% horizontal difference is the step length.
foot_st_0   = P_st(x0(1:p.nq));
foot_sw_end = P_sw(out.x_end(1:p.nq));
out.L_step  = foot_sw_end(1) - foot_st_0(1);

[out.x_next, out.impulse] = ch3_impact(out.x_end, p);

end

% ---------------------------------------------------------------------------
function [t_all, X_all, fired] = integrate_zoh(x0, alpha, p, T_cap, opts)
%INTEGRATE_ZOH  Sampled-data integration: one control solve per period.
%
% The controller is evaluated ONCE at the start of each period of length
% p.control_dt and then held, so inside a period the integrand is the smooth
% f(x) + g(x) u with u a fixed vector. ode45 therefore sees no active-set
% kinks and takes full steps; the guard is still checked continuously, since
% opts carries the same Events function.
%
% DENSE OUTPUT.  Each period produces its own ode45 solution, so there is no
% single struct deval can interpolate. ch3_step returns out.sol = [] here.
% The only consumer of out.sol is ch3_col_seed, which seeds from a PD rollout
% (p.control_dt = 0) and so never takes this branch -- but a caller that needs
% dense output under sampling should resample out.t / out.x instead.

dt    = p.control_dt;
t_all = 0;
X_all = x0(:);
tk    = 0;
xk    = x0(:);
fired = false;

while tk < T_cap - eps(T_cap)
    u  = ch3_control(xk, alpha, p);          % one QP solve for this period
    te = min(tk + dt, T_cap);

    s = ode45(@(~, x) zoh_rhs(x, u, p), [tk te], xk, opts);

    % ode45 returns at least the endpoints; drop the duplicated start point.
    t_all = [t_all, s.x(2:end)];             %#ok<AGROW>
    X_all = [X_all, s.y(:, 2:end)];          %#ok<AGROW>

    tk = s.x(end);
    xk = s.y(:, end);

    if ~isempty(s.ie)
        fired = true;
        return;
    end
end
end

% ---------------------------------------------------------------------------
function xdot = zoh_rhs(x, u, p)
%ZOH_RHS  Plant only: the control is already decided for this period.
[f, g] = ch3_control_affine(x, p);
xdot   = f + g*u;
end

% ---------------------------------------------------------------------------
function [value, isterminal, direction] = guard_event(t, x, p)
if t < p.guard_min_time
    value = 1;                       % hold positive; foot has not lifted yet
else
    value = ch3_guard(x, p);         % true swing-foot height
end
isterminal = 1;
direction  = -1;                     % fire only on a DOWNWARD crossing
end
