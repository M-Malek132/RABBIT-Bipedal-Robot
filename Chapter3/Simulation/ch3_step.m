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

sol = ode45(@(t,x) ch3_ode_rhs(t, x, alpha, p), [0 T_cap], x0(:), opts);

t = sol.x;
X = sol.y;

out = struct();
out.t     = t;
out.x     = X;
out.x_end = X(:, end);
out.T     = t(end);
out.ok    = ~isempty(sol.ie);
out.sol   = sol;

% stance-foot position at the start, swing-foot position at strike: their
% horizontal difference is the step length.
foot_st_0   = P_st(x0(1:p.nq));
foot_sw_end = P_sw(out.x_end(1:p.nq));
out.L_step  = foot_sw_end(1) - foot_st_0(1);

[out.x_next, out.impulse] = ch3_impact(out.x_end, p);

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
