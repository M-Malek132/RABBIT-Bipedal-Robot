function out = ch4_step(x0, xi0, alpha, p)
%CH4_STEP  Integrate ONE step of the hybrid system under model uncertainty.
%
%   out = ch4_step(x0, xi0, alpha, p)
%
% The Chapter-3 step integrator, extended in the three ways Chapter 4 needs:
%
%   1. the plant is the TRUE model (ch4_control_affine with p.uncertainty)
%      while the controller is built from the nominal one;
%   2. the controller may carry state, which is advanced across the step and
%      handed back so the next step continues it;
%   3. the reset map is ch4_impact, on the true model.
%
% xi0 = [] initializes the controller state from the start-of-step transverse
% state, so a caller that does not care about controller state can ignore it
% entirely.
%
% ---------------------------------------------------------- the two integrators
% p.control_dt = 0   ode45 on the augmented [x; xi]. Exact, and appropriate for
%                    the smooth laws, but see the warning below.
% p.control_dt > 0   sampled data: one control decision per period, held; the
%                    plant integrates under the held torque and the controller
%                    state advances by RK4 (ch4_l1_advance).
%
% CHAPTER 4 DEFAULTS TO SAMPLED DATA and that default should not be casually
% overridden. Two independent reasons:
%
%   * the constrained QPs are only piecewise smooth in x, so an adaptive
%     explicit solver stalls on the active-set kinks -- the Chapter-3 lesson,
%     measured in the p.control_dt note in ch3_params;
%   * with Gamma = 1e4 the adaptation is far stiffer than the robot, so a
%     continuous run makes ode45 resolve the ESTIMATOR at every step of the
%     PLANT. That is not just slow, it also silently reports a controller that
%     no digital implementation could run.
%
% Inputs
%   x0    : 14x1 start-of-step state
%   xi0   : controller state, or [] to initialize
%   alpha : ny x n_ctrl coefficients
%   p     : parameter struct
%
% Outputs
%   out : struct with
%           .t .x .x_end .x_next .T .impulse .L_step .ok .sol   as ch3_step
%           .xi       controller state on the SAMPLE grid (empty if stateless)
%           .t_xi     the sample times those correspond to
%           .xi_end   controller state at the end of the step (pre-reset)
%           .xi_next  controller state to start the next step (post-reset)
%
% See also CH3_STEP, CH4_ODE_RHS, CH4_IMPACT, CH4_L1_ADVANCE, CH4_L1_STATE.

stateful = ch4_is_stateful(p);

if stateful && isempty(xi0)
    [~, ~, o0] = ch3_outputs(x0, alpha, p);
    xi0 = ch4_l1_state('init', p, o0.eta);
elseif ~stateful
    xi0 = zeros(0,1);
end

opts = odeset('RelTol',  p.ode_reltol, ...
              'AbsTol',  p.ode_abstol, ...
              'MaxStep', p.ode_maxstep, ...
              'Events',  @(t,z) guard_event(t, z, p));

T_cap = p.T_max * 2;

if p.control_dt > 0
    [t, X, XI, t_xi, fired] = integrate_zoh(x0, xi0, alpha, p, T_cap, opts);
    sol = [];
else
    z0  = [x0(:); xi0(:)];
    sol = ode45(@(t,z) ch4_ode_rhs(t, z, alpha, p), [0 T_cap], z0, opts);
    t   = sol.x;
    X   = sol.y(1:p.nx, :);
    XI  = sol.y(p.nx+1:end, :);
    t_xi = t;
    fired = ~isempty(sol.ie);
end

out = struct();
out.t      = t;
out.x      = X;
out.x_end  = X(:, end);
out.T      = t(end);
out.ok     = fired;
out.sol    = sol;
out.xi     = XI;
out.t_xi   = t_xi;
out.xi_end = XI(:, end);

foot_st_0   = P_st(x0(1:p.nq));
foot_sw_end = P_sw(out.x_end(1:p.nq));
out.L_step  = foot_sw_end(1) - foot_st_0(1);

% --- reset map on the TRUE model -----------------------------------------
[out.x_next, out.impulse] = ch4_impact(out.x_end, p);

% --- carry the controller state across the impact ------------------------
% eta jumps here. p.l1.reset_predictor decides whether the predictor is told;
% the parameter estimates are carried regardless (see ch4_l1_state).
if stateful
    [~, ~, o_plus] = ch3_outputs(out.x_next, alpha, p);
    out.xi_next    = ch4_l1_state('reset', p, out.xi_end, o_plus.eta);
else
    out.xi_next = zeros(0,1);
end

end

% ---------------------------------------------------------------------------
function [t_all, X_all, XI_all, t_xi, fired] = integrate_zoh(x0, xi0, alpha, p, T_cap, opts)
%INTEGRATE_ZOH  One control decision per period; plant and controller advance.

dt = p.control_dt;

stateful = ch4_is_stateful(p);
clf      = ch3_res_clf(p);

t_all  = 0;
X_all  = x0(:);
XI_all = xi0(:);
t_xi   = 0;

tk    = 0;
xk    = x0(:);
xik   = xi0(:);
fired = false;

while tk < T_cap - eps(T_cap)

    % ---- bail out rather than spin --------------------------------------
    % A diverging controller drives the state to Inf/NaN, and ode45 then
    % returns without advancing time. Without this check the loop would make no
    % progress and never terminate -- a hang, reported as nothing at all,
    % instead of a failed run reported as a failed run. Chapter 4 deliberately
    % runs controllers past the point where they cope, so this path is reached
    % in normal use, not only when something is broken.
    if ~all(isfinite(xk)) || ~all(isfinite(xik))
        break;
    end

    % Catch the divergence BEFORE it reaches Inf. Once the state is merely
    % large the dynamics stiffen, and ode45 then spends thousands of internal
    % steps crawling through a single 1 ms period -- the run still fails, it
    % just takes minutes to do it. p.dq_max is the collocation's own box on
    % joint velocity, so ten times it is far outside any trajectory the gait
    % could produce and unambiguously means "this controller has lost the
    % robot". Reported as a failed step, which is the honest answer.
    if norm(xk(p.nq+1:end), inf) > 10 * p.dq_max
        break;
    end

    % ---- one control decision for this period ---------------------------
    [u, ~, ci] = ch4_control(xk, xik, alpha, p);

    if ~all(isfinite(u))
        break;
    end

    te = min(tk + dt, T_cap);
    s  = ode45(@(~, x) zoh_rhs(x, u, p), [tk te], xk, opts);

    t_all = [t_all, s.x(2:end)];        %#ok<AGROW>
    X_all = [X_all, s.y(:, 2:end)];     %#ok<AGROW>

    dt_actual = s.x(end) - tk;          % shorter than dt if the guard fired
    if dt_actual <= 0
        break;                          % solver could not advance: give up
    end
    tk = s.x(end);
    xk = s.y(:, end);

    % ---- advance the controller state over the SAME interval ------------
    if stateful
        xik = ch4_l1_advance(xik, ci.eta, ci.l1.mu1, ci.l1.mu1_hat, ...
                             clf, p, dt_actual);
        XI_all = [XI_all, xik];         %#ok<AGROW>
        t_xi   = [t_xi, tk];            %#ok<AGROW>
    end

    if ~isempty(s.ie)
        fired = true;
        break;
    end
end

if ~stateful
    XI_all = zeros(0, 1);
    t_xi   = 0;
end
end

% ---------------------------------------------------------------------------
function xdot = zoh_rhs(x, u, p)
%ZOH_RHS  TRUE plant only: the control is already decided for this period.
[f, g] = ch4_control_affine(x, p);
xdot   = f + g*u;
end

% ---------------------------------------------------------------------------
function [value, isterminal, direction] = guard_event(t, z, p)
% z may be the plain 14x1 state (ZOH inner solve) or the augmented [x; xi].
if t < p.guard_min_time
    value = 1;                       % hold positive; the foot has not lifted
else
    value = ch3_guard(z(1:p.nx), p);
end
isterminal = 1;
direction  = -1;
end
