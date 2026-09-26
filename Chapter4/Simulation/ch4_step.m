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
%           .u        nu x m torque held over each control period (sampled
%                     control only; empty for continuous control)
%           .t_u      1 x m start time of each period
%           .lambda   2 x nt TRUE contact force [Fx; Fz] at every point of .t,
%                     from the true model under the torque held there. The
%                     simulation pins the stance foot, so nothing here stops
%                     Fz < 0 or |Fx| > mu Fz; this is what lets a run be checked
%                     for physical validity afterwards (ch4_validity) at the
%                     full solver resolution, for the price of one contact
%                     solve per point.
%           .contact_invalid, .invalid   as ch3_step: with p.stop_on_invalid
%                     the step ends AT its first sample the ground could not
%                     supply (Fz <= 0 or |Fx|/Fz > p.limits.mu_s) and is not ok
%           .u_last   the last COMMANDED torque, which a one-sample actuation
%                     delay applies at the start of the next step
%
% IMPLEMENTATION EFFECTS of a structured uncertainty (ch4_uncertainty_set),
% sampled control only:
%   uncertainty.noise   .q .dq std of the measurement noise; the controller
%                       (and the L1 predictor) sees x + noise, the plant and
%                       the guard the true x. One draw per sample, from
%                       .noise.stream (ch4_simulate creates it from .seed).
%   uncertainty.delay   1: the robot receives each command one sample late
%                       (the controller does not know); 0: none.
%   uncertainty.u_prev  the command pending from the previous step.
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
    [t, X, XI, t_xi, fired, U, t_u, LAM, inv, u_last] = ...
        integrate_zoh(x0, xi0, alpha, p, T_cap, opts);
    sol = [];
else
    U = zeros(p.nu, 0); t_u = zeros(1, 0); LAM = zeros(2, 0);
    inv = []; u_last = [];
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
out.u      = U;
out.t_u    = t_u;
out.lambda = LAM;
out.contact_invalid = ~isempty(inv);
out.invalid = inv;
out.u_last  = u_last;

foot_st_0   = P_st(x0(1:p.nq));
foot_sw_end = P_sw(out.x_end(1:p.nq));
out.L_step  = foot_sw_end(1) - foot_st_0(1);

% --- reset map on the TRUE model -----------------------------------------
[out.x_next, out.impulse] = ch4_impact(out.x_end, p);

% --- carry the controller state across the impact ------------------------
% eta jumps here. p.l1.reset_predictor decides whether the predictor is told;
% the parameter estimates are carried regardless (see ch4_l1_state).
if stateful
    [~, ~, o_minus] = ch3_outputs(out.x_end,  alpha, p);
    [~, ~, o_plus]  = ch3_outputs(out.x_next, alpha, p);
    out.xi_next     = ch4_l1_state('reset', p, out.xi_end, o_plus.eta, o_minus.eta);
else
    out.xi_next = zeros(0,1);
end

end

% ---------------------------------------------------------------------------
function [t_all, X_all, XI_all, t_xi, fired, U_all, t_u, LAM_all, inv, u_last] = integrate_zoh(x0, xi0, alpha, p, T_cap, opts)
%INTEGRATE_ZOH  One control decision per period; plant and controller advance.
% Also records the torque the robot received in each period and the TRUE
% contact force at every solver point under it (see .lambda in the header).

dt = p.control_dt;

stateful = ch4_is_stateful(p);
clf      = ch3_res_clf(p);

plant_predictor = false;
if stateful
    l1o = ch4_l1_opts(p);
    plant_predictor = strcmp(l1o.predictor, 'plant');
end

% measurement noise and actuation delay of a structured uncertainty
[sig_q, sig_dq, rs, delay, u_hold] = implementation(p);
measure = @(x) x;
if sig_q > 0 || sig_dq > 0
    measure = @(x) x + [sig_q * randn(rs, p.nq, 1); sig_dq * randn(rs, p.nq, 1)];
end

stop_inv = isfield(p, 'stop_on_invalid') && ~isempty(p.stop_on_invalid) && p.stop_on_invalid;
n_seen   = 0;
inv      = [];

t_all  = 0;
X_all  = x0(:);
XI_all = xi0(:);
t_xi   = 0;
U_all   = zeros(p.nu, 0);
t_u     = zeros(1, 0);
LAM_all = nan(2, 1);                    % t = 0 is filled by the first period

tk    = 0;
xk    = x0(:);
xik   = xi0(:);
fired = false;
x_meas = measure(xk);                   % what the controller sees at t_0

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
    % from the MEASURED state (the true one when there is no noise)
    [u, ~, ci] = ch4_control(x_meas, xik, alpha, p);

    if ~all(isfinite(u))
        break;
    end

    % A one-sample actuation delay: the robot receives the previous command.
    % The first period of a run has no previous command and applies its own.
    if delay > 0
        if isempty(u_hold), u_hold = u; end
        u_apply = u_hold;
        u_hold  = u;
    else
        u_apply = u;
    end

    te = min(tk + dt, T_cap);
    s  = ode45(@(~, x) zoh_rhs(x, u_apply, p), [tk te], xk, opts);

    t_all = [t_all, s.x(2:end)];        %#ok<AGROW>
    X_all = [X_all, s.y(:, 2:end)];     %#ok<AGROW>

    U_all = [U_all, u_apply(:)];        %#ok<AGROW>
    t_u   = [t_u, tk];                  %#ok<AGROW>
    if tk == 0
        LAM_all(:, 1) = true_lambda(xk, u_apply, p);
    end
    lam_s = zeros(2, numel(s.x) - 1);
    for j = 2:numel(s.x)
        lam_s(:, j-1) = true_lambda(s.y(:, j), u_apply, p);
    end
    LAM_all = [LAM_all, lam_s];         %#ok<AGROW>

    % End the step at its first contact-invalid sample (p.stop_on_invalid).
    if stop_inv
        [j, inv] = first_invalid(t_all, LAM_all, n_seen, p.limits.mu_s);
        if ~isempty(j)
            t_all   = t_all(1:j);
            X_all   = X_all(:, 1:j);
            LAM_all = LAM_all(:, 1:j);
            fired   = false;
            break;
        end
        n_seen = numel(t_all);
    end

    dt_actual = s.x(end) - tk;          % shorter than dt if the guard fired
    if dt_actual <= 0
        break;                          % solver could not advance: give up
    end
    tk = s.x(end);
    xk = s.y(:, end);
    x_meas = measure(xk);               % the sample at t_{k+1}, used twice below

    % ---- advance the controller state over the SAME interval ------------
    % The plant-input predictor is also told where eta ended up: xk has
    % already advanced, so sampling it here is what a digital controller
    % does at t_{k+1}, before it computes the next control (ch4_l1_advance).
    % It reads the same MEASUREMENT the next control decision will.
    if stateful
        smp = struct('eta', ci.eta, 'eta_next', [], 'mu', ci.mu, ...
                     'mu1_hat', ci.l1.mu1_hat);
        if plant_predictor
            [~, ~, o_next] = ch3_outputs(x_meas, alpha, p);
            smp.eta_next   = o_next.eta;
        end
        xik = ch4_l1_advance(xik, smp, clf, p, dt_actual);
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
u_last = u_hold;
end

% ---------------------------------------------------------------------------
function [sig_q, sig_dq, rs, delay, u_hold] = implementation(p)
%IMPLEMENTATION  Measurement noise, its random stream, and the actuation delay.
sig_q = 0; sig_dq = 0; rs = []; delay = 0; u_hold = [];
u = p.uncertainty;
if isempty(u) || ~isstruct(u), return; end
if isfield(u, 'noise') && isstruct(u.noise) && ~isempty(u.noise)
    if isfield(u.noise, 'q')  && ~isempty(u.noise.q),  sig_q  = u.noise.q;  end
    if isfield(u.noise, 'dq') && ~isempty(u.noise.dq), sig_dq = u.noise.dq; end
    if isfield(u.noise, 'stream') && ~isempty(u.noise.stream)
        rs = u.noise.stream;
    elseif sig_q > 0 || sig_dq > 0
        % called outside ch4_simulate: a stream from the seed, so the run is
        % still reproducible (each such call restarts it)
        sd = 0;
        if isfield(u.noise, 'seed') && ~isempty(u.noise.seed), sd = u.noise.seed; end
        rs = RandStream('mt19937ar', 'Seed', sd);
    end
end
if isfield(u, 'delay') && ~isempty(u.delay), delay = u.delay; end
if ~(delay == 0 || delay == 1)
    error('ch4_step:delay', ...
          'uncertainty.delay must be 0 or 1 control sample (got %g).', delay);
end
if isfield(u, 'u_prev'), u_hold = u.u_prev; end
end

% ---------------------------------------------------------------------------
function [j, inv] = first_invalid(t_all, LAM, n_seen, mu_s)
%FIRST_INVALID  First sample after n_seen the ground could not have supplied.
% The test is ch3_validity's: Fz <= 0 (lift-off) or |Fx|/Fz > mu_s (slip).
j = []; inv = [];
cols = n_seen+1 : size(LAM, 2);
Fx = LAM(1, cols);  Fz = LAM(2, cols);
ok   = isfinite(Fx) & isfinite(Fz);
lift = ok & (Fz <= 0);
slip = ok & (Fz > 0) & (abs(Fx) ./ max(Fz, realmin) > mu_s);
k = find(lift | slip, 1);
if isempty(k), return; end
j = cols(k);
if lift(k)
    kind = 'lift-off'; mu = NaN;
else
    kind = 'slip';     mu = abs(Fx(k)) / Fz(k);
end
inv = struct('t', t_all(j), 'kind', kind, 'Fz', Fz(k), 'mu', mu);
end

% ---------------------------------------------------------------------------
function lam = true_lambda(x, u, p)
%TRUE_LAMBDA  Stance contact force of the TRUE model at x under held torque u.
[~, ~, aux] = ch4_control_affine(x, p);
lam = aux.lam_drift + aux.lam_in * u;
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
