function F = ch4_forces(t, X, alpha, p, XI, t_xi, max_samples)
%CH4_FORCES  Recover torques, forces and estimator signals under uncertainty.
%
%   F = ch4_forces(t, X, alpha, p)
%   F = ch4_forces(t, X, alpha, p, XI, t_xi)
%   F = ch4_forces(t, X, alpha, p, XI, t_xi, max_samples)
%
% The Chapter-3 analogue with one correction that matters a great deal and one
% addition.
%
% ---------------------------------------------- the correction: WHICH GRF?
% ch3_forces recovers the ground reaction force as
%
%       lambda = ci.aux.lam_drift + ci.aux.lam_in * u
%
% where ci.aux came from the controller's own dynamics evaluation. In Chapter 3
% that is correct, because there is only one model. Here it would report the
% force the CONTROLLER PREDICTS, not the force the robot produces -- and under
% a mass scale of 1.5 those differ by 50%. Reporting the nominal prediction as
% though it were the measured force would make every friction and unilateral
% contact check meaningless in exactly the situation the chapter is about.
%
% So this function computes lambda from the TRUE model and returns the nominal
% prediction alongside it as F.lambda_nom. The gap between them is the concrete
% content of Remark 4.4: constraints written on the nominal model are not
% robust, and F.grf_pred_error measures by how much.
%
% -------------------------------------------------- the addition: estimator
% For the L1 laws, the controller state is a signal worth reading: theta_hat is
% the controller's running estimate of the model error, and eta_tilde is the
% prediction error that drives it. Pass XI and t_xi (from ch4_step /
% ch4_simulate) and both are recovered here.
%
% ZERO-ORDER HOLD LOOKUP. The controller state is defined on the SAMPLE grid,
% not the solver grid, and between samples it is held -- that is what a digital
% controller does. So xi is looked up with 'previous' interpolation, never
% linearly. Interpolating it linearly would invent a state the controller never
% had and would make theta_hat look smoother than it is.
%
% ------------------------------------------------------------- COST, AND WHY
% This function re-solves the controller at every sample it is given -- two QPs
% and several KKT factorizations apiece -- so it is FAR more expensive per point
% than the simulation that produced the points. A 3-step run at 1 kHz emits
% ~11000 solver points, and analysing all of them costs minutes per run while
% telling you nothing extra: the control was HELD at 1 kHz, so between samples
% every quantity here is either constant or a smooth interpolation of its
% endpoints.
%
% max_samples therefore decimates the grid (endpoints always kept) and F.t
% reports the grid actually used. Leave it Inf only when you genuinely need
% every solver point.
%
% Inputs
%   t     : 1 x nt time grid
%   X     : 14 x nt state trajectory
%   alpha : ny x n_ctrl coefficients
%   p     : parameter struct
%   XI    : 5ny x m controller state samples (optional; required for L1)
%   t_xi  : 1 x m sample times for XI
%   max_samples : cap on analysed points (default 2000; Inf for all)
%
% Outputs
%   F : struct with columns aligned to F.t (NOT to the input t unless no
%       decimation happened) --
%         .t          1 x nt the grid actually analysed
%         .u          nu x nt commanded joint torque [Nm]
%         .lambda     2 x nt  TRUE ground reaction [Fx; Fz] [N]
%         .lambda_nom 2 x nt  the controller's nominal prediction of it
%         .mu_req     1 x nt  |Fx|/Fz demanded of the TRUE contact
%         .y .ydot    ny x nt tracking error and rate
%         .V          1 x nt  CLF value
%         .delta      1 x nt  CLF slack (constrained laws)
%         .s          1 x nt  phase
%         .sw_height  1 x nt  swing foot height [m]
%         .theta_hat  ny x nt L1 uncertainty estimate  (NaN if not L1)
%         .eta_tilde  2ny x nt L1 prediction error     (NaN if not L1)
%         .mu2        ny x nt L1 adaptive component    (NaN if not L1)
%         .Delta1 .Delta2_n  the TRUE induced uncertainty (4.4) along the run
%         summary scalars: .int_u2 .torque_max .Fz_min .mu_max .qp_infeasible
%                          .delta_max .grf_pred_error .theta_err_max
%
% See also CH3_FORCES, CH4_UNCERTAINTY, CH4_CONTROL, CH4_SIMULATE.

if nargin < 5, XI = []; t_xi = []; end
if nargin < 7 || isempty(max_samples), max_samples = 2000; end

% --- decimate before doing any work --------------------------------------
if numel(t) > max_samples
    keep = unique([round(linspace(1, numel(t), max_samples)), numel(t)]);
    t = t(keep);
    X = X(:, keep);
end

nt = numel(t);
nu = p.nu;
ny = p.ny;

stateful = ch4_is_stateful(p);
if stateful && isempty(XI)
    error('ch4_forces:missingState', ...
          ['p.controller = "%s" carries state, so XI and t_xi are required. ' ...
           'Pass sim.xi and sim.t_xi from ch4_simulate.'], p.controller);
end

F = struct();
F.t          = t;
F.x          = X;      % the decimated states, so callers can stay aligned
F.u          = zeros(nu, nt);
F.lambda     = zeros(2,  nt);
F.lambda_nom = zeros(2,  nt);
F.mu_req     = zeros(1,  nt);
F.y          = zeros(ny, nt);
F.ydot       = zeros(ny, nt);
F.V          = nan(1,   nt);
F.delta      = zeros(1, nt);
F.s          = zeros(1, nt);
F.sw_height  = zeros(1, nt);
F.theta_hat  = nan(ny, nt);
F.eta_tilde  = nan(2*ny, nt);
F.mu2        = nan(ny, nt);
F.Delta1     = zeros(ny, nt);
F.Delta2_n   = zeros(1, nt);
F.theta_true = zeros(ny, nt);

n_infeas = 0;

for k = 1:nt
    x = X(:, k);

    xi = zeros(0,1);
    if stateful
        xi = zoh_lookup(XI, t_xi, t(k));
    end

    [u, ~, ci] = ch4_control(x, xi, alpha, p);

    F.u(:, k) = u;

    % --- TRUE contact force, not the controller's prediction --------------
    [~, ~, aux_true] = ch4_control_affine(x, p);
    lam     = aux_true.lam_drift + aux_true.lam_in * u;
    lam_nom = ci.aux.lam_drift   + ci.aux.lam_in   * u;

    F.lambda(:, k)     = lam;
    F.lambda_nom(:, k) = lam_nom;

    if lam(2) > 1e-6
        F.mu_req(k) = abs(lam(1)) / lam(2);
    else
        F.mu_req(k) = Inf;
    end

    F.y(:, k)    = ci.y;
    F.ydot(:, k) = ci.ydot;
    F.V(k)       = ci.V;
    F.delta(k)   = ci.delta;
    F.s(k)       = ci.o.s;
    if ~ci.qp_feasible, n_infeas = n_infeas + 1; end

    if stateful && ~isempty(ci.l1)
        F.theta_hat(:, k) = ci.l1.theta_hat;
        F.eta_tilde(:, k) = ci.l1.eta_tilde;
        F.mu2(:, k)       = ci.l1.mu2;
    end

    % --- the uncertainty the controller was actually facing ---------------
    D = ch4_uncertainty(x, alpha, p);
    F.Delta1(:, k) = D.Delta1;
    F.Delta2_n(k)  = D.n2;

    % The quantity the L1 estimator is actually chasing, eq (4.15):
    %   theta = Delta1 + Delta2*mu,
    % evaluated at the mu the controller really applied. Comparing theta_hat
    % against Delta1 alone would credit the estimator with an error it does
    % not have (or hide one it does) whenever mu is large.
    F.theta_true(:, k) = D.Delta1 + D.Delta2 * ci.mu;

    Psw = P_sw(x(1:p.nq));
    F.sw_height(k) = Psw(2);
end

u2 = sum(F.u.^2, 1);
if nt > 1
    F.int_u2 = trapz(t, u2);
else
    F.int_u2 = 0;
end

F.torque_max    = max(abs(F.u(:)));
F.Fz_min        = min(F.lambda(2, :));
F.Fz_max        = max(F.lambda(2, :));
F.mu_max        = max(F.mu_req(isfinite(F.mu_req)));
if isempty(F.mu_max), F.mu_max = Inf; end
F.qp_infeasible = n_infeas;
F.delta_max     = max(F.delta);

% How wrong the controller's own force prediction was -- the size of the gap
% Remark 4.4 warns about, in newtons.
F.grf_pred_error = max(vecnorm(F.lambda - F.lambda_nom, 2, 1));

% How well the L1 estimate tracked the uncertainty it was estimating, against
% the full theta = Delta1 + Delta2*mu rather than against Delta1 alone.
%
% Do not expect this to be small pointwise, and do not read a large value as a
% broken estimator. theta_hat is what the ADAPTATION holds; what reaches the
% joints is the filtered -C(s)*theta_hat, and the L1 guarantee (4.35) is a
% bound on the state and parameter errors, not a claim that theta_hat converges
% to theta pointwise. A fast transient at each footstrike is expected and is
% exactly what the filter exists to keep out of the torque.
if stateful
    F.theta_err_max = max(vecnorm(F.theta_hat - F.theta_true, 2, 1));
else
    F.theta_err_max = NaN;
end

end

% ---------------------------------------------------------------------------
function xi = zoh_lookup(XI, t_xi, tq)
%ZOH_LOOKUP  The controller state in force at time tq: the most recent sample.
idx = find(t_xi <= tq + 1e-12, 1, 'last');
if isempty(idx), idx = 1; end
xi = XI(:, min(idx, size(XI, 2)));
end
