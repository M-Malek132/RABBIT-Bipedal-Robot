function [mu, u, qp] = ch3_ctrl_clf_qp(Lf2y, LgLfy, u_ff, info, p, constrained)
%CH3_CTRL_CLF_QP  Stages 7 and 8: the CLF quadratic program.
%
%   [mu, u, qp] = ch3_ctrl_clf_qp(Lf2y, LgLfy, u_ff, info, p, constrained)
%
% STAGE 7 -- unconstrained CLF-QP (constrained = false).
% Because Vdot is affine in mu, the stability requirement is a single linear
% inequality, so we can ask for the SMALLEST control that satisfies it:
%
%       min_mu  mu'mu   s.t.   LgV mu <= -LfV - (c3/eps) V
%
% A one-inequality least-norm QP has a closed form -- no solver needed:
%
%       mu* = 0                                if the constraint already holds
%       mu* = -(psi / ||LgV||^2) LgV'          otherwise,
%       psi = LfV + (c3/eps) V
%
% This is the pointwise minimum-norm RES-CLF control.  It is used here rather
% than calling quadprog because this function sits inside an ODE right-hand
% side that runs thousands of times per step; the closed form is exact and
% roughly two orders of magnitude cheaper.
%
% STAGE 8 -- constrained CLF-QP (constrained = true).  THE ACTUAL PAYOFF.
% Now that stability is a linear inequality, anything else affine in u can be
% bolted on:
%
%       min_{u,delta}  mu'mu + p_slack delta^2
%       s.t.   LgV mu <= -LfV - (c3/eps) V + delta        (relaxed CLF)
%              u_min <= u <= u_max                        (torque saturation)
%              |Fx| <= mu_s Fz                            (friction cone)
%              Fz   >= Fz_min                             (stay loaded)
%
% with mu tied to u by  mu = LgLfy (u - u_ff), so every constraint above is
% LINEAR IN u:  the GRF is lambda = lam_drift + lam_in u, straight out of the
% same KKT solve that produced f and g, so the friction cone and the minimum
% normal force are exact affine constraints, not approximations.
%
% THE SLACK.  delta, penalized by a large p_slack, lets the CLF condition bend
% when the torque box or the friction cone would otherwise make the QP
% infeasible.  This is precisely the thing PD control cannot do: PD has no
% mechanism to trade convergence rate against a saturating actuator; it just
% saturates and voids its own guarantee.
%
% THE TRADEOFF (Remark 3.2).  Relaxing the CLF condition can cost stability --
% the exponential bound holds only while delta stays at zero, and beyond that
% only under extra sufficient conditions.  qp.delta is returned so callers can
% report how often, and how hard, the relaxation was actually used.
%
% Inputs
%   Lf2y, LgLfy, u_ff, info : from ch3_io_lin
%   p           : parameter struct
%   constrained : logical, stage 8 if true
%
% Outputs
%   mu : ny x 1 virtual input actually applied
%   u  : nu x 1 joint torque
%   qp : struct .V .LfV .LgV .psi .delta .active .exitflag .feasible
%
% See also CH3_RES_CLF, CH3_CLF_EVAL, CH3_CTRL_PD.

if nargin < 6, constrained = false; end

clf = ch3_res_clf(p);
eta = info.eta;

[V, LfV, LgV] = ch3_clf_eval(eta, clf, p.eps);
psi = LfV + (clf.c3 / p.eps) * V;        % the RES-CLF residual to be made <= 0

qp = struct('V', V, 'LfV', LfV, 'LgV', LgV, 'psi', psi, ...
            'delta', 0, 'active', false, 'exitflag', 1, 'feasible', true);

%% ---------------------------------------------------------------- stage 7
if ~constrained
    nrm2 = LgV * LgV.';
    if psi <= 0 || nrm2 < eps_tol()
        mu = zeros(p.ny, 1);             % already decreasing fast enough
    else
        mu = -(psi / nrm2) * LgV.';      % least-norm point on the halfspace
        qp.active = true;
    end
    if nargout > 1
        u = u_ff + solve_decoupling(LgLfy, mu, info);
    end
    return;
end

%% ---------------------------------------------------------------- stage 8
% Decision vector z = [u; delta].  Working in u rather than mu makes the
% torque box a plain bound and the friction cone a plain linear row.
nu = p.nu;
A  = LgLfy;
AtA = A.' * A;

% cost: ||A(u - u_ff)||^2 + p_slack delta^2
%
% SCALING.  p.clf_slack_penalty is a RELATIVE weight, not an absolute one. The
% decoupling matrix carries the inverse inertia, so AtA's entries land around
% 1e-3 here; pairing that with a literal 1e6 on delta would hand quadprog a
% Hessian with a condition number near 1e9 and make the returned delta
% numerically meaningless. Scaling the penalty by the mu-block's own magnitude
% makes the ratio exactly clf_slack_penalty, independent of the robot's
% inertia and of the units u happens to be in.
mu_scale = mean(diag(AtA));
if ~isfinite(mu_scale) || mu_scale <= 0, mu_scale = 1; end
p_slack = p.clf_slack_penalty * mu_scale;

H = 2 * blkdiag(AtA, p_slack);
H = (H + H.')/2;
fq = [-2 * (AtA * u_ff); 0];

% --- relaxed CLF row:  (LgV A) u - delta <= -psi + LgV A u_ff -------------
Aineq = [LgV * A, -1];
bineq = -psi + LgV * A * u_ff;

% --- friction cone and minimum normal force ------------------------------
% lambda = lam_drift + lam_in u = [Fx; Fz], world frame, z up-positive.
lam_d = info.aux.lam_drift;
lam_i = info.aux.lam_in;

if p.limits.enable.friction
    ms = p.limits.mu_s;
    %  Fx - ms Fz <= 0   and   -Fx - ms Fz <= 0
    Aineq = [Aineq; ...
             [ lam_i(1,:) - ms*lam_i(2,:), 0]; ...
             [-lam_i(1,:) - ms*lam_i(2,:), 0]];
    bineq = [bineq; ...
             -( lam_d(1) - ms*lam_d(2)); ...
             -(-lam_d(1) - ms*lam_d(2))];
end

if p.limits.enable.grf
    % Fz >= Fz_min   ->   -lam_i(2,:) u <= lam_d(2) - Fz_min
    Aineq = [Aineq; [-lam_i(2,:), 0]];
    bineq = [bineq;  lam_d(2) - p.limits.Fz_min];
end

% --- torque box ----------------------------------------------------------
if p.limits.enable.torque
    lb = [-p.limits.u_max * ones(nu,1); 0];
    ub = [ p.limits.u_max * ones(nu,1); inf];
else
    lb = [-inf(nu,1); 0];
    ub = [ inf(nu,1); inf];
end

opts = qp_options();

z0 = [u_ff; 0];
[z, ~, exitflag] = quadprog(H, fq, Aineq, bineq, [], [], lb, ub, z0, opts);

if exitflag <= 0 || isempty(z)
    % Infeasible or failed: fall back to the saturated feedforward so the ODE
    % keeps integrating, and record it. With the slack variable present this
    % should essentially never fire -- delta can always absorb the CLF row --
    % so a nonzero count here means a genuinely infeasible torque/friction set.
    u  = u_ff;
    if p.limits.enable.torque
        u = min(max(u, -p.limits.u_max), p.limits.u_max);
    end
    mu = A * (u - u_ff);
    qp.exitflag = exitflag;
    qp.feasible = false;
    return;
end

u  = z(1:nu);
mu = A * (u - u_ff);

qp.delta    = z(end);
qp.exitflag = exitflag;
qp.active   = psi > 0;

end

% ---------------------------------------------------------------------------
function v = solve_decoupling(A, mu, info)
if isfield(info,'rcond') && (~isfinite(info.rcond) || info.rcond < 1e-12)
    v = pinv(A) * mu;
else
    v = A \ mu;
end
end

function t = eps_tol()
t = 1e-12;
end

function opts = qp_options()
%QP_OPTIONS  The quadprog options object, built once and reused.
%
% This is not premature micro-optimization. Measured on this machine, building
% the options object costs 2.80 ms while the 5-variable QP it configures solves
% in 0.98 ms -- so constructing it per call spent 74% of the controller's time
% on bookkeeping. Because this function sits inside an ODE right-hand side that
% runs thousands of times per step, hoisting it into a persistent turns a
% ~3.8 ms control step into a ~1.0 ms one, and it was the difference between
% ch3_compare_controllers finishing in ~1 min per step versus 20+.
%
% Keep this a plain persistent (not a struct field on p): the options object is
% immutable configuration, identical for every call, and threading it through
% the ODE would change the signature of every controller in the chain.
persistent o
if isempty(o)
    o = optimoptions('quadprog', 'Display', 'off', ...
                     'Algorithm', 'interior-point-convex');
end
opts = o;
end
