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
% THE ONE STATE THE CLOSED FORM CANNOT ANSWER.  If LgV = 0 while psi > 0, the
% constraint reads 0 <= -psi < 0: no mu satisfies it, and the least-norm
% formula divides by zero. That is NOT "the constraint already holds", and
% returning mu = 0 as though it were -- which this used to do -- reports a
% satisfied CLF condition at a state where it is unmeetable. Flagged instead,
% through qp.feasible.
%
% DEFENSIVE, NOT A LIVE PATH. The check is exact: psi is the quadratic form
% eta'(M_eps + rate P_eps)eta and LgV vanishes precisely on Z = null(G'P_eps),
% so "can psi > 0 while LgV = 0" is decided by eig(Z'(M_eps + rate P_eps)Z) --
% MEASURED negative for both constructions at every eps in [0.05, 2], so
% neither can reach the state on this robot. The guard costs one comparison and
% covers 'lyap', which unlike 'care' has no theorem behind it (ch3_res_clf).
% Same spirit as the pinv fallback in ch3_io_lin.
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
%   qp : struct .V .LfV .LgV .psi .delta .res .active .exitflag .feasible
%        .res is the ACHIEVED CLF residual, psi + LgV mu - delta: <= 0 means
%        the rapid-exponential condition was met at this state.
%
% See also CH3_RES_CLF, CH3_CLF_EVAL, CH3_CTRL_PD.

if nargin < 6, constrained = false; end

clf = ch3_res_clf(p);
eta = info.eta;

[V, LfV, LgV] = ch3_clf_eval(eta, clf, p.eps);
psi = LfV + clf.rate * V;                % the RES-CLF residual to be made <= 0

qp = struct('V', V, 'LfV', LfV, 'LgV', LgV, 'psi', psi, ...
            'delta', 0, 'res', psi, 'active', false, ...
            'exitflag', 1, 'feasible', true);

%% ---------------------------------------------------------------- stage 7
if ~constrained
    nrm2 = LgV * LgV.';

    if psi <= 0
        mu = zeros(p.ny, 1);             % already decreasing fast enough
    elseif nrm2 < lgv_tol()
        % Unmeetable: the constraint is 0 <= -psi with psi > 0. Apply nothing
        % and say so -- see the header note.
        mu = zeros(p.ny, 1);
        qp.feasible = false;
        qp.exitflag = -2;                % quadprog's "no feasible point"
    else
        mu = -(psi / nrm2) * LgV.';      % least-norm point on the halfspace
        qp.active = true;
    end

    qp.res = psi + LgV * mu;
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
% SCALING.  p.clf_slack_penalty is a RELATIVE weight, not an absolute one.
% AtA carries the SQUARE of the decoupling matrix, so its magnitude is set by
% the robot's inertia and by whatever units u is in -- MEASURED on the
% reference gait, mean(diag(AtA)) = 35.7 and cond(AtA) = 7.5e3. Pairing an
% inertia-dependent mu block with a literal 1e6 on delta hands quadprog a
% Hessian whose two blocks are four to five orders apart and makes the returned
% delta numerically meaningless. Scaling the penalty by the mu block's own
% magnitude makes the ratio exactly clf_slack_penalty, independent of the
% robot's inertia and of the units u happens to be in.
mu_scale = mean(diag(AtA));
if ~isfinite(mu_scale) || mu_scale <= 0, mu_scale = 1; end
p_slack = p.clf_slack_penalty * mu_scale;

H = 2 * blkdiag(AtA, p_slack);
H = (H + H.')/2;
fq = [-2 * (AtA * u_ff); 0];

% --- relaxed CLF row:  (LgV A) u - delta <= -psi + LgV A u_ff -------------
LgVA  = LgV * A;
Aineq = [LgVA, -1];
bineq = -psi + LgVA * u_ff;

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

% --- starting point ------------------------------------------------------
% NO WARM START, AND THAT IS A MEASURED CHOICE. Chapter 5 starts its QPs from
% the stage-7 closed form and documents the ~1% of solver failures that fixes
% (ch5_min_norm_mu). Tried here, it does nothing: over 400 QPs at perturbed
% states with the torque box swept 3..300 Nm, cold and warm gave the same
% failure count (6 of 400), agreed to 3.8e-07, and cold was fractionally faster
% (154.6 vs 157.7 us) since the warm start costs an extra 4x4 solve.
%
% The difference is what the two QPs contain: Chapter 5's barrier row drags mu
% far from the CLF's own answer (psi ~ 1e5, ||mu|| ~ 1e3), so the starting
% point is most of the problem. Stage 8 only trims the CLF answer to fit a box
% and a cone, so either start is already close.
z0 = [u_ff; 0];

[z, ~, exitflag] = quadprog(H, fq, Aineq, bineq, [], [], lb, ub, z0, qp_options());

if exitflag <= 0 || isempty(z)
    % Infeasible or failed: fall back to the saturated feedforward so the ODE
    % keeps integrating, and record it. With the slack variable present this
    % should essentially never fire -- delta can always absorb the CLF row --
    % so a nonzero count here means a genuinely infeasible torque/friction set,
    % i.e. no torque inside the box keeps the foot inside the cone. Note the
    % fallback does not pretend otherwise: clamping u_ff satisfies the box but
    % need not satisfy the cone, and qp.feasible = false is the honest report.
    u  = u_ff;
    if p.limits.enable.torque
        u = min(max(u, -p.limits.u_max), p.limits.u_max);
    end
    mu = A * (u - u_ff);
    qp.res      = psi + LgV * mu;
    qp.exitflag = exitflag;
    qp.feasible = false;
    return;
end

u  = z(1:nu);
mu = A * (u - u_ff);

qp.delta    = z(end);
qp.exitflag = exitflag;
qp.res      = psi + LgV * mu - qp.delta;
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

function t = lgv_tol()
%LGV_TOL  Floor on ||LgV||^2 below which the CLF row has no usable gradient.
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
%
% ---------------------------------------------- why 'active-set', not the default
% Five variables and at most four rows is precisely the dense small-scale
% regime quadprog's 'active-set' algorithm is written for. MEASURED against the
% interior-point-convex this used to run:
%
%       time on the reference state    234.6 us  ->  154.6 us   (1.5x)
%       failures over 400 hard QPs       8/400  ->    6/400
%
% (the 400 being perturbed states with the torque box swept 3..300 Nm, so
% deliberately including boxes no torque can satisfy). Faster and no less
% robust. Chapter 5 hit the stronger version: interior-point-convex returned
% exitflag -3 ("unbounded") on a strictly convex QP for 222 of 30001 samples,
% all of which active-set solved. See ch5_solve_qp.
persistent o
if isempty(o)
    o = optimoptions('quadprog', 'Display', 'off', ...
                     'Algorithm', 'active-set');
end
opts = o;
end
