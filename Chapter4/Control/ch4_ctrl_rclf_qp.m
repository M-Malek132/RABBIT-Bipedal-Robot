function [mu, u, qp] = ch4_ctrl_rclf_qp(Lf2y, LgLfy, u_ff, info, p, constrained)
%CH4_CTRL_RCLF_QP  Section 4.1: the robust CLF-QP, eq (4.12) and (4.13).
%
%   [mu, u, qp] = ch4_ctrl_rclf_qp(Lf2y, LgLfy, u_ff, info, p, constrained)
%
% All inputs come from the NOMINAL model (ch4_io_lin with unc = []). This
% controller never sees the true plant; that is the whole point.
%
% ---------------------------------------------------------------- the setup
% Under uncertainty the transverse dynamics (4.5) are no longer etadot =
% F eta + G mu but
%
%       etadot = F eta + [0; Delta1] + (G + [0; Delta2]) mu
%
% so from (4.8), and using G = [0; I] so that [0; Delta1] = G*Delta1,
%
%       Vdot = LfV + LgV*Delta1 + LgV*(I + Delta2)*mu
%
% with LfV = eta'(F'Peps + Peps F)eta and LgV = 2 eta' Peps G exactly as in
% Chapter 3. Writing psi = LfV + (c3/eps)V for the Chapter-3 residual, the
% robust RES condition (4.11) is
%
%       max          psi + LgV*Delta1 + LgV*mu + LgV*Delta2*mu  <=  0.
%    ||Delta_i||<=..
%
% ------------------------------------------------- taking the max in closed form
% The Delta1 term is a linear functional over a ball, so its max is immediate
% and is a CONSTANT -- it does not involve mu at all:
%
%       max_{||Delta1|| <= D1} LgV*Delta1 = D1*||LgV||.
%
% Fold it into psi and define the ROBUSTIFIED RESIDUAL
%
%       a := psi + D1*||LgV||.
%
% The Delta2 term does involve mu, and how its max behaves depends on what
% "||Delta2|| <= D2" is taken to mean. p.rclf.delta2_model selects:
%
%   'scalar'   Delta2 = d2*I with |d2| <= D2.  Then LgV*Delta2*mu = d2*(LgV*mu)
%              and the max over the interval is D2*|LgV*mu|. A bound on an
%              absolute value is TWO LINEAR INEQUALITIES, one per sign:
%
%                   a + (1 + D2)*(LgV*mu) <= 0
%                   a + (1 - D2)*(LgV*mu) <= 0
%
%              This is the reduction that keeps (4.12) an honest QP (Remark 4.5).
%
%   'matrix'   the unstructured ball ||Delta2||_2 <= D2, the literal reading of
%              (4.10). The max is D2*||LgV||*||mu||_2, a SECOND-ORDER CONE
%              constraint:
%
%                   a + LgV*mu + D2*||LgV||*||mu||_2 <= 0.
%
% ------------------------------------------- the two models agree where it counts
% For the UNCONSTRAINED problem (4.12) both models give the SAME control, and
% it is available in closed form. By symmetry the least-norm feasible point
% lies along -LgV', so put mu = -t*LgV'/||LgV|| with t >= 0:
%
%   scalar:  the (1-D2) row binds (it is the weaker coefficient), giving
%            a - (1-D2)*t*||LgV|| <= 0
%   matrix:  a - t*||LgV|| + D2*||LgV||*t <= 0, the same inequality.
%
% Either way t = a / ((1-D2)*||LgV||) and
%
%       mu* = - a*LgV' / ((1 - D2)*||LgV||^2)         when a > 0
%       mu* = 0                                       when a <= 0.
%
% Setting D1 = D2 = 0 recovers Chapter 3's mu* = -(psi/||LgV||^2)*LgV' exactly,
% so the robust controller is a strict generalization rather than a variant.
% The models differ only in the CONSTRAINED problem, where mu is no longer free
% to point along -LgV'.
%
% --------------------------------------------------------------- feasibility
% The closed form needs 1 - D2 > 0. That is not an artifact of the algebra: at
% D2 = 1 the worst-case model inside the bound can cancel the control's entire
% effect on V, and past it, reverse the sign -- the controller would be pushing
% the wrong way. So
%
%       Delta2max < 1
%
% is the precise form of Remark 4.3's "chosen within this region", and
% ch4_delta_bounds reports it. Violating it is reported through qp.feasible
% rather than silently producing a number.
%
% ------------------------------------------------------------ the price paid
% Note a >= psi always, and strictly greater whenever D1 > 0 and eta ~= 0. The
% robust controller therefore demands MORE decrease than Chapter 3 does even
% when the model is perfect -- it is defending against a worst case that is not
% happening. That is the limitation Section 4.1.4 closes on and the entire
% motivation for the L1 controller in Section 4.2.
%
% ------------------------------------- the exact law is a sliding-mode law
% Near the orbit psi is O(||eta||^2) while D1*||LgV|| is O(||eta||), so a is
% dominated by the robust term and the closed form tends to
%
%       mu*  ->  - M * LgV' / ||LgV||,        M = D1 / (1 - D2),
%
% a correction of FIXED magnitude M (574 rad/s^2 on posture_195) along a unit
% vector that flips whenever the state crosses the hyperplane LgV = 0, however
% small eta is. That is unit-vector sliding-mode control. In continuous time
% it slides; held for a sample period T it cannot, because one sample of it
% moves LgV by
%
%       phi_T = ||2 P22|| * M * T,          P22 = the ydot-ydot block of Peps
%
% (0.029 on posture_195 at 1 kHz), and whatever was left of LgV is overshot. So
% the sampled law chatters at the sample rate; the numbers are below.
%
% ---------------------------------------------------------- the boundary layer
% The standard remedy is to replace the unit vector by a saturation over a
% layer of thickness phi in ||LgV||:
%
%       D1*||LgV||   ->   D1*||LgV|| * min(1, ||LgV||/phi)
%
% Outside the layer nothing changes. Inside it the robust term becomes
% D1||LgV||^2/phi, so the robust part of mu* is the LINEAR feedback
% -(M/phi) LgV', which vanishes continuously at LgV = 0 instead of flipping.
%
% THE LAYER IS SIZED IN SAMPLES. One held sample of that linear feedback, on a
% plant whose true input gain is (1 + d2) with d2 in [-D2, D2], multiplies LgV by
%
%       1 - (1 + d2) * phi_T / phi,
%
% so the thickness at which the WORST gain in the bound, d2 = D2, lands exactly
% on LgV = 0 in one sample is (1 + D2)*phi_T. The layer is set in those units,
%
%       phi = kappa * (1 + D2) * phi_T,         kappa = p.rclf.boundary_layer,
%
% and the multiplier is at worst 1 - 1/kappa: the sampled loop along LgV
% converges for every d2 in the bound iff kappa > 1/2, and does so WITHOUT
% overshoot -- without chatter -- iff kappa >= 1. Sizing phi in these units
% keeps that statement true when T, the bounds or the CLF change; a fixed phi
% would silently re-open the chatter at a longer sample period.
%
% MEASURED on posture_195 at 1 kHz and eps = 0.35, the first 3 steps of
% 'rclfqp_con' in Cases I-III (scales 1 / 1.5 / 0.7; the multiplier itself
% does not depend on eps), as the median change in the held torque from one
% sample to the next: 253 / 573 / 481 Nm with no layer. At kappa = 0.33 the
% multiplier predicts -1.0 / -0.33 / -1.86 for the three true gains, and the
% x0.7 case still chattered at 514 Nm while the other two fell to a 1-2 Nm
% median. From kappa = 0.66 up, 1-2 Nm in every case.
%
% Two consequences worth knowing:
%   * Inside the layer the robust correction -LgV'/(kappa (1+D2) ||2 P22|| T)
%     does not depend on D1. Near the orbit a larger Delta1 bound no longer
%     buys a larger control; only the layer's edge moves out.
%   * kappa = 0, or control_dt = 0 (continuous control, nothing to overshoot),
%     gives phi = 0: the exact law of (4.12)/(4.13), chatter included.
%
% WHAT THE GUARANTEE BECOMES. The QP now enforces the smoothed row, so the
% robust RES condition (4.11) holds EXACTLY wherever ||LgV|| >= phi, and inside
% the layer the worst case can exceed it by at most
%
%       gap = D1*||LgV||*(1 - ||LgV||/phi)  <=  D1*phi/4,
%
% which qp.bl_gap reports per call. Vdot <= -(c3/eps)V + D1*phi/4 (plus the
% slack, in (4.13)) then makes the tracking error uniformly ultimately bounded
% rather than convergent -- the honest version of "errors go to zero" for any
% implementation that holds its control for a sample.
%
% -------------------------------------------------- constrained form (4.13)
% Identical structure to Chapter 3 stage 8: decision vector [u; delta], the CLF
% row relaxed by a penalized slack, and torque / friction / GRF rows added.
%
% REMARK 4.4 APPLIES AND IS NOT COSMETIC. The added constraints are evaluated
% on the NOMINAL model:
%   * the torque box IS invariant to model uncertainty -- u is computed from
%     the nominal model and applied verbatim, so a bound on it is exact;
%   * the friction cone and the minimum normal force are NOT. They use the
%     nominal lam_drift, lam_in, and the true contact force differs. Enabling
%     them here constrains a PREDICTION, not the real force. Chapter 8 is where
%     that gets fixed; until then, treat those two as advisory.
% qp.robust_constraints records which of the active rows are actually robust.
%
% Inputs
%   Lf2y, LgLfy, u_ff, info : from ch4_io_lin on the NOMINAL model
%   p                       : parameter struct (uses p.rclf, p.eps, p.limits,
%                             p.control_dt)
%   constrained             : logical, eq (4.13) if true
%
% Outputs
%   mu : ny x 1 virtual input
%   u  : nu x 1 joint torque
%   qp : struct .V .LfV .LgV .psi .a .delta .active .exitflag .feasible
%               .margin (how much slack the robust row had) .robust_constraints
%               .phi (boundary-layer thickness, 0 = exact law)
%               .bl_gap (how far the worst case may exceed (4.11) at this
%               state because of the layer; 0 outside it)
%
% See also CH3_CTRL_CLF_QP, CH4_UNCERTAINTY, CH4_DELTA_BOUNDS.

if nargin < 6, constrained = false; end

clf = ch3_res_clf(p);
eta = info.eta;

[V, LfV, LgV] = ch3_clf_eval(eta, clf, p.eps);
psi = LfV + (clf.c3 / p.eps) * V;

D1 = p.rclf.delta1_max;
D2 = p.rclf.delta2_max;

nrm_LgV = norm(LgV, 2);

% Boundary layer, sized in samples -- see the header. D2 >= 1 has no finite M
% and is reported infeasible below, so it gets no layer. A parameter struct
% saved before the layer existed (an old ch4_result .mat) has no field, and
% re-analysing it must reproduce the exact law that run actually used.
kappa = 0;
if isfield(p.rclf, 'boundary_layer'), kappa = p.rclf.boundary_layer; end

phi = 0;
if kappa > 0 && D2 < 1 && p.control_dt > 0
    M_rob = D1 / (1 - D2);
    phi_T = norm(clf.PG2_eps(p.ny+1:end, :), 2) * M_rob * p.control_dt;
    phi   = kappa * (1 + D2) * phi_T;
end

rob = D1 * nrm_LgV;                    % max over the Delta1 ball
if nrm_LgV < phi
    rob = rob * (nrm_LgV / phi);       % saturated inside the layer
end
a = psi + rob;                         % robustified residual

qp = struct('V', V, 'LfV', LfV, 'LgV', LgV, 'psi', psi, 'a', a, ...
            'delta', 0, 'active', false, 'exitflag', 1, 'feasible', true, ...
            'margin', 0, 'robust_constraints', true, ...
            'phi', phi, 'bl_gap', D1 * nrm_LgV - rob);

if D2 >= 1
    % Report rather than divide by a non-positive number. The caller sees
    % feasible = false and the fallback below keeps the simulation running.
    qp.feasible = false;
    qp.exitflag = -3;
    mu = zeros(p.ny, 1);
    if nargout > 1, u = u_ff; end
    return;
end

%% ------------------------------------------------------ unconstrained (4.12)
if ~constrained
    if a <= 0 || nrm_LgV < 1e-12
        mu = zeros(p.ny, 1);           % worst case already decreasing enough
    else
        mu = -(a / ((1 - D2) * nrm_LgV^2)) * LgV.';
        qp.active = true;
    end
    qp.margin = -(a + LgV*mu + D2*worst_case_term(LgV, mu, p));
    if nargout > 1
        u = u_ff + solve_decoupling(LgLfy, mu, info);
    end
    return;
end

%% -------------------------------------------------------- constrained (4.13)
nu  = p.nu;
A   = LgLfy;
AtA = A.' * A;

% Same relative-scaling argument as ch3_ctrl_clf_qp: the decoupling matrix
% carries the inverse inertia, so pairing AtA with a literal 1e6 would hand
% quadprog a hopelessly conditioned Hessian and make delta meaningless.
mu_scale = mean(diag(AtA));
if ~isfinite(mu_scale) || mu_scale <= 0, mu_scale = 1; end
p_slack = p.clf_slack_penalty * mu_scale;

bA = LgV * A;                          % 1 x nu, the row acting on u

switch lower(p.rclf.delta2_model)

    case 'scalar'
        % z = [u; delta].  Two robust CLF rows, one per sign of d2.
        n_extra = 0;
        H  = 2 * blkdiag(AtA, p_slack);
        fq = [-2 * (AtA * u_ff); 0];

        Aineq = [ (1 + D2) * bA, -1 ; ...
                  (1 - D2) * bA, -1 ];
        bineq = [ -a + (1 + D2) * (bA * u_ff) ; ...
                  -a + (1 - D2) * (bA * u_ff) ];

    case 'matrix'
        % z = [u; delta; w], w in R^ny bounding |mu| elementwise so that
        % sum(w) >= ||mu||_1 >= ||mu||_2.
        %
        % CONSERVATISM, STATED PLAINLY. This replaces the exact SOC term
        % D2*||LgV||*||mu||_2 with D2*||LgV||*sum(w), an upper bound that can
        % overshoot by up to sqrt(ny) = 2. The resulting control is therefore
        % SAFE but may be more aggressive than the exact min-max solution. The
        % 'scalar' model has no such gap, which is the other reason it is the
        % default. The UNCONSTRAINED path above is exact for both models.
        n_extra = p.ny;
        H  = 2 * blkdiag(AtA, p_slack, zeros(n_extra));
        fq = [-2 * (AtA * u_ff); 0; zeros(n_extra,1)];

        cone = D2 * nrm_LgV;

        %   a + LgV*A*(u - u_ff) + cone*sum(w) - delta <= 0
        Aineq = [ bA, -1, cone * ones(1, n_extra) ];
        bineq = -a + bA * u_ff;

        %   w_i >= +(A(u-u_ff))_i    ->    A(i,:)u - w_i <=  A(i,:)u_ff
        %   w_i >= -(A(u-u_ff))_i    -> -A(i,:)u - w_i <= -A(i,:)u_ff
        Aineq = [ Aineq ; ...
                  A,  zeros(p.ny,1), -eye(n_extra) ; ...
                 -A,  zeros(p.ny,1), -eye(n_extra) ];
        bineq = [ bineq ; A * u_ff ; -A * u_ff ];

    otherwise
        error('ch4_ctrl_rclf_qp:delta2Model', ...
              'Unknown p.rclf.delta2_model "%s" (expected scalar|matrix).', ...
              p.rclf.delta2_model);
end

H = (H + H.')/2;

% --- friction cone and minimum normal force ------------------------------
% NOMINAL quantities: see the Remark 4.4 note in the header.
lam_d = info.aux.lam_drift;
lam_i = info.aux.lam_in;
pad   = zeros(1, n_extra);

if p.limits.enable.friction
    ms = p.limits.mu_s;
    Aineq = [Aineq; ...
             [ lam_i(1,:) - ms*lam_i(2,:), 0, pad]; ...
             [-lam_i(1,:) - ms*lam_i(2,:), 0, pad]];
    bineq = [bineq; ...
             -( lam_d(1) - ms*lam_d(2)); ...
             -(-lam_d(1) - ms*lam_d(2))];
    qp.robust_constraints = false;
end

if p.limits.enable.grf
    Aineq = [Aineq; [-lam_i(2,:), 0, pad]];
    bineq = [bineq;  lam_d(2) - p.limits.Fz_min];
    qp.robust_constraints = false;
end

% --- torque box (genuinely invariant to the uncertainty) -----------------
if p.limits.enable.torque
    lb = [-p.limits.u_max * ones(nu,1); 0; zeros(n_extra,1)];
    ub = [ p.limits.u_max * ones(nu,1); inf; inf(n_extra,1)];
else
    lb = [-inf(nu,1); 0; zeros(n_extra,1)];
    ub = [ inf(nu,1); inf; inf(n_extra,1)];
end

z0 = [u_ff; 0; zeros(n_extra,1)];
[z, ~, exitflag] = quadprog(H, fq, Aineq, bineq, [], [], lb, ub, z0, qp_options());

if exitflag <= 0 || isempty(z)
    u = u_ff;
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

qp.delta    = z(nu+1);
qp.exitflag = exitflag;
qp.active   = a > 0;
qp.margin   = qp.delta - (a + LgV*mu + D2*worst_case_term(LgV, mu, p));

end

% ---------------------------------------------------------------------------
function w = worst_case_term(LgV, mu, p)
%WORST_CASE_TERM  max over the Delta2 ball of LgV*Delta2*mu, per the model.
%
% Used only to report the realized margin, so the caller can check that the
% robust inequality it asked for is the one that actually holds.
switch lower(p.rclf.delta2_model)
    case 'scalar', w = abs(LgV * mu);
    otherwise,     w = norm(LgV, 2) * norm(mu, 2);
end
end

function v = solve_decoupling(A, mu, info)
if isfield(info,'rcond') && (~isfinite(info.rcond) || info.rcond < 1e-12)
    v = pinv(A) * mu;
else
    v = A \ mu;
end
end

function opts = qp_options()
%QP_OPTIONS  Built once and reused -- see the note in ch3_ctrl_clf_qp.
persistent o
if isempty(o)
    o = optimoptions('quadprog', 'Display', 'off', ...
                     'Algorithm', 'interior-point-convex');
end
opts = o;
end
