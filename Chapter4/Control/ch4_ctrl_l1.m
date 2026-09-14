function [mu, u, xidot, l1] = ch4_ctrl_l1(Lf2y, LgLfy, u_ff, info, xi, p, constrained)
%CH4_CTRL_L1  Section 4.2: L1 adaptive control on a CLF-QP reference model.
%
%   [mu, u, xidot, l1] = ch4_ctrl_l1(Lf2y, LgLfy, u_ff, info, xi, p, constrained)
%
% All model quantities are NOMINAL (ch4_io_lin with unc = []). The controller
% never sees the true plant -- it infers the gap from the prediction error.
%
% ------------------------------------------------------------- the structure
% The applied virtual input is split (Section 4.2.2):
%
%       mu = mu1 + mu2
%
%   mu1  makes the system follow a REFERENCE MODEL. Chapter 3's min-norm
%        RES-CLF controller is that reference model, so mu1 is exactly the
%        stage-7 QP (4.21) -- or the torque-saturated QP (4.36) when
%        constrained. Note what this means: the reference model here is
%        NONLINEAR and has no closed-form expression, unlike the usual linear
%        reference model of textbook MRAC. That is the contribution of Section
%        4.2.2, and it is why the guarantee inherited is the RES rate rather
%        than a pole placement.
%
%   mu2  cancels the estimated uncertainty, mu2 = -C(s) theta_hat (4.23).
%        It is a STATE of the controller (the filter output), not an algebraic
%        term, so it is read from xi rather than recomputed here.
%
% If the estimate were perfect and the filter instantaneous, mu2 would exactly
% remove theta from (4.16) and leave eta_dot = F eta + G mu1, which is (4.22):
% the Chapter-3 closed loop, restored. Everything else in the method is about
% how close to that we actually get, and (4.35) is the answer.
%
% ----------------------------------------------- one QP solve, or two
% mu1 is the QP on the REAL eta. Under p.l1.predictor = 'thesis' a second QP,
% mu1_hat, is solved on the PREDICTOR state eta_hat, eq (4.21) / (4.37),
% because the thesis predictor is driven by its own copy of the reference
% model. Under 'plant' the predictor is driven by the input the robot actually
% received, so mu1_hat does not exist and each call solves one QP, not two. See
% ch4_l1_deriv for why the thesis form cannot deliver the error dynamics (4.28)
% claims for it.
%
% ---------------------------------------------- torque saturation (Sect 4.2.3)
% With constrained = true, mu1 comes from Chapter 3's stage-8 QP with the torque
% box at p.l1.u_max. What the box -- and any contact rows -- bound depends on
% p.l1.constrain_applied.
%
%   false  Section 4.2.3 as written: THE SATURATION IS ON mu1 ONLY, and only the
%          torque box is imposed. The realized torque is
%          u_ff + LgLfy^-1 (mu1 + mu2), and only the mu1 part was inside the
%          box; the adaptive term can push it back out.
%
%   true   the rows bound the torque the robot actually receives. mu2 enters
%          the QP as a known offset: the QP is handed u_ff + LgLfy^-1 mu2 as
%          its feedforward, so its cost and its CLF row still act on mu1 alone,
%          while the torque box and the friction and normal-force rows (as
%          p.limits.enable has them) are written on the total torque. That is
%          what the realizability problem needs. Measured at eps = 0.20 over
%          25 steps, 'l1' and 'l1_con' in the thesis form put the true normal
%          force below zero on 11-14% of the samples under either
%          perturbation, because mu2 acted outside every constraint. The rows
%          are on the NOMINAL model, so they are advisory in the sense of
%          Remark 4.4 -- but they are no longer bypassed.
%
% l1.u_box_excess reports how far the realized torque left the box, per call:
% whatever mu2 added under false, zero by construction under true.
%
% Inputs
%   Lf2y, LgLfy, u_ff, info : from ch4_io_lin on the NOMINAL model
%   xi                      : 5ny x 1 controller state (ch4_l1_state)
%   p                       : parameter struct
%   constrained             : logical, eq (4.36)-(4.37) if true
%
% Outputs
%   mu    : ny x 1 applied virtual input, mu1 + mu2
%   u     : nu x 1 joint torque
%   xidot : 5ny x 1 controller state derivative
%   l1    : struct .mu1 .mu1_hat .mu2 .theta_hat .eta_tilde .V .V_hat
%                  .delta .qp_feasible .u_box_excess .alpha_hat .beta_hat
%
% See also CH4_L1_DERIV, CH4_L1_ADVANCE, CH4_L1_STATE, CH3_CTRL_CLF_QP.

if nargin < 7, constrained = false; end

o   = ch4_l1_opts(p);
clf = ch3_res_clf(p);
st  = ch4_l1_state('unpack', p, xi);

eta = info.eta;
mu2 = st.mu2;

% --- the inner CLF-QP(s) --------------------------------------------------
% u_ref is the torque mu1 is measured from. It is u_ff itself unless the rows
% are to bound the applied torque, in which case mu2 is folded into it.
pq    = p;
u_ref = u_ff;
if constrained
    pq.limits.u_max         = p.l1.u_max;
    pq.limits.enable.torque = true;
    if o.constrain_applied
        u_ref = u_ff + solve_decoupling(LgLfy, mu2, info);
    else
        pq.limits.enable.friction = false;  % (4.36) specifies the box only
        pq.limits.enable.grf      = false;
    end
end

% mu1: reference model driven by the REAL transverse state
[mu1, ~, qp1] = ch3_ctrl_clf_qp(Lf2y, LgLfy, u_ref, info, pq, constrained);

% mu1_hat: the same QP driven by the PREDICTOR state -- thesis predictor only.
% Same Lf2y / LgLfy / u_ref: eq (4.37) uses the same decoupling matrix and
% feedforward, since those are properties of the nominal model at the current
% x, not of eta.
mu1_hat = nan(p.ny, 1);
qp1h    = struct('V', NaN, 'feasible', true);
if strcmp(o.predictor, 'thesis')
    info_hat     = info;
    info_hat.eta = st.eta_hat;
    [mu1_hat, ~, qp1h] = ch3_ctrl_clf_qp(Lf2y, LgLfy, u_ref, info_hat, ...
                                          pq, constrained);
end

% --- applied control ------------------------------------------------------
mu = mu1 + mu2;
u  = u_ff + solve_decoupling(LgLfy, mu, info);

% --- controller state derivative -----------------------------------------
% sig.mu is the input the plant receives at this instant, which is what the
% plant-input predictor is driven by; the thesis predictor reads mu1_hat.
sig = struct('eta', eta, 'mu', mu, 'mu1_hat', mu1_hat);
[xidot, d] = ch4_l1_deriv(xi, sig, clf, p);

if nargout > 3
    if constrained
        excess = max(max(abs(u)) - p.l1.u_max, 0);
    else
        excess = 0;
    end
    l1 = struct('mu1', mu1, 'mu1_hat', mu1_hat, 'mu2', mu2, ...
                'theta_hat', d.theta_hat, 'eta_tilde', d.eta_tilde, ...
                'V', qp1.V, 'V_hat', qp1h.V, ...
                'delta', qp1.delta, ...
                'qp_feasible', qp1.feasible && qp1h.feasible, ...
                'u_box_excess', excess, ...
                'alpha_hat', st.alpha_hat, 'beta_hat', st.beta_hat);
end

end

% ---------------------------------------------------------------------------
function v = solve_decoupling(A, mu, info)
if isfield(info,'rcond') && (~isfinite(info.rcond) || info.rcond < 1e-12)
    v = pinv(A) * mu;
else
    v = A \ mu;
end
end
