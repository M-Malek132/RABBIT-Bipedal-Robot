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
% -------------------------------------------------- why there are TWO QP solves
% mu1 is the QP on the REAL eta. mu1_hat is the same QP on the PREDICTOR state
% eta_hat, eq (4.21) / (4.37). Both are needed and they are not interchangeable:
% the predictor has to be driven by its own control for eta_tilde to isolate
% the uncertainty. Substituting mu1 for mu1_hat would fold the reference
% model's own tracking behaviour into the prediction error and the adaptation
% would chase it.
%
% ---------------------------------------------- torque saturation (Sect 4.2.3)
% With constrained = true, mu1 and mu1_hat come from the CLF-QP with the torque
% box (4.36)-(4.37) at p.l1.u_max, via Chapter 3's stage-8 QP.
%
% THE SATURATION IS ON mu1 ONLY. Section 4.2.3 says so explicitly, and it has a
% visible consequence: the realized torque is u_ff + LgLfy^-1 (mu1 + mu2), and
% only the mu1 part was inside the box. The adaptive term can push it back out.
% l1.u_box_excess reports exactly how far, per call, instead of leaving the
% reader to discover it from a plot.
%
% Only the torque box is imposed on the inner QPs -- friction and GRF rows are
% disabled even if p.limits.enable has them on -- because (4.36) specifies the
% torque box and nothing else.
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

clf = ch3_res_clf(p);
st  = ch4_l1_state('unpack', p, xi);

eta = info.eta;

% --- the inner CLF-QPs ----------------------------------------------------
pq = p;
if constrained
    pq.limits.u_max          = p.l1.u_max;
    pq.limits.enable.torque  = true;
    pq.limits.enable.friction = false;      % (4.36) specifies the box only
    pq.limits.enable.grf      = false;
end

% mu1: reference model driven by the REAL transverse state
[mu1, ~, qp1] = ch3_ctrl_clf_qp(Lf2y, LgLfy, u_ff, info, pq, constrained);

% mu1_hat: the same QP driven by the PREDICTOR state. Same Lf2y / LgLfy / u_ff
% -- eq (4.37) uses the same decoupling matrix and feedforward, since those are
% properties of the nominal model at the current x, not of eta.
info_hat     = info;
info_hat.eta = st.eta_hat;
[mu1_hat, ~, qp1h] = ch3_ctrl_clf_qp(Lf2y, LgLfy, u_ff, info_hat, pq, constrained);

% --- applied control ------------------------------------------------------
mu2 = st.mu2;
mu  = mu1 + mu2;

u = u_ff + solve_decoupling(LgLfy, mu, info);

% --- controller state derivative -----------------------------------------
[xidot, d] = ch4_l1_deriv(xi, eta, mu1, mu1_hat, clf, p);

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
