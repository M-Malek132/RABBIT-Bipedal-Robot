function [xidot, d] = ch4_l1_deriv(xi, eta, mu1, mu1_hat, clf, p)
%CH4_L1_DERIV  Time derivative of the L1 controller state.
%
%   [xidot, d] = ch4_l1_deriv(xi, eta, mu1, mu1_hat, clf, p)
%
% The four coupled pieces of Section 4.2.2, in one place so that the control
% law (ch4_ctrl_l1) and the sampled-data advance (ch4_l1_advance) cannot drift
% apart.
%
% 1. UNCERTAINTY PARAMETRIZATION (4.17), (4.20)
%
%       theta_hat = alpha_hat*||eta|| + beta_hat
%
%    Eq (4.17) is an existence statement, not an assumption: for each instant
%    SOME (alpha, beta) reproduces the true theta in this form. So the estimator
%    is not restricting what it can represent; it is choosing coordinates in
%    which the adaptation is linear and the projection bounds are meaningful.
%
% 2. STATE PREDICTOR (4.19)
%
%       eta_hat_dot = F eta_hat + G mu1_hat + G (mu2 + theta_hat)
%
%    Compare the true transverse dynamics (4.16), eta_dot = F eta + G(mu + theta)
%    with mu = mu1 + mu2. The predictor is the SAME system with the estimated
%    uncertainty in place of the real one and its own CLF-QP output mu1_hat in
%    place of mu1. Subtracting gives (4.24) exactly:
%
%       eta_tilde_dot = F eta_tilde + G mu1_tilde + G(alpha_tilde||eta|| + beta_tilde)
%
%    which is why the prediction error is driven by, and only by, the
%    estimation error. That is the whole reason the predictor exists: eta_tilde
%    is measurable, alpha_tilde and beta_tilde are not.
%
% 3. ADAPTATION LAWS (4.26) with the functions (4.30)
%
%       y_alpha = -G' Peps eta_tilde ||eta||,   y_beta = -G' Peps eta_tilde
%       alpha_hat_dot = Gamma Proj(alpha_hat, y_alpha)
%       beta_hat_dot  = Gamma Proj(beta_hat,  y_beta)
%
%    These are not free choices. Differentiating the composite Lyapunov
%    function (4.27) produces the cross term 2 eta_tilde' Peps G (alpha_tilde
%    ||eta|| + beta_tilde), and (4.30) is exactly the y that cancels it against
%    the parameter terms, leaving (4.31). G'Peps eta_tilde is recovered here as
%    LgV(eta_tilde)'/2 from ch3_clf_eval rather than rebuilding Peps, so there
%    is one definition of Peps in the repository and not two.
%
% 4. LOW-PASS FILTER (4.23)
%
%       mu2_dot = omega_c (-theta_hat - mu2)      i.e.  mu2 = -C(s) theta_hat
%
%    with C(s) = omega_c/(s + omega_c), unit DC gain. In steady state
%    mu2 -> -theta_hat and the uncertainty is cancelled, which is (4.22).
%
%    This filter is the difference between L1 adaptive control and ordinary
%    fast adaptation. Gamma is deliberately huge, so theta_hat is fast and
%    ragged; feeding it straight to the joints would put that content into the
%    ground reaction force and lift the foot. The filter separates HOW FAST WE
%    ESTIMATE from HOW FAST WE ACT, and only the second has to respect the
%    contact.
%
% Inputs
%   xi       : 5ny x 1 controller state (see ch4_l1_state)
%   eta      : 2ny x 1 REAL transverse state
%   mu1      : ny x 1 CLF-QP output for the real system
%   mu1_hat  : ny x 1 CLF-QP output for the predictor
%   clf      : struct from ch3_res_clf
%   p        : parameter struct
%
% Outputs
%   xidot : 5ny x 1
%   d     : struct .theta_hat .eta_tilde .y_alpha .y_beta, for analysis
%
% See also CH4_CTRL_L1, CH4_L1_ADVANCE, CH4_PROJ, CH4_L1_STATE.

ny = p.ny;
s  = ch4_l1_state('unpack', p, xi);

nrm_eta = norm(eta, 2);

% --- 1. estimated uncertainty --------------------------------------------
theta_hat = s.alpha_hat * nrm_eta + s.beta_hat;

% --- 2. state predictor ---------------------------------------------------
eta_hat_dot = clf.F * s.eta_hat + clf.G * mu1_hat + clf.G * (s.mu2 + theta_hat);

% --- 3. adaptation --------------------------------------------------------
eta_tilde = s.eta_hat - eta;

[~, ~, LgV_t] = ch3_clf_eval(eta_tilde, clf, p.eps);
GP_eta_tilde  = LgV_t.' / 2;                 % = G' Peps eta_tilde   (ny x 1)

y_alpha = -GP_eta_tilde * nrm_eta;
y_beta  = -GP_eta_tilde;

alpha_hat_dot = p.l1.Gamma * ch4_proj(s.alpha_hat, y_alpha, ...
                                      p.l1.alpha_max, p.l1.proj_eps);
beta_hat_dot  = p.l1.Gamma * ch4_proj(s.beta_hat,  y_beta,  ...
                                      p.l1.beta_max,  p.l1.proj_eps);

% --- 4. low-pass filter ---------------------------------------------------
mu2_dot = p.l1.omega_c * (-theta_hat - s.mu2);

xidot = [eta_hat_dot; alpha_hat_dot; beta_hat_dot; mu2_dot];

if nargout > 1
    d = struct('theta_hat', theta_hat, 'eta_tilde', eta_tilde, ...
               'y_alpha', y_alpha, 'y_beta', y_beta);
end

end
