function [xidot, d] = ch4_l1_deriv(xi, sig, clf, p)
%CH4_L1_DERIV  Time derivative of the L1 controller state.
%
%   [xidot, d] = ch4_l1_deriv(xi, sig, clf, p)
%
% The four coupled pieces of Section 4.2.2, in one place so that the control
% law (ch4_ctrl_l1) and the sampled-data advance (ch4_l1_advance) cannot drift
% apart. p.l1.predictor selects between two state predictors (piece 2); the
% other three pieces are shared.
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
% 2a. THE THESIS PREDICTOR (4.19), p.l1.predictor = 'thesis'
%
%       eta_hat_dot = F eta_hat + G mu1_hat + G (mu2 + theta_hat)
%
%    with mu1_hat the CLF-QP solved on eta_hat. Subtracting the true dynamics
%    (4.16), eta_dot = F eta + G(mu1 + mu2 + theta), gives (4.24) exactly:
%
%       eta_tilde_dot = F eta_tilde + G mu1_tilde + G(alpha_tilde||eta|| + beta_tilde)
%
%    THE STEP FROM (4.24) TO (4.28) DOES NOT FOLLOW. (4.28) needs
%    F eta_tilde + G mu1_tilde to decrease Peps along eta_tilde. What the two
%    CLF-QPs guarantee is a decrease along eta and along eta_hat separately, and
%    the min-norm law is nonlinear, so mu1(eta_hat) - mu1(eta) is not the law
%    evaluated at eta_tilde. The difference forces eta_tilde exactly as a model
%    error would, and the adaptation absorbs it. Measured on posture_195 at
%    eps = 0.35: theta_hat overshot the true theta 3-5x under perturbation, and
%    matched its size only when a faster reference model (eps = 0.25) shrank
%    the mismatch.
%
% 2b. THE PLANT-INPUT PREDICTOR, p.l1.predictor = 'plant'
%
%       eta_hat_dot = F eta + G (mu + theta_hat) - a (eta_hat - eta)
%
%    with mu = mu1 + mu2 the input the plant actually received and
%    a = p.l1.predictor_rate. This is the standard L1 state predictor. Against
%    the true dynamics the F eta and G mu terms cancel identically, leaving
%
%       eta_tilde_dot = -a eta_tilde + G (alpha_tilde||eta|| + beta_tilde)
%
%    a stable linear system driven by the estimation error and by nothing else.
%    The reference model is no longer inside the predictor at all: it acts
%    through mu1 on the real robot, and the predictor only has to explain what
%    the plant did with the input it was given.
%
% 3. ADAPTATION LAWS (4.26)
%
%       alpha_hat_dot = Gamma_alpha Proj(alpha_hat, y_beta ||eta||)
%       beta_hat_dot  = Gamma       Proj(beta_hat,  y_beta)
%
%    y is not a free choice: it is the one that cancels the cross term
%    2 eta_tilde' P G (alpha_tilde||eta|| + beta_tilde) in the composite
%    Lyapunov function (4.27), with P the matrix that certifies the predictor's
%    error dynamics:
%
%       thesis:  y_beta = -G' Peps eta_tilde    (4.30), Peps from the CLF --
%                recovered as LgV(eta_tilde)'/2 from ch3_clf_eval, so there is
%                one definition of Peps in the repository, not two
%       plant:   y_beta = -G' eta_tilde         P = I, since -aI is its own
%                Lyapunov certificate: (-aI)'I + I(-aI) = -2aI
%
%    Under 'plant' this is the velocity prediction error, and the beta channel
%    closes a second-order loop s^2 + a s + Gamma: natural frequency
%    sqrt(Gamma), damping a/(2 sqrt(Gamma)).
%
%    ALPHA HAS ITS OWN GAIN, p.l1.Gamma_alpha (empty: Gamma). Its regressor is
%    ||eta||, so its loop gain is Gamma_alpha*||eta||^2: the alpha channel
%    stiffens exactly when tracking degrades, which is when it can least afford
%    to. Gamma_alpha = 0 freezes alpha_hat at its initial zero and leaves beta
%    to carry theta, which (4.17) says it can.
%
%    THE PARAMETRIZATION (4.17) IS REDUNDANT, and the adaptation cannot see the
%    redundancy: while ||eta|| stays near a value m, every pair with
%    alpha m + beta = theta estimates theta equally well, so nothing holds the
%    split. Over long runs alpha_hat and beta_hat drift into near-opposite
%    directions (median cosine -0.9 at 0.7x and 1.5x) and alpha_hat presses its
%    projection bound. That drift is real but it is a symptom, not what makes
%    L1 fall: leaking alpha_hat back to zero (2-50 /s), or keeping theta_hat
%    continuous across footstrikes, made the long runs fall sooner.
%
%    AN OPTIONAL LEAK, p.l1.alpha_leak = lambda [1/s], a sigma-modification
%    aimed at exactly that unconstrained direction:
%
%       alpha_hat_dot = Gamma_alpha/m^2 Proj(alpha_hat, y_alpha) - lambda alpha_hat
%
%    It pulls alpha_hat back to zero along the line the data cannot see, so
%    beta_hat carries the steady part of theta, and it points inward, so the
%    projection ball stays invariant. The price is the usual one: in (4.31) the
%    leak adds -2 lambda alpha_tilde' Gamma^-1 alpha_hat, which bounds
%    alpha_tilde but no longer lets it converge, so (4.35) gains a term in
%    lambda ||alpha||^2. Off by default (ch4_params has why).
%
%    AND AN OPTIONAL CAP ON ITS REGRESSOR, p.l1.alpha_regressor_rate = kappa:
%    in pieces 1 and 3, ||eta|| is replaced by phi = min(||eta||, phi_max), the
%    same phi in both, so (4.24)-(4.31) go through unchanged and (4.17) still
%    holds as an existence statement (any bounded phi admits some alpha, beta).
%
%    What does make L1 fall is the estimator loop's speed against the sample
%    rate. Under 'plant' the coupled prediction error and estimates oscillate at
%    about sqrt(Gamma + Gamma_alpha phi^2) rad/s, so alpha's share grows with
%    the tracking error. After a bad footstrike ||eta|| reaches 12-14: about
%    4100 rad/s at the defaults, 4.1 rad per 1 ms sample, past the RK4 advance's
%    stability limit (about 2.8) and past what a 1 kHz loop can realize. In the
%    three long-run L1 falls re-simulated sample by sample, theta_hat reached
%    ~2700-2900 within 1-3 samples of such an impact, against a true theta of
%    30-300, and the step collapsed. At a 0.5 ms control period, with no cap,
%    most of the long-run falls do not happen. ch4_l1_opts turns kappa into
%    phi_max = sqrt((kappa/dt)^2 - Gamma) / sqrt(Gamma_alpha); capping removes
%    most falls at 1 kHz too, at a price in tracking (see ch4_params).
%
%    OR NORMALIZE BOTH LAWS, p.l1.normalized_rate = kappa_n (the default, at
%    0.75 rad per sample -- ch4_params has the measurements): divide them by
%
%       m^2 = max(1, (Gamma + Gamma_alpha phi^2) / (kappa_n/dt)^2)
%
%    The loop then runs at no more than kappa_n/dt at any tracking error, the
%    law is untouched wherever it was already slower, and theta_hat keeps
%    alpha_hat*||eta|| in full -- which also means that, unlike the cap, it
%    does nothing about the jump alpha_hat*(change in ||eta||) a footstrike
%    puts straight into theta_hat. Under 'plant', the cancellation behind
%    (4.27) survives with the prediction error weighted by 1/m^2:
%    V = eta_tilde'eta_tilde/m^2 + alpha_tilde'alpha_tilde/Gamma_alpha
%    + beta_tilde'beta_tilde/Gamma has
%    Vdot <= -(2a + d(ln m^2)/dt) ||eta_tilde||^2 / m^2, so the bound on the
%    prediction error loosens by the factor m, and holds only while m^2 decays
%    slower than exp(-2at), which a footstrike's jump in ||eta|| can break.
%
%    OR REPLACE THE LAW, p.l1.adaptation = 'pwc': the piecewise-constant
%    adaptation of the sampled-data L1 literature, which has no estimator
%    loop at all (see ch4_params and ch4_l1_advance, where the update lives).
%    Here it only holds the estimate over the period.
%
% 4. LOW-PASS FILTER (4.23)
%
%       mu2_dot = omega_c (-theta_hat - mu2)      i.e.  mu2 = -C(s) theta_hat
%
%    with C(s) = omega_c/(s + omega_c), unit DC gain. In steady state
%    mu2 -> -theta_hat and the uncertainty is cancelled, which is (4.22).
%
%    This filter is the difference between L1 adaptive control and ordinary
%    fast adaptation. Gamma is deliberately large, so theta_hat is fast and
%    ragged; feeding it straight to the joints would put that content into the
%    ground reaction force and lift the foot. The filter separates HOW FAST WE
%    ESTIMATE from HOW FAST WE ACT, and only the second has to respect the
%    contact.
%
% Inputs
%   xi  : 5ny x 1 controller state (see ch4_l1_state)
%   sig : struct of what the controller sees at this instant --
%           .eta      2ny x 1 REAL transverse state
%           .mu       ny x 1  input the plant received, mu1 + mu2 ('plant')
%           .mu1_hat  ny x 1  CLF-QP output for the predictor ('thesis')
%   clf : struct from ch3_res_clf
%   p   : parameter struct
%
% Outputs
%   xidot : 5ny x 1
%   d     : struct .theta_hat .eta_tilde .y_alpha .y_beta .m2, for analysis
%
% See also CH4_CTRL_L1, CH4_L1_ADVANCE, CH4_L1_OPTS, CH4_PROJ, CH4_L1_STATE.

o  = ch4_l1_opts(p);
s  = ch4_l1_state('unpack', p, xi);

eta     = sig.eta;
nrm_eta = min(norm(eta, 2), o.phi_max);     % alpha's regressor, capped

% --- the piecewise-constant law (p.l1.adaptation = 'pwc') -----------------
% No adaptation dynamics: theta_hat is held in the beta_hat slot for the whole
% period (alpha_hat stays 0) and is reset at the END of each period by
% ch4_l1_advance. Inside the period only the predictor and the filter move,
% the predictor with its own rate a_s = p.l1.pwc_rate, the one the update
% inverts:  eta_tilde_dot = -a_s eta_tilde + G (theta_hat - theta).
if strcmp(o.adaptation, 'pwc')
    theta_hat   = s.beta_hat;
    eta_tilde   = s.eta_hat - eta;
    eta_hat_dot = clf.F * eta + clf.G * (sig.mu + theta_hat) - o.pwc_rate * eta_tilde;
    mu2_dot     = p.l1.omega_c * (-theta_hat - s.mu2);
    ny = numel(s.beta_hat);
    xidot = [eta_hat_dot; zeros(ny, 1); zeros(ny, 1); mu2_dot];
    if nargout > 1
        d = struct('theta_hat', theta_hat, 'eta_tilde', eta_tilde, ...
                   'y_alpha', zeros(ny, 1), 'y_beta', zeros(ny, 1), 'm2', 1);
    end
    return;
end

% --- 1. estimated uncertainty --------------------------------------------
theta_hat = s.alpha_hat * nrm_eta + s.beta_hat;
eta_tilde = s.eta_hat - eta;

% --- 2. state predictor, and the P that certifies its error ---------------
switch o.predictor
    case 'thesis'
        eta_hat_dot = clf.F * s.eta_hat + clf.G * sig.mu1_hat ...
                      + clf.G * (s.mu2 + theta_hat);
        [~, ~, LgV_t] = ch3_clf_eval(eta_tilde, clf, p.eps);
        GP_eta_tilde  = LgV_t.' / 2;          % = G' Peps eta_tilde
    case 'plant'
        eta_hat_dot = clf.F * eta + clf.G * (sig.mu + theta_hat) ...
                      - o.predictor_rate * eta_tilde;
        GP_eta_tilde  = clf.G.' * eta_tilde;  % = G' eta_tilde   (P = I)
end

% --- 3. adaptation --------------------------------------------------------
y_beta  = -GP_eta_tilde;
y_alpha =  y_beta * nrm_eta;

% normalization: exactly 1 when off or when the loop is under its ceiling
m2 = max(1, (o.Gamma + o.Gamma_alpha * nrm_eta^2) / o.loop_gain_max);

% the leak is a rate on alpha_hat itself, outside the normalization: it acts
% where the data exerts no force, so there is no loop gain for m^2 to limit
alpha_hat_dot = o.Gamma_alpha / m2 * ch4_proj(s.alpha_hat, y_alpha, ...
                                              p.l1.alpha_max, p.l1.proj_eps) ...
                - o.alpha_leak * s.alpha_hat;
beta_hat_dot  = o.Gamma / m2       * ch4_proj(s.beta_hat,  y_beta,  ...
                                              p.l1.beta_max,  p.l1.proj_eps);

% --- 4. low-pass filter ---------------------------------------------------
mu2_dot = p.l1.omega_c * (-theta_hat - s.mu2);

xidot = [eta_hat_dot; alpha_hat_dot; beta_hat_dot; mu2_dot];

if nargout > 1
    d = struct('theta_hat', theta_hat, 'eta_tilde', eta_tilde, ...
               'y_alpha', y_alpha, 'y_beta', y_beta, 'm2', m2);
end

end
