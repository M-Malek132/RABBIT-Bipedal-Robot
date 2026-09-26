function xi = ch4_l1_advance(xi, smp, clf, p, dt)
%CH4_L1_ADVANCE  Advance the L1 controller state over one control period.
%
%   xi = ch4_l1_advance(xi, smp, clf, p, dt)
%
% Chapter 4 runs sampled-data (p.control_dt > 0), so the controller reads x
% once per period, decides u, and holds it. Its INTERNAL state still has to
% advance over that period, and this is where.
%
% WHY NOT JUST FORWARD EULER.  The adaptation gain is deliberately large --
% Gamma = 1e4 by default, because the error bound (4.35) shrinks like
% 1/||Gamma||. The coupled estimator loop therefore runs at a frequency not far
% below the 1 kHz sample rate. Euler at that ratio does not merely lose
% accuracy, it can add energy and make the estimator diverge -- and the
% divergence would look exactly like "L1 is unstable here" rather than like an
% integration artifact. RK4 over the period costs four evaluations of cheap
% algebra and removes the ambiguity.
%
% WHAT THE PREDICTOR IS FED OVER THE PERIOD depends on p.l1.predictor, and it is
% the whole difference between the two.
%
%   'thesis'  eta, mu1 and mu1_hat are held at their sampled values, while mu2
%             moves with the filter inside the RK4 stages. That is two
%             inconsistencies with the plant at once. The plant received mu2 AT
%             THE SAMPLE and held it; the predictor sees it change. And the
%             plant's eta moves across the period while the predictor's copy
%             of it stands still. The adaptation reads both as model error.
%             Kept verbatim so the Section 4.2 results reproduce.
%
%   'plant'   the predictor gets exactly what the plant got: the applied input
%             mu = mu1 + mu2, held, and eta read at BOTH ends of the period,
%             interpolated linearly between them. This is causal -- the
%             advance runs at t_{k+1}, after eta has been sampled there and
%             before the next control is computed from the new state. Where
%             the output acceleration is constant across the period the
%             transverse velocity is affine in time, so the interpolation
%             reproduces it exactly. On the robot the acceleration drifts
%             within a period even under a held torque, and the error is
%             O(dt^2) rather than zero. The velocity prediction error is the
%             only part of eta_tilde the adaptation reads (ch4_l1_deriv,
%             P = I). Freezing eta instead injects a*ydd*tau into exactly
%             that channel: a bias of order a*ydd*dt/2 on theta_hat whenever
%             the robot accelerates. ch4_test_l1 check 9 measures the
%             difference.
%
% Inputs
%   xi  : 5ny x 1 controller state at the start of the period
%   smp : struct of the period's signals --
%           .eta       2ny x 1 transverse state sampled at the start
%           .eta_next  2ny x 1 transverse state sampled at the end ('plant')
%           .mu        ny x 1  applied virtual input, held ('plant')
%           .mu1_hat   ny x 1  predictor's CLF-QP output, held ('thesis')
%   clf : struct from ch3_res_clf
%   p   : parameter struct
%   dt  : length of the period actually integrated [s] (shorter than
%         p.control_dt on the period the guard cuts off)
%
% Output
%   xi : 5ny x 1 controller state at the end of the period
%
% See also CH4_L1_DERIV, CH4_CTRL_L1, CH4_STEP, CH4_L1_OPTS.

o = ch4_l1_opts(p);

sig0 = struct('eta', smp.eta, 'mu', smp.mu, 'mu1_hat', smp.mu1_hat);

switch o.predictor
    case 'thesis'
        sig_mid = sig0;
        sig_end = sig0;
    case 'plant'
        sig_mid = sig0;  sig_mid.eta = 0.5 * (smp.eta + smp.eta_next);
        sig_end = sig0;  sig_end.eta = smp.eta_next;
end

k1 = ch4_l1_deriv(xi,             sig0,    clf, p);
k2 = ch4_l1_deriv(xi + dt/2 * k1, sig_mid, clf, p);
k3 = ch4_l1_deriv(xi + dt/2 * k2, sig_mid, clf, p);
k4 = ch4_l1_deriv(xi + dt   * k3, sig_end, clf, p);

xi = xi + dt/6 * (k1 + 2*k2 + 2*k3 + k4);

% THE PIECEWISE-CONSTANT LAW (p.l1.adaptation = 'pwc') updates HERE, once per
% period, from the prediction error the period left behind. The RK4 stages
% above held theta_hat (ch4_l1_deriv gives it zero rate), so eta_tilde at the
% end of the period is exactly what the held estimate failed to explain.
if strcmp(o.adaptation, 'pwc')
    xi = pwc_update(xi, smp, clf, p, o, dt);
    return;
end

% THE PROJECTION BALLS HOLD IN CONTINUOUS TIME, NOT AFTER A STEP. Proj removes
% the outward component of the update, which keeps the estimate on the ball
% only in the limit of small steps; a finite step along the tangent lands
% outside it (ch4_test_l1 check 1 measures the drift). Where the loop is soft
% that is negligible. Where it is stiff it is not: with both estimates at
% Gamma = 1e4 and the box scaled for the 1.5x robot, 'l1_con' let theta_hat
% reach 1e7 within 25 steps. So under the plant predictor the estimates are
% returned radially to the ball each period, the discrete-time counterpart of
% the continuous guarantee (4.32) rests on. The thesis path is left verbatim.
if strcmp(o.predictor, 'plant')
    s = ch4_l1_state('unpack', p, xi);
    s.alpha_hat = to_ball(s.alpha_hat, p.l1.alpha_max * sqrt(1 + p.l1.proj_eps));
    s.beta_hat  = to_ball(s.beta_hat,  p.l1.beta_max  * sqrt(1 + p.l1.proj_eps));
    xi = ch4_l1_state('pack', p, s);
end

end

% ---------------------------------------------------------------------------
function xi = pwc_update(xi, smp, clf, p, o, dt)
%PWC_UPDATE  theta_hat(k+1) = -Phi(T)^-1 e^(A_s T) G' eta_tilde(t_{k+1}).
%
% With A_s = -a_s I the error dynamics over one period are
%   eta_tilde_v(t_{k+1}) = e^(-a_s T) eta_tilde_v(t_k) + Phi (theta_hat_k - theta_bar),
%   Phi = (1 - e^(-a_s T)) / a_s,  theta_bar the (weighted) mean uncertainty,
% and choosing theta_hat_{k+1} = -e^(-a_s T)/Phi eta_tilde_v(t_{k+1}) makes the
% next period's error depend on the uncertainty alone -- no memory, no loop.
% At a_s = 0: theta_hat_{k+1} = -eta_tilde_v(t_{k+1})/T = theta_bar over the
% last period, exactly, one sample late.
%
% A period the guard cut to a sliver (under a tenth of the control period)
% carries too little signal for 1/T to amplify, so the estimate is kept.
if dt < 0.1 * p.control_dt
    return;
end
s  = ch4_l1_state('unpack', p, xi);
ev = clf.G.' * (s.eta_hat - smp.eta_next);        % velocity prediction error
a  = o.pwc_rate;
if a > 0
    gain = a / expm1(a * dt);                      % e^(-aT) / Phi(T)
else
    gain = 1 / dt;
end
th = -gain * ev;
if isfinite(o.pwc_max), th = to_ball(th, o.pwc_max); end
s.alpha_hat = zeros(size(s.alpha_hat));
s.beta_hat  = th;
xi = ch4_l1_state('pack', p, s);
end

% ---------------------------------------------------------------------------
function v = to_ball(v, r)
n = norm(v);
if n > r, v = v * (r / n); end
end
