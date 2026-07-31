function xi = ch4_l1_advance(xi, eta, mu1, mu1_hat, clf, p, dt)
%CH4_L1_ADVANCE  Advance the L1 controller state over one control period.
%
%   xi = ch4_l1_advance(xi, eta, mu1, mu1_hat, clf, p, dt)
%
% Chapter 4 runs sampled-data (p.control_dt > 0), so the controller reads x
% once per period, decides u, and holds it. Its INTERNAL state still has to
% advance over that period, and this is where.
%
% WHY NOT JUST FORWARD EULER.  The adaptation gain is deliberately large --
% Gamma = 1e4 by default, because the error bound (4.35) shrinks like
% 1/||Gamma||. The coupled loop (eta_tilde -> alpha_hat, beta_hat -> mu2 ->
% eta_tilde) therefore runs at a natural frequency of order sqrt(Gamma*||Peps||
% *||eta||), which at the defaults is comparable to the 1 kHz sample rate
% itself. Euler at that ratio does not merely lose accuracy, it can add energy
% and make the estimator diverge -- and the divergence would look exactly like
% "L1 is unstable here" rather than like an integration artifact. RK4 over the
% period costs four evaluations of cheap algebra and removes the ambiguity.
%
% WHAT IS FROZEN, AND WHY IT IS ALLOWED.  Over the period we hold
%
%   eta       the plant state the controller sampled -- this is what a digital
%             controller genuinely does, not an approximation of it;
%   mu1       already committed: it is inside the u being held;
%   mu1_hat   a forcing term in the predictor equation.
%
% Freezing mu1_hat is the one real approximation: eta_hat moves within the
% period, so the predictor's own QP output would move a little too. It is
% first order in dt and, more to the point, mu1_hat does not participate in the
% stiff loop above -- it drives the predictor, it is not driven by it. So RK4
% still resolves the part that needed resolving. The alternative, a QP solve
% per RK4 stage, would quadruple the dominant cost of the whole simulation to
% refine a term that is already O(dt) correct.
%
% Inputs
%   xi              : 5ny x 1 controller state at the start of the period
%   eta             : 2ny x 1 sampled transverse state
%   mu1, mu1_hat    : ny x 1 held CLF-QP outputs
%   clf             : struct from ch3_res_clf
%   p               : parameter struct
%   dt              : control period [s]
%
% Output
%   xi : 5ny x 1 controller state at the end of the period
%
% See also CH4_L1_DERIV, CH4_CTRL_L1, CH4_STEP.

k1 = ch4_l1_deriv(xi,            eta, mu1, mu1_hat, clf, p);
k2 = ch4_l1_deriv(xi + dt/2*k1,  eta, mu1, mu1_hat, clf, p);
k3 = ch4_l1_deriv(xi + dt/2*k2,  eta, mu1, mu1_hat, clf, p);
k4 = ch4_l1_deriv(xi + dt*k3,    eta, mu1, mu1_hat, clf, p);

xi = xi + dt/6 * (k1 + 2*k2 + 2*k3 + k4);

end
