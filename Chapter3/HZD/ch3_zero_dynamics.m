function Z = ch3_zero_dynamics(alpha, p, n_grid)
%CH3_ZERO_DYNAMICS  The hybrid zero dynamics, and the (5.79)/(5.80) conditions.
%
%   Z = ch3_zero_dynamics(alpha, p)
%   Z = ch3_zero_dynamics(alpha, p, n_grid)
%
% This is the machinery NEC4 and NEC5 of Section 6.3.4 are stated in.  It takes
% the virtual-constraint coefficients alone -- no trajectory, no simulation --
% and returns the closed-form fixed point of the restricted Poincare map
% together with the two conditions for that fixed point to EXIST and to be
% STABLE.
%
% -------------------------------------------------------------------------
% FROM THE RESTRICTED ODE TO (kappa1, kappa2)
%
% ch3_zd_point reduces the dynamics on Z to a scalar second-order ODE
%
%       a(theta) thetaddot + b(theta) thetadot^2 + c(theta) = 0.                (1)
%
% The book writes the zero dynamics in the first-order form
%
%       thetadot = kappa1(theta) sigma,      sigmadot = kappa2(theta),          (2)
%
% where sigma is the momentum conjugate to theta -- for this robot, the angular
% momentum of the whole biped about the stance-foot contact point.  The two
% forms are related by an INTEGRATING FACTOR.  Setting sigma := m(theta)
% thetadot and asking that sigmadot depend on theta alone:
%
%       sigmadot = m' thetadot^2 + m thetaddot
%                = (m' - m b/a) thetadot^2 - m c/a,
%
% so the thetadot^2 term cancels exactly when m' = m b/a, i.e.
%
%       m(theta) = exp( int_{theta_0}^{theta} b/a dxi ),   m(theta_0) = 1,      (3)
%
% and then kappa1 = 1/m and kappa2 = -m c / a.
%
% This is why ch3_zd_point may normalize its projection row w however it likes:
% rescaling w scales a, b and c together, only b/a and c/a reach this file, and
% both are invariant.  The one convention that IS fixed is m(theta_0) = 1, which
% pins the scale of sigma -- and delta_zero below is a RATIO of sigmas across
% the impact, so it is invariant to that too.
%
% -------------------------------------------------------------------------
% THE CONSERVED-ENERGY COORDINATE AND THE POINCARE MAP
%
% With zeta := (1/2) sigma^2, the chain rule along (2) collapses the two states
% to a quadrature in theta alone:
%
%       dzeta/dtheta = sigma sigmadot / thetadot = kappa2/kappa1 = -m^2 c / a
%
%       V_zero(theta) := -int_{theta_0}^{theta} (kappa2/kappa1) dxi
%                      =  int_{theta_0}^{theta} m^2 c / a dxi                   (4)
%
%       zeta(theta) = zeta^+ - V_zero(theta).                                   (5)
%
% The impact contracts sigma by a constant factor, sigma^+ = delta_zero sigma^-,
% because Delta is LINEAR in velocity and both endpoints lie on Z:
%
%       delta_zero = sigma^+ / sigma^-
%                  = [c_theta * Delta_qdot(q_z(theta_f)) * v_z(theta_f)] / m(theta_f),
%
% computed here by running the actual ch3_impact on the unit-rate point of
% S ∩ Z, so the reset map used is the same code the simulator and the
% collocation use rather than a second derivation that can drift out of sync.
%
% Composing (5) with the impact gives a scalar affine Poincare map on zeta
% taken at the END of the step,
%
%       zeta^-_{k+1} = delta_zero^2 zeta^-_k - V_zero(theta_f),                 (6)
%
% whose fixed point and slope are immediate:
%
%       zeta*_2 = -V_zero(theta_f) / (1 - delta_zero^2),   d/dzeta = delta_zero^2.
%
% -------------------------------------------------------------------------
% THE TWO CONDITIONS
%
%   NEC5 / (5.80)   0 < delta_zero^2 < 1.
%       The slope of a scalar affine map IS its eigenvalue, so this is exactly
%       exponential stability of the fixed point -- and it also makes zeta*_2
%       positive, since V_zero(theta_f) < 0 for a gait that gains kinetic
%       energy over the step to pay for the impact loss.
%
%   NEC4 / (5.79)   zeta*_2 > V_zero^MAX / delta_zero^2.
%       EXISTENCE.  Multiply through by delta_zero^2: the left side becomes
%       zeta^+, the kinetic coordinate at the START of the step, and (5) then
%       says zeta(theta) > 0 for every theta in the step.  zeta = sigma^2/2, so
%       this is the statement that the robot never stalls mid-step -- thetadot
%       does not reach zero before the swing foot strikes.  A gait violating it
%       has a perfectly good fixed point of (6) that no actual solution reaches.
%
% NOTE ON "NEC".  The book files these under nonlinear EQUALITY constraints, but
% both are written -- there and here -- as strict INEQUALITIES on the fixed
% point, and that is how ch3_col_constraints imposes them.  NEC1, NEC2 and NEC3
% are the ones that are genuinely equalities or one-sided bounds on the
% trajectory itself.
%
% -------------------------------------------------------------------------
% Inputs
%   alpha  : ny x n_ctrl Bezier coefficients
%   p      : parameter struct (uses theta_minus/plus, c_theta, hzd_grid, tolerances)
%   n_grid : optional override of p.hzd_grid
%
% Output
%   Z : struct with
%         .theta .a .b .c .m .kappa1 .kappa2 .V_zero    (1 x n_grid)
%         .q_z .v_z                                     (nq x n_grid)
%         .V_max      max of V_zero over the step
%         .V_end      V_zero(theta_f)
%         .delta_zero .delta2
%         .zeta_star  the fixed point zeta*_2  (PRE-impact, end of step)
%         .zeta_plus  delta2 * zeta_star, the post-impact value
%         .sigma_star .thetadot_star   the same fixed point in readable units
%         .exists .stable   logicals for NEC4 / NEC5
%         .nec4 .nec5       signed residuals, <= 0 means satisfied
%         .ok               exists && stable
%
% COST.  n_grid points, each about four 7x7 solves plus the mass matrix -- so
% the default 81 is roughly 320 dynamics evaluations.  That is cheap for a
% report and expensive inside fmincon's finite differencing, which is why
% p.limits.enable.hzd defaults to false and the report computes it anyway.
%
% See also CH3_ZD_POINT, CH3_COL_CONSTRAINTS, CH3_POINCARE.

if nargin < 3 || isempty(n_grid)
    if isfield(p, 'hzd_grid'), n_grid = p.hzd_grid; else, n_grid = 81; end
end
n_grid = max(11, round(n_grid));

nq = p.nq;
th = linspace(p.theta_minus, p.theta_plus, n_grid);

a = zeros(1, n_grid);
b = zeros(1, n_grid);
c = zeros(1, n_grid);
q_z = zeros(nq, n_grid);
v_z = zeros(nq, n_grid);
res_foot = 0;

for k = 1:n_grid
    Zp = ch3_zd_point(th(k), alpha, p);
    a(k) = Zp.a;  b(k) = Zp.b;  c(k) = Zp.c;
    q_z(:,k) = Zp.q;  v_z(:,k) = Zp.v;
    res_foot = max(res_foot, Zp.res_foot);
end

% --- integrating factor and the book's (kappa1, kappa2) -------------------
m      = exp(cumtrapz(th, b ./ a));      % m(1) = 1 by construction
kappa1 = 1 ./ m;
kappa2 = -m .* c ./ a;

% --- V_zero, eq. (4) ------------------------------------------------------
V_zero = cumtrapz(th, (m.^2) .* c ./ a);
V_max  = max(V_zero);
V_end  = V_zero(end);

% --- delta_zero from the real reset map ----------------------------------
% Unit rate at the end of the step: thetadot^- = 1, so sigma^- = m(theta_f).
x_end  = [q_z(:,end); v_z(:,end)];
x_post = ch3_impact(x_end, p);
d_theta = p.c_theta(:).' * x_post(nq+1:2*nq);     % thetadot^+ for thetadot^- = 1

delta_zero = d_theta / m(end);
delta2     = delta_zero^2;

% --- fixed point ----------------------------------------------------------
if delta2 < 1 - 1e-12
    zeta_star = -V_end / (1 - delta2);
else
    zeta_star = NaN;          % map is non-contracting; no attracting fixed point
end
zeta_plus = delta2 * zeta_star;

sigma_star = sign(zeta_star) * sqrt(max(zeta_star, 0) * 2);
thetadot_star = sigma_star / m(end);      % the pre-impact rate at the fixed point

% --- NEC5 / (5.80):  0 < delta_zero^2 < 1 --------------------------------
% Written as a signed residual so ch3_col_constraints can use it directly.
% The lower guard is separate: delta2 = 0 would mean the impact annihilates the
% momentum entirely, which makes the map degenerate rather than merely unstable.
tol = hzd_tol(p);
nec5 = max(delta2 - (1 - tol.delta_margin), tol.delta_min - delta2);
stable = (delta2 > tol.delta_min) && (delta2 < 1 - tol.delta_margin);

% --- NEC4 / (5.79):  zeta*_2 > V_zero^MAX / delta_zero^2 -----------------
% When the map is not contracting the fixed point is meaningless, so the
% residual is reported as +Inf-in-spirit (a large positive number) rather than
% NaN: fmincon must see a finite violated constraint it can push against, and a
% NaN would abort the solve instead.
if stable && isfinite(zeta_star)
    nec4 = V_max / delta2 - zeta_star + tol.zeta_margin;
    exists = nec4 <= 0;
else
    nec4 = 1 + abs(V_max) + abs(V_end);
    exists = false;
end

Z = struct('theta', th, 'a', a, 'b', b, 'c', c, ...
           'm', m, 'kappa1', kappa1, 'kappa2', kappa2, ...
           'q_z', q_z, 'v_z', v_z, 'res_foot', res_foot, ...
           'V_zero', V_zero, 'V_max', V_max, 'V_end', V_end, ...
           'delta_zero', delta_zero, 'delta2', delta2, ...
           'zeta_star', zeta_star, 'zeta_plus', zeta_plus, ...
           'sigma_star', sigma_star, 'thetadot_star', thetadot_star, ...
           'exists', exists, 'stable', stable, ...
           'nec4', nec4, 'nec5', nec5, ...
           'ok', exists && stable, 'n_grid', n_grid);

end

% ---------------------------------------------------------------------------
function tol = hzd_tol(p)
%HZD_TOL  Margins for the two strict inequalities.
%
% Both conditions are STRICT in the book, and a numerical optimizer cannot meet
% a strict inequality -- it will happily park on the boundary, where delta2 = 1
% exactly and the fixed point is neutrally stable. The margins turn each into a
% closed constraint with a documented amount of daylight.
tol = struct('delta_margin', 1e-3, ...   % delta2 <= 1 - this
             'delta_min',    1e-8, ...   % delta2 >= this
             'zeta_margin',  1e-6);      % zeta* exceeds V_max/delta2 by this
if isfield(p, 'hzd_tol')
    fn = fieldnames(p.hzd_tol);
    for i = 1:numel(fn), tol.(fn{i}) = p.hzd_tol.(fn{i}); end
end
end
