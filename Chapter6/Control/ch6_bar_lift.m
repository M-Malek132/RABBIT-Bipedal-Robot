function B = ch6_bar_lift(G, k, label)
%CH6_BAR_LIFT  Lift a geometric constraint g(r,t) to the barrier row's terms.
%
%   B = ch6_bar_lift(G, k, label)
%
% G says what the constraint is in the PLANE (value, gradient, curvature, time
% partials -- from ch6_bar_circle or ch6_bar_affine). k says how the plane point
% moves with the ROBOT (from ch6_kin). This is the only place the two meet, and
% it produces the four numbers every Chapter-6 barrier row is built from:
%
%       g, gdot, and   gddot = Lf2g + LgLfg u.
%
% ------------------------------------------------------------------ the chain
% With r = r(q(t)) and g = G(r, t),
%
%   gdot  = dG/dt + dG/dr rdot
%   gddot = d2G/dt2 + 2 (d2G/dtdr) rdot + rdot' (d2G/dr2) rdot + dG/dr rddot
%
% and rddot = Lf2r + LgLfr u from ch6_kin, so
%
%   Lf2g  = d2G/dt2 + 2 (d2G/dtdr) rdot + rdot' (d2G/dr2) rdot + dG/dr Lf2r
%   LgLfg = dG/dr LgLfr                                       (1 x nu)
%
% THE FACTOR OF 2 is the pair of mixed partials d/dt(dG/dr) and d/dr(dG/dt),
% which are equal, not one term. Dropping it costs nothing on a static stone
% (both are zero) and silently halves the feed-forward on a moving one, which
% is precisely the case Section 6.2 exists to handle -- so it is the kind of
% error that passes every test written against Section 6.1.
%
% ---------------------------------------------------- relative degree, checked
% LgLfg = dG/dr * (Jr ddq_in) is the first place u appears, so rb = 2 -- IF it
% is nonzero. It can vanish: dG/dr is a direction in the (l_f, h_f) plane, and
% at a pose where the swing-foot Jacobian cannot accelerate the foot along that
% direction with any torque, the row loses its grip on u. B.controllable
% reports it rather than letting the QP silently drop a constraint it appears
% to be enforcing. This is the same honesty ch5_ctrl_ecbf_clf_qp applies to its
% Lb, and for the same reason: the guarantee genuinely has a gap there.
%
% Inputs
%   G     : from ch6_bar_circle / ch6_bar_affine
%   k     : from ch6_kin
%   label : short string used in reports and plots
%
% Output
%   B : struct
%         .g .gdot           the barrier and its rate
%         .eta_b   2x1       [g; gdot]                            (5.10), rb=2
%         .Lf2g              scalar, the u-free part of gddot
%         .LgLfg   1 x nu    the coefficient of u in gddot
%         .controllable      false when LgLfg is numerically zero
%         .singular          G was evaluated at a circle centre
%         .label
%
% See also CH6_KIN, CH6_CBF_ROW, CH6_BARRIER.

if nargin < 3, label = ''; end

rdot = k.rdot;

B = struct();
B.label = label;
B.g     = G.g;
B.gdot  = G.gt + G.dg * rdot;
B.eta_b = [B.g; B.gdot];

B.Lf2g  = G.gtt ...
        + 2 * (G.gtr * rdot) ...
        + rdot.' * G.Hg * rdot ...
        + G.dg * k.Lf2r;

B.LgLfg = G.dg * k.LgLfr;              % 1 x nu

B.controllable = norm(B.LgLfg, inf) > 1e-10;
B.singular     = G.singular;

end
