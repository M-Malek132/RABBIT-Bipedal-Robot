function [B, geom] = ch6_bar_width(k, wp)
%CH6_BAR_WIDTH  The step-width barriers of Section 6.3.2, (6.19)-(6.20).
%
%   [B, geom] = ch6_bar_width(k, wp)
%
% Fig. 6.12, a TOP view: the swing foot lives in the (l_f, w_f) plane and the
% step width at touchdown must land in [w_min, w_max]. Two circles do it, by
% exactly the containment argument of Section 6.3.2:
%
%   O3  tangent to  w = w_max  from below, through the initial foot position,
%       foot INSIDE   =>  w_f <= w_max  for the whole step
%   O4  tangent to  w = w_min  from above, through the initial foot position,
%       foot INSIDE   =>  w_f >= w_min  for the whole step
%
% ------------------------------------------------------------ the radius, (6.19)
% Tangency to w = w_max from below puts the centre at w_c = w_max - R3. Passing
% through the initial point (l_f0, w0), at horizontal offset l3 from the centre:
%
%       l3^2 + (w0 - w_c)^2 = R3^2
%       l3^2 + (R3 - (w_max - w0))^2 = R3^2
%   =>  R3 = ( l3^2 + (w_max - w0)^2 ) / ( 2 (w_max - w0) )              (6.19)
%
% and the same algebra with (w0 - w_min) gives R4. l3 is the only free
% parameter and it is a CONSERVATISM dial: R grows like l3^2, so a large l3
% gives a nearly-flat arc that constrains w_f and barely touches l_f, while a
% small l3 gives a tight circle that couples the two. Case 3 of Section 6.3.3
% varies step length and width together, so l3 must be large enough that O3 and
% O4 are not also fighting over l_f.
%
% ==================================================================== ON O4's SENSE
% THE THESIS WRITES THE LOWER BOUND AS "O4F >= R4", i.e. the foot OUTSIDE the
% disc. That cannot bound w_f from below, and the counterexample is immediate.
% Take the outside condition at the tangency abscissa l_f = l0':
%
%       (w_f - (w_min + R4))^2 >= R4^2   =>   w_f >= w_min + 2 R4   or
%                                             w_f <= w_min
%
% -- the second branch is the WRONG SIDE of the line, and it is the branch a
% foot drifting downward is on. More generally, the complement of a finite disc
% contains points arbitrarily far below w = w_min in every direction, so no
% "outside a disc" condition can imply a lower bound on w_f.
%
% What DOES imply it is the same principle O3 uses -- containment. A disc
% tangent to w = w_min from above lies entirely in {w >= w_min}, so being inside
% it gives the bound directly, and the R4 of (6.19) is exactly the radius that
% makes such a disc pass through the initial foot position. So "the same
% principle can be applied" is right and only the inequality sense is not.
%
% It is implemented as containment. wp.o4_sense = 'outside' selects the literal
% reading so the failure is reproducible rather than asserted, and
% ch6_test_barrier checks BOTH directions: containment implies w_f >= w_min on
% sampled points, and the outside form admits points with w_f < w_min.
% =============================================================================
%
% ---------------------------------------------------- what this is validated on
% Section 6.3 runs on DURUS, a 23-DoF 3D humanoid. There is no DURUS model in
% this repository and fabricating one would make every number in Section 6.3.3
% unfalsifiable. What IS checked here is what does not need the robot: the
% geometry above (exactly, by sampling), the derivative chain (against finite
% differences), and the closed loop on ch6_foot3d, a 3-DoF double-integrator
% swing foot -- which is the subsystem the constraint is actually written in.
% See ch6_foot3d and docs/CH6_STEPPING_STONES.md for what that does and does
% not establish.
%
% Inputs
%   k  : kinematics struct in the (l_f, w_f) plane -- same fields as ch6_kin
%   wp : p.width struct -- .l3 .w_min .w_max .w0, optionally .l_f0 (default 0),
%        .o4_sense ('inside', default | 'outside', the literal thesis form) and
%        .dw_floor (default 0.01 m; see the note on (6.19) below)
%
% Outputs
%   B    : 1x2 struct array, B(1) = w <= w_max (O3), B(2) = w >= w_min (O4)
%   geom : struct .R3 .R4 .c3 .c4 .l_f .w_f .w_min .w_max .recovering
%
% See also CH6_BAR_CIRCLE, CH6_BAR_LIFT, CH6_FOOT3D.

l_f0 = 0;
if isfield(wp, 'l_f0') && ~isempty(wp.l_f0), l_f0 = wp.l_f0; end

o4_sense = 'inside';
if isfield(wp, 'o4_sense') && ~isempty(wp.o4_sense), o4_sense = wp.o4_sense; end

% THE FLOOR IS A DIVISION GUARD, NOT A DESIGN PARAMETER, so it is microns. A
% floor big enough to "look safe" is worse than none: at dw_floor = 0.01 a
% perfectly legal gap of 0.008 gets replaced by 0.010, the circle no longer
% passes through the start point, and g(0) comes out at -2e-3 -- the barrier
% reports a violated initial condition for a case that was fine. Measured, and
% it looked like a controller failure.
%
% A genuinely small POSITIVE gap needs no floor at all: R = (l3^2 + dw^2)/(2 dw)
% simply grows, and a huge circle tangent to the boundary IS the half-plane
% locally. That is the least conservative correct answer, and the construction
% degenerates to it gracefully.
dw_floor = 1e-6;
if isfield(wp, 'dw_floor') && ~isempty(wp.dw_floor), dw_floor = wp.dw_floor; end

dw_hi = wp.w_max - wp.w0;
dw_lo = wp.w0    - wp.w_min;

% ================================== (6.19) IS UNDEFINED WHEN w0 LEAVES THE WINDOW
% R = (l3^2 + dw^2) / (2 dw) divides by the gap between the previous step width
% and the boundary, so it needs w_min < w0 < w_max. That is not a numerical
% edge case, it is what the construction MEANS: O3 is "tangent to the maximum
% boundary AND containing the initial swing foot position", and there is no
% circle tangent to w = w_max from below through a point that is already above
% it.
%
% So the containment argument of Section 6.3.2 only covers step-width changes
% small enough that the NEW window still contains the PREVIOUS width. Section
% 6.3.3 draws w_d randomly from [12, 33] cm with +-2.5 cm windows, and
% consecutive draws routinely violate that -- so the geometry as written cannot
% cover its own numerical study, and this is worth knowing rather than dividing
% by a small number and carrying on.
%
% WHAT HAPPENS INSTEAD. The gap is floored at dw_floor and the circle is built
% anyway. It is then tangent to the boundary but does NOT contain the starting
% point, so g starts NEGATIVE and the barrier is being asked to recover rather
% than to maintain. That is a legitimate thing to ask of an ECBF and not of a
% reciprocal CBF -- Chapter 5 makes exactly this point: a polynomial condition
% on the derivatives of h stays well defined at h < 0, while B = 1/h does not
% (ch5_barrier, .rec_defined). The forward-invariance guarantee does not apply
% until g comes back up, and geom.recovering says so on every call.
recovering = (dw_hi <= 0) || (dw_lo <= 0);

dw_hi = max(dw_hi, dw_floor);
dw_lo = max(dw_lo, dw_floor);

%% ------------------------------------------------------ O3: w_f <= w_max
R3 = (wp.l3^2 + dw_hi^2) / (2 * dw_hi);
c3 = [l_f0 + wp.l3; wp.w_max - R3];
G3 = ch6_bar_circle(k.r, c3, R3, 'inside');

%% ------------------------------------------------------ O4: w_f >= w_min
R4 = (wp.l3^2 + dw_lo^2) / (2 * dw_lo);
c4 = [l_f0 + wp.l3; wp.w_min + R4];
G4 = ch6_bar_circle(k.r, c4, R4, o4_sense);

B = [ch6_bar_lift(G3, k, 'ws <= wmax'), ...
     ch6_bar_lift(G4, k, 'ws >= wmin')];

geom = struct('R3', R3, 'R4', R4, 'c3', c3, 'c4', c4, ...
              'l_f', k.r(1), 'w_f', k.r(2), ...
              'w_min', wp.w_min, 'w_max', wp.w_max, 'w0', wp.w0, ...
              'o4_sense', o4_sense, 'recovering', recovering);

end
