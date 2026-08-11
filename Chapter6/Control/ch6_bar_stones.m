function [B, geom] = ch6_bar_stones(k, stone, t)
%CH6_BAR_STONES  The two footstep-placement barriers, (6.9) and (6.16).
%
%   [B, geom] = ch6_bar_stones(k, stone, t)
%
% Fig. 6.3, in swing-foot coordinates r = (l_f, h_f) measured from the stance
% foot. Two circles bracket the foot's trajectory so that the step length AT
% IMPACT lands inside the foothold window:
%
%   g_ST1 = R1 + l_max - |r - (-R1, 0)|                      >= 0     inside O1
%   g_ST2 = |r - (l_min/2, -R2)| - sqrt(R2^2 + (l_min/2)^2)  >= 0     outside O2
%
% ---------------------------------------------- why this implies (6.7) at impact
% Set h_f = 0, the definition of impact, and both collapse to statements about
% l_s alone:
%
%   from O1:  |R1 + l_s| <= R1 + l_max                  =>  l_s <= l_max
%   from O2:  (l_s - l_min/2)^2 >= (l_min/2)^2          =>  l_s <= 0 or >= l_min
%
% and the swing foot is in front of the stance foot at strike, so the second
% branch is the live one. That is (6.8). It is CONDITIONAL on l_s > 0 -- the
% barrier alone does not forbid stepping backwards past the stance foot -- and
% ch6_test_barrier states it that way rather than claiming more.
%
% Note what this construction gives for free. Nothing above is a constraint at
% the impact instant; both are constraints on the WHOLE STEP, enforced
% continuously. The impact-time bound is then a corollary of a set the foot
% never leaves, which is the entire reason a barrier is the right tool here:
% an equality constraint at touchdown would be a boundary-value problem needing
% a re-plan every step, and this needs neither.
%
% ------------------------------------------------------ O2 does two jobs (Rmk 6.2)
% O2 passes through (0,0) and (l_min,0) and arches above the ground between
% them, so "outside O2" is the lower step-length bound AND the swing-foot
% clearance. R2 sets the arch height; see ch6_R2_from_clearance. This is why
% Chapter 6 needs no separate scuffing constraint while the Chapter-3
% optimization needed p.limits.sw_clear_min.
%
% ------------------------------------------------------------ moving stones (6.2)
% l_min and l_max come from ch6_stone_level, which supplies their first and
% second time derivatives too. Both circles then move:
%
%   O1: centre fixed, radius rho1 = R1 + l_max(t)
%   O2: centre (l_min(t)/2, -R2) AND radius sqrt(R2^2 + (l_min(t)/2)^2)
%
% so O2 has a moving centre and a changing radius simultaneously. R2 itself is
% held FIXED across the step even when it was derived from l_min: it encodes
% the desired foot clearance, which is a design choice about the swing
% trajectory and has no business chasing a stone. ch6_simulate resolves it once
% per step and stores it in stone.R2.
%
% Inputs
%   k     : swing-foot kinematics from ch6_kin(x, aux, p, 'swing')
%   stone : struct with .l_min0 .l_max0 .motion .R1 .R2 (+ motion parameters)
%   t     : time since the start of the current step [s]
%
% Outputs
%   B    : 1x2 struct array of barrier terms from ch6_bar_lift
%          B(1) = upper bound (O1), B(2) = lower bound + scuffing (O2)
%   geom : struct .l_min .l_max .l_f .h_f .R1 .R2, for plots and reports
%
% See also CH6_BAR_LIFT, CH6_BAR_CIRCLE, CH6_STONE_LEVEL, CH6_R2_FROM_CLEARANCE.

if nargin < 3, t = 0; end

lv = ch6_stone_level(stone, t);

R1 = stone.R1;
R2 = stone.R2;

%% ------------------------------------------------------- O1: upper bound
c1     = [-R1; 0];
rho1   = R1 + lv.l_max;
drho1  = lv.dl_max;
ddrho1 = lv.ddl_max;

G1 = ch6_bar_circle(k.r, c1, rho1, 'inside', [], [], drho1, ddrho1);

%% --------------------------------------- O2: lower bound + scuffing avoidance
% a = l_min/2 is both the centre abscissa and the half-chord, which is what
% makes the circle pass through the two ground points; rho2 follows.
a    = lv.l_min   / 2;
da   = lv.dl_min  / 2;
dda  = lv.ddl_min / 2;

rho2   = sqrt(R2^2 + a^2);
drho2  = a*da / rho2;
ddrho2 = (da^2 + a*dda) / rho2 - (a*da)^2 / rho2^3;

c2   = [a;   -R2];
dc2  = [da;   0];
ddc2 = [dda;  0];

G2 = ch6_bar_circle(k.r, c2, rho2, 'outside', dc2, ddc2, drho2, ddrho2);

%% ------------------------------------------------------------------ lift
B = [ch6_bar_lift(G1, k, 'ls <= lmax'), ...
     ch6_bar_lift(G2, k, 'ls >= lmin')];

geom = struct('l_min', lv.l_min, 'l_max', lv.l_max, ...
              'l_f', k.r(1), 'h_f', k.r(2), 'R1', R1, 'R2', R2);

end
