function [R2, c_used] = ch6_R2_from_clearance(l_min, c_des)
%CH6_R2_FROM_CLEARANCE  Design R2 in Fig. 6.3 from the swing-foot clearance.
%
%   [R2, c_used] = ch6_R2_from_clearance(l_min, c_des)
%
% Remark 6.2: "The choice of R2 in Fig. 6.3 can be designed based on the
% desired mid-step swing foot clearance." This is that design, inverted.
%
% Circle O2 is centred at (l_min/2, -R2) with radius sqrt(R2^2 + a^2),
% a = l_min/2, so it passes through (0,0) and (l_min,0) and bulges upward in
% between. Keeping the swing foot OUTSIDE it therefore does two jobs at once:
% it is the lower step-length bound l_s >= l_min AND a foot-scuffing
% constraint. The bulge height at mid-step, l_f = a, is
%
%       c = sqrt(R2^2 + a^2) - R2      ->      R2 = (a^2 - c^2) / (2c)
%
% so R2 is not a free tuning constant.
%
% ------------------------------------------------------- what limits c, and why
% c -> a gives R2 -> 0: the circle degenerates to the semicircle on the segment
% [0, l_min], whose apex is a = l_min/2. There is no circle through both ground
% points with a higher apex, so a SHORT STONE BUYS LESS CLEARANCE -- a 10 cm
% step cannot be given a 6 cm arch by this construction, whatever R2 is. That is
% a real geometric limit of the barrier in (6.9), not a numerical one, so it is
% clamped and reported (c_used) rather than warned about on every call: on a
% random-terrain run with l_min varying step to step it would fire constantly
% and mean nothing.
%
% The clamp is at 0.9a rather than a because R2 -> 0 makes O2's curvature blow
% up (Hg ~ 1/rho), and a barrier whose second derivative is 1e3 has a gddot the
% torque box cannot supply. 0.9a leaves R2 = 0.105 a, small but finite.
%
% Inputs
%   l_min : lower edge of the foothold window [m]
%   c_des : desired mid-step clearance [m]
%
% Outputs
%   R2     : the radius parameter of circle O2 [m]
%   c_used : the clearance actually achievable, = min(c_des, 0.9 l_min/2)
%
% See also CH6_BAR_STONES, CH6_PARAMS.

a = l_min / 2;

if a <= 0
    error('ch6_R2_from_clearance:l_min', ...
          'l_min must be positive to define circle O2 (got %.4f).', l_min);
end

c_used = min(c_des, 0.9 * a);
c_used = max(c_used, 1e-6);

R2 = (a^2 - c_used^2) / (2 * c_used);

end
