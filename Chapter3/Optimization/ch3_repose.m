function [z_new, info] = ch3_repose(z, p, delta)
%CH3_REPOSE  Pitch the torso by delta, holding both legs fixed in the world.
%
%   [z_new, info] = ch3_repose(z, p, delta)
%
% Rotates the torso forward by delta (rad, positive = forward) at every node
% and moves the Bezier coefficients to match, WITHOUT moving either leg.
% Intended as the warm start for a change of p.qt_range: it lands the gait
% inside the new box while leaving almost every constraint exactly satisfied.
%
% WHY THIS IS NEEDED.  The obvious way to move the torso -- raise the qt_range
% bound and let fmincon clamp the old solution into it -- does not work, and
% the reason is theta.  The phase clock is
%
%       theta = qt + q1 + q2/2
%
% so shifting qt alone shifts theta at every node.  The virtual constraints
% are evaluated at s(theta), so every yd(s) is now read at the wrong phase,
% y = H q - yd(s) is suddenly nonzero, the I/O linearization answers with a
% large corrective torque, and the Hermite-Simpson defects blow up.  Measured
% on the N = 61 upright gait, clamping qt by a mere 0.1 rad took the equality
% residual from 7e-06 to 1.25e+03, and SQP then stalled at ~2e+02 with the
% step length collapsing -- a worse starting point than the seed.
%
% THE FIX.  q1 and q3 are the hip angles RELATIVE to the torso, so the legs'
% absolute orientations are qt+q1 (stance thigh) and qt+q3 (swing thigh).
% Subtracting delta from both hips while adding it to qt therefore rotates the
% torso and leaves both legs exactly where they were:
%
%       qt -> qt + delta      q1 -> q1 - delta      q3 -> q3 - delta
%
% and because yd must follow the joints it drives, the same shift applies to
% the Bezier rows for those two outputs.  A Bernstein basis is a partition of
% unity, so subtracting delta from every control point in a row shifts that
% curve by exactly delta at every s -- no reparametrization, no refit.
%
% WHAT IS EXACTLY PRESERVED.  theta (qt+q1 is unchanged, q2 untouched), hence
% the phase clock and both guard conditions; the absolute angle of every link,
% hence both foot positions, the stance contact, the step length and the hip
% height; y and ydot, hence the zero dynamics surface; and all velocities,
% which are not touched at all.
%
% WHAT IS NOT.  The torso's own orientation is the point of the exercise, so
% its gravitational moment and its inertial coupling to the legs DO change.
% The accelerations therefore change and the Hermite-Simpson defects pick up a
% residual -- this is a warm start, not a symmetry of the dynamics. That
% residual is the honest cost of the move, and it is orders of magnitude
% smaller than the one clamping produces.
%
% Inputs
%   z     : decision vector
%   p     : parameter struct
%   delta : torso pitch shift [rad], positive leans FORWARD
%
% Outputs
%   z_new : reposed decision vector
%   info  : struct .qt_before .qt_after [min max] of torso pitch
%
% See also CH3_POSTURE_MARCH, CH3_COL_PACK, CH3_PARAMS.

% Coordinate layout q = [px pz qt q1 q2 q3 q4]: the same hardcoded indices
% ch3_col_bounds uses for the qt box. Outputs are q(p.iact) = [q1 q2 q3 q4],
% so the two hips are alpha rows 1 and 3.
IQT = 3;  IQ1 = 4;  IQ3 = 6;
ROW_Q1 = 1;  ROW_Q3 = 3;

[X, T, alpha] = ch3_col_unpack(z, p);

info.qt_before = [min(X(IQT,:)), max(X(IQT,:))];

X(IQT,:) = X(IQT,:) + delta;
X(IQ1,:) = X(IQ1,:) - delta;
X(IQ3,:) = X(IQ3,:) - delta;

alpha(ROW_Q1,:) = alpha(ROW_Q1,:) - delta;
alpha(ROW_Q3,:) = alpha(ROW_Q3,:) - delta;

info.qt_after = [min(X(IQT,:)), max(X(IQT,:))];

z_new = ch3_col_pack(X, T, alpha, p);

end
