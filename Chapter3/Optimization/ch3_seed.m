function [alpha, x0] = ch3_seed(p, x0)
%CH3_SEED  Build a starting alpha (and start state) for the optimizer.
%
%   [alpha, x0] = ch3_seed(p)
%   [alpha, x0] = ch3_seed(p, x0)
%
% Produces Bezier coefficients that are consistent with a given start state
% rather than merely near it.  Two of the columns are pinned analytically,
% exploiting the Bezier endpoint identities:
%
%   yd(0)     = alpha_0                 ->  alpha_0 = H q0          gives y(0) = 0
%   dyd/ds(0) = M (alpha_1 - alpha_0)   ->  alpha_1 = alpha_0 + (H dq0 / sdot0)/M
%                                                                  gives ydot(0) = 0
%
% so the seed starts exactly ON the zero dynamics surface Z = {eta = 0}, with
% no transient to be absorbed before the gait is even evaluated.  The final
% column is set to the RELABELED start pose -- the guess that after a periodic
% step the legs have swapped into the configuration they began in -- and the
% interior columns interpolate linearly between.
%
% This is a SEED, not a gait: it is not periodic and not optimal. It exists to
% give the collocation solve something feasible-shaped to start from.
%
% Inputs
%   p  : parameter struct
%   x0 : optional 14x1 start state. Default is a hand-tuned post-impact pose:
%        stance foot on the ground, swing foot 0.28 m BEHIND it with clearance,
%        both knees bent, theta = theta_minus so the phase sweeps forward.
%
% Outputs
%   alpha : ny x n_ctrl coefficients
%   x0    : the start state actually used
%
% See also CH3_BEZIER, CH3_COL_SOLVE.

if nargin < 2 || isempty(x0)
    x0 = [ +0.0000; -0.9310; +0.2999; -0.7934; +0.6869; -0.6318; +0.9810; ...
           +0.3952; -0.0419; +0.1847; +0.1847; +0.1045; +0.0000; +0.0000];
end
x0 = x0(:);

nq = p.nq;

% --- put the start state exactly ON the surface theta = theta_minus --------
% theta_minus is a DESIGN CONSTANT of the gait (s must be 0 at the start of a
% step, and the collocation solve enforces theta(x_1) = theta_minus as an
% equality). A hand-tuned pose will sit a hair off it -- the default x0 above
% is at theta = -0.150055 against theta_minus = -0.15 -- which puts s slightly
% negative at t = 0.
%
% theta = qt + q1 + q2/2 is linear in q, so shifting the torso pitch qt by the
% deficit moves theta by exactly that amount. Rotating the torso also moves
% the stance foot, so the foot is re-planted on the ground afterwards.
d_theta = p.theta_minus - p.c_theta(:).' * x0(1:nq);
x0(3)   = x0(3) + d_theta;
foot    = P_st(x0(1:nq));
x0(2)   = x0(2) + foot(2);          % foot_z = -pz - C, so pz += foot_z zeroes it

% --- put the stance foot at the ORIGIN, not merely on the ground ----------
% ch3_col_constraints pins P_st(q_1) = [0;0] to fix the translation gauge (px
% is otherwise free, and leaving it free would make the whole problem
% invariant under a horizontal shift). Planting the foot only at z = 0 leaves
% its x wherever the hand-tuned pose put it, which shows up as a large node-1
% residual that has nothing to do with the gait. px enters the kinematics as a
% pure translation, so this shift changes no angle, velocity or force.
foot  = P_st(x0(1:nq));
x0(1) = x0(1) - foot(1);

q0  = x0(1:nq);
dq0 = x0(nq+1:2*nq);

M_deg  = p.n_ctrl - 1;
n_ctrl = p.n_ctrl;

y_start = p.H * q0;                       % actuated joints now
y_end   = [q0(6); q0(7); q0(4); q0(5)];   % relabeled: legs have swapped

alpha = zeros(p.ny, n_ctrl);
for i = 1:n_ctrl
    lam = (i-1) / (n_ctrl-1);
    alpha(:, i) = (1-lam)*y_start + lam*y_end;
end

% Pin the first column so y(0) = 0 exactly.
alpha(:, 1) = y_start;

% Pin the second column so ydot(0) = 0 exactly, provided the phase is actually
% advancing at t = 0 (sdot0 ~ 0 would mean the gait clock is not running, and
% there is then no slope to match).
[~, ds_dq] = ch3_phase(x0, p);
sdot0 = ds_dq * dq0;
if abs(sdot0) > 1e-9 && M_deg >= 1
    alpha(:, 2) = alpha(:, 1) + (p.H * dq0 / sdot0) / M_deg;
end

% Keep the last column as the relabeled target.
alpha(:, end) = y_end;

% --- mid-step swing-knee flexion, so the foot actually clears the ground ---
% Without this the seed does not walk AT ALL. A monotone interpolation from
% the start pose to the relabeled target sweeps the swing leg forward while
% its foot descends monotonically, so the foot scuffs down onto the stance
% foot: measured L_step = 0.002 m, i.e. a step in place, and the step
% terminates at s = 0.52 having never completed the phase.
%
% The swing knee is the natural clearance mechanism. Its sensitivity is
%     d(swing foot height)/d(swing knee) = sin(qt + q3 + q4)/2 = +0.30
% at the seed pose, so a POSITIVE bump raises the foot. The bump is applied to
% the interior control points only: columns 1 and 2 are pinned by the y(0) = 0
% and ydot(0) = 0 conditions above, and the last column is the target pose, so
% perturbing any of them would undo work already done.
%
% Bezier control points are not interpolated, so the realized lift is smaller
% than the coefficient bump -- roughly 60% of it for a degree-5 curve.
if p.seed_knee_lift ~= 0 && n_ctrl >= 4
    knee_row = 4;                       % outputs are [st hip, st knee, sw hip, sw knee]
    for i = 3:n_ctrl-1
        lam = (i-1)/(n_ctrl-1);
        alpha(knee_row, i) = alpha(knee_row, i) + p.seed_knee_lift * sin(pi*lam);
    end
end

end
