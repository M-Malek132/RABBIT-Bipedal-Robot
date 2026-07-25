function [x_plus, impulse] = ch3_impact(x_minus, p)
%CH3_IMPACT  Stage 1: the reset map Delta.   x+ = Delta(x-)  for x- in S.
%
%   [x_plus, impulse] = ch3_impact(x_minus, p)
%
% Delta is the composition of three things, applied in THIS ORDER:
%
%   1. RIGID PLASTIC IMPACT.  The swing foot strikes and sticks (no rebound,
%      no slip), imposing J_sw*dq+ = 0.  Configuration q is CONTINUOUS -- only
%      velocities jump -- so the momentum balance is
%
%          [ M   -J_sw' ] [ dq+     ]   [ M dq- ]
%          [ J_sw   0   ] [ impulse ] = [   0   ]
%
%      Both M and J_sw are evaluated at the PRE-IMPACT configuration.
%
%   2. RELABEL.  The landed foot becomes stance, so leg indices swap.
%
%   3. RE-PLANT.  pz is shifted so the NEW stance foot sits exactly at z = 0,
%      absorbing the small residual foot height left by the event tolerance.
%
% ORDER MATTERS.  The impact uses the PRE-relabel swing Jacobian J_sw.
% Relabel first and the impulse is applied to the wrong foot.
%
% px is deliberately NOT reset -- it is the translation gauge and must
% accumulate so the robot advances in the world.  Every periodicity check
% therefore excludes index 1.
%
% Inputs
%   x_minus : 14x1 pre-impact state
%   p       : parameter struct (uses p.nq)
%
% Outputs
%   x_plus  : 14x1 start-of-next-step state
%   impulse : 2x1 contact impulse [Ix; Iz] at the impacting foot [Ns], world
%             frame, z up-positive.  This is the Table 3.1 impulse quantity.
%
% See also CH3_RELABEL, CH3_GUARD.

nq = p.nq;
q  = x_minus(1:nq);
dq = x_minus(nq+1:2*nq);

% --- 1. plastic impact ----------------------------------------------------
M_mat = M(q);
Jsw   = J_sw(q);                 % PRE-impact swing-foot Jacobian
nc    = size(Jsw, 1);

A = [M_mat, -Jsw.'; ...
     Jsw,    zeros(nc)];
b = [M_mat * dq; ...
     zeros(nc, 1)];

sol     = A \ b;
dq_plus = sol(1:nq);
impulse = sol(nq+1:end);

x_post = [q; dq_plus];           % q unchanged through impact

% --- 2. relabel -----------------------------------------------------------
x_rel = ch3_relabel(x_post, p);

% --- 3. re-plant the new stance foot on the ground ------------------------
% foot_z = -pz - (joint terms), so driving foot_z to zero needs
%   pz_new = pz + foot_z.
foot   = P_st(x_rel(1:nq));
x_plus = x_rel;
x_plus(2) = x_plus(2) + foot(2);

end
