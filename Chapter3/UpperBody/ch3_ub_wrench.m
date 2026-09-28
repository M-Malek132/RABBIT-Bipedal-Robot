function [Q, W, r_sh, J_sh] = ch3_ub_wrench(q, s, p)
%CH3_UB_WRENCH  Generalized force of the upper-body wrench on RABBIT.
%
%   [Q, W, r_sh, J_sh] = ch3_ub_wrench(q, s, p)
%
% Maps the torso wrench W = [Fx; Fz; M] (see CH3_UB_DEFAULTS) into the 7
% generalized coordinates q = [px pz qt q1 q2 q3 q4] by virtual work:
%
%       Q = J_sh(q)' * [Fx; Fz]  +  e_3 * M
%
%   r_sh(q) = world [x; z] of the shoulder point, Tt(q) * [0; -l_sh; 0; 1]
%             (same torso-axis convention as ch3_body_points' torso_top).
%   J_sh    = d r_sh / d q, 2x7.
%
% J_sh is taken by COMPLEX STEP through the generated Tt.m rather than being
% hand-derived, so it inherits the repository's own sign conventions (pz is
% measured downward in q; world z is up) with no chance of a transcription
% slip, and it is exact to machine precision -- complex step has no
% subtractive cancellation, unlike a finite difference.  Only columns 1..3 can
% be non-zero (Tt depends on px, pz, qt alone); columns 4..7 are set to zero
% exactly rather than computed.
%
% Inputs
%   q : 7x1 configuration
%   s : phase in [0,1] (clamped here: the wrench profile is not extrapolated)
%   p : parameter struct carrying p.ub
%
% Outputs
%   Q    : 7x1 generalized force (zeros when the wrench is inactive)
%   W    : 3x1 wrench [Fx; Fz; M] at this phase
%   r_sh : 2x1 shoulder position, world frame
%   J_sh : 2x7 shoulder Jacobian
%
% See also CH3_UB_DEFAULTS, CH3_CONTROL_AFFINE, CH3_ZD_POINT.

nq = numel(q);
Q  = zeros(nq, 1);
W  = zeros(3, 1);

active = ch3_ub_info(p);
if ~active && nargout < 3
    return;
end

s = min(max(s, 0), 1);

if active
    switch p.ub.mode
        case 'profile'
            W = p.ub.fun(s);
        case 'bezier'
            W = ch3_bezier(p.ub.beta, s);
    end
    W = W(:);
end

l_sh = p.ub.l_sh;
pt   = [0; -l_sh; 0; 1];

T0   = Tt(q);
r4   = T0 * pt;
r_sh = [r4(1); r4(3)];

h    = 1e-20;
J_sh = zeros(2, nq);
for k = 1:3
    qc      = complex(q(:));
    qc(k)   = qc(k) + 1i*h;
    rc      = Tt(qc) * pt;
    J_sh(:, k) = imag([rc(1); rc(3)]) / h;
end

if active
    Q    = J_sh.' * W(1:2);
    Q(3) = Q(3) + W(3);
end
end
