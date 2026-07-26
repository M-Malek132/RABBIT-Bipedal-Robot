function [h, dh_dq] = ch3_guard(x, p)
%CH3_GUARD  Stage 1: the switching surface S (swing-foot strike).
%
%   [h, dh_dq] = ch3_guard(x, p)
%
% Chapter 3's hybrid model switches when x reaches
%
%       S = { x : swing-foot height = 0,  swing foot moving DOWN }
%
% h(x) is the swing-foot height in the WORLD frame, z up-positive, ground at
% z = 0.  h > 0 above ground, h < 0 penetrating.  The impact fires on the
% DOWNWARD crossing h: + -> -.
%
% Sign trap worth restating: pz in the generalized coordinates is stored
% DOWN-positive, so world height is -pz.  P_sw already returns world-frame
% up-positive coordinates, so read the height off P_sw and never off pz.
%
% Inputs
%   x : 14x1 state
%   p : parameter struct (uses p.nq)
%
% Outputs
%   h     : scalar swing-foot height [m]
%   dh_dq : 1x7 gradient (row 2 of the swing-foot Jacobian)
%
% See also CH3_IMPACT, CH3_STEP.

q   = x(1:p.nq);
Psw = P_sw(q);
h   = Psw(2);

if nargout > 1
    Jsw   = J_sw(q);
    dh_dq = Jsw(2, :);
end

end
