function [y, ydot, o] = ch3_outputs(x, alpha, p)
%CH3_OUTPUTS  Stage 2: the virtual constraints y = y0(q) - yd(s(q), alpha).
%
%   [y, ydot, o] = ch3_outputs(x, alpha, p)
%
% Chapter 3 designs the gait by picking things to control and driving them to
% a desired profile.  For RABBIT the controlled variables are the four
% actuated joint angles, so y0(q) = H*q, and the profile is indexed by the
% phase s(q) rather than by time:
%
%       y    = H q - yd(s(q), alpha)
%       ydot = H dq - (dyd/ds) sdot         with sdot = (ds/dq) dq
%            = Jy(q) dq,      Jy = H - (dyd/ds)(ds/dq)
%
% Differentiating once more (this is what stage 4 needs) and using the fact
% that ds/dq is CONSTANT -- the payoff of a phase variable linear in q --
%
%       ydd = Jy ddq - (d2yd/ds2) sdot^2
%
% so the only curvature term is the profile's own, and the whole q-dependence
% of ydd sits in Jy.  Substituting ddq = ddq_drift + ddq_in*u makes ydd affine
% in u, which is exactly the structure stages 4, 7 and 8 exploit.
%
% Inputs
%   x     : 14x1 state
%   alpha : ny x n_ctrl coefficients
%   p     : parameter struct
%
% Outputs
%   y    : ny x 1 output error
%   ydot : ny x 1 output error rate
%   o    : struct of intermediates reused downstream --
%             .Jy      ny x nq   output Jacobian, ydot = Jy*dq
%             .curv    ny x 1    (d2yd/ds2)*sdot^2, the drift curvature term
%             .s .sdot .ds_dq .theta
%             .yd .dyd .d2yd
%             .eta     2*ny x 1  [y; ydot], the transverse variables
%
% See also CH3_YD, CH3_PHASE, CH3_IO_LIN.

nq = p.nq;
q  = x(1:nq);
dq = x(nq+1:2*nq);

[s, ds_dq, theta] = ch3_phase(x, p);
sdot = ds_dq * dq;

[yd, dyd, d2yd] = ch3_yd(alpha, s, p);

y    = p.H * q - yd;
Jy   = p.H - dyd * ds_dq;              % ny x nq
ydot = Jy * dq;

if nargout > 2
    o = struct('Jy',   Jy, ...
               'curv', d2yd * sdot^2, ...
               's',    s, 'sdot', sdot, 'ds_dq', ds_dq, 'theta', theta, ...
               'yd',   yd, 'dyd', dyd, 'd2yd', d2yd, ...
               'eta',  [y; ydot]);
end

end
