function [V, LfV, LgV] = ch5_clf_eval(eta, clf, eps_scale)
%CH5_CLF_EVAL  Evaluate V and its Lie derivatives along etadot = F eta + G mu.
%
%   [V, LfV, LgV] = ch5_clf_eval(eta, clf, eps_scale)
%
% With P_eps = I_eps P I_eps and, at relative degree r,
%
%       I_eps = blkdiag( I/eps^(r-1), I/eps^(r-2), ..., I/eps, I )
%
% (each block ny x ny, so at r = 2 this is ch3's blkdiag(I/eps, I)):
%
%       V    = eta' P_eps eta
%       Vdot = eta' (F' P_eps + P_eps F) eta  +  2 eta' P_eps G mu
%              \__________ LfV ____________/     \____ LgV ____/
%
% Vdot IS AFFINE IN mu at every relative degree -- that fact does not weaken as
% r grows, and it is the reason the stability requirement stays a single linear
% inequality in (5.31) rather than turning into something a QP cannot carry.
% The barrier row of Section 5.2 is built to have exactly this property, which
% is what lets the two sit in the same program.
%
% Inputs
%   eta       : (r*ny) x 1 transverse variables
%   clf       : struct from ch5_res_clf
%   eps_scale : the eps knob (1 = plain CLF)
%
% Outputs
%   V   : scalar
%   LfV : scalar
%   LgV : 1 x ny row
%
% See also CH5_RES_CLF, CH5_CTRL_CLF_QP.

r  = clf.r;
ny = clf.ny;

if eps_scale == 1
    P_eps = clf.P;
else
    scales = eps_scale .^ -(r-1:-1:0);          % 1/eps^(r-1) ... 1/eps^0
    I_eps  = kron(diag(scales), eye(ny));
    P_eps  = I_eps * clf.P * I_eps;
end

eta = eta(:);

V   = eta.' * P_eps * eta;
LfV = eta.' * (clf.F.'*P_eps + P_eps*clf.F) * eta;
LgV = 2 * eta.' * P_eps * clf.G;

end
