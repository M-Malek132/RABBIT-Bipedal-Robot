function [V, LfV, LgV] = ch3_clf_eval(eta, clf, eps)
%CH3_CLF_EVAL  Evaluate V_eps and its Lie derivatives along etadot = F eta + G mu.
%
%   [V, LfV, LgV] = ch3_clf_eval(eta, clf, eps)
%
% With P_eps = I_eps P I_eps and I_eps = blkdiag((1/eps) I, I),
%
%       V_eps = eta' P_eps eta
%       Vdot  = eta' (F' P_eps + P_eps F) eta  +  2 eta' P_eps G mu
%               \_________ LfV ____________/     \____ LgV ____/
%
% The decisive structural fact for stages 7 and 8: Vdot is AFFINE IN mu, so
% the rapid-exponential-stability requirement
%
%       LfV + LgV mu + (c3/eps) V <= 0
%
% is a single LINEAR INEQUALITY in mu.  That is what makes the CLF condition
% something a quadratic program can carry as a constraint.
%
% Inputs
%   eta : 2ny x 1 transverse variables [y; ydot]
%   clf : struct from ch3_res_clf
%   eps : the epsilon knob (smaller = faster required convergence)
%
% Outputs
%   V   : scalar V_eps(eta)
%   LfV : scalar
%   LgV : 1 x ny row
%
% See also CH3_RES_CLF, CH3_CTRL_CLF_QP.

ny = size(clf.G, 2);

I_eps = blkdiag(eye(ny)/eps, eye(ny));
P_eps = I_eps * clf.P * I_eps;

V   = eta.' * P_eps * eta;
LfV = eta.' * (clf.F.'*P_eps + P_eps*clf.F) * eta;
LgV = 2 * eta.' * P_eps * clf.G;

end
