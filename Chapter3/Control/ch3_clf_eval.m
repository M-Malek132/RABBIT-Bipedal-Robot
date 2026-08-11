function [V, LfV, LgV] = ch3_clf_eval(eta, clf, eps_clf)
%CH3_CLF_EVAL  Evaluate V_eps and its Lie derivatives along etadot = F eta + G mu.
%
%   [V, LfV, LgV] = ch3_clf_eval(eta, clf, eps_clf)
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
% ------------------------------------------------------------ WHERE THE WORK IS
% None of P_eps, F'P_eps + P_eps F, or 2 P_eps G depends on eta. They depend
% only on P and eps, both fixed for a whole rollout -- so ch3_res_clf builds
% them once and this function is three quadratic forms and nothing else.
% MEASURED: rebuilding them per call cost 9.5 us, which was 31% of a 30.4 us
% stage-7 control call, inside an ODE right-hand side evaluated thousands of
% times per step.
%
% The eps argument is kept in the signature -- rather than read off clf -- so
% that callers can evaluate ONE certificate at a different eps (Chapter 4's L1
% adaptation and ch3_test_control's cross-certificate check both do). When the
% eps handed in is the one the struct was built for, the precomputed block is
% used; otherwise the scaling is rebuilt here for that call alone. So the fast
% path is automatic and the general path still works.
%
% The argument is named eps_clf, not eps: MATLAB's eps is machine epsilon, and
% shadowing it inside a numerical routine is how a tolerance silently becomes a
% control gain.
%
% Inputs
%   eta     : 2ny x 1 transverse variables [y; ydot]
%   clf     : struct from ch3_res_clf
%   eps_clf : the epsilon knob (smaller = faster required convergence)
%
% Outputs
%   V   : scalar V_eps(eta)
%   LfV : scalar
%   LgV : 1 x ny row
%
% See also CH3_RES_CLF, CH3_CTRL_CLF_QP.

if isfield(clf, 'P_eps') && clf.eps == eps_clf
    P_eps   = clf.P_eps;
    M_eps   = clf.M_eps;
    PG2_eps = clf.PG2_eps;
else
    % A certificate evaluated at an eps it was not built for -- rebuild the
    % scaling for this call only, leaving the cached block untouched.
    ny      = size(clf.G, 2);
    I_eps   = blkdiag(eye(ny)/eps_clf, eye(ny));
    P_eps   = I_eps * clf.P * I_eps;
    M_eps   = clf.F.'*P_eps + P_eps*clf.F;
    PG2_eps = 2 * P_eps * clf.G;
end

V = eta.' * P_eps * eta;

if nargout > 1
    LfV = eta.' * M_eps * eta;
end
if nargout > 2
    LgV = eta.' * PG2_eps;
end

end
