function [z, exitflag] = ch5_solve_qp(H, f, Aineq, bineq, lb, ub, z0)
%CH5_SOLVE_QP  quadprog with the options built once.
%
%   [z, exitflag] = ch5_solve_qp(H, f, Aineq, bineq, lb, ub, z0)
%
% Every QP in this chapter is tiny -- at most four variables and a handful of
% rows -- and is solved tens of thousands of times per study. At that size the
% dominant cost is not the factorization, it is optimoptions() rebuilding an
% options object on every call, so it is built once and cached. Same reasoning
% as the qp_options helper in ch3_ctrl_clf_qp.
%
% ------------------------------------------------- why 'active-set', not the default
% Chapters 3 and 4 use 'interior-point-convex'. It is the wrong choice here,
% and the difference is not marginal.
%
% When the barrier row is pushing hard, the ECBF constraint and the CLF
% objective disagree by a wide margin -- the CLF alone would use ||mu|| ~ 8
% while the barrier row demands ~93. Measured over a pendulum run,
% interior-point-convex failed on 222 of 30001 samples with exitflag -3
% ("unbounded") on a problem whose Hessian is positive definite and which
% therefore cannot be unbounded. Re-solving those same 222 QPs:
%
%       interior-point-convex     0 / 222
%       active-set              222 / 222,  worst residual 1.4e-11
%
% 'active-set' is also the better structural fit: these problems have at most
% four variables and four rows, which is exactly the dense small-scale regime
% it is written for, and unlike interior-point it USES the starting point --
% which is why ch5_min_norm_mu's warm start is worth passing.
%
% H is symmetrized here rather than at each call site: quadprog warns about
% asymmetry, and a warning printed 30,000 times inside a control loop drowns
% out anything worth reading.
%
% See also CH5_MIN_NORM_MU, CH5_SCALE_ROW, CH5_CTRL_ECBF_CLF_QP.

persistent opts
if isempty(opts)
    opts = optimoptions('quadprog', 'Display', 'off', ...
                        'Algorithm', 'active-set');
end

H = (H + H.') / 2;

[z, ~, exitflag] = quadprog(H, f, Aineq, bineq, [], [], lb, ub, z0, opts);

end
