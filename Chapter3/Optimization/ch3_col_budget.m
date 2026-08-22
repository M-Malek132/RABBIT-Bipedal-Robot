function p = ch3_col_budget(p, iters, z)
%CH3_COL_BUDGET  Make MaxIterations the binding cap, not MaxFunctionEvaluations.
%
%   p = ch3_col_budget(p)              size the eval cap to p.max_iter
%   p = ch3_col_budget(p, iters)       set p.max_iter = iters first
%   p = ch3_col_budget(p, iters, z)    take the mesh from z, not p.N_nodes
%
% WHY THIS IS NOT OPTIONAL FOR A MARCH.  ch3_col_solve finite-differences
% centrally, so one gradient costs about 2*n_vars evaluations.  Against the
% default p.max_fun_evals = 3e5 that imposes an IMPLICIT iteration limit which
% is well below p.max_iter at any mesh worth solving on:
%
%       N_nodes   n_vars   implied cap   vs p.max_iter = 300
%          21       319        470       ok
%          41       599        250       CAPS
%          61       879        170       CAPS
%          81      1159        129       CAPS
%
% A driver that sets p.max_iter = 300 and stops there therefore gets ~170
% iterations at N = 61, and fmincon returns exitflag 0 -- which reads exactly
% like "MaxIterations reached", i.e. like the posture/speed/ratio step was too
% big, when the truth is that the run was never given the budget it asked for.
% That misdiagnosis costs a re-solve with a smaller step that was never the
% problem.  ch3_col_solve warns (ch3_col_solve:evalCapBinds), but a march
% prints hundreds of lines and one warning is easy to lose.
%
% The 2.2 factor is 2 for the central differences plus 10% headroom for the
% line-search evaluations SQP makes on top of each gradient.
%
% THE EVAL CAP IS SET, NOT RAISED.  If p.max_fun_evals was already larger than
% the computed budget it comes DOWN to it.  That is deliberate: the point is to
% make the two caps agree so that exitflag 0 has one unambiguous meaning, which
% a cap left far above the iteration budget would not achieve.
%
% PASS z WHEN YOU HAVE IT.  n_vars is read from p.N_nodes by default, which is
% correct only while p and the decision vector agree about the mesh -- true for
% a fresh ch3_col_seed and preserved by ch3_col_remesh, but NOT for a warm start
% loaded from a .mat whose gait was refined after that p was written.  Passing z
% takes the mesh from the vector actually being solved, which is the one
% fmincon will differentiate.
%
% Inputs
%   p     : parameter struct
%   iters : optional; sets p.max_iter when non-empty
%   z     : optional decision vector; its length fixes n_vars
%
% See also CH3_COL_SOLVE, CH3_COL_REMESH.

if nargin >= 2 && ~isempty(iters)
    p.max_iter = iters;
end

if nargin >= 3 && ~isempty(z)
    n_vars = numel(z);
else
    n_vars = p.nx * p.N_nodes + 1 + p.ny * p.n_ctrl;
end

p.max_fun_evals = ceil(2.2 * n_vars * p.max_iter);

end
