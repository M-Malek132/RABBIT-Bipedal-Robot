function [lb, ub] = ch3_col_bounds(p, N)
%CH3_COL_BOUNDS  Simple box bounds on the collocation decision vector.
%
%   [lb, ub] = ch3_col_bounds(p, N)
%
% These are sanity rails, not physics -- the physics lives in
% ch3_col_constraints.  Their job is to stop fmincon wandering into poses
% where the generated trig is meaningless (a torso rotated past vertical, a
% knee folded through itself) or velocities large enough that the KKT solve
% loses conditioning, since a single bad evaluation early on can poison the
% whole solve.
%
% See also CH3_COL_CONSTRAINTS, CH3_COL_SOLVE.

nq = p.nq;

q_lo  = [-3.0; -1.30; -1.0; -2.5; -2.5; -2.5; -2.5];
q_hi  = [ 3.0; -0.40;  1.0;  2.5;  2.5;  2.5;  2.5];
dq_lo = -p.dq_max * ones(nq, 1);
dq_hi =  p.dq_max * ones(nq, 1);

x_lo = [q_lo; dq_lo];
x_hi = [q_hi; dq_hi];

lb = [repmat(x_lo, N, 1); p.T_min; -3*ones(p.ny*p.n_ctrl, 1)];
ub = [repmat(x_hi, N, 1); p.T_max;  3*ones(p.ny*p.n_ctrl, 1)];

end
