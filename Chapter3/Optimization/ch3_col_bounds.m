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

% The torso pitch box is the one entry here that is a DESIGN CHOICE rather
% than a rail, so it comes from p -- see the p.qt_range comment in ch3_params.
% qt is unactuated and appears nowhere else in the constraints except inside
% theta, so this bound is the only thing standing between the optimizer and a
% gait that leans 46 degrees to save torque.
q_lo(3) = p.qt_range(1);
q_hi(3) = p.qt_range(2);
dq_lo = -p.dq_max * ones(nq, 1);
dq_hi =  p.dq_max * ones(nq, 1);

x_lo = [q_lo; dq_lo];
x_hi = [q_hi; dq_hi];

lb = [repmat(x_lo, N, 1); p.T_min; -3*ones(p.ny*p.n_ctrl, 1)];
ub = [repmat(x_hi, N, 1); p.T_max;  3*ones(p.ny*p.n_ctrl, 1)];

% The phase endpoints, when they are decision variables (p.free_theta). The
% two ranges are disjoint, so theta_plus > theta_minus holds at every iterate
% SQP visits (it keeps bounds satisfied) and s(q) is always defined.
if isfield(p, 'free_theta') && ~isempty(p.free_theta) && p.free_theta
    tb = p.theta_bounds;
    if ~(isequal(size(tb), [2 2]) && tb(1,1) < tb(1,2) && tb(2,1) < tb(2,2) ...
            && tb(1,2) < tb(2,1))
        error('ch3_col_bounds:thetaBounds', ...
              ['p.theta_bounds must be [lo hi] for theta_minus over [lo hi] for ' ...
               'theta_plus, with the theta_minus range entirely below the ' ...
               'theta_plus one (got %s).'], mat2str(tb));
    end
    lb = [lb; tb(1,1); tb(2,1)];
    ub = [ub; tb(1,2); tb(2,2)];
end

end
