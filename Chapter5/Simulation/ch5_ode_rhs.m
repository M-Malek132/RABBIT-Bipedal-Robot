function dx = ch5_ode_rhs(~, x, u, p)
%CH5_ODE_RHS  The plant, with the control held.
%
%   dx = ch5_ode_rhs(t, x, u, p)
%
%       xdot = f(x) + g(x) u
%
% u is a fixed vector, not a function: Chapter 5 runs sampled-data, so the
% control is computed once per period by ch5_step and held across the interval.
% Nothing in this file solves a QP, which is the point -- see the note in
% ch5_params on why the QP does not belong inside an adaptive ODE solver.
%
% Unlike Chapter 3 and 4 there is no guard and no impact map: both validation
% systems are smooth and continuous for all time, so this is the entire
% simulation model.
%
% See also CH5_STEP, CH5_SIMULATE, CH5_CONTROL_AFFINE.

[f, g] = ch5_control_affine(x, p);
dx = f + g * u(:);

end
