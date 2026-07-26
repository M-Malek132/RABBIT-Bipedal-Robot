function [xdot, u, aux] = ch3_col_dynamics(x, alpha, p)
%CH3_COL_DYNAMICS  Closed-loop RHS used inside the collocation constraints.
%
%   [xdot, u, aux] = ch3_col_dynamics(x, alpha, p)
%
% Inside the collocation solve the input is the FEEDFORWARD alone,
%
%       u = u_ff(x, alpha) = -(LgLf y)^-1 Lf^2 y,
%
% i.e. exactly the torque that holds ydd = 0.  There is no feedback term, and
% that is deliberate: the solve constrains the trajectory to start ON the zero
% dynamics surface Z = {eta = 0}, and u_ff renders Z INVARIANT, so eta stays
% zero along the whole step and any PD or CLF term would be multiplying zero.
% Adding feedback here would contribute nothing but extra nonlinearity for
% fmincon to differentiate through.
%
% This is also why the stage-3 cost is the honest one: the torque being
% integrated is the torque the gait actually requires, not a tracking
% transient.
%
% One ch3_io_lin call supplies both the input and the dynamics split, so a
% node costs a single KKT factorization.
%
% See also CH3_IO_LIN, CH3_COL_CONSTRAINTS, CH3_COL_COST.

[~, ~, u, info] = ch3_io_lin(x, alpha, p);

nq   = p.nq;
dq   = x(nq+1:2*nq);
ddq  = info.aux.ddq_drift + info.aux.ddq_in * u;
xdot = [dq; ddq];

if nargout > 2
    aux = info.aux;
end

end
