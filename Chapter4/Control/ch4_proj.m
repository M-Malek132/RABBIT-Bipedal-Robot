function pr = ch4_proj(theta, y, theta_max, eps_p)
%CH4_PROJ  The projection operator used by the L1 adaptation laws (4.26).
%
%   pr = ch4_proj(theta, y, theta_max)
%   pr = ch4_proj(theta, y, theta_max, eps_p)
%
% Plain integration thetadot = Gamma*y is unbounded: nothing stops the estimate
% from drifting, and the error bound (4.35) is stated in terms of ||alpha_tilde||
% and ||beta_tilde|| being bounded, so an unbounded estimate does not merely
% look untidy -- it voids the result. The projection operator confines the
% estimate to a ball while preserving the ONE property the Lyapunov argument
% actually uses, eq (4.29):
%
%       (theta - theta*)' (Proj(theta, y) - y)  <=  0
%
% for every theta* inside the ball. In (4.31) that inequality is what cancels
% the cross terms and leaves Vdot_tilde + (c3/eps)Vtilde <= (c3/eps)*delta_V.
%
% CONSTRUCTION.  With the convex indicator
%
%       f(theta) = (||theta||^2 - theta_max^2) / (eps_p * theta_max^2)
%
% f < 0 strictly inside the ball, f = 0 on it, f = 1 on the outer edge of a
% smoothing shell of relative width eps_p. Then
%
%       Proj(theta,y) = y - (grad_f grad_f' / ||grad_f||^2) y f(theta)
%                                       if f > 0 and grad_f'y > 0
%                     = y               otherwise.
%
% In words: only when the estimate is outside the ball AND the update would
% push it further out is the OUTWARD RADIAL COMPONENT removed, and even then it
% is faded in over the shell rather than switched. The fade matters -- a hard
% switch would make thetadot discontinuous, and this sits inside an ODE right
% hand side where that costs the integrator dearly (the same lesson as
% p.control_dt in ch3_params).
%
% The tangential component always survives, which is why the estimate can still
% travel around the boundary and is not stuck once it arrives there.
%
% RESULTING BOUND.  The invariant set is {f <= 1}, so
%
%       ||theta|| <= theta_max * sqrt(1 + eps_p)
%
% which is the ||alpha_tilde|| <= alpha_b, ||beta_tilde|| <= beta_b of (4.32),
% up to the true parameter's own bound. ch4_test_l1 checks both the invariance
% and inequality (4.29) numerically rather than trusting the derivation.
%
% THAT INVARIANCE IS A CONTINUOUS-TIME STATEMENT, AND THE DISTINCTION BITES.
% Outside the ball the projection removes the entire outward radial component,
% so thetadot is exactly TANGENTIAL and d/dt||theta||^2 = 0. But a finite step
% along a tangent lands outside the circle it was tangent to: every explicit
% integrator inflates the radius by O(step^2) per step and will eventually
% spiral out, however correct this operator is. Measured here, Euler at
% Gamma*dt = 0.1 against theta_max = 3 walks the estimate out to ||theta|| = 10.
%
% This is why ch4_l1_advance integrates the controller state with RK4 rather
% than Euler. It is also why the bound should be treated as a guarantee held in
% reserve rather than an operating point -- in the 1.5x mass case alpha_hat
% peaks near 150 against alpha_max = 200 and never engages the boundary at all.
% If you find the estimates riding the bound, raise the bound or lower Gamma;
% do not rely on the projection to hold a hard limit against a coarse step.
%
% Inputs
%   theta     : n x 1 current estimate
%   y         : n x 1 unprojected update direction
%   theta_max : radius of the ball
%   eps_p     : relative width of the smoothing shell, in (0,1]. Default 0.1
%
% Outputs
%   pr : n x 1 projected update
%
% See also CH4_CTRL_L1, CH4_L1_STATE.

if nargin < 4 || isempty(eps_p), eps_p = 0.1; end

pr = y;

if theta_max <= 0 || ~isfinite(theta_max)
    return;                       % no ball to project onto
end

nrm2 = theta.' * theta;
f    = (nrm2 - theta_max^2) / (eps_p * theta_max^2);

if f <= 0
    return;                       % strictly inside: nothing to do
end

grad = 2 * theta / (eps_p * theta_max^2);
gy   = grad.' * y;

if gy <= 0
    return;                       % update already points back inward
end

gg = grad.' * grad;
if gg < 1e-300
    return;                       % theta = 0 cannot be outside a positive ball
end

pr = y - (grad * (gy / gg)) * min(f, 1);

end
