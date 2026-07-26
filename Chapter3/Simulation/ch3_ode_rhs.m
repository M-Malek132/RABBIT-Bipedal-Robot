function xdot = ch3_ode_rhs(~, x, alpha, p)
%CH3_ODE_RHS  Closed-loop swing-phase right-hand side.
%
%   xdot = ch3_ode_rhs(t, x, alpha, p)
%
% Assembles the two halves of the pipeline: stage 1 supplies f and g, the
% selected controller supplies u, and the closed loop is
%
%       xdot = f(x) + g(x) u(x)
%
% Time is passed but unused -- the whole point of the phase variable is that
% the feedback is time-invariant.  It stays in the signature only because
% ode45 requires it.
%
% See also CH3_CONTROL, CH3_STEP.

u        = ch3_control(x, alpha, p);
[f, g]   = ch3_control_affine(x, p);
xdot     = f + g*u;

end
