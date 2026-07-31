function zdot = ch4_ode_rhs(~, z, alpha, p)
%CH4_ODE_RHS  Closed-loop right-hand side: TRUE plant, NOMINAL controller.
%
%   zdot = ch4_ode_rhs(t, z, alpha, p)
%
% The augmented state is
%
%       z = [x; xi]     x  14x1 robot state
%                       xi 5ny x 1 controller state (absent for stateless laws)
%
% and the two halves are deliberately built from DIFFERENT MODELS:
%
%       xdot  = f(x) + g(x) u        <- TRUE model, ch4_control_affine
%       u     = u(x, xi)             <- NOMINAL model, ch4_control
%       xidot = xidot(x, xi)         <- NOMINAL model, ch4_control
%
% This one function is where the entire premise of Chapter 4 lives. Chapter 3's
% ch3_ode_rhs called ch3_control and ch3_control_affine, which are the same
% model twice. Here the plant is p.uncertainty away from what the controller
% believes, and every result in the chapter is a statement about that gap.
%
% Time is passed but unused: the virtual constraints are indexed by phase, and
% the L1 controller is driven by the prediction error rather than by a clock.
% It stays in the signature because ode45 requires it.
%
% See also CH4_CONTROL, CH4_CONTROL_AFFINE, CH4_STEP.

nx = p.nx;

x  = z(1:nx);
xi = z(nx+1:end);

[u, xidot] = ch4_control(x, xi, alpha, p);
[f, g]     = ch4_control_affine(x, p);        % TRUE model (p.uncertainty)

zdot = [f + g*u; xidot];

end
