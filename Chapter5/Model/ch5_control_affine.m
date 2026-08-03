function [f, g] = ch5_control_affine(x, p)
%CH5_CONTROL_AFFINE  The control-affine form (3.6) of whichever plant is selected.
%
%   [f, g] = ch5_control_affine(x, p)
%
%       xdot = f(x) + g(x) u
%
% This is the single interface everything above the model layer talks to. The
% two plants of Chapter 5 could hardly be more different -- one linear with one
% input, one nonlinear with two -- and the point of the chapter is that the
% ECBF construction does not care, because it only ever asks for Lie
% derivatives of a scalar constraint along this pair.
%
% Inputs
%   x : nx x 1 state
%   p : parameter struct
%
% Outputs
%   f : nx x 1
%   g : nx x nu
%
% See also CH5_SPRINGMASS, CH5_PENDULUM, CH5_ODE_RHS.

switch lower(p.system)

    case 'springmass'
        sm = ch5_springmass(p);
        f  = sm.A * x(:);
        g  = sm.B;

    case 'pendulum'
        d = ch5_pendulum(x(:), p.plant.pv);
        f = d.f;
        g = d.g;

    otherwise
        error('ch5_control_affine:unknownSystem', ...
              'Unknown p.system "%s".', p.system);
end

end
