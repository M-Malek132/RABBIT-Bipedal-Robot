function [rho, A, info] = ch3_poincare(x_star, alpha, p)
%CH3_POINCARE  Orbital stability of a periodic gait: spectral radius of dP/dx.
%
%   [rho, A, info] = ch3_poincare(x_star, alpha, p)
%
% PERIODIC IS NOT STABLE, AND THE DIFFERENCE IS THE WHOLE BALL GAME.
% The stage-3 periodicity equality makes x_star a FIXED POINT of the
% step-to-step (Poincare) map
%
%       P : x  ->  Delta( flow(x) up to the guard )
%
% but a fixed point can be repelling. A gait can satisfy periodicity to 1e-8
% and still fall over in three steps, because each step's small error is
% multiplied by the map's Jacobian. Only rho < 1 means the orbit ATTRACTS,
% i.e. that the robot recovers from the disturbance every real step supplies.
%
% The Jacobian is built by central-differencing P about x_star in the 13
% directions that matter -- px is excluded, since it is the translation gauge
% and advances by design, and including it would contribute a guaranteed unit
% eigenvalue that says nothing about stability.
%
% COST: 26 step simulations. That is why this is a post-hoc DIAGNOSTIC rather
% than a constraint inside the solve -- as an inequality, fmincon would finite
% difference it and pay 26 sims per variable per gradient.
%
% Inputs
%   x_star : 14x1 fixed point (start-of-step state)
%   alpha  : ny x n_ctrl coefficients
%   p      : parameter struct
%
% Outputs
%   rho  : spectral radius. rho < 1 => stable/walkable, rho >= 1 => falls.
%   A    : 13x13 Poincare Jacobian
%   info : struct .eigs .fixed_point_residual .ok
%
% See also CH3_SIMULATE, CH3_REPORT.

idx = 2:p.nx;                 % everything except px
n   = numel(idx);
h   = 1e-5;

% residual of the fixed-point condition itself
s0 = ch3_step(x_star, alpha, p);
fp_res = norm(s0.x_next(idx) - x_star(idx), inf);
ok = s0.ok;

A = zeros(n, n);
for j = 1:n
    xp = x_star; xp(idx(j)) = xp(idx(j)) + h;
    xm = x_star; xm(idx(j)) = xm(idx(j)) - h;

    sp = ch3_step(xp, alpha, p);
    sm = ch3_step(xm, alpha, p);

    if ~sp.ok || ~sm.ok
        ok = false;
        A(:, j) = NaN;
        continue;
    end

    A(:, j) = (sp.x_next(idx) - sm.x_next(idx)) / (2*h);
end

if any(~isfinite(A(:)))
    rho = Inf;
    ev  = NaN;
else
    ev  = eig(A);
    rho = max(abs(ev));
end

info = struct('eigs', ev, 'fixed_point_residual', fp_res, 'ok', ok);

end
