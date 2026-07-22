function [rho, eigvals, A] = poincare_stability(coeffs, x_start, p, h)
%POINCARE_STABILITY  Orbital stability of a walking gait via the step-to-step
% (Poincare) map linearization.
%
%   [rho, eigvals, A] = poincare_stability(coeffs, x_start, p, h)
%
% One walking step defines the Poincare (step-to-step) map
%       x_next = P(x_start) = reset( impact( simulate_one_step(x_start) ) ).
% A periodic gait is a FIXED POINT of P (that is what the periodicity
% equality enforces). The gait is orbitally STABLE iff every eigenvalue of
% the map's Jacobian A = dP/dx at the fixed point lies strictly inside the
% unit circle, i.e. the spectral radius
%       rho = max |eig(A)| < 1.
% This is Westervelt's stability requirement (NEC5): a periodic orbit that is
% not attracting is a limit cycle the robot falls off after a few steps.
%
% A is built by CENTRAL finite differences on the 13 periodic states
% (index 2:14; px is the translation gauge and is excluded). Each column
% costs two step simulations, so a call is ~2*13 = 26 simulations -- this is
% expensive; see hzd_constraints / p.enforce_stability for how it is used.
%
% Inputs
%   coeffs, x_start, p : as elsewhere in the pipeline
%   h                  : finite-difference step (default 1e-5)
% Outputs
%   rho     : spectral radius max|eig(A)|   (< 1 => stable)
%   eigvals : the 13 eigenvalues of A
%   A       : 13x13 Poincare-map Jacobian

    if nargin < 4 || isempty(h), h = 1e-5; end

    idx = 2:14;                 % periodic states (exclude px gauge)
    n   = numel(idx);

    A = zeros(n);
    for j = 1:n
        xp = x_start;  xp(idx(j)) = xp(idx(j)) + h;
        xm = x_start;  xm(idx(j)) = xm(idx(j)) - h;

        xnp = step_map(coeffs, xp, p);
        xnm = step_map(coeffs, xm, p);

        A(:,j) = (xnp(idx) - xnm(idx)) / (2*h);
    end

    eigvals = eig(A);
    rho     = max(abs(eigvals));
end

% ----------------------------------------------------------------------
function x_next = step_map(coeffs, x_start, p)
% One step of the Poincare map. A failed step returns a large state so the
% resulting Jacobian column (and hence rho) is large -> flagged unstable.
    [x_end, ~, ~, status] = simulate_hzd_gait(coeffs, x_start, p);
    if status < 0
        x_next = 1e3 * ones(2*p.nq, 1);
        return;
    end
    x_next = rabbit_reset_map(rabbit_impact_map(x_end));
end
