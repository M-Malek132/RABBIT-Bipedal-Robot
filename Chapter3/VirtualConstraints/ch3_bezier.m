function [yd, dyd, d2yd] = ch3_bezier(alpha, s)
%CH3_BEZIER  Bezier (Bernstein) curve and its first two derivatives in s.
%
%   [yd, dyd, d2yd] = ch3_bezier(alpha, s)
%
% Chapter 3 writes the desired output profile as "a vector of Bezier
% polynomials with coefficients alpha":
%
%       yd(s) = sum_{k=0}^{M} alpha_k * B_{k,M}(s),
%       B_{k,M}(s) = nchoosek(M,k) * s^k * (1-s)^(M-k)
%
% Both derivatives are ANALYTIC via the standard difference identities
%
%       d/ds  yd = M       * sum_{k=0}^{M-1} (a_{k+1} - a_k)             B_{k,M-1}
%       d2/ds2 yd = M(M-1) * sum_{k=0}^{M-2} (a_{k+2} - 2a_{k+1} + a_k)  B_{k,M-2}
%
% Analytic (not finite-differenced) matters here: d2yd/ds2 appears directly in
% Lf^2 y, so a noisy second derivative would inject noise straight into the
% feedforward torque and into every constraint gradient fmincon estimates.
%
% Two properties the pipeline uses:
%   * ENDPOINT INTERPOLATION.  yd(0) = alpha_0 and yd(1) = alpha_M exactly, so
%     the boundary values of the profile are decision variables in their own
%     right -- the seed can match a pose exactly with no fitting.
%   * ENDPOINT SLOPE.  dyd/ds(0) = M(alpha_1 - alpha_0), likewise at s = 1, so
%     the periodicity condition on ydot is linear in alpha.
%
% Inputs
%   alpha : ny x (M+1) coefficient matrix (one ROW per output)
%   s     : scalar phase
%
% Outputs
%   yd, dyd, d2yd : ny x 1 value, d/ds, d2/ds2
%
% See also CH3_YD, CH3_OUTPUTS.

n_ctrl = size(alpha, 2);
M      = n_ctrl - 1;

yd = alpha * bern(M, s).';

if nargout > 1
    if M >= 1
        d1  = diff(alpha, 1, 2);                 % ny x M
        dyd = M * (d1 * bern(M-1, s).');
    else
        dyd = zeros(size(alpha,1), 1);
    end
end

if nargout > 2
    if M >= 2
        d2   = diff(alpha, 2, 2);                % ny x (M-1)
        d2yd = M*(M-1) * (d2 * bern(M-2, s).');
    else
        d2yd = zeros(size(alpha,1), 1);
    end
end

end

% ---------------------------------------------------------------------------
function B = bern(M, s)
%BERN  Row vector of the M+1 Bernstein basis polynomials of degree M at s.
%
% Binomials are built by the multiplicative recurrence rather than nchoosek so
% this stays exact and allocation-free for the small degrees used here.
B = zeros(1, M+1);
if M < 0
    B = zeros(1,0);
    return;
end
c = 1;                                   % nchoosek(M,k), updated in place
for k = 0:M
    B(k+1) = c * s^k * (1-s)^(M-k);
    c = c * (M-k) / (k+1);
end
end
