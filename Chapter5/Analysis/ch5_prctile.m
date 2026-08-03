function v = ch5_prctile(x, q)
%CH5_PRCTILE  Percentile, without the Statistics Toolbox.
%
%   v = ch5_prctile(x, q)     q in [0, 100]
%
% Linear interpolation between the order statistics at (i - 0.5)/n, which is
% MATLAB's own prctile convention, so results match where both are available.
%
% Written out rather than called because Chapter 5's only hard dependencies are
% the Optimization Toolbox (quadprog) and, for REGENERATING the pendulum
% derivatives only, the Symbolic Math Toolbox -- matching what the repository
% README promises. Reaching for prctile would quietly add a third that is
% needed to print a table.
%
% See also CH5_REPORT, CH5_ROBUST_ULIM.

x = sort(x(isfinite(x(:))));
n = numel(x);

if n == 0, v = NaN; return; end
if n == 1, v = x;   return; end

p = 100 * ((1:n) - 0.5) / n;

v = interp1(p, x, q, 'linear');

% Outside the interpolation range, clamp to the extremes -- interp1 returns NaN
% there, and a NaN in a summary table reads as a failure rather than as "q is
% beyond what n samples can resolve".
if isnan(v)
    if q < p(1), v = x(1); else, v = x(end); end
end

end
