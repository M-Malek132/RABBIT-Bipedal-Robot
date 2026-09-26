function [lo, hi, phat] = ch4_wilson(k, n, conf)
%CH4_WILSON  Wilson score interval for a binomial proportion.
%
%   [lo, hi]       = ch4_wilson(k, n)          95% interval for k of n
%   [lo, hi, phat] = ch4_wilson(k, n, conf)
%
% The fall rates in Chapter 4 are counts over a handful of runs (0 of 6, 4 of
% 6), where the textbook Wald interval p +/- z sqrt(p(1-p)/n) collapses to a
% point at 0 and 1. The Wilson interval does not, and it is what the report
% quotes: 0 of 6 is [0, 0.39], 4 of 6 is [0.30, 0.90].
%
% Standard library only (no Statistics Toolbox): z from erfcinv.
%
% See also CH4_FISHER_EXACT, CH4_SEED_STUDY.

if nargin < 3 || isempty(conf), conf = 0.95; end
if n <= 0
    lo = NaN; hi = NaN; phat = NaN;
    return;
end
z    = sqrt(2) * erfcinv(1 - conf);          % 1.959964 at 95%
phat = k / n;
den  = 1 + z^2 / n;
ctr  = (phat + z^2 / (2*n)) / den;
half = z * sqrt(phat * (1 - phat) / n + z^2 / (4*n^2)) / den;
lo   = max(0, ctr - half);
hi   = min(1, ctr + half);
end
