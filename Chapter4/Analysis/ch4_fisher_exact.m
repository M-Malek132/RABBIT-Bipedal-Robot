function p2 = ch4_fisher_exact(a, n1, b, n2)
%CH4_FISHER_EXACT  Two-sided Fisher exact test for two proportions.
%
%   p = ch4_fisher_exact(a, n1, b, n2)     a of n1 against b of n2
%
% The 2x2 table [a, n1-a; b, n2-b] with its margins fixed: the probability of
% every table at least as unlikely as the observed one, under the hypergeometric
% law. It is what "0 of 6 against 4 of 6" needs -- p = 0.061, suggestive and not
% significant -- and what the seed study reports for every pair of variants.
%
% Standard library only: log-factorials via gammaln, and a 1e-7 relative
% tolerance on "as unlikely" (added in log space), the usual guard against
% rounding ties such as the symmetric 0-of-6 / 4-of-6 table.
%
% See also CH4_WILSON, CH4_SEED_STUDY.

N  = n1 + n2;
K  = a + b;                                   % total "successes" (falls)
lo = max(0, K - n2);
hi = min(K, n1);
lp = @(x) lchoose(K, x) + lchoose(N - K, n1 - x) - lchoose(N, n1);
p_obs = lp(a);
p2 = 0;
for x = lo:hi
    lx = lp(x);
    if lx <= p_obs + 1e-7
        p2 = p2 + exp(lx);
    end
end
p2 = min(1, p2);
end

function v = lchoose(n, k)
v = gammaln(n + 1) - gammaln(k + 1) - gammaln(n - k + 1);
end
