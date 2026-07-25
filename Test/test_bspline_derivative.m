function test_bspline_derivative()
%TEST_BSPLINE_DERIVATIVE  Validate BSpline_derivative across degrees 2..5.
%
%   test_bspline_derivative
%
% Regression test for a bug that made BSpline_derivative correct at degree 3 --
% the only degree the HZD pipeline calls -- but wrong at degree 5, where the
% recursion's loop range shrank with degree and dropped the top basis
% functions. The error was invisible at degree 3 because the shrinking range
% still happened to reach the knot span there.
%
% Three independent checks:
%
%   1. FINITE DIFFERENCE.  The analytic derivative of a random spline must
%      match a central difference of BSpline's own curve, at every degree
%      2..5 and across s in (0,1). This is the check that would have caught
%      the original bug: it compares the derivative against the VALUE routine,
%      which was never broken.
%
%   2. BERNSTEIN IDENTITY.  A clamped B-spline of degree M with exactly M+1
%      control points IS the degree-M Bezier curve, so
%          dyd/ds = M * sum_k (a_{k+1} - a_k) B_{k,M-1}(s)
%      must hold exactly. This is a fully independent oracle -- a different
%      formula, implemented in a different file -- rather than another way of
%      running the same recurrence.
%
%   3. PARTITION OF UNITY.  The basis functions sum to 1 for all s, so their
%      derivatives must sum to 0. Catches a whole basis function being
%      dropped, which is precisely the original failure mode.
%
% Endpoints s = 0 and s = 1 are included: they used to be handled by a
% one-sided finite difference with step 1e-8 and are now analytic.

fprintf('\n=== test_bspline_derivative ===\n');
rng(11);
pass = true;

s_grid = [0, 0.001, 0.05:0.05:0.95, 0.999, 1];

%% 1. analytic vs central difference of BSpline, degrees 2..5
fprintf('\n  [1] analytic derivative vs central difference of BSpline\n');
for deg = 2:5
    n = 5;                                   % 6 control points
    a = randn(1, n+1);
    worst = 0; worst_s = NaN;
    for s = s_grid
        d_ana = a * BSpline_derivative(n, deg, s).';

        % Step in from the ends so the difference stays inside [0,1].
        h  = 1e-6;
        sm = max(s-h, 0); sp = min(s+h, 1);
        d_fd = (a*BSpline(n,deg,sp).' - a*BSpline(n,deg,sm).') / (sp - sm);

        err = abs(d_ana - d_fd);
        if err > worst, worst = err; worst_s = s; end
    end
    ok = worst < 1e-4;
    pass = pass && ok;
    fprintf('      [%s] degree %d : worst err = %.3e at s = %.3f\n', ...
            tf(ok), deg, worst, worst_s);
end

%% 2. Bezier equivalence: degree M with M+1 control points
fprintf('\n  [2] vs the independent Bernstein derivative identity\n');
for M = 2:5
    n = M;                                   % M+1 control points
    a = randn(1, n+1);
    worst = 0;
    for s = s_grid
        d_bs = a * BSpline_derivative(n, M, s).';
        [~, d_bez] = ch3_bezier(a, s);       % different formula, different file
        worst = max(worst, abs(d_bs - d_bez));
    end
    ok = worst < 1e-9;
    pass = pass && ok;
    fprintf('      [%s] degree %d, %d ctrl pts : worst err = %.3e\n', ...
            tf(ok), M, n+1, worst);
end

%% 3. partition of unity => derivatives sum to zero
fprintf('\n  [3] sum of basis derivatives = 0 (partition of unity)\n');
for deg = 2:5
    n = 5;
    worst_sum = 0; worst_val = 0;
    for s = s_grid
        worst_sum = max(worst_sum, abs(sum(BSpline_derivative(n, deg, s))));
        worst_val = max(worst_val, abs(sum(BSpline(n, deg, s)) - 1));
    end
    ok = worst_sum < 1e-9 && worst_val < 1e-9;
    pass = pass && ok;
    fprintf('      [%s] degree %d : sum(dN) = %.3e , sum(N)-1 = %.3e\n', ...
            tf(ok), deg, worst_sum, worst_val);
end

%% 4. degree 3 is unchanged -- the pipeline's degree must not move
% The HZD pipeline runs bsp_deg = 3 and was already correct there, so the fix
% must be a no-op at that degree. Checked against the finite difference to the
% same tolerance the original satisfied.
fprintf('\n  [4] degree 3 (the degree the HZD pipeline uses) still exact\n');
n = 5; a = randn(1, n+1); worst = 0;
for s = 0.05:0.01:0.95
    d_ana = a * BSpline_derivative(n, 3, s).';
    h = 1e-7;
    d_fd = (a*BSpline(n,3,s+h).' - a*BSpline(n,3,s-h).') / (2*h);
    worst = max(worst, abs(d_ana - d_fd));
end
ok = worst < 1e-6;
pass = pass && ok;
fprintf('      [%s] degree 3 : worst err = %.3e\n', tf(ok), worst);

fprintf('\n--- test_bspline_derivative: %s ---\n\n', tf(pass));

if ~pass
    error('test_bspline_derivative:failed', 'BSpline_derivative validation failed.');
end

end

function s = tf(b)
if b, s = 'PASS'; else, s = 'FAIL'; end
end
