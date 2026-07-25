function [dN] = BSpline_derivative(n, p, u)
% Compute derivative of B-spline basis functions analytically
%
% d/ds N_{i,p}(s) = p * [ N_{i,p-1}(s) / (u_{i+p} - u_i)
%                        - N_{i+1,p-1}(s) / (u_{i+p+1} - u_{i+1}) ]
%
% Inputs:
%   n: number of data points (control points = n+1)
%   p: degree of B-spline
%   u: parameter value [0,1]
%
% Output:
%   dN: 1 x (n+1) vector of basis function derivatives
%
% Knot vector matches BSpline.m exactly (clamped, uniform interior):
%   U_0..U_p = 0,  U_p..U_{n+1} = linspace(0,1,n-p+2),  U_{n+1}..U_m = 1
% with m = n+p+1. The two files MUST agree or the derivative describes a
% different curve than the value.
%
% ---------------------------------------------------------------------------
% FIXED: this routine used to be correct at degree 3 but WRONG at degree 5.
%
% Two defects, both now repaired:
%
% 1. THE RECURSION LOOP RANGE SHRANK WITH DEGREE.  The inner loop ran
%    `for i = 0:(n+2-d)`, so the range of basis-function indices it updated got
%    SMALLER as the recursion climbed in degree -- exactly backwards. The
%    indices it skipped were left at zero.
%
%    For a clamped knot vector the span index satisfies p <= span <= n, and
%    N_{i,d} is nonzero precisely for span-d <= i <= span. So the loop must
%    reach i = span. With n = 5:
%      p = 3: final degree d = 2 covers i = 0..5, and span <= 5 -- just barely
%             covered, which is why degree 3 always looked correct;
%      p = 5: final degree d = 4 covers only i = 0..3, while span = 5 and the
%             nonzero range is i = 1..5. N_{4,4} and N_{5,4} were never
%             written, so the derivative silently lost two terms.
%    Measured before the fix (6 control points, degree 5): the analytic result
%    drifted from a central difference of BSpline's own curve by 1.9e-4 at
%    s = 0.05 rising monotonically to 0.40 at s = 0.95.
%
%    The loop now covers the full index range 0..n+1 at every degree. That is
%    safe for both arrays: nonzero values live at 1 <= i <= n, and every knot
%    lookup stays inside U because i+d+2 <= (n+1)+(p-1)+2 = m+1 = numel(U).
%
% 2. THE ENDPOINTS WERE FINITE-DIFFERENCED.  u = 0 and u = 1 were handled by a
%    one-sided difference with step 1e-8, because the span search
%    `u >= U(i) && u < U(i+1)` finds nothing at u = 1 and the repeated knots of
%    a clamped vector give 0/0. That capped endpoint accuracy at ~1e-8 in a
%    routine advertised as analytic. Clamping the span to [p, n] -- the
%    standard treatment, NURBS Book algorithm A2.1 -- makes the recurrence
%    produce the correct one-sided derivative at both ends with no special
%    case, and the existing zero-denominator guards already absorb the
%    repeated knots.
%
% Validated by TEST_BSPLINE_DERIVATIVE across degrees 2..5, against central
% differences of BSpline and, in the degree-M/M+1-control-point case, against
% the independent Bernstein identity in Chapter3/VirtualConstraints.
% ---------------------------------------------------------------------------

    if p == 0
        dN = zeros(1, n+1);
        return;
    end

    m = n + p + 1;

    U = zeros(1, m+1);
    U(:, p+1:n+2) = linspace(0, 1, m+1-2*p);
    U(:, n+2:m+1) = 1;

    u = min(max(u, U(1)), U(end));

    % --- knot span containing u, clamped to [p, n] -----------------------
    % Clamping is what lets u = 0 and u = 1 go through the same code path as
    % the interior: at u = 1 the half-open test U_i <= u < U_{i+1} matches no
    % span, and below U_p every candidate span is degenerate.
    span = p;
    for i = p:n
        if u >= U(i+1) && u < U(i+2)
            span = i;
            break;
        end
        if i == n
            span = n;                 % u at (or past) the right end
        end
    end

    % --- degree 0 ---------------------------------------------------------
    N_lower = zeros(1, n+2);          % holds N_{i,d} for i = 0..n+1
    N_lower(span+1) = 1;

    % --- recurse up to degree p-1 ----------------------------------------
    for d = 1:p-1
        N_temp = zeros(1, n+2);
        for i = 0:(n+1)               % full range: see note 1 above
            idx = i + 1;

            term1 = 0; term2 = 0;

            if N_lower(idx) ~= 0
                denom1 = U(i+d+1) - U(i+1);          % U_{i+d} - U_i
                if denom1 ~= 0
                    term1 = ((u - U(i+1)) / denom1) * N_lower(idx);
                end
            end

            if idx+1 <= numel(N_lower) && N_lower(idx+1) ~= 0
                denom2 = U(i+d+2) - U(i+2);          % U_{i+d+1} - U_{i+1}
                if denom2 ~= 0
                    term2 = ((U(i+d+2) - u) / denom2) * N_lower(idx+1);
                end
            end

            N_temp(idx) = term1 + term2;
        end
        N_lower = N_temp;
    end

    % --- derivative formula ----------------------------------------------
    % A zero denominator means a repeated knot, where that term drops out.
    dN = zeros(1, n+1);
    for i = 0:n
        idx = i + 1;
        term1 = 0; term2 = 0;

        if N_lower(idx) ~= 0
            denom1 = U(i+p+1) - U(i+1);              % U_{i+p} - U_i
            if denom1 ~= 0
                term1 = p * N_lower(idx) / denom1;
            end
        end

        if idx+1 <= numel(N_lower) && N_lower(idx+1) ~= 0
            denom2 = U(i+p+2) - U(i+2);              % U_{i+p+1} - U_{i+1}
            if denom2 ~= 0
                term2 = p * N_lower(idx+1) / denom2;
            end
        end

        dN(idx) = term1 - term2;
    end

end
