function [yd, dyd, d2yd] = ch3_yd(alpha, s, p)
%CH3_YD  Desired output profile yd(s, alpha) and its first two s-derivatives.
%
%   [yd, dyd, d2yd] = ch3_yd(alpha, s, p)
%
% Thin dispatcher over the curve basis, so that every other file in the
% pipeline is basis-agnostic and p.basis is the single switch:
%
%   'bezier'   Chapter-3 form. Analytic derivatives (ch3_bezier).
%   'bspline'  Clamped B-spline through the repo's existing BSpline.m /
%              BSpline_derivative.m, kept as an independent cross-check.
%              Its second derivative is central-differenced, because
%              BSpline_derivative offers no analytic second derivative.
%
% These two bases are not merely similar -- a clamped B-spline of degree M
% with exactly M+1 control points IS the degree-M Bezier curve.  Setting
% p.basis='bspline' with p.bsp_deg = size(alpha,2)-1 therefore reproduces the
% Bezier VALUES to machine precision, and ch3_test_vc asserts that.
%
% HISTORY.  BSpline_derivative.m used to be correct at degree 3 -- the degree
% the existing pipeline runs -- but wrong at degree 5, disagreeing with a
% finite difference of its own curve by up to ~0.4 near s = 1, because its
% recursion loop range shrank with degree and dropped the top basis functions.
% That is fixed (see Test/test_bspline_derivative.m), so both bases now agree
% on values AND derivatives at every degree, which ch3_test_vc asserts.
%
% 'bezier' remains the default anyway: its derivatives are closed-form rather
% than recursive, and it supplies an ANALYTIC second derivative. The B-spline
% path still central-differences d2yd/ds2, and Lf^2 y depends on it directly.
%
% Inputs
%   alpha : ny x n_ctrl coefficient matrix (one ROW per output)
%   s     : scalar phase in [0,1]
%   p     : parameter struct (uses p.basis, p.bsp_deg)
%
% Outputs
%   yd, dyd, d2yd : ny x 1
%
% See also CH3_BEZIER, CH3_OUTPUTS, CH3_PHASE.

switch lower(p.basis)

    case 'bezier'
        switch nargout
            case {0,1}, yd = ch3_bezier(alpha, s);
            case 2,     [yd, dyd] = ch3_bezier(alpha, s);
            otherwise,  [yd, dyd, d2yd] = ch3_bezier(alpha, s);
        end

    case 'bspline'
        n   = size(alpha, 2) - 1;      % BSpline.m's "n": n+1 control points
        deg = p.bsp_deg;
        ny  = size(alpha, 1);

        % A clamped knot vector defines the curve only on [0,1], so s HAS to be
        % clamped -- unlike the Bezier branch, which is a polynomial and needs
        % no clamp at all.  But once s is clamped the curve being reported is
        % CONSTANT, so the derivatives reported alongside it must be zero to
        % match.  Returning the frozen boundary slope instead (as this used to)
        % hands the caller a dyd that is not the derivative of the yd it came
        % with, and ch3_outputs then builds Jy = H - dyd*ds_dq from it and adds
        % a term that should not be there.  ch3_phase zeroes ds_dq under its own
        % clamp for exactly this reason.
        clamped = (s < 0) || (s > 1);
        s       = min(max(s, 0), 1);

        yd = alpha * BSpline(n, deg, s).';

        if nargout > 1
            if clamped
                dyd = zeros(ny, 1);
            else
                dyd = alpha * BSpline_derivative(n, deg, s).';
            end
        end
        if nargout > 2
            if clamped
                d2yd = zeros(ny, 1);
            else
                d2yd = d2_bspline(alpha, n, deg, s);
            end
        end

    otherwise
        error('ch3_yd:basis', 'Unknown basis "%s" (expected bezier|bspline).', p.basis);
end

end

% ---------------------------------------------------------------------------
function d2 = d2_bspline(alpha, n, deg, s)
%D2_BSPLINE  d2yd/ds2 by differencing the analytic first derivative.
%
% BSpline_derivative offers no analytic second derivative, so this differences
% it.  The stencil is kept ENTIRELY inside [0,1] rather than clipped against
% the endpoint: clipping silently turns the central difference into a
% half-width one-sided formula, which drops it from O(h^2) to O(h) and costs
% five orders of accuracy at s = 1 -- measured 1.76e-04 there against 4e-09 in
% the interior.  That is touchdown on every step, and it is where |d2yd/ds2| is
% largest, so it feeds straight into Lf^2 y at the worst moment.  Within h of
% an end, a second-order ONE-SIDED formula is used instead.

h  = 1e-5;
d1 = @(t) alpha * BSpline_derivative(n, deg, t).';

if s - h >= 0 && s + h <= 1
    d2 = (d1(s+h) - d1(s-h)) / (2*h);                    % central,  O(h^2)
elseif s - h < 0
    d2 = (-3*d1(s) + 4*d1(s+h) - d1(s+2*h)) / (2*h);     % forward,  O(h^2)
else
    d2 = ( 3*d1(s) - 4*d1(s-h) + d1(s-2*h)) / (2*h);     % backward, O(h^2)
end

end
