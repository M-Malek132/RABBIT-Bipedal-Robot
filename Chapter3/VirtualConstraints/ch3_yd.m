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
        s   = min(max(s, 0), 1);

        yd = alpha * BSpline(n, deg, s).';

        if nargout > 1
            dyd = alpha * BSpline_derivative(n, deg, s).';
        end
        if nargout > 2
            % Central difference of the analytic first derivative, stepped in
            % from the ends so the clamped knot vector's repeated knots are
            % never straddled.
            h  = 1e-5;
            sm = min(max(s-h, 0), 1);
            sp = min(max(s+h, 0), 1);
            if sp == sm
                d2yd = zeros(size(alpha,1), 1);
            else
                d2yd = (alpha * BSpline_derivative(n, deg, sp).' - ...
                        alpha * BSpline_derivative(n, deg, sm).') / (sp - sm);
            end
        end

    otherwise
        error('ch3_yd:basis', 'Unknown basis "%s" (expected bezier|bspline).', p.basis);
end

end
