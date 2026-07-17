%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% B-SPLINE EVAL
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [b, db] = bspline_eval(c, s, p)
    n = p.n_coeffs - 1;
    s = max(0, min(1, s));

    try
        N  = BSpline(n, p.bs_degree, s);
        dN = BSpline_derivative(n, p.bs_degree, s);
    catch
        N  = zeros(n+1, 1);
        dN = zeros(n+1, 1);
    end

    b  = c(:).' * N(:);
    db = c(:).' * dN(:);
end