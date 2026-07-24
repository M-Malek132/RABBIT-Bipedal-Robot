%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% B-SPLINE EVAL
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [b, db, ddb] = bspline_eval(c, s, p)
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

    % Second derivative d^2b/ds^2 (only when requested). BSpline_derivative has
    % no analytic second-derivative sibling, so central-difference the first
    % derivative. Needed by the CLF-QP controller's L_f^2 y term
    % (rabbit_clf_qp_controller.m); ordinary [b,db] callers are unaffected.
    if nargout > 2
        h  = 1e-3;
        sp = min(1, s + h);
        sm = max(0, s - h);
        try
            dNp = BSpline_derivative(n, p.bs_degree, sp);
            dNm = BSpline_derivative(n, p.bs_degree, sm);
            ddb = c(:).' * ((dNp(:) - dNm(:)) / (sp - sm));
        catch
            ddb = 0;
        end
    end
end