function [Aineq, bineq] = ch5_box_rows(Aineq, bineq, io, p, n_tail)
%CH5_BOX_ROWS  Append |u| <= u_max as rows in mu, the Ac(x)u <= bc(x) of (5.8).
%
%   [Aineq, bineq] = ch5_box_rows(Aineq, bineq, io, p, n_tail)
%
% Every QP in this chapter optimizes over mu, not u, because that is how (5.31)
% is written and because mu'mu is the objective. The input box therefore has to
% be pushed through the (IO) row,
%
%       u = u_ff + Ainv mu,        Ainv = (L_g L_f^(r-1) y)^-1
%
% which keeps it linear:
%
%       Ainv mu <=  u_max - u_ff
%      -Ainv mu <=  u_max + u_ff
%
% Note this is an EXACT restatement, not a relaxation: the map from mu to u is
% an invertible affine bijection, so the feasible sets correspond one to one.
% There is no conservatism introduced by working in mu rather than u.
%
% n_tail is the number of decision variables after mu (delta, and mu_b in the
% VIOL form), which get zero columns.
%
% See also CH5_SOLVE_QP, CH5_CTRL_CLF_QP.

if isempty(p.u_max)
    return;
end

pad = zeros(p.sys.nu, n_tail);
umx = p.u_max(:);
if isscalar(umx), umx = umx * ones(p.sys.nu, 1); end

[Abox, bbox] = ch5_scale_row([ io.Ainv, pad ;  ...
                              -io.Ainv, pad ], ...
                             [ umx - io.u_ff ; ...
                               umx + io.u_ff ]);

Aineq = [Aineq; Abox];
bineq = [bineq; bbox];

end
