function [z, exitflag] = ch5_solve_qp_eq(H, f, Aineq, bineq, Aeq, beq, lb, ub, z0)
%CH5_SOLVE_QP_EQ  quadprog with equality rows, for the VIOL programs.
%
%   [z, exitflag] = ch5_solve_qp_eq(H, f, Aineq, bineq, Aeq, beq, lb, ub, z0)
%
% Same cached-options reasoning as ch5_solve_qp; split out because the VIOL
% formulations of (5.15) and (5.31) carry the definition of mu_b as an EQUALITY
% row rather than substituting it away.
%
% Keeping it an equality is the whole point of Virtual Input-Output
% Linearization: mu_b is a decision variable that the QP is free to choose,
% tied to mu by the barrier dynamics, and the safety condition is then a
% statement about mu_b alone. Substituting it out recovers the direct form and
% loses the structure that generalizes to high relative degree.
%
% 'active-set' for the reason given in ch5_solve_qp: interior-point fails
% outright on the samples where the barrier row is pushing hard against the
% CLF. H is positive SEMI-definite here (mu_b carries no cost), which
% active-set handles.
%
% See also CH5_SOLVE_QP, CH5_CTRL_CBF_CLF_QP, CH5_CTRL_ECBF_CLF_QP.

persistent opts
if isempty(opts)
    % A raised iteration cap because the VIOL form carries an extra variable
    % and an equality row, so it needs more pivots than the direct form for the
    % same answer -- and at badly scaled states it was hitting the default cap
    % and returning the origin, which is feasible for nothing.
    opts = optimoptions('quadprog', 'Display', 'off', ...
                        'Algorithm', 'active-set', ...
                        'MaxIterations', 400);
end

H = (H + H.') / 2;

[z, ~, exitflag] = quadprog(H, f, Aineq, bineq, Aeq, beq, lb, ub, z0, opts);

end
