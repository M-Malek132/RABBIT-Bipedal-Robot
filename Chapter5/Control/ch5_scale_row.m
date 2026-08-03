function [A, b] = ch5_scale_row(A, b, slack_col)
%CH5_SCALE_ROW  Rescale QP inequality rows for conditioning.
%
%   [A, b] = ch5_scale_row(A, b)              hard rows: exactly a no-op
%   [A, b] = ch5_scale_row(A, b, true)        a row whose LAST column is -delta
%
% WHY ANY OF THIS. The CLF row carries LgV and psi, both of which scale with
% the CLF's Q. On the pendulum with Q = 1000 I they reach ~1e5 while the
% objective's mu block is O(1) -- five orders of magnitude inside a
% three-variable program. quadprog's interior-point method reported exitflag
% -3 ("unbounded") on ~0.3% of samples of a problem with a positive definite
% Hessian, which cannot be unbounded. After scaling, 0%.
%
% ============================================================================
%  HARD ROWS: divide by ||coefficients||, NOT by ||[coefficients, rhs]||.
% ============================================================================
% Either is an exact no-op -- dividing an inequality by a positive constant
% does not move its feasible set -- but only one leaves a usable gradient.
%
% The ECBF row is  -L_b mu <= b0 + Kb eta_b, and when the barrier is pushing
% hard the right-hand side dominates: measured on the pendulum, |L_b| ~ 1.8
% against a right-hand side of -167. Including the right-hand side in the norm
% divides the row by 167 and leaves  -0.011 mu <= -1 -- an almost flat
% constraint that still demands ||mu|| >= 93. quadprog's interior-point method
% returned exitflag -3 ("unbounded") on every one of the 222 samples where this
% happened. Normalizing by the coefficients alone leaves -1.0 mu <= -93, which
% says the same thing with a gradient of unit length, and solves.
%
% So: the right-hand side never sets the units, for either kind of row.
% ============================================================================
%
% ============================================================================
%  THE SLACKED ROW: divide by ||coefficients|| ONLY, and leave -delta alone.
% ============================================================================
% This one is not a no-op and the choice took two wrong turns to get right, so
% both are recorded here -- each looked correct and each broke something a long
% way from this file.
%
%   ATTEMPT 1, divide the mu coefficients and the RHS by s = ||[LgV, psi]||,
%   keep delta's coefficient at -1. That redefines the variable as
%   delta' = delta/s, so the penalty p*delta'^2 is p/s^2 in the original units.
%   With psi ~ 1e5 dominating s, the slack became nearly FREE and the QP
%   started relaxing a CLF condition that was perfectly satisfiable: the
%   ECBF-CLF-QP with its barrier row inactive stopped reproducing the plain
%   CLF-QP's closed form, drifting 5e-5 for no visible reason.
%
%   ATTEMPT 2, divide the entire row including delta's coefficient. Now delta
%   keeps its absolute meaning -- but its coefficient becomes -1/s, so
%   relaxing the row by Delta costs p*s^2*Delta^2. At s ~ 1e5 that is ~1e10*p,
%   and the slack became numerically UNUSABLE. On the pendulum, where the CLF
%   row and the barrier row genuinely do conflict (the CLF wants theta2 -> 0
%   while the barrier needs the arm folded), the solver could no longer pay for
%   the relaxation, fell back to the feedforward, and the constraint was
%   violated by a full metre.
%
%   WHAT WORKS: scale by s = ||coefficients||_inf, EXCLUDING the right-hand
%   side, and keep delta's coefficient at -1. Then
%
%     * the mu coefficients come out O(1), which is the conditioning fix;
%     * a unit of delta still buys a unit of relaxation, so the slack stays
%       usable when the rows really do conflict;
%     * the effective penalty is p/||LgV||^2, a RELATIVE penalty -- the same
%       choice, for the same reason, as the mu_scale factor in ch3_ctrl_clf_qp
%       and ch4_ctrl_rclf_qp;
%     * the gap to the hard-constrained closed form is 1/p, independent of psi,
%       because the normalized coefficient vector has unit norm by construction.
%
% The right-hand side is deliberately left out of s: psi is how far the CLF
% condition currently is from being satisfied, and letting that quantity set
% the units is what made both wrong attempts wrong.
%
% See also CH5_SOLVE_QP, CH5_BOX_ROWS, CH5_CTRL_CLF_QP, CH5_CTRL_ECBF_CLF_QP.

if nargin < 3, slack_col = false; end
if isempty(A), return; end

b = b(:);

if slack_col
    % delta's column is left at -1 so a unit of slack still buys a unit of
    % relaxation; only the mu coefficients and the right-hand side are scaled.
    s = max(1, max(abs(A(:, 1:end-1)), [], 2));
    A(:, 1:end-1) = A(:, 1:end-1) ./ s;
    b = b ./ s;
else
    s = max(1, max(abs(A), [], 2));
    A = A ./ s;
    b = b ./ s;
end

end
