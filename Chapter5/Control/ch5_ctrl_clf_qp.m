function [mu, u, qp] = ch5_ctrl_clf_qp(io, p)
%CH5_CTRL_CLF_QP  The stability-only baseline. No barrier row at all.
%
%   [mu, u, qp] = ch5_ctrl_clf_qp(io, p)
%
%       min_{mu,delta}  mu'mu + p delta^2
%       s.t.  LfV + LgV mu + lambda V <= delta            (CLF)
%             |u| <= u_max                                (optional)
%             u = u_ff + (L_g L_f^(r-1) y)^-1 mu          (IO)
%
% THIS CONTROLLER IS IN THE CHAPTER TO FAIL. It is Fig. 5.3a and Fig. 5.4a: a
% perfectly good stabilizing controller that walks straight through the safety
% constraint, because nothing ever told it the constraint exists. Everything
% Section 5.2 adds is measured against this line.
%
% Its failure is also not a tuning problem, which is worth being clear about.
% Slowing the CLF down would reduce the overshoot in Fig. 5.3a, but the
% constraint would still be enforced only by luck, and the same controller
% would violate it again from a different initial condition. The barrier row
% is a different KIND of guarantee, not a better-tuned version of this one.
%
% ------------------------------------------------------- the two solve paths
% With no input box the problem is a least-norm point on a single halfspace and
% has the closed form ch3_ctrl_clf_qp uses -- exact, and far cheaper than a
% solver inside a sampled-data loop. The slack is then structurally
% unnecessary: with 'care' the CLF inequality is pointwise feasible, so the
% optimal delta is zero and the closed form IS the p -> inf limit. With a box,
% quadprog, and the slack starts earning its penalty.
%
% Inputs
%   io : from ch5_io_lin
%   p  : parameter struct
%
% Outputs
%   mu : ny x 1
%   u  : nu x 1
%   qp : struct .V .LfV .LgV .psi .lambda .delta .active .exitflag .feasible
%
% See also CH5_RES_CLF, CH5_CLF_EVAL, CH5_CTRL_ECBF_CLF_QP.

clf = ch5_res_clf(p);
ny  = p.sys.ny;

[V, LfV, LgV] = ch5_clf_eval(io.eta, clf, p.clf.eps);
psi = LfV + clf.lambda * V;

% cbf_active and barrier_controllable are always false here and are carried
% anyway, so that every controller in the chapter returns the same record and
% ch5_simulate can log a baseline run and a barrier run through one code path.
% A baseline whose diagnostics are shaped differently from the thing it is
% being compared against is a baseline you cannot plot on the same axes.
qp = struct('V', V, 'LfV', LfV, 'LgV', LgV, 'psi', psi, ...
            'lambda', clf.lambda, 'delta', 0, 'active', false, ...
            'exitflag', 1, 'feasible', true, ...
            'cbf_active', false, 'barrier_controllable', false);

%% -------------------------------------------------- closed form, no input box
mu_ls = ch5_min_norm_mu(psi, LgV, ny);

if isempty(p.u_max)
    mu = mu_ls;
    qp.active = any(mu ~= 0);
    u = io.u_ff + io.Ainv * mu;
    return;
end

%% ------------------------------------------------------- with the input box
% z = [mu; delta].
H  = 2 * blkdiag(eye(ny), p.clf.slack_penalty);
fq = zeros(ny+1, 1);

[Aineq, bineq] = ch5_scale_row([LgV, -1], -psi, true);

[Aineq, bineq] = ch5_box_rows(Aineq, bineq, io, p, 1);

lb = [-inf(ny,1); 0];
ub =  inf(ny+1, 1);

[z, exitflag] = ch5_solve_qp(H, fq, Aineq, bineq, lb, ub, [mu_ls; 0]);

if exitflag <= 0 || isempty(z)
    % The box alone made it infeasible, which for a slacked CLF row means the
    % box is empty -- report rather than invent a control.
    mu = zeros(ny, 1);
    u  = min(max(io.u_ff, -p.u_max), p.u_max);
    qp.exitflag = exitflag;
    qp.feasible = false;
    return;
end

mu = z(1:ny);
u  = io.u_ff + io.Ainv * mu;

qp.delta    = z(end);
qp.exitflag = exitflag;
qp.active   = psi > 0;

end
