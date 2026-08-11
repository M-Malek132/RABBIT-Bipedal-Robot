function [mu, u, qp] = ch5_ctrl_cbf_clf_qp(io, b, p, use_viol)
%CH5_CTRL_CBF_CLF_QP  Section 5.1: the RECIPROCAL CBF-CLF-QP, (5.7)/(5.8)/(5.15).
%
%   [mu, u, qp] = ch5_ctrl_cbf_clf_qp(io, b, p)             direct form (5.7)
%   [mu, u, qp] = ch5_ctrl_cbf_clf_qp(io, b, p, true)       VIOL form  (5.15)
%
% The reciprocal candidate (5.6) is B = 1/h, and the CBF condition (5.3) is
%
%       Bdot(x,u) <= gamma / B(x).
%
% Substituting B = 1/h, so Bdot = -hdot/h^2 and gamma/B = gamma h, gives the
% form actually implemented:
%
%       hdot >= -gamma h^3,        hdot = L_f h + L_g h u.            (*)
%
% Both sides are written out in ch5_barrier so this file never re-derives them.
%
% ============================================================================
%  THE POINT OF THIS FILE IS THE CONDITION IT CANNOT ENFORCE
% ============================================================================
% (*) constrains u through L_g h. For a constraint of relative degree rb > 1,
% L_g h IS IDENTICALLY ZERO -- that is the definition of relative degree -- and
% the barrier row degenerates to
%
%       0 = L_g h u >= -gamma h^3 - L_f h
%
% which contains no control. Two outcomes, and this function distinguishes them
% rather than lumping both into "it failed":
%
%   * right-hand side <= 0: the row is VACUOUSLY TRUE. The QP solves fine and
%     returns the plain CLF-QP control. The barrier is doing nothing at all,
%     and it will keep doing nothing right up until the state crosses the
%     boundary. This is the quiet failure and it is the dangerous one.
%
%   * right-hand side > 0: the row is UNSATISFIABLE and the QP IS INFEASIBLE --
%     the chapter's "the controller will require high control inputs or become
%     infeasible", reproduced exactly.
%
% qp.barrier_controllable and qp.feasible report which. Both Chapter-5 plants
% put their real constraint at relative degree 6 and 4 respectively, so this is
% not an edge case here, it is the generic case -- and it is the entire reason
% Section 5.2 exists. ch5_test_qp asserts it.
%
% For a genuine relative-degree-1 constraint (p.constraint.type = 'v1_max')
% this controller works exactly as advertised, which is how the tests establish
% that the failure above is about relative degree and not about this code.
%
% ---------------------------------------------------------- the two forms
% Remark 5.4 claims (5.7) and the VIOL program (5.15) have IDENTICAL solutions.
% They are both implemented so the claim can be checked rather than believed:
%
%   direct (5.7)   z = [mu; delta],       one inequality row in mu
%   VIOL   (5.15)  z = [mu; mu_b; delta], an EQUALITY defining mu_b and an
%                  inequality on mu_b alone
%
% The VIOL form looks like a detour -- it adds a variable to say something the
% direct form already said. Its value is entirely structural: writing the
% barrier dynamics as "mu_b, with mu_b tied to mu by an equality" is what
% generalizes to rb > 1, where the direct form has nothing to write down. This
% file is where that reformulation is shown to cost nothing.
%
% Inputs
%   io       : from ch5_io_lin
%   b        : from ch5_barrier
%   p        : parameter struct (uses p.cbf.gamma)
%   use_viol : logical, use (5.15) instead of (5.7)
%
% Outputs
%   mu, u, qp  -- qp adds .h .B_rec .mu_b .cbf_active .barrier_controllable
%                 .barrier_defined .rhs
%
% See also CH5_BARRIER, CH5_CTRL_ECBF_CLF_QP, CH5_CTRL_CLF_QP.

if nargin < 4, use_viol = false; end

clf = ch5_res_clf(p);
ny  = p.sys.ny;
gam = p.cbf.gamma;

[V, LfV, LgV] = ch5_clf_eval(io.eta, clf, p.clf.eps);
psi = LfV + clf.lambda * V;

% Barrier row in u:   L_g h u >= -gamma h^3 - L_f h
% and in mu, via u = u_ff + Ainv mu:
%       (L_g h Ainv) mu >= -gamma h^3 - L_f h - L_g h u_ff
Lb  = b.Lgh * io.Ainv;                       % 1 x ny
rhs = -gam * b.h^3 - b.Lfh - b.Lgh * io.u_ff;

qp = struct('V', V, 'LfV', LfV, 'LgV', LgV, 'psi', psi, ...
            'lambda', clf.lambda, 'delta', 0, 'active', psi > 0, ...
            'exitflag', 1, 'feasible', true, ...
            'h', b.h, 'B_rec', b.B_rec, 'mu_b', NaN, 'cbf_active', false, ...
            'barrier_controllable', b.rec_controllable, ...
            'barrier_defined', b.rec_defined, 'rhs', rhs);

%% ----------------------------------- outside the set, B = 1/h is not a barrier
% (5.2) sandwiches B between class-K functions of h for x in Int(C) = {h > 0},
% and that is not a technicality about the domain -- at h <= 0 the reciprocal
% is negative or infinite, so "B large near the boundary" stops being true and
% the condition Bdot <= gamma/B stops meaning anything. ch5_barrier returns
% NaN there rather than a usable number, and passing NaN to quadprog raises
% rather than returning a bad answer, so the state is reported here instead.
%
% This is the OTHER half of the chapter's complaint about reciprocal CBFs: not
% only can they not be written for high relative degree, they also have no
% recovery behaviour if the state ever leaves the set. An ECBF does -- its
% condition is a polynomial in the derivatives of h and stays perfectly well
% defined at h < 0, which is why ch5_ctrl_ecbf_clf_qp needs no equivalent guard.
%
% The row is DROPPED rather than the controller stopping. Stopping would leave
% a simulation frozen on a feedforward the moment it first left the set, and
% the resulting divergence would be an artifact of this decision rather than a
% property of Section 5.1. Falling back to the plain CLF-QP shows what a
% reciprocal-CBF controller actually does past the boundary -- nothing -- and
% qp.barrier_defined records that the barrier had stopped participating.
if ~b.rec_defined
    Lb  = zeros(1, ny);
    rhs = -inf;
end

%% ------------------------------------ the relative-degree > 1 degeneracy
if ~b.rec_controllable
    if rhs > 0
        % Unsatisfiable: no u exists. Report it as infeasible and hold the
        % feedforward, exactly as the constrained solvers below would.
        qp.feasible = false;
        qp.exitflag = -2;
        mu = zeros(ny, 1);
        u  = io.u_ff;
        return;
    end
    % Vacuous: drop the row and fall through to a plain CLF-QP. The barrier is
    % present in name only, which is what the plots then show.
    Lb  = zeros(1, ny);
    rhs = -inf;
end

%% ------------------------------------------------------------- direct (5.7)
% Warm start at the min-norm CLF control -- see ch5_min_norm_mu.
mu_ls = ch5_min_norm_mu(psi, LgV, ny);

if ~use_viol
    % z = [mu; delta].  Barrier row negated into <= form.
    H  = 2 * blkdiag(eye(ny), p.clf.slack_penalty);
    fq = zeros(ny+1, 1);

    [Aineq, bineq] = ch5_scale_row([LgV, -1], -psi, true);   % CLF, slacked

    if isfinite(rhs)
        [ra, rbs] = ch5_scale_row([-Lb, 0], -rhs);           % CBF, NOT slacked
        Aineq = [Aineq; ra];
        bineq = [bineq; rbs];
    end

    [Aineq, bineq] = ch5_box_rows(Aineq, bineq, io, p, 1);

    lb = [-inf(ny,1); 0];
    ub =  inf(ny+1, 1);
    z0 = [mu_ls; 0];
    n_mu = ny;

%% -------------------------------------------------------------- VIOL (5.15)
else
    % z = [mu; mu_b; delta], with mu_b = Bdot(x,mu) imposed as an equality and
    % the CBF condition (5.14) written on mu_b alone.
    %
    % mu_b is the derivative of the RECIPROCAL barrier here, matching (5.13):
    %       mu_b = L_f B_rec + L_g B_rec u,     mu_b <= gamma / B_rec = gamma h.
    if b.rec_defined
        LbB = b.LgB_rec * io.Ainv;                   % 1 x ny
        b0B = b.LfB_rec + b.LgB_rec * io.u_ff;       % scalar
    else
        % Bdot is not defined at h <= 0, so there is nothing for the (VIOL)
        % equality to tie mu_b to. Zero it out; the barrier row is already
        % dropped above, so mu_b is a free variable with no cost and no
        % constraint, and the program reduces to the CLF-QP as intended.
        LbB = zeros(1, ny);
        b0B = 0;
    end

    % w = mu_b / sb; same units argument as in ch5_ctrl_ecbf_clf_qp.
    sb = max(1, max(abs(b0B), abs(gam * b.h)));

    % Token weight on w so H is definite -- see the note in
    % ch5_ctrl_ecbf_clf_qp; the (VIOL) equality pins w to mu regardless.
    H  = 2 * blkdiag(eye(ny), 1e-10, p.clf.slack_penalty);
    fq = zeros(ny+2, 1);

    seq = max(1, norm([LbB, sb, b0B], inf));         % (VIOL) sb w - LbB mu = b0B
    Aeq = [-LbB, sb, 0] / seq;
    beq = b0B / seq;

    [Aineq, bineq] = ch5_scale_row([LgV, 0, -1], -psi, true);   % CLF, slacked

    if isfinite(rhs)
        [ra, rbs] = ch5_scale_row([0*LbB, sb, 0], gam * b.h);   % (5.14)
        Aineq = [Aineq; ra];
        bineq = [bineq; rbs];
    end

    [Aineq, bineq] = ch5_box_rows(Aineq, bineq, io, p, 2);

    lb = [-inf(ny,1); -inf; 0];
    ub =  inf(ny+2, 1);
    z0 = [mu_ls; (b0B + LbB*mu_ls)/sb; 0];
    n_mu = ny;

    % mu_b is a free auxiliary with no cost, so H is only positive
    % SEMI-definite. quadprog's interior-point method handles that, but the
    % equality row is what pins mu_b, so pass it rather than folding it in --
    % folding it in would be the direct form again and defeat the check.
    [z, exitflag] = ch5_solve_qp_eq(H, fq, Aineq, bineq, Aeq, beq, lb, ub, z0);

    if exitflag <= 0 || isempty(z)
        qp.exitflag = exitflag;
        qp.feasible = false;
        mu = zeros(ny, 1);
        u  = io.u_ff;
        return;
    end

    mu = z(1:n_mu);
    u  = io.u_ff + io.Ainv * mu;

    % Recomputed from mu, not read back from the scaled auxiliary, so both
    % solve paths report the identical quantity -- which is what makes the
    % Remark 5.4 comparison in ch5_test_qp a test of the FORMULATIONS rather
    % than of how each one happened to store its answer.
    if b.rec_defined
        qp.mu_b = b.LfB_rec + b.LgB_rec * u;
    end
    qp.delta      = z(end);
    qp.exitflag   = exitflag;
    qp.cbf_active = isfinite(rhs) && (qp.mu_b >= gam*b.h - 1e-8);
    return;
end

%% ------------------------------------------------------------------ solve
[z, exitflag] = ch5_solve_qp(H, fq, Aineq, bineq, lb, ub, z0);

if exitflag <= 0 || isempty(z)
    qp.exitflag = exitflag;
    qp.feasible = false;
    mu = zeros(ny, 1);
    u  = io.u_ff;
    return;
end

mu = z(1:n_mu);
u  = io.u_ff + io.Ainv * mu;

qp.delta      = z(end);
qp.exitflag   = exitflag;
qp.mu_b       = b.LfB_rec + b.LgB_rec * u;
qp.cbf_active = isfinite(rhs) && (Lb*mu - rhs <= 1e-8);

end
