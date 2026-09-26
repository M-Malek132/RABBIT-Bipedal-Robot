function [mu, u, qp] = ch5_ctrl_ecbf_clf_qp(io, b, e, p, use_viol, X)
%CH5_CTRL_ECBF_CLF_QP  Section 5.2: the ECBF-CLF-QP, eq (5.31). The chapter's result.
%
%   [mu, u, qp] = ch5_ctrl_ecbf_clf_qp(io, b, e, p)
%   [mu, u, qp] = ch5_ctrl_ecbf_clf_qp(io, b, e, p, true)   explicit mu_b variable
%   [mu, u, qp] = ch5_ctrl_ecbf_clf_qp(io, b, e, p, false, X)
%                 with the EXTRA barriers of ch5_extra_barriers: one more row
%                 per barrier, of the same form, none of them slacked
%
%       min_{mu,delta}  mu'mu + p delta^2                              (5.31)
%       s.t.  Vdot(eta,mu) + lambda V(eta) <= delta          (CLF)
%             Ac(x) mu <= bc(x)                              (Constraints)
%             mu_b >= -Kb eta_b                              (Exponential CBF)
%             h^(rb)(x,mu) = mu_b                            (VIOL)
%
% ------------------------------------------------ what makes this row possible
% Section 5.1 could not write a barrier row for these plants because Bdot did
% not contain u. VIOL fixes that by refusing to differentiate only once.
% Differentiate h all the way to rb, which is by definition where u first
% appears, and apply (5.21):
%
%       h^(rb) = L_f^rb h + L_g L_f^(rb-1) h * u
%              = L_f^rb h + L_g L_f^(rb-1) h * (u_ff + Ainv mu)
%              =: b0(x) + Lb(x) mu   =:  mu_b                       (VIOL)
%
% mu_b is AFFINE IN mu -- one scalar row, exactly like the CLF row -- so the
% safety condition (5.25),
%
%       mu_b >= -Kb eta_b,
%
% is a single linear inequality and the program stays a quadratic program at
% ANY relative degree. That is the whole technical content of Section 5.2.
%
% -------------------------------------------------------- Remark 5.6, in code
% The row above is identical to y_rb(x) >= 0 for the family (5.28), because Kb
% was built from those same poles:
%
%       h^(rb) + a1 h^(rb-1) + ... + a_rb h  =  (d/dt+p_1)o...o(d/dt+p_rb) h.
%
% qp.y_rb returns that value, so a run can be checked against Theorem 5.1's
% chain of sets rather than only against h >= 0. ch5_test_ecbf asserts the
% identity.
%
% -------------------------------------------------- what is hard and what bends
% THE BARRIER ROW IS NEVER SLACKED. Only the CLF row gets delta. When the two
% conflict -- and on the pendulum they genuinely do, since the CLF wants
% theta2 -> 0 while the barrier needs the arm folded -- safety wins and
% tracking waits. Reversing that would produce a controller that is
% occasionally unsafe and always convergent, which is the wrong trade and not
% what (5.31) says.
%
% ------------------------------------------------------------- degeneracy
% Lb = L_g L_f^(rb-1) h * Ainv can vanish at isolated states: on the pendulum
% it equals dpy/dtheta, which is zero exactly when the arm is straight up or
% straight down. There the barrier row loses its grip on mu for that instant.
% It is reported (qp.barrier_controllable) rather than smoothed over, because
% the honest statement is that the guarantee has a gap at those configurations
% -- and the run then shows whether it mattered. For Fig. 5.4 it does not: the
% straight-up configuration is where h is at its largest, so the row is far
% from active there.
%
% Inputs
%   io       : from ch5_io_lin
%   b        : from ch5_barrier
%   e        : from ch5_ecbf_gain
%   p        : parameter struct
%   use_viol : keep mu_b as an explicit decision variable (default false).
%              Mathematically identical, per Remark 5.4; ch5_test_qp checks it.
%
% Outputs
%   mu, u, qp  -- qp adds .h .eta_b .mu_b .y_rb .Kb_eta .cbf_active
%                 .barrier_controllable .margin .Lb_norm
%                 .extra  1 x n struct .h .y_rb .active .controllable, one per
%                         extra barrier (empty without them)
%
% See also CH5_ECBF_GAIN, CH5_BARRIER, CH5_ECBF_ADMISSIBLE, CH5_CTRL_CBF_CLF_QP,
%          CH5_EXTRA_BARRIERS.

if nargin < 5 || isempty(use_viol), use_viol = false; end
if nargin < 6, X = []; end
if use_viol && ~isempty(X)
    error('ch5_ctrl_ecbf_clf_qp:extraViol', ...
          'Extra barriers are implemented in the direct form only (use_viol = false).');
end

clf = ch5_res_clf(p);
ny  = p.sys.ny;

[V, LfV, LgV] = ch5_clf_eval(io.eta, clf, p.clf.eps);
psi = LfV + clf.lambda * V;

%% -------------------------------------------------------------- VIOL terms
Lb = b.LgLfrb1 * io.Ainv;                    % 1 x ny, coefficient of mu
b0 = b.Lfrb + b.LgLfrb1 * io.u_ff;           % scalar, the mu-free part

Kb_eta = e.Kb * b.eta_b;                     % scalar

% ECBF row:  mu_b >= -Kb eta_b   <=>   b0 + Lb mu >= -Kb eta_b
%                                 <=>  -Lb mu <= b0 + Kb eta_b
row_A = -Lb;
row_b = b0 + Kb_eta;

controllable = norm(Lb, inf) > 1e-12;

qp = struct('V', V, 'LfV', LfV, 'LgV', LgV, 'psi', psi, ...
            'lambda', clf.lambda, 'delta', 0, 'active', psi > 0, ...
            'exitflag', 1, 'feasible', true, ...
            'h', b.h, 'eta_b', b.eta_b, 'mu_b', NaN, 'y_rb', NaN, ...
            'Kb_eta', Kb_eta, 'cbf_active', false, ...
            'barrier_controllable', controllable, 'margin', NaN, ...
            'Lb_norm', norm(Lb));

%% ------------------------------------------------------- extra barriers
% Same construction as the main row, one per barrier: b0x + Lbx mu >= -Kbx eta_bx.
nx_ = numel(X);
XA = zeros(0, ny); Xb = zeros(0, 1);
ext = struct('h', cell(1, nx_), 'y_rb', NaN, 'active', false, ...
             'controllable', true, 'Lb', [], 'b0', NaN, 'Kb_eta', NaN);
for k = 1:nx_
    bx  = X(k).b;
    Lbx = bx.LgLfrb1 * io.Ainv;
    b0x = bx.Lfrb + bx.LgLfrb1 * io.u_ff;
    Kx  = X(k).e.Kb * bx.eta_b;
    ext(k).h = bx.h; ext(k).Lb = Lbx; ext(k).b0 = b0x; ext(k).Kb_eta = Kx;
    ext(k).controllable = norm(Lbx, inf) > 1e-12;
    if ext(k).controllable
        XA = [XA; -Lbx];            %#ok<AGROW>
        Xb = [Xb; b0x + Kx];        %#ok<AGROW>
    elseif b0x + Kx < 0
        qp.feasible = false;        % an extra row with no grip, violated
    end
end
qp.extra = ext;
if ~qp.feasible
    qp.exitflag = -2;
    mu = zeros(ny, 1);
    u  = io.u_ff;
    return;
end

%% ------------------------------------------ the row has no grip on mu here
if ~controllable
    if row_b < 0
        qp.feasible = false;
        qp.exitflag = -2;
        mu = zeros(ny, 1);
        u  = io.u_ff;
        qp.mu_b = b0;
        qp.y_rb = b0 + Kb_eta;
        return;
    end
    row_A = [];                              % vacuously satisfied; drop it
    row_b = [];
end

%% --------------------------------------------------------------- assemble
% Warm start at the min-norm CLF control. See ch5_min_norm_mu for why this is
% a correctness fix and not a speed one.
mu_ls = ch5_min_norm_mu(psi, LgV, ny);

if ~use_viol
    % z = [mu; delta]
    H  = 2 * blkdiag(eye(ny), p.clf.slack_penalty);
    fq = zeros(ny+1, 1);

    [Aineq, bineq] = ch5_scale_row([LgV, -1], -psi, true);

    if ~isempty(row_A)
        [ra, rbs] = ch5_scale_row([row_A, 0], row_b);
        Aineq = [Aineq; ra];
        bineq = [bineq; rbs];
    end
    for k = 1:size(XA, 1)
        [ra, rbs] = ch5_scale_row([XA(k, :), 0], Xb(k));
        Aineq = [Aineq; ra];         %#ok<AGROW>
        bineq = [bineq; rbs];        %#ok<AGROW>
    end

    [Aineq, bineq] = ch5_box_rows(Aineq, bineq, io, p, 1);

    lb = [-inf(ny,1); 0];
    ub =  inf(ny+1, 1);

    [z, exitflag] = ch5_solve_qp(H, fq, Aineq, bineq, lb, ub, [mu_ls; 0]);

else
    % z = [mu; w; delta], the literal reading of (5.31) with both the (VIOL)
    % equality and the (Exponential CBF) inequality present as written.
    %
    % w = mu_b / sb, not mu_b itself. On the pendulum mu_b = h^(4) runs to 1e5
    % while mu is O(10), and a zero-cost variable five orders of magnitude
    % larger than the costed ones is what makes an interior-point solver report
    % a positive-definite problem unbounded. sb is a pure change of units and
    % is undone on the way out.
    sb = max(1, max(abs(b0), abs(Kb_eta)));

    % w carries NO cost in (5.31), which leaves H positive SEMI-definite. The
    % active-set solver needs a definite Hessian and returns a different point
    % without one, so w gets a token weight. It is not a design choice and it
    % does not bend the formulation: the (VIOL) equality determines w uniquely
    % from mu, so the feasible set is a graph over mu and this perturbs the
    % optimum by O(reg) -- with w scaled to O(1) by sb, that is 1e-10.
    % ch5_test_qp checks the two forms still agree.
    reg = 1e-10;

    H  = 2 * blkdiag(eye(ny), reg, p.clf.slack_penalty);
    fq = zeros(ny+2, 1);

    % (VIOL)  sb w - Lb mu = b0
    seq = max(1, norm([Lb, sb, b0], inf));
    Aeq = [-Lb, sb, 0] / seq;
    beq = b0 / seq;

    [Aineq, bineq] = ch5_scale_row([LgV, 0, -1], -psi, true);

    % (Exponential CBF)  mu_b >= -Kb eta_b   ->   -sb w <= Kb eta_b
    [ra, rbs] = ch5_scale_row([zeros(1,ny), -sb, 0], Kb_eta);
    Aineq = [Aineq; ra];
    bineq = [bineq; rbs];

    [Aineq, bineq] = ch5_box_rows(Aineq, bineq, io, p, 2);

    lb = [-inf(ny,1); -inf; 0];
    ub =  inf(ny+2, 1);

    [z, exitflag] = ch5_solve_qp_eq(H, fq, Aineq, bineq, Aeq, beq, lb, ub, ...
                                    [mu_ls; (b0 + Lb*mu_ls)/sb; 0]);
end

%% ------------------------------------------------------------------ report
if exitflag <= 0 || isempty(z)
    % Infeasible means the barrier row and the input box cannot both hold: the
    % CLF row is slacked and the barrier row alone is always satisfiable in mu
    % when Lb ~= 0, so with no box this branch is unreachable. Hold u_ff.
    qp.exitflag = exitflag;
    qp.feasible = false;
    mu = zeros(ny, 1);
    u  = io.u_ff;
    qp.mu_b = b0;
    qp.y_rb = b0 + Kb_eta;
    return;
end

mu = z(1:ny);
u  = io.u_ff + io.Ainv * mu;

qp.delta    = z(end);
qp.exitflag = exitflag;

% Recompute mu_b from mu rather than reading the auxiliary variable, so the
% two solve paths report the same quantity and the VIOL equality is checked
% implicitly on every call.
qp.mu_b   = b0 + Lb * mu;
qp.y_rb   = qp.mu_b + Kb_eta;                % Remark 5.6: this is y_rb(x)
qp.margin = qp.y_rb;                         % >= 0 is the guarantee holding
qp.cbf_active = controllable && (qp.y_rb <= 1e-7 * max(1, abs(Kb_eta)));

for k = 1:nx_
    ext(k).y_rb   = ext(k).b0 + ext(k).Lb * mu + ext(k).Kb_eta;
    ext(k).active = ext(k).controllable && ...
                    (ext(k).y_rb <= 1e-7 * max(1, abs(ext(k).Kb_eta)));
end
qp.extra = ext;

end
