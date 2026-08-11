function clf = ch3_res_clf(p)
%CH3_RES_CLF  Stage 6: build the rapidly exponentially stabilizing CLF.
%
%   clf = ch3_res_clf(p)
%
% After stage 4 the transverse dynamics are a pure double integrator in the
% transverse variables eta = [y; ydot]:
%
%       etadot = F eta + G mu,     F = [0 I; 0 0],   G = [0; I]
%
% Chapter 3 replaces the fixed PD gain with a Lyapunov function plus a
% REQUIRED CONVERGENCE RATE.  Take
%
%       V_eps(eta) = eta' P_eps eta,   P_eps = I_eps P I_eps,
%       I_eps      = blkdiag((1/eps) I, I)
%
% and demand rapid exponential stability
%
%       Vdot_eps + (c3/eps) V_eps <= 0     with c3 = lam_min(Q)/lam_max(P)
%
% which yields  ||eta(t)|| <= (1/eps) sqrt(c2/c1) exp(-c3 t / 2eps) ||eta(0)||.
% The point of the eps-scaling is that the bound's RATE is 1/eps: eps is the
% knob that buys convergence fast enough to beat the expansion the impact map
% applies to eta at every step.
%
% Any mu satisfying the inequality is admissible -- that is the whole payoff,
% because it turns "which control law?" into "which point in a feasible set?",
% and stages 7-8 then pick the cheapest admissible point.
%
% TWO CONSTRUCTIONS (p.clf_construction):
%
%   'care'  F'P + PF - P G G' P + Q = 0.
%           The standard RES-CLF.  Because it is the CONTROL ARE, the
%           inequality is pointwise feasible for every eta -- so the stage-7 QP
%           is always solvable without slack.  Solved here from the stable
%           invariant subspace of the Hamiltonian matrix, which needs only
%           schur/ordschur (no Control System Toolbox).
%
%   'lyap'  A'P + PA = -Q with A = [0 I; -Kp -Kd], the literal equation
%           written in Chapter 3.  Also a valid CLF -- the stage-5 PD mu always
%           satisfies the decrease condition -- but it is tied to the PD gains
%           rather than being gain-free.  Solved by the Kronecker form, again
%           toolbox-free.
%
%           NOTE THE ASYMMETRY IN WHAT IS GUARANTEED. 'care' is pointwise
%           feasible BY CONSTRUCTION -- the ARE forces LgV = 0 to imply
%           LfV + (c3/eps)V <= 0, for any Q and any eps -- which is why
%           Chapter 5 defaults to it. 'lyap' has no such theorem.
%
%           It nonetheless IS pointwise feasible here, by an exact check rather
%           than a sample: the question is the sign of eig(Z'(M_eps + (c3/eps)
%           P_eps)Z) on Z = null(G'P_eps), MEASURED negative for every eps in
%           [0.05, 2]. But its margin thins with eps -- -0.274 at eps = 0.5
%           against -4.9e-03 at eps = 2, where 'care' holds -3.122 and -0.105.
%           So the default is safe, resting on a measurement rather than a
%           proof; stage 7 carries a guard for the state rather than dividing
%           by zero in it.
%
% Output
%   clf : struct with
%           .P     2ny x 2ny  solution of the (C)ARE / Lyapunov equation
%           .F .G             transverse dynamics matrices
%           .c3               guaranteed rate coefficient
%           .lam_min_P .lam_max_P   the c1, c2 of the bound above
%           .residual         equation residual, for the test to assert on
%           .construction     which construction was used
%         plus the eps-dependent block consumed by ch3_clf_eval:
%           .eps .P_eps .M_eps .PG2_eps .rate
%
% ---------------------------------------------------------------- THE CACHE
% This sits behind an ODE right-hand side that runs thousands of times per
% step, so the construction is cached. Both properties of the cache were
% MEASURED to matter:
%
%   FOUR ENTRIES, NOT ONE. ch3_test_control and ch3_compare_controllers both
%   ALTERNATE constructions, so a one-entry cache missed on every call in
%   exactly the loops the chapter exists to produce: 220.9 us per care/lyap
%   pair against 8.8 us now.
%
%   A CHEAP KEY. The old key was a CELL of matrices compared with isequaln,
%   costing 10.3 us on a HIT -- a third of a 30.4 us stage-7 call, spent
%   proving nothing had changed. Comparing scalars first takes a hit to 3.7 us
%   against a 0.11 us empty-call floor; the rest is the three isequal calls on
%   Q, Kp and Kd, which a hit cannot skip without giving up the guarantee that
%   the cached P belongs to the parameters in hand. The scalars short-circuit
%   MISSES, which is where the alternating-loop win comes from.
%
% See also CH3_CLF_EVAL, CH3_CTRL_CLF_QP.

persistent CACHE

ny = p.ny;
n  = 2*ny;

construction = lower(p.clf_construction);
Q  = p.Q_clf;
Kp = p.Kp;
Kd = p.Kd;

% eps enters only the derived block, not P -- but it is part of the key
% because that block is cached alongside P.
if isfield(p, 'eps') && ~isempty(p.eps)
    eps_clf = p.eps;
else
    eps_clf = 1;
end

%% ------------------------------------------------------------- cache lookup
% Cheap fields first: ny and eps are scalars and settle almost every query
% before Q, Kp or Kd are ever compared.
%
% Read fields through the index -- CACHE(k).ny -- rather than pulling the entry
% out with c = CACHE(k) first. MEASURED: worth nothing on a hit (2.80 vs
% 2.82 us -- copy-on-write means the assignment duplicates nothing) but 3.4x on
% a MISS (1.50 vs 5.03 us over a four-entry scan), because it lets the scalar
% short-circuit fire before the entry is materialized.
for k = 1:numel(CACHE)
    if CACHE(k).ny == ny && CACHE(k).eps == eps_clf ...
            && strcmp(CACHE(k).construction, construction) ...
            && isequal(CACHE(k).Q,  Q)  ...
            && isequal(CACHE(k).Kp, Kp) ...
            && isequal(CACHE(k).Kd, Kd)
        clf = CACHE(k).clf;
        return;
    end
end

F = [zeros(ny), eye(ny); zeros(ny), zeros(ny)];
G = [zeros(ny); eye(ny)];

switch construction

    case 'care'
        % Stable invariant subspace of the Hamiltonian
        %   H = [ F  -G G' ; -Q  -F' ]
        % Columns of [X1; X2] spanning the stable subspace give P = X2/X1.
        %
        % WHY schur/ordschur AND NOT eig. Both find the same subspace here and
        % agree to 1e-14 on this problem, so this is not a bug fix -- it is the
        % numerically supported route, and the one MATLAB's own care() takes.
        % eig has to sort the eigenvalues and gather the selected eigenVECTORS,
        % which is exactly the step that loses meaning when the Hamiltonian is
        % defective or has eigenvalues near the imaginary axis; the ordered real
        % Schur form delivers an orthonormal basis for the stable subspace
        % without ever forming an eigenvector. The residual check below is what
        % actually certifies the answer either way.
        Ham = [F, -(G*G.'); -Q, -F.'];
        [U, T] = schur(Ham, 'real');
        [U, ~] = ordschur(U, T, 'lhp');     % stable block to the top-left
        X1 = U(1:n,     1:n);
        X2 = U(n+1:end, 1:n);
        rc1 = rcond(X1);
        if ~isfinite(rc1) || rc1 < 1e-12
            error('ch3_res_clf:singularSubspace', ...
                  ['The stable invariant subspace of the Hamiltonian is not ' ...
                   'graphable (rcond(X1) = %.3e), so P = X2/X1 is not ' ...
                   'defined. Check p.Q_clf > 0.'], rc1);
        end
        P = real(X2 / X1);
        P = (P + P.')/2;                    % symmetrize
        residual = norm(F.'*P + P*F - P*(G*G.')*P + Q, 'fro');

    case 'lyap'
        A = [zeros(ny), eye(ny); -Kp, -Kd];
        if max(real(eig(A))) >= 0
            error('ch3_res_clf:notHurwitz', ...
                  ['A = [0 I; -Kp -Kd] is not Hurwitz, so the Lyapunov ' ...
                   'construction has no positive definite solution. ' ...
                   'Check p.Kp / p.Kd.']);
        end
        % vec(A'P + PA) = (I kron A' + A' kron I) vec(P) = -vec(Q)
        I_n = eye(n);
        Kmat = kron(I_n, A.') + kron(A.', I_n);
        P = reshape(-(Kmat \ Q(:)), n, n);
        P = (P + P.')/2;
        residual = norm(A.'*P + P*A + Q, 'fro');

    otherwise
        error('ch3_res_clf:construction', ...
              'Unknown p.clf_construction "%s" (expected care|lyap).', ...
              p.clf_construction);
end

%% ------------------------------------------------------------- certify it
% The residual was already computed and returned for ch3_test_control to
% assert on. Check it HERE as well: the test runs on the default parameters,
% while this function is called with whatever Q, Kp and Kd a study happens to
% sweep, and a silently unsolved equation would come back as a V that is not a
% CLF at all -- a wrong answer downstream rather than an error here. MEASURED
% on the defaults: 1.2e-14 (care) and 0 (lyap), so the gate is many orders
% clear of anything the good path produces.
res_tol = 1e-8 * max(1, norm(Q, 'fro'));
if ~isfinite(residual) || residual > res_tol
    error('ch3_res_clf:unsolved', ...
          ['The %s equation was not solved: residual %.3e exceeds %.3e. ' ...
           'P is not a valid CLF certificate.'], construction, residual, res_tol);
end

ev = eig(P);
if min(ev) <= 0
    error('ch3_res_clf:notPositiveDefinite', ...
          'P is not positive definite (min eig %.3e).', min(ev));
end

clf = struct('P', P, 'F', F, 'G', G, ...
             'lam_min_P', min(ev), 'lam_max_P', max(ev), ...
             'c3', min(eig(Q)) / max(ev), ...
             'residual', residual, ...
             'construction', construction);

%% ------------------------------------------- the eps-scaled block, once
% ch3_clf_eval used to rebuild all of this on EVERY evaluation, inside the ODE
% right-hand side, at 9.5 us a call -- for quantities that depend only on P and
% eps and so never change during a rollout. Same argument, and same size of
% win, as hoisting the quadprog options in ch3_ctrl_clf_qp.
%
% F'P_eps + P_eps F is NOT zero despite F being nilpotent: it is the
% off-diagonal coupling that makes LfV nonzero, i.e. the reason the CLF
% condition has content at all.
I_eps = blkdiag(eye(ny)/eps_clf, eye(ny));

clf.eps     = eps_clf;
clf.P_eps   = I_eps * P * I_eps;
clf.P_eps   = (clf.P_eps + clf.P_eps.')/2;
clf.M_eps   = F.'*clf.P_eps + clf.P_eps*F;      % LfV = eta' M_eps eta
clf.PG2_eps = 2 * clf.P_eps * G;                % LgV = eta' PG2_eps  (1 x ny)
clf.rate    = clf.c3 / eps_clf;                 % the (c3/eps) of the condition

%% ------------------------------------------------------------- cache store
entry = struct('ny', ny, 'eps', eps_clf, 'construction', construction, ...
               'Q', Q, 'Kp', Kp, 'Kd', Kd, 'clf', clf);

if isempty(CACHE)
    CACHE = entry;
else
    CACHE = [entry, CACHE];                 % most recent first
    if numel(CACHE) > 4
        CACHE = CACHE(1:4);
    end
end

end
