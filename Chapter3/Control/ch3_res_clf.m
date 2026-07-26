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
%           invariant subspace of the Hamiltonian matrix, which needs only eig
%           (no Control System Toolbox).
%
%   'lyap'  A'P + PA = -Q with A = [0 I; -Kp -Kd], the literal equation
%           written in Chapter 3.  Also a valid CLF -- the stage-5 PD mu always
%           satisfies the decrease condition -- but it is tied to the PD gains
%           rather than being gain-free.  Solved by the Kronecker form, again
%           toolbox-free.
%
% Output
%   clf : struct with
%           .P     2ny x 2ny  solution of the (C)ARE / Lyapunov equation
%           .F .G             transverse dynamics matrices
%           .c3               guaranteed rate coefficient
%           .lam_min_P .lam_max_P   the c1, c2 of the bound above
%           .residual         equation residual, for the test to assert on
%           .construction     which construction was used
%
% Results are cached across calls (the construction depends only on Q, the
% gains and the choice of equation), because this sits behind an ODE
% right-hand side that runs thousands of times per step.
%
% See also CH3_CLF_EVAL, CH3_CTRL_CLF_QP.

persistent cache_key cache_val

ny = p.ny;
n  = 2*ny;

F = [zeros(ny), eye(ny); zeros(ny), zeros(ny)];
G = [zeros(ny); eye(ny)];
Q = p.Q_clf;

key = {lower(p.clf_construction), Q, p.Kp, p.Kd, ny};
if ~isempty(cache_key) && isequaln(cache_key, key)
    clf = cache_val;
    return;
end

switch lower(p.clf_construction)

    case 'care'
        % Stable invariant subspace of the Hamiltonian
        %   H = [ F  -G G' ; -Q  -F' ]
        % Columns of [X1; X2] spanning the stable subspace give P = X2/X1.
        Ham = [F, -(G*G.'); -Q, -F.'];
        [Vec, D] = eig(Ham);
        lam = diag(D);
        [~, order] = sort(real(lam), 'ascend');
        sel = order(1:n);                       % n most stable eigenvalues
        Xs  = Vec(:, sel);
        X1  = Xs(1:n, :);
        X2  = Xs(n+1:end, :);
        P   = real(X2 / X1);
        P   = (P + P.')/2;                      % symmetrize
        residual = norm(F.'*P + P*F - P*(G*G.')*P + Q, 'fro');

    case 'lyap'
        A = [zeros(ny), eye(ny); -p.Kp, -p.Kd];
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

ev = eig(P);
if min(ev) <= 0
    error('ch3_res_clf:notPositiveDefinite', ...
          'P is not positive definite (min eig %.3e).', min(ev));
end

clf = struct('P', P, 'F', F, 'G', G, ...
             'lam_min_P', min(ev), 'lam_max_P', max(ev), ...
             'c3', min(eig(Q)) / max(ev), ...
             'residual', residual, ...
             'construction', lower(p.clf_construction));

cache_key = key;
cache_val = clf;

end
