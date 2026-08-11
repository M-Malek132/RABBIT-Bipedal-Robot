function clf = ch5_res_clf(p)
%CH5_RES_CLF  The CLF for the stability row, at ARBITRARY relative degree.
%
%   clf = ch5_res_clf(p)
%
% Chapter 3's ch3_res_clf is hard-wired to relative degree 2, because that is
% what a rigid robot with direct joint torques gives you. Chapter 5's whole
% subject is what happens when that assumption fails, so this is the same
% construction written for general r:
%
%       eta = [y; ydot; ...; y^(r-1)]  in R^(r*ny)
%       etadot = F eta + G mu,   F = kron(shift_r, I_ny),  G = kron(e_r, I_ny)
%
% i.e. ny parallel chains of r integrators, which is exactly what input-output
% linearization delivers once the decoupling matrix is inverted.
%
% ------------------------------------------------------------- the two knobs
% p.clf.eps    the eta-scaling I_eps = blkdiag(I/eps^(r-1), ..., I/eps, I).
%              At r = 2 this is ch3's blkdiag(I/eps, I). eps = 1 is the plain
%              CLF, which is the form (5.31) is written in -- Chapter 5 has no
%              impact map to out-run, so the rapid-exponential machinery is
%              available but not needed, and the default is eps = 1.
%
% p.clf.lambda the rate demanded in  Vdot + lambda V <= delta.  Empty means
%              "use the rate the construction actually guarantees",
%              c3/eps with c3 = lam_min(Q)/lam_max(P).
%
% ---------------------------------------------------- why 'care' is the default
% For the CONTROL algebraic Riccati equation the CLF inequality is pointwise
% feasible at every eta, so the stability row can never be the thing that makes
% the QP infeasible. In this chapter that matters more than it did in Chapter
% 3, because the barrier row is HARD (never slacked). If both rows could go
% infeasible there would be no way to tell which one did; with 'care', an
% infeasible QP is unambiguously the safety constraint talking.
%
% Output
%   clf : struct .P .F .G .r .ny .n .c3 .lam_min_P .lam_max_P
%                .residual .construction .lambda
%
% Cached, because this sits behind a sampled-data control loop.
%
% See also CH5_CLF_EVAL, CH5_CTRL_CLF_QP, CH3_RES_CLF.

persistent cache_key cache_val

r  = p.sys.r;
ny = p.sys.ny;
n  = r * ny;

Q = p.clf.Q;
if isempty(Q), Q = eye(n); end
if ~isequal(size(Q), [n n])
    error('ch5_res_clf:Qsize', ...
          'p.clf.Q must be %dx%d for r=%d, ny=%d (got %dx%d).', ...
          n, n, r, ny, size(Q,1), size(Q,2));
end

key = {lower(p.clf.construction), Q, r, ny, p.clf.eps, p.clf.lambda, p.clf.poles};
if ~isempty(cache_key) && isequaln(cache_key, key)
    clf = cache_val;
    return;
end

% shift_r has ones on the superdiagonal: block i of etadot is block i+1 of eta.
shift = diag(ones(r-1,1), 1);
F = kron(shift, eye(ny));
G = kron([zeros(r-1,1); 1], eye(ny));

switch lower(p.clf.construction)

    case 'care'
        % Stable invariant subspace of the Hamiltonian, same toolbox-free
        % route as ch3_res_clf -- no Control System Toolbox needed.
        Ham = [F, -(G*G.'); -Q, -F.'];
        [Vec, D] = eig(Ham);
        [~, order] = sort(real(diag(D)), 'ascend');
        Xs = Vec(:, order(1:n));
        P  = real(Xs(n+1:end, :) / Xs(1:n, :));
        P  = (P + P.')/2;
        residual = norm(F.'*P + P*F - P*(G*G.')*P + Q, 'fro');

    case 'lyap'
        % A'P + PA = -Q with A = F - G*K pole-placed at p.clf.poles. Kept for
        % parity with ch3_res_clf; the poles must be supplied because unlike
        % r = 2 there is no natural Kp/Kd pair to fall back on.
        if isempty(p.clf.poles)
            error('ch5_res_clf:noPoles', ...
                  ['p.clf.construction = ''lyap'' needs p.clf.poles ' ...
                   '(%d positive values, used as -poles).'], r);
        end
        K = ch5_pole_gain(p.clf.poles, r);        % 1 x r, per output chain
        A = F - G * kron(K, eye(ny));
        if max(real(eig(A))) >= 0
            error('ch5_res_clf:notHurwitz', 'A = F - G K is not Hurwitz.');
        end
        Kmat = kron(eye(n), A.') + kron(A.', eye(n));
        P = reshape(-(Kmat \ Q(:)), n, n);
        P = (P + P.')/2;
        residual = norm(A.'*P + P*A + Q, 'fro');

    otherwise
        error('ch5_res_clf:construction', ...
              'Unknown p.clf.construction "%s" (expected care|lyap).', ...
              p.clf.construction);
end

ev = eig(P);
if min(ev) <= 0
    error('ch5_res_clf:notPositiveDefinite', ...
          'P is not positive definite (min eig %.3e).', min(ev));
end

c3 = min(eig(Q)) / max(ev);

lambda = p.clf.lambda;
if isempty(lambda)
    lambda = c3 / p.clf.eps;
end

clf = struct('P', P, 'F', F, 'G', G, 'r', r, 'ny', ny, 'n', n, ...
             'lam_min_P', min(ev), 'lam_max_P', max(ev), 'c3', c3, ...
             'lambda', lambda, 'residual', residual, ...
             'construction', lower(p.clf.construction));

cache_key = key;
cache_val = clf;

end
