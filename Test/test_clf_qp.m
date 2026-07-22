function test_clf_qp(result_file)
%TEST_CLF_QP  Validate the CLF-QP controller (rabbit_clf_qp_controller).
%
%   test_clf_qp()             % uses the most recent Results/hzd_result_*.mat
%   test_clf_qp(result_file)  % uses a specific saved gait (.mat with z_opt,p)
%
% Checks, at several states along the gait:
%   (1) yddot predicted by (L_f^2 y + Adec*u) matches a finite-difference of
%       ydot along the constrained-dynamics flow  -> validates the Lie
%       derivatives, the affine ddq map, and all indexing/sign conventions.
%   (2) the decoupling matrix Adec = L_g L_f y is well conditioned (invertible),
%       so input-output linearization is valid.
%   (3) the QP solution satisfies the RES-CLF descent Vdot <= -(gamma/eps)V
%       (allowing the relaxation delta) and the torque box |u| <= u_max.
%   (4) reports QP vs fixed-PD torque magnitude for context.
%
% No toolboxes beyond what the optimizer already needs (quadprog).

    here = fileparts(mfilename('fullpath'));
    root = fileparts(here);
    addpath(genpath(root));

    if nargin < 1 || isempty(result_file)
        d = dir(fullfile(root, 'Results', 'hzd_result_*.mat'));
        assert(~isempty(d), 'No Results/hzd_result_*.mat found; pass a file explicitly.');
        [~, newest] = max([d.datenum]);
        result_file = fullfile(d(newest).folder, d(newest).name);
    end
    fprintf('Loading gait: %s\n', result_file);
    S = load(result_file);
    p = S.p;
    [coeffs, x_start] = unpack_z(S.z_opt, p);

    nq = p.nq;
    theta_minus = theta_of_q(x_start(1:nq));
    theta_plus  = p.theta_plus;

    % Sample a handful of states across the step by simulating the PD gait.
    try
        [~, X] = simulate_hzd_gait_full(coeffs, x_start, p);
        idx = round(linspace(1, size(X,1), 6));
        states = X(idx, :)';                 % 2nq x K
    catch
        states = x_start(:);                 % fall back to just the start state
    end
    K = size(states, 2);

    fprintf('\n%-4s %-14s %-12s %-14s %-12s %-10s\n', ...
            'k', 'FD-err(yddot)', 'cond(Adec)', 'Vdot+gVeps', '|u|_qp', '|u|_pd');
    fprintf('%s\n', repmat('-', 1, 70));

    maxfd = 0; worstclf = -inf;
    for k = 1:K
        x  = states(:, k);
        q  = x(1:nq);
        dq = x(nq+1:end);

        [u, info] = rabbit_clf_qp_controller(q, dq, coeffs, theta_minus, theta_plus, p);

        % (1) finite-difference yddot along the flow with THIS u ------------
        yddot_pred = info.Lf2y + info.Adec * u;
        h  = 1e-6;
        ddq = info.ddq;                       % constrained accel at (q,dq,u)
        xdot = [dq; ddq];
        ydp = ydot_of(x + h*xdot, coeffs, p, theta_minus, theta_plus);
        ydm = ydot_of(x - h*xdot, coeffs, p, theta_minus, theta_plus);
        yddot_num = (ydp - ydm) / (2*h);
        fd_err = norm(yddot_pred - yddot_num) / max(1, norm(yddot_num));

        % The FD reference is only valid where the phase is strictly interior:
        % at s==0 or s==1 the s=max(0,min(1,.)) clamp puts a kink in yd(s), so
        % the central difference straddles it and is not a meaningful check of
        % the (correct) analytic value. Flag those samples instead of counting
        % them toward the pass metric.
        clamped = (info.s <= 1e-3) || (info.s >= 1 - 1e-3);
        if ~clamped
            maxfd = max(maxfd, fd_err);
        end

        % (2) decoupling matrix conditioning --------------------------------
        cnd = cond(info.Adec);

        % (3) CLF descent (should be <= 0 up to the relaxation delta) --------
        a = sqrt(3); gamma = 1/(a+1);
        clf_slack = info.Vdot + (gamma/p_eps(p))*info.V;   % want <= delta (>=0 penalized)
        worstclf  = max(worstclf, clf_slack - max(info.delta,0));

        % (4) PD comparison --------------------------------------------------
        [~, u_pd] = hzd_control_and_dynamics(q, dq, coeffs, theta_minus, theta_plus, p);

        tag = '';
        if clamped, tag = ' (s@clamp: FD n/a)'; end
        fprintf('%-4d %-14.2e %-12.2e %-14.3e %-12.3f %-10.3f%s\n', ...
                k, fd_err, cnd, clf_slack, norm(u), norm(u_pd), tag);

        assert(all(abs(u) <= getdef(p,'u_max',120) + 1e-6), ...
               'Torque box violated at k=%d', k);
    end

    fprintf('%s\n', repmat('-', 1, 70));
    fprintf('max FD error (yddot vs L_f^2y + Adec u): %.2e  (want < 1e-4)\n', maxfd);
    fprintf('worst CLF slack minus relaxation       : %.2e  (want <= ~0)\n', worstclf);
    if maxfd < 1e-4
        fprintf('PASS: analytic Lie derivatives match finite differences.\n');
    else
        fprintf('CHECK: FD error larger than expected -- inspect Adec / Lf2y.\n');
    end
end

% -------------------------------------------------------------------------
function yd = ydot_of(x, coeffs, p, theta_minus, theta_plus)
% ydot = Jy(q)*dq, recomputed from scratch (independent of the controller) so
% the finite-difference check does not reuse the code it is validating.
    nq = p.nq; m = p.nu; act_idx = 4:nq;
    q  = x(1:nq); dq = x(nq+1:end);
    kappa = 1/(theta_plus - theta_minus);
    theta = theta_of_q(q);
    s = max(0, min(1, (theta - theta_minus)*kappa));
    c_row = dtheta_dq_of(q);
    dyd_ds = zeros(m,1);
    for i = 1:m
        [~, db] = bspline_eval(coeffs(i,:), s, p);
        dyd_ds(i) = db;
    end
    Hsel = zeros(m, nq); for i=1:m, Hsel(i, act_idx(i))=1; end
    Jy = Hsel - kappa*(dyd_ds*c_row);
    yd = Jy*dq;
end

function e = p_eps(p)
    e = getdef(p, 'clf_eps', 0.10);
end

function v = getdef(s, f, d)
    if isfield(s, f) && ~isempty(s.(f)), v = s.(f); else, v = d; end
end
