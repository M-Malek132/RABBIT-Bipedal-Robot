function pass = ch6_test_qp()
%CH6_TEST_QP  The Chapter-6 quadratic program, (6.12) / (6.26).
%
% Checks:
%   1. with no barriers the QP reduces EXACTLY to Chapter 3's stage-8 QP
%   2. the returned u satisfies every row it was given -- barrier rows are hard,
%      so a "solution" that violates one is worse than an infeasible report
%   3. the barrier row is what it claims: recomputing gddot from the returned u
%      reproduces the margin the QP reports
%   4. the CLF row is the only one that slacks -- forcing a conflict with a tight
%      torque box shows delta > 0 while the barrier margin stays >= 0
%   5. a position constraint has L_g g == 0, so the Section 5.1 reciprocal CBF
%      cannot be written for it -- the motivation for (6.1) is measured here,
%      not asserted
%   6. an infeasible barrier/torque pair is REPORTED, not silently relaxed
%
% See also CH6_CTRL_CBF_CLF_QP, CH6_CBF_ROW, CH3_CTRL_CLF_QP.

fprintf('\n=== ch6_test_qp ===\n');
pass = true;
rng(6);

p = ch6_params();
p.limits.enable.torque = false;

%% ------------------------------------- 1. nb = 0 reduces to ch3 stage 8
q = p;  q.cbf.problem = 'none';
err = 0;
for i = 1:5
    x = sample_state();
    alpha = 0.2*randn(p.ny, p.n_ctrl);
    [Lf2y, LgLfy, u_ff, info] = ch3_io_lin(x, alpha, p);

    u6      = ch6_ctrl_cbf_clf_qp(Lf2y, LgLfy, u_ff, info, ...
                                  ch6_barrier(x, info.aux, q, 0), q);
    [~, u3] = ch3_ctrl_clf_qp(Lf2y, LgLfy, u_ff, info, p, true);
    err = max(err, norm(u6 - u3) / max(1, norm(u3)));
end
ok = err < 1e-6;
fprintf('  [%s] %-42s rel diff %.2e\n', tf(ok), ...
        'nb = 0 == ch3_ctrl_clf_qp(constrained)', err);
pass = pass && ok;

%% ---------------- 1b. reference 'pd' is a SAFETY FILTER on the PD law
% The property the default rests on: where NO ROW IS BINDING AT u_pd, the QP
% must return u_pd exactly. If it does not, it is not a filter on PD, it is a
% different controller that happens to carry barrier rows -- and the 10-step
% stability that motivated the default would be an accident.
%
% "No row binding" has to include the CLF row, not just the barriers. PD does
% not inherit the RES-CLF rate (Chapter 3 says so explicitly and measures PD
% exceeding the CARE bound by 3.5x transiently), so at a state where the CLF row
% is violated at u_pd the QP is SUPPOSED to move away from PD -- that is the row
% doing its job. Testing only the barriers would therefore fail for a correct
% controller, which is what the first version of this check did.
q = ch6_params();
q.cbf.problem = 'stones';
q.stone = ch6_resolve_stone(q.stones, 0.05, 1.20);   % a window nothing can bind
clf = ch3_res_clf(q);
err = 0;  n_slack = 0;
for i = 1:40
    x = sample_state();
    alpha = 0.2*randn(q.ny, q.n_ctrl);

    [Lf2y, LgLfy, u_ff, info] = ch3_io_lin(x, alpha, q);
    [~, u_pd] = ch3_ctrl_pd(Lf2y, LgLfy, u_ff, info, q);

    [V, LfV, LgV] = ch3_clf_eval(info.eta, clf, q.eps);  %#ok<ASGLU>
    psi = LfV + (clf.c3 / q.eps) * V;
    clf_slack_at_pd = -psi - LgV * (LgLfy * (u_pd - u_ff));

    B = ch6_barrier(x, info.aux, q, 0);
    bar_slack = min(arrayfun(@(b) ch6_cbf_row(b, q).b - ...
                                  ch6_cbf_row(b, q).A * u_pd, B));

    if clf_slack_at_pd > 1e-6 && bar_slack > 1e-6 && ...
            max(abs(u_pd)) < q.limits.u_max
        n_slack = n_slack + 1;
        u_qp = ch6_control(0, x, alpha, q);
        err  = max(err, norm(u_qp - u_pd) / max(1, norm(u_pd)));
    end
end
ok = n_slack > 0 && err < 1e-6;
fprintf('  [%s] %-42s rel diff %.2e over %d fully slack states\n', tf(ok), ...
        'reference pd, no row binding => u == u_PD', err, n_slack);
pass = pass && ok;

%% ------------------------------------ 2 & 3. the rows the QP was given hold
worst_row = inf; worst_recompute = 0; n_feas = 0;
for i = 1:12
    x = sample_state();
    alpha = 0.2*randn(p.ny, p.n_ctrl);
    [Lf2y, LgLfy, u_ff, info] = ch3_io_lin(x, alpha, p);
    B = ch6_barrier(x, info.aux, p, 0);
    [u, qp] = ch6_ctrl_cbf_clf_qp(Lf2y, LgLfy, u_ff, info, B, p);

    if ~qp.feasible, continue; end
    n_feas = n_feas + 1;

    for j = 1:numel(B)
        row = ch6_cbf_row(B(j), p);
        if ~row.live, continue; end
        slack = row.b - row.A*u;
        worst_row = min(worst_row, slack / max(1, abs(row.b)));

        % Recompute the ECBF statement directly from the applied u:
        %    gddot + Kb eta_b >= 0
        gdd = B(j).Lf2g + B(j).LgLfg * u;
        ref = gdd + row.Kb * B(j).eta_b;
        worst_recompute = max(worst_recompute, abs(ref - qp.margin(j)));
    end
end
ok = n_feas > 0 && worst_row > -1e-7;
fprintf('  [%s] %-42s worst slack %.2e over %d solves\n', tf(ok), ...
        'returned u satisfies every barrier row', worst_row, n_feas);
pass = pass && ok;

ok = worst_recompute < 1e-8;
fprintf('  [%s] %-42s max diff %.2e\n', tf(ok), ...
        'reported margin == gddot + Kb eta_b', worst_recompute);
pass = pass && ok;

%% -------------------------------- 4. only the CLF row bends, and it does bend
% Squeeze the torque box until the CLF row cannot be met, and check that the
% barrier margin stays non-negative while delta goes positive. This is the
% behaviour (6.26) specifies and the thing that would be silently wrong if the
% slack were attached to the wrong row.
q = ch6_params();
q.limits.enable.torque = true;
q.limits.u_max = 12;                 % well below what the gait needs
saw_delta = false; min_margin = inf;
for i = 1:20
    x = sample_state();
    alpha = 0.2*randn(q.ny, q.n_ctrl);
    [Lf2y, LgLfy, u_ff, info] = ch3_io_lin(x, alpha, q);
    B = ch6_barrier(x, info.aux, q, 0);
    [~, qp] = ch6_ctrl_cbf_clf_qp(Lf2y, LgLfy, u_ff, info, B, q);
    if ~qp.feasible, continue; end
    saw_delta = saw_delta || qp.delta > 1e-6;
    m = qp.margin(~isnan(qp.margin));
    if ~isempty(m), min_margin = min(min_margin, min(m)); end
end
ok = saw_delta && min_margin > -1e-6;
fprintf('  [%s] %-42s delta>0 seen: %d, min margin %.2e\n', tf(ok), ...
        'CLF slacks, barriers do not', saw_delta, min_margin);
pass = pass && ok;

%% ---------------------- 5. L_g g == 0: why (6.1) is needed at all
% The Section 5.1 reciprocal barrier's row is  L_f g + L_g g u >= -gamma g^3.
% For a POSITION constraint L_g g is identically zero -- that is what relative
% degree 2 means -- so the row contains no control, and Section 6.1.1 exists to
% get around it. That is measured here rather than asserted, the same evidence
% ch5_test_qp collects for its two plants.
%
% L_g g = (dg/dx) g(x), and g(x) = [0; ddq_in] has a zero configuration block,
% so L_g g vanishes exactly when g does not depend on dq. Perturbing every
% velocity coordinate is a stronger and more direct check than forming the
% product, because it would also catch a barrier that had picked up a velocity
% dependence by accident.
dev = 0;
for i = 1:6
    x  = sample_state();
    g0 = stones_g(x, p);
    for j = p.nq+1:2*p.nq
        xj = x;  xj(j) = xj(j) + 0.5;
        dev = max(dev, max(abs(stones_g(xj, p) - g0)));
    end
end
ok = dev < 1e-14;
fprintf('  [%s] %-42s dg/d(dq) = %.1e  => L_g g = 0\n', tf(ok), ...
        'position constraints have rel. degree 2', dev);
pass = pass && ok;

%% ----------------------------- 6. infeasibility is reported, not relaxed
% A foothold window the gait cannot reach, plus a torque box that leaves no
% authority. The QP must come back infeasible rather than returning something
% that violates a hard row.
r = ch6_params();
r.limits.enable.torque = true;
r.limits.u_max = 0.5;                       % essentially no actuation
r.stones.l_min = 0.55;  r.stones.l_max = 0.56;
r.stone = ch6_resolve_stone(r.stones, 0.55, 0.56);
n_inf = 0; viol = 0;
for i = 1:25
    x = sample_state();
    alpha = 0.2*randn(r.ny, r.n_ctrl);
    [Lf2y, LgLfy, u_ff, info] = ch3_io_lin(x, alpha, r);
    B = ch6_barrier(x, info.aux, r, 0);
    [u, qp] = ch6_ctrl_cbf_clf_qp(Lf2y, LgLfy, u_ff, info, B, r);
    if ~qp.feasible
        n_inf = n_inf + 1;
    else
        for j = 1:numel(B)
            row = ch6_cbf_row(B(j), r);
            if row.live && (row.b - row.A*u) < -1e-6*max(1,abs(row.b))
                viol = viol + 1;
            end
        end
    end
end
ok = viol == 0;
fprintf('  [%s] %-42s %d infeasible, %d silent violations\n', tf(ok), ...
        'infeasible is reported, never relaxed', n_inf, viol);
pass = pass && ok;

fprintf('  --> %s\n', upper(tf(pass)));

end

% ---------------------------------------------------------------------------
function g = stones_g(x, p)
[~, ~, aux] = ch3_control_affine(x, p);
k = ch6_kin(x, aux, p, 'swing');
B = ch6_bar_stones(k, p.stone, 0);
g = [B.g];
end

function x = sample_state()
q  = [0.1*randn; -0.9 + 0.05*randn; 0.2*randn; ...
      0.3*randn; -0.4 - 0.3*rand;  0.3*randn; -0.4 - 0.3*rand];
dq = 0.8*randn(7,1);
x  = [q; dq];
end

function s = tf(b)
if b, s = 'pass'; else, s = 'FAIL'; end
end
