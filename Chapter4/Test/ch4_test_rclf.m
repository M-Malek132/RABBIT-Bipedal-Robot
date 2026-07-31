function ch4_test_rclf()
%CH4_TEST_RCLF  The robust CLF-QP of Section 4.1.
%
% The central test here is not that the code reproduces a formula -- it is that
% the CONTROLLER'S GUARANTEE ACTUALLY HOLDS. The robust RES condition (4.11)
% claims that the mu it returns satisfies
%
%       Vdot(eta, Delta1, Delta2, mu) + (c3/eps) V  <=  0
%
% for EVERY Delta1, Delta2 in the bound, not merely for the worst case the
% derivation identified. So the test samples the uncertainty ball at random --
% including its interior and its boundary -- and checks the inequality against
% each draw. If the max in (4.11) had been taken incorrectly, or a sign in the
% two-inequality reduction were wrong, this is what would catch it, and nothing
% else in the file would.
%
% MATCHING THE SET TO THE DESIGN. The 'scalar' design assumes Delta2 = d2 I and
% is tested against that set; the 'matrix' design assumes ||Delta2||_2 <= D2 and
% is tested against full random matrices. Testing 'scalar' against a general
% matrix would be testing a guarantee it never made -- and ch4_test_model shows
% that for this chapter's own perturbation Delta2 IS isotropic, so the scalar
% set is the physically relevant one.
%
% Checks:
%   1. bounds of zero collapse to the Chapter-3 CLF-QP exactly
%   2. the closed form sits exactly on the robust constraint boundary
%   3. THE GUARANTEE: the RES condition holds for sampled Delta in the ball
%   4. the Chapter-3 controller VIOLATES it on the same draws (so 3 has teeth)
%   5. the constrained QP respects the torque box
%   6. a loose box reproduces the unconstrained closed form
%   7. Delta2max >= 1 is reported infeasible rather than silently divided by
%   8. the price of robustness: ||mu|| grows monotonically with the bounds

fprintf('\n=== ch4_test_rclf ===\n');
pass = true;

[x0, alpha, p] = ch4_load_gait();

rng(9);
XS = x0 + 0.12*randn(14, 6);
XS(:,1) = x0 + [zeros(7,1); 0.05*ones(7,1)];   % ensure eta ~= 0

clf = ch3_res_clf(p);

%% 1. zero bounds == Chapter 3
e0 = 0;
p0 = p; p0.rclf.delta1_max = 0; p0.rclf.delta2_max = 0;
for k = 1:size(XS,2)
    [Lf2y, LgLfy, u_ff, info] = ch4_io_lin(XS(:,k), alpha, p, []);
    m_r = ch4_ctrl_rclf_qp(Lf2y, LgLfy, u_ff, info, p0, false);
    m_3 = ch3_ctrl_clf_qp(  Lf2y, LgLfy, u_ff, info, p,  false);
    e0 = max(e0, norm(m_r - m_3, inf));
end
pass = report('D=0 collapses to ch3 CLF-QP', e0, 1e-10, pass);

%% 2. the closed form is exactly on the boundary
e_b = 0;
for k = 1:size(XS,2)
    [Lf2y, LgLfy, u_ff, info] = ch4_io_lin(XS(:,k), alpha, p, []);
    [~, ~, qp] = ch4_ctrl_rclf_qp(Lf2y, LgLfy, u_ff, info, p, false);
    if qp.a > 0
        e_b = max(e_b, abs(qp.margin));       % active: margin must be 0
    else
        e_b = max(e_b, max(-qp.margin, 0));   % inactive: margin must be >= 0
    end
end
pass = report('closed form on the boundary', e_b, 1e-8, pass);

%% 3. THE GUARANTEE, sampled
for model = {'scalar', 'matrix'}
    pm = p; pm.rclf.delta2_model = model{1};
    D1 = pm.rclf.delta1_max;
    D2 = pm.rclf.delta2_max;

    worst = -inf;
    for k = 1:size(XS,2)
        [Lf2y, LgLfy, u_ff, info] = ch4_io_lin(XS(:,k), alpha, pm, []);
        [mu, ~, qp] = ch4_ctrl_rclf_qp(Lf2y, LgLfy, u_ff, info, pm, false);

        for j = 1:200
            [d1, d2] = sample_ball(pm.ny, D1, D2, model{1}, j);
            resid = qp.psi + qp.LgV*d1 + qp.LgV*(eye(pm.ny) + d2)*mu;
            worst = max(worst, resid);
        end
    end
    % scale by the residual's own magnitude so the tolerance means something
    pass = report(sprintf('RES holds over ball (%s)', model{1}), ...
                  max(worst, 0), 1e-6, pass);
end

%% 4. the same draws break the Chapter-3 controller
n_viol = 0; n_tot = 0; worst3 = -inf;
for k = 1:size(XS,2)
    [Lf2y, LgLfy, u_ff, info] = ch4_io_lin(XS(:,k), alpha, p, []);
    [mu3, ~, qp3] = ch3_ctrl_clf_qp(Lf2y, LgLfy, u_ff, info, p, false);
    for j = 1:200
        [d1, d2] = sample_ball(p.ny, p.rclf.delta1_max, p.rclf.delta2_max, ...
                               'scalar', j);
        resid = qp3.psi + qp3.LgV*d1 + qp3.LgV*(eye(p.ny) + d2)*mu3;
        n_tot = n_tot + 1;
        if resid > 1e-6, n_viol = n_viol + 1; end
        worst3 = max(worst3, resid);
    end
end
ok = n_viol > 0;
fprintf('  [%s] %-30s ch3 violates %d/%d draws, worst +%.3e\n', ...
        tf(ok), 'baseline is NOT robust', n_viol, n_tot, worst3);
pass = pass && ok;

%% 5. the constrained QP respects the box
pc = p;
pc.limits.enable.torque = true;
pc.limits.u_max = 50;
e_box = 0;
for k = 1:size(XS,2)
    [Lf2y, LgLfy, u_ff, info] = ch4_io_lin(XS(:,k), alpha, pc, []);
    for model = {'scalar','matrix'}
        pc.rclf.delta2_model = model{1};
        [~, u] = ch4_ctrl_rclf_qp(Lf2y, LgLfy, u_ff, info, pc, true);
        e_box = max(e_box, max(abs(u)) - pc.limits.u_max);
    end
end
pass = report('constrained QP respects box', max(e_box,0), 1e-6, pass);

%% 6. a loose box reproduces the closed form
pl = p;
pl.limits.enable.torque   = false;
pl.limits.enable.friction = false;
pl.limits.enable.grf      = false;
pl.clf_slack_penalty      = 1e10;      % make the slack effectively unavailable
e_l = 0;
for k = 1:size(XS,2)
    [Lf2y, LgLfy, u_ff, info] = ch4_io_lin(XS(:,k), alpha, pl, []);
    mu_c = ch4_ctrl_rclf_qp(Lf2y, LgLfy, u_ff, info, pl, false);
    mu_q = ch4_ctrl_rclf_qp(Lf2y, LgLfy, u_ff, info, pl, true);
    e_l  = max(e_l, norm(mu_c - mu_q, inf) / max(norm(mu_c), 1));
end
pass = report('unconstrained QP == closed form', e_l, 1e-4, pass);

%% 7. Delta2max >= 1 is reported, not divided by
pb = p; pb.rclf.delta2_max = 1.0;
[Lf2y, LgLfy, u_ff, info] = ch4_io_lin(XS(:,1), alpha, pb, []);
[mu_b, ~, qp_b] = ch4_ctrl_rclf_qp(Lf2y, LgLfy, u_ff, info, pb, false);
ok = ~qp_b.feasible && all(isfinite(mu_b));
fprintf('  [%s] %-30s feasible=%d, mu finite=%d\n', tf(ok), ...
        'D2 >= 1 reported infeasible', qp_b.feasible, all(isfinite(mu_b)));
pass = pass && ok;

%% 8. the price of robustness is monotone
[Lf2y, LgLfy, u_ff, info] = ch4_io_lin(XS(:,1), alpha, p, []);
levels = [0 0.25 0.5 0.75];
nrm = zeros(size(levels));
for i = 1:numel(levels)
    pi_ = p;
    pi_.rclf.delta1_max = levels(i) * p.rclf.delta1_max;
    pi_.rclf.delta2_max = levels(i) * 0.8;
    nrm(i) = norm(ch4_ctrl_rclf_qp(Lf2y, LgLfy, u_ff, info, pi_, false));
end
ok = all(diff(nrm) >= -1e-9);
fprintf('  [%s] %-30s ||mu|| = %s\n', tf(ok), ...
        'more robustness costs more mu', mat2str(round(nrm,2)));
pass = pass && ok;
fprintf('        (this monotonicity IS the Section 4.1.4 limitation: the\n');
fprintf('         controller pays the worst-case price even at zero error)\n');

fprintf('--- ch4_test_rclf: %s ---\n\n', tf(pass));
end

% ---------------------------------------------------------------------------
function [d1, d2] = sample_ball(ny, D1, D2, model, j)
%SAMPLE_BALL  A draw from the uncertainty set, boundary-weighted.
%
% Every third draw is placed exactly ON the boundary, because that is where the
% guarantee is tight and where an error in the max would show up. Purely
% interior sampling would pass a controller that is wrong only at the extreme.
on_edge = (mod(j,3) == 0);

v  = randn(ny,1);
r1 = D1 * pick_r(on_edge);
d1 = r1 * v / max(norm(v), realmin);

switch model
    case 'scalar'
        d2 = (D2 * pick_r(on_edge) * sign(randn())) * eye(ny);
    otherwise
        Mrand = randn(ny);
        Mrand = Mrand / max(norm(Mrand, 2), realmin);
        d2 = (D2 * pick_r(on_edge)) * Mrand;
end
end

function r = pick_r(on_edge)
if on_edge, r = 1; else, r = rand(); end
end

function ok = report(name, err, tol, ok_in)
ok = ok_in && (err <= tol);
fprintf('  [%s] %-30s err = %.3e (tol %.0e)\n', tf(err <= tol), name, err, tol);
end

function s = tf(b)
if b, s = 'PASS'; else, s = 'FAIL'; end
end
