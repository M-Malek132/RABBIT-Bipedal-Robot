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
% THE BOUNDARY LAYER CHANGES WHAT IS GUARANTEED, so checks 3 and 8 pin
% p.rclf.boundary_layer = 0: they are statements about the exact min-max law of
% (4.12). Checks 9-12 are the corresponding statements about the layer, which
% is what the pipeline actually runs.
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
%   9. with the layer, (4.11) holds exactly outside it and is exceeded by at
%      most qp.bl_gap <= D1*phi/4 inside it, over the same sampled ball
%  10. at the orbit the exact law keeps a correction of norm D1/(1-D2); the
%      layer's vanishes linearly with eta
%  11. under sample-and-hold the exact law chatters and the layer does not,
%      and the layer's error settles inside its ultimate bound
%  12. inside the layer the robust correction no longer grows with D1

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
    pm.rclf.boundary_layer = 0;              % the exact law; the layer is check 9
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
    pi_.rclf.boundary_layer = 0;             % the exact law's price; see 12
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

%% 9. the layer's guarantee, sampled over the same ball
% Each state's eta is also scaled up and down so that both sides of the layer
% are exercised; a check that only ever landed outside it would pass a layer
% that did nothing, and one only inside it would never test the exact region.
for model = {'scalar', 'matrix'}
    pm = p; pm.rclf.delta2_model = model{1};
    D1 = pm.rclf.delta1_max;
    D2 = pm.rclf.delta2_max;

    excess = -inf; gap_over_cap = -inf; n_in = 0; n_out = 0;
    for k = 1:size(XS,2)
        [Lf2y, LgLfy, u_ff, info] = ch4_io_lin(XS(:,k), alpha, pm, []);
        for sc = [10 1 1e-1 1e-2]
            info_s = info; info_s.eta = sc * info.eta;
            [mu, ~, qp] = ch4_ctrl_rclf_qp(Lf2y, LgLfy, u_ff, info_s, pm, false);
            if qp.bl_gap > 0, n_in = n_in + 1; else, n_out = n_out + 1; end
            gap_over_cap = max(gap_over_cap, qp.bl_gap - D1 * qp.phi / 4);
            for j = 1:100
                [d1, d2] = sample_ball(pm.ny, D1, D2, model{1}, j);
                resid = qp.psi + qp.LgV*d1 + qp.LgV*(eye(pm.ny) + d2)*mu;
                excess = max(excess, resid - qp.bl_gap);
            end
        end
    end
    ok = excess <= 1e-6 && gap_over_cap <= 1e-9 && n_in > 0 && n_out > 0;
    fprintf('  [%s] %-30s excess %.1e, gap-cap %.1e, in/out %d/%d\n', tf(ok), ...
            sprintf('layer guarantee (%s)', model{1}), max(excess, 0), ...
            max(gap_over_cap, 0), n_in, n_out);
    pass = pass && ok;
end

%% 10. what the layer removes: the fixed-magnitude term at the orbit
% Along a ray into eta = 0, psi shrinks like t^2 and ||LgV|| like t. The exact
% law's ||mu|| therefore tends to M = D1/(1-D2) -- a finite correction applied
% at an arbitrarily small error, whose direction flips with the sign of eta.
% Inside the layer a is O(t^2) too, so the layer's ||mu|| is exactly linear in t.
[Lf2y, LgLfy, u_ff, info] = ch4_io_lin(XS(:,1), alpha, p, []);
M  = p.rclf.delta1_max / (1 - p.rclf.delta2_max);
pe = p; pe.rclf.boundary_layer = 0;
ts = [1e-2 1e-4 1e-6];
n_exact = zeros(size(ts)); n_layer = zeros(size(ts));
for i = 1:numel(ts)
    info_t = info; info_t.eta = ts(i) * info.eta;
    n_exact(i) = norm(ch4_ctrl_rclf_qp(Lf2y, LgLfy, u_ff, info_t, pe, false));
    n_layer(i) = norm(ch4_ctrl_rclf_qp(Lf2y, LgLfy, u_ff, info_t, p,  false));
end
slope = n_layer ./ ts;
ok = abs(n_exact(end) / M - 1) < 1e-3 ...
     && n_layer(end) < 1e-6 * M ...
     && (max(slope) - min(slope)) <= 1e-6 * max(max(slope), realmin);
fprintf('  [%s] %-30s exact %s, layer %s (M = %.1f)\n', tf(ok), ...
        'no fixed term at the orbit', mat2str(n_exact, 4), ...
        mat2str(n_layer, 3), M);
pass = pass && ok;

%% 11. sample-and-hold: the exact law chatters, the layer does not
% The transverse dynamics are exactly etadot = F eta + G((1 + d2) mu + Delta1),
% so the sampled loop is simulated here without the robot: zero-order hold at
% p.control_dt, d2 at 0.9 D2 -- the high input gain, where overshoot is worst --
% and two disturbances, judged over the second half of a 2 s run by the change
% in mu from one sample to the next:
%
%   Delta1 = 0                  Case I. The exact law still flips a correction
%                               of size M every sample; the layer settles to
%                               eta = 0 with mu constant.
%   Delta1 = 0.2 D1 sin(3 Hz)   a disturbance inside the ball that moves on the
%                               timescale of a step. The layer's error must stay
%                               inside its ultimate bound eps*D1*phi/(4*c3).
%
% Chatter is a sustained train of large jumps, so it is counted: samples whose
% jump exceeds 10% of M. A max would also flag the min-norm term's isolated kick
% where the sinusoid changes sign (measured: one of 22 rad/s^2 at eps = 0.20,
% against a median of 0.34), which is not chatter and not the robust term.
%
% WHAT THIS DOES NOT CLAIM. A CONSTANT Delta1 pushed against the high gain still
% chatters with the layer (measured at eps = 0.35: 0.8 D1 gives a median |dmu|
% of 277 against 815 without it). That is not the robust term either. Such a
% disturbance settles where ydot = 0 and y ~= 0, and there the Chapter-3
% min-norm term's own gain along LgV, psi/((1-D2)||LgV||^2), is of order 1e4 --
% far past what a 1 ms hold carries. The layer does not touch that term.
T  = p.control_dt;
ny = p.ny;
Ad = [eye(ny), T*eye(ny); zeros(ny), eye(ny)];
Bd = [T^2/2 * eye(ny); T * eye(ny)];
rng(4);
v    = randn(ny, 1);  v = v / norm(v);
gain = 1 + 0.9 * p.rclf.delta2_max;
K    = 2000;
dist = {@(k) zeros(ny, 1), ...
        @(k) 0.2 * p.rclf.delta1_max * sin(2*pi*3*k*T) * v};

kaps  = [0, p.rclf.boundary_layer];
dmu_med = zeros(2, 2); dmu_max = zeros(2, 2); n_big = zeros(2, 2);
V_end = zeros(2, 2); V_hi = zeros(2, 2); V_0 = 0; phi_l = 0;
for id = 1:2
    for ik = 1:2
        pk  = p; pk.rclf.boundary_layer = kaps(ik);
        eta = [0.02 * ones(ny,1); zeros(ny,1)];
        mu_prev = zeros(ny,1); dm = zeros(1, K/2);
        for k = 1:K
            [mu, ~, qp] = ch4_ctrl_rclf_qp(Lf2y, LgLfy, u_ff, ...
                                           struct('eta', eta), pk, false);
            if k == 1, V_0 = qp.V; end
            if k > K/2
                dm(k - K/2) = norm(mu - mu_prev, inf);
                V_hi(id, ik) = max(V_hi(id, ik), qp.V);
            end
            eta = Ad * eta + Bd * (gain * mu + dist{id}(k));
            mu_prev = mu;
        end
        dmu_med(id, ik) = median(dm);
        dmu_max(id, ik) = max(dm);
        n_big(id, ik)   = sum(dm > 0.1 * M);
        V_end(id, ik)   = qp.V;
        if kaps(ik) > 0, phi_l = qp.phi; end
    end
end
V_ult = p.eps * p.rclf.delta1_max * phi_l / (4 * clf.c3);
ok = all(n_big(:, 1) >= K/4) ...                           % exact law chatters
     && dmu_max(1, 2) <= 1e-3 * M && V_end(1, 2) <= 1e-3 * V_0 ...
     && n_big(2, 2) == 0 && dmu_med(2, 2) <= 5e-3 * M ...
     && V_hi(2, 2) <= V_ult;
fprintf(['  [%s] %-30s jumps > 0.1 M: %d/%d -> %d/%d of %d, ' ...
         'V %.1e <= %.1e\n'], tf(ok), 'sampled loop: no chatter', ...
        n_big(1,1), n_big(2,1), n_big(1,2), n_big(2,2), K/2, V_hi(2,2), V_ult);
pass = pass && ok;

%% 12. inside the layer a larger bound does not buy a larger control
% phi grows in proportion to D1, so D1*||LgV||^2/phi -- and with it mu -- is
% independent of D1 while ||LgV|| stays inside the layer. The exact law's
% correction instead grows with D1: that is check 8's price, which the layer
% caps near the orbit.
[Lf2y, LgLfy, u_ff, info] = ch4_io_lin(XS(:,1), alpha, p, []);
info_s = info; info_s.eta = 1e-3 * info.eta;
d1_levels = [0.25 0.5 1] * p.rclf.delta1_max;
n_l = zeros(size(d1_levels)); n_e = zeros(size(d1_levels)); inside = true;
for i = 1:numel(d1_levels)
    pl_ = p;  pl_.rclf.delta1_max = d1_levels(i);
    pe_ = pl_; pe_.rclf.boundary_layer = 0;
    [m_l, ~, q_l] = ch4_ctrl_rclf_qp(Lf2y, LgLfy, u_ff, info_s, pl_, false);
    n_l(i) = norm(m_l);
    n_e(i) = norm(ch4_ctrl_rclf_qp(Lf2y, LgLfy, u_ff, info_s, pe_, false));
    inside = inside && q_l.bl_gap > 0;
end
ok = inside && all(n_l > 0) && (max(n_l) - min(n_l)) <= 1e-9 * max(n_l) ...
     && all(diff(n_e) > 0);
fprintf('  [%s] %-30s layer %s, exact %s\n', tf(ok), ...
        'layer caps the price', mat2str(n_l, 4), mat2str(n_e, 4));
pass = pass && ok;

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
