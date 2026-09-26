function ch3_test_model()
%CH3_TEST_MODEL  Stage-1 verification: control-affine split, guard, Delta.
%
% Checks, against independent references:
%   1. f(x) + g(x)u reproduces the KKT constrained dynamics exactly, for
%      several random u  (reference: Dynamics/rabbit_constrained_dynamics).
%   2. lambda = lam_drift + lam_in*u reproduces the KKT contact force.
%   3. The stance contact constraint J_st*ddq + Jdot*dq = 0 holds.
%   4. Impact enforces J_sw*dq+ = 0 (foot sticks) and leaves q untouched.
%   5. Delta leaves the NEW stance foot exactly on the ground, and preserves
%      px (the translation gauge).
%   6. The guard equals the swing-foot world height.
%   7. ch3_mvg is M.m/V.m/G.m bit for bit by default; p.model_blend reaches
%      the 30 kg model (restored from git) at 0 and blends linearly between;
%      ch3_model_signature / ch3_model_check tell match, mismatch and
%      unrecorded apart.

fprintf('\n=== ch3_test_model ===\n');
p = ch3_params();
tol = 1e-9;
pass = true;

% A physically sensible start-of-step pose (stance foot down, swing behind).
x0 = [ +0.0000; -0.9310; +0.2999; -0.7934; +0.6869; -0.6318; +0.9810; ...
       +0.3952; -0.0419; +0.1847; +0.1847; +0.1045; +0.0000; +0.0000];

rng(0);
states = [x0, x0 + [zeros(7,1); 0.3*randn(7,1)], x0 + 0.05*randn(14,1)];

%% 1-3. control-affine split vs the KKT reference
worst_ddq = 0; worst_lam = 0; worst_con = 0;
for k = 1:size(states,2)
    x = states(:,k);
    q = x(1:7); dq = x(8:14);
    [f, g, aux] = ch3_control_affine(x, p);

    for trial = 1:4
        u = 50*randn(4,1);

        % ours
        xdot_ours = f + g*u;
        ddq_ours  = xdot_ours(8:14);
        lam_ours  = aux.lam_drift + aux.lam_in*u;

        % reference
        [ddq_ref, lam_ref] = rabbit_constrained_dynamics(q, dq, u);

        worst_ddq = max(worst_ddq, norm(ddq_ours - ddq_ref, inf));
        worst_lam = max(worst_lam, norm(lam_ours - lam_ref, inf));

        % the qdot block of f must be dq itself
        worst_ddq = max(worst_ddq, norm(xdot_ours(1:7) - dq, inf));

        % holonomic stance constraint at the acceleration level
        con = J_st(q)*ddq_ours + Jdotdq_st(q, dq);
        worst_con = max(worst_con, norm(con, inf));
    end
end
pass = report('f+g*u == KKT ddq',      worst_ddq, 1e-8, pass);
pass = report('lambda affine in u',    worst_lam, 1e-8, pass);
pass = report('stance constraint = 0', worst_con, 1e-8, pass);

%% 4. impact
x_pre = states(:,2);
[x_plus, impulse] = ch3_impact(x_pre, p);

% foot sticks: evaluate BEFORE relabel, so redo the raw impact here
q_pre = x_pre(1:7);
Mm = M(q_pre); Jsw = J_sw(q_pre);
sol = [Mm, -Jsw.'; Jsw, zeros(2)] \ [Mm*x_pre(8:14); zeros(2,1)];
dqp = sol(1:7);
pass = report('impact: J_sw*dq+ = 0', norm(Jsw*dqp, inf), 1e-9, pass);
pass = report('impulse matches KKT',  norm(impulse - sol(8:end), inf), tol, pass);

%% 5. Delta: new stance foot on the ground, px preserved
foot_new = P_st(x_plus(1:7));
pass = report('Delta: new stance foot z = 0', abs(foot_new(2)), 1e-12, pass);
pass = report('Delta: px preserved',          abs(x_plus(1) - x_pre(1)), tol, pass);

% q continuity through impact (before the re-plant shifts y)
x_relabel_check = ch3_relabel([q_pre; dqp], p);
pass = report('relabel is an involution', ...
              norm(ch3_relabel(x_relabel_check, p) - [q_pre; dqp], inf), tol, pass);

%% 6. guard
h    = ch3_guard(x0, p);
Psw0 = P_sw(x0(1:7));
pass = report('guard == P_sw height', abs(h - Psw0(2)), tol, pass);

% guard gradient vs finite difference
[~, dh] = ch3_guard(x0, p);
fd = zeros(1,7);
for i = 1:7
    dxp = x0; dxp(i) = dxp(i) + 1e-6;
    dxm = x0; dxm(i) = dxm(i) - 1e-6;
    fd(i) = (ch3_guard(dxp,p) - ch3_guard(dxm,p)) / 2e-6;
end
pass = report('guard gradient vs FD', norm(dh - fd, inf), 1e-6, pass);

%% 7. the model hook (ch3_mvg) and the model signature
% blend [] and blend 1 must BE M.m / V.m / G.m, bit for bit: every existing
% result depends on the default path being untouched.
xs = states(:,2); qs = xs(1:7); dqs = xs(8:14);
[M0, V0, G0] = ch3_mvg(qs, dqs, p);
err = norm(M0 - M(qs), inf) + norm(V0 - V(xs), inf) + norm(G0 - G(qs), inf);
pass = report('ch3_mvg default == M/V/G', err, 0, pass);
p1 = p; p1.model_blend = 1;
[M1, V1, G1] = ch3_mvg(qs, dqs, p1);
err = norm(M1 - M0, inf) + norm(V1 - V0, inf) + norm(G1 - G0, inf);
pass = report('blend 1 == today''s model', err, 0, pass);

% blend 0 is the 30 kg model restored from git; the implied mass says so.
pz = p; pz.model_blend = 0;
[Mz, ~, Gz] = ch3_mvg(qs, dqs, pz);
pass = report('blend 0 == M_m30', norm(Mz - M_m30(qs), inf), 1e-13, pass);
pass = report('blend 0 weighs 30 kg', abs(-Gz(2)/p.g0 - 30), 1e-9, pass);

% halfway: the parameters are blended, so everything is the average -- and
% the result is still a mechanical system (M symmetric positive definite).
ph = p; ph.model_blend = 0.5;
[Mh, Vh, Gh] = ch3_mvg(qs, dqs, ph);
Vz = V_m30(xs);
err = norm(Mh - (Mz + M0)/2, inf) + norm(Vh - (Vz + V0)/2, inf) + ...
      norm(Gh - (Gz + G0)/2, inf);
pass = report('blend 0.5 == average', err, 1e-12, pass);
pass = report('blend 0.5 weighs 52 kg', abs(-Gh(2)/p.g0 - 52), 1e-9, pass);
err = max(norm(Mh - Mh.', inf), -min(eig((Mh + Mh.')/2)));
pass = report('blend 0.5: M sym. pos. def.', max(err, 0), 1e-10, pass);

% and the blend reaches the dynamics the simulation integrates
[fh] = ch3_control_affine(xs, ph);
[f0] = ch3_control_affine(xs, p);
ok = norm(fh - f0, inf) > 1e-3;
fprintf('  [%s] %-30s |f(0.5) - f(1)| = %.3e\n', tf(ok), 'blend changes the dynamics', ...
        norm(fh - f0, inf));
pass = pass && ok;

% the signature: same model -> match; a model whose numbers moved -> mismatch;
% a file with none -> unrecorded.
sig = ch3_model_signature(p);
R1  = ch3_model_check(struct('p', p, 'model_sig', sig), 'quiet');
ok  = strcmp(R1.status, 'match') && R1.rel_err == 0;
fprintf('  [%s] %-30s %s (%s, %.1f kg)\n', tf(ok), 'signature matches itself', ...
        R1.status, sig.digest, sig.mass);
pass = pass && ok;
bad = sig; bad.vals = bad.vals * (1 + 1e-6);
R2  = ch3_model_check(struct('p', p, 'model_sig', bad), 'quiet');
ok  = strcmp(R2.status, 'mismatch');
fprintf('  [%s] %-30s %s (rel %.1e)\n', tf(ok), 'perturbed model -> mismatch', ...
        R2.status, R2.rel_err);
pass = pass && ok;
R3  = ch3_model_check(struct('p', p), 'quiet');
ok  = strcmp(R3.status, 'unrecorded');
fprintf('  [%s] %-30s %s\n', tf(ok), 'no signature -> unrecorded', R3.status);
pass = pass && ok;
sz  = ch3_model_signature(pz);
ok  = abs(sz.mass - 30) < 1e-9 && ~strcmp(sz.digest, sig.digest);
fprintf('  [%s] %-30s %s (%.1f kg) vs %s\n', tf(ok), 'blend 0 has its own signature', ...
        sz.digest, sz.mass, sig.digest);
pass = pass && ok;

fprintf('--- ch3_test_model: %s ---\n\n', tf(pass));
end

function ok = report(name, err, tol, ok_in)
ok = ok_in && (err <= tol);
fprintf('  [%s] %-30s err = %.3e (tol %.0e)\n', tf(err <= tol), name, err, tol);
end

function s = tf(b)
if b, s = 'PASS'; else, s = 'FAIL'; end
end
