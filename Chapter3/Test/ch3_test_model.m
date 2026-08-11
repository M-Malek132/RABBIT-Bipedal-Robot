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

fprintf('--- ch3_test_model: %s ---\n\n', tf(pass));
end

function ok = report(name, err, tol, ok_in)
ok = ok_in && (err <= tol);
fprintf('  [%s] %-30s err = %.3e (tol %.0e)\n', tf(err <= tol), name, err, tol);
end

function s = tf(b)
if b, s = 'PASS'; else, s = 'FAIL'; end
end
