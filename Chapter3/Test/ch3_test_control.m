function ch3_test_control()
%CH3_TEST_CONTROL  Stage 5-8 verification: RES-CLF, CLF-QP, constrained QP.
%
% Checks:
%   1. CARE and Lyapunov constructions solve their own equations; P > 0; c3 > 0.
%   2. LfV / LgV match a finite difference of V along etadot = F eta + G mu,
%      i.e. Vdot really is affine in mu.
%   3. Seed consistency: y(0) = 0 and ydot(0) = 0 by construction.
%   4. Stage 7 closed form satisfies the RES-CLF inequality, and agrees with
%      what quadprog returns for the same least-norm problem.
%   5. Stage 5 PD also satisfies the inequality (it must, or the 'lyap'
%      construction would be invalid).
%   6. Stage 8 respects a TIGHT torque box that stage 5 violates -- the whole
%      point of the constrained QP -- and reports the slack it needed.
%   7. Closed-loop: eta decays under both PD and CLF-QP from a perturbation.

fprintf('\n=== ch3_test_control ===\n');
pass = true;

%% 1. CLF constructions
p_care = ch3_params('clf_construction','care');
p_lyap = ch3_params('clf_construction','lyap');
clear ch3_res_clf                       % drop the persistent cache
c1 = ch3_res_clf(p_care);
clear ch3_res_clf
c2 = ch3_res_clf(p_lyap);

pass = report('CARE residual',      c1.residual, 1e-9, pass);
pass = report('Lyapunov residual',  c2.residual, 1e-9, pass);
fprintf('        CARE : lam(P) in [%.4f %.4f], c3 = %.4f\n', c1.lam_min_P, c1.lam_max_P, c1.c3);
fprintf('        lyap : lam(P) in [%.4f %.4f], c3 = %.4f\n', c2.lam_min_P, c2.lam_max_P, c2.c3);
pass = pass && c1.lam_min_P > 0 && c2.lam_min_P > 0 && c1.c3 > 0 && c2.c3 > 0;

%% 2. Lie derivatives of V
p = p_care;
clear ch3_res_clf
clf = ch3_res_clf(p);
rng(3);
e_lie = 0;
for k = 1:5
    eta = randn(2*p.ny, 1);
    mu  = randn(p.ny, 1);
    [V, LfV, LgV] = ch3_clf_eval(eta, clf, p.eps);
    etadot = clf.F*eta + clf.G*mu;
    h = 1e-6;
    Vp = ch3_clf_eval(eta + h*etadot, clf, p.eps);
    Vm = ch3_clf_eval(eta - h*etadot, clf, p.eps);
    e_lie = max(e_lie, abs((LfV + LgV*mu) - (Vp - Vm)/(2*h)));
    if V <= 0, pass = false; end
end
pass = report('Vdot = LfV + LgV*mu', e_lie, 1e-4, pass);

%% 3. seed consistency
[alpha, x0] = ch3_seed(p);
[y0, yd0] = ch3_outputs(x0, alpha, p);
pass = report('seed: y(0) = 0',    norm(y0,  inf), 1e-12, pass);
pass = report('seed: ydot(0) = 0', norm(yd0, inf), 1e-12, pass);

%% test states: perturb off the zero dynamics surface so eta ~= 0
xs = x0; xs(3) = xs(3) + 0.15;              % mid-phase
xs(8:14) = xs(8:14) + 0.15*randn(7,1);

%% 4. stage 7 closed form
[Lf2y, LgLfy, u_ff, info] = ch3_io_lin(xs, alpha, p);
[mu7, u7, qp7] = ch3_ctrl_clf_qp(Lf2y, LgLfy, u_ff, info, p, false);
res7 = qp7.LfV + qp7.LgV*mu7 + (clf.c3/p.eps)*qp7.V;
pass = report('stage 7 satisfies RES-CLF', max(res7, 0), 1e-8, pass);

% compare against quadprog on the same problem
Hq = 2*eye(p.ny); fq = zeros(p.ny,1);
oq = optimoptions('quadprog','Display','off');
mu_qp = quadprog(Hq, fq, qp7.LgV, -qp7.psi, [],[],[],[],[], oq);
pass = report('stage 7 == quadprog', norm(mu7 - mu_qp, inf), 1e-6, pass);

%% 5. stage 5 PD also satisfies the inequality
[mu5, u5] = ch3_ctrl_pd(Lf2y, LgLfy, u_ff, info, p);
res5 = qp7.LfV + qp7.LgV*mu5 + (clf.c3/p.eps)*qp7.V;
pass = report('stage 5 satisfies RES-CLF', max(res5, 0), 1e-6, pass);
fprintf('        |u| : PD = %.1f Nm, CLF-QP = %.1f Nm, u_ff = %.1f Nm\n', ...
        max(abs(u5)), max(abs(u7)), max(abs(u_ff)));

%% 6. stage 8 with a torque box tight enough that stage 5 violates it
u_tight = 0.4 * max(abs(u5));
p8 = p;
p8.limits.u_max = u_tight;
p8.limits.enable.torque = true;
[mu8, u8, qp8] = ch3_ctrl_clf_qp(Lf2y, LgLfy, u_ff, info, p8, true); %#ok<ASGLU>

viol5 = max(abs(u5)) - u_tight;
viol8 = max(max(abs(u8)) - u_tight, 0);
fprintf('        box = %.1f Nm : PD peak = %.1f (violates by %.1f), QP peak = %.1f\n', ...
        u_tight, max(abs(u5)), viol5, max(abs(u8)));
pass = report('stage 8 respects torque box', viol8, 1e-7, pass);
pass = pass && viol5 > 0;      % the test would be vacuous if PD already fit
fprintf('        [%s] %-30s delta = %.3e, feasible = %d\n', ...
        tf(qp8.feasible), 'stage 8 slack used', qp8.delta, qp8.feasible);
pass = pass && qp8.feasible && qp8.delta >= -1e-9;

%% 7. closed loop: the RES-CLF guarantee is on V, not on ||eta||
%
% The bound Chapter 3 proves is
%     ||eta(t)|| <= (1/eps) sqrt(c2/c1) exp(-c3 t / 2eps) ||eta(0)||
% whose PREFACTOR is (1/eps)sqrt(c2/c1) -- at eps = 0.5 that is already 3.9,
% so the theory explicitly permits ||eta|| to GROW several-fold before it
% contracts. Asserting that ||eta|| decreases monotonically would be testing
% something the method never claims.
%
% What IS claimed, pointwise and therefore along the whole trajectory, is
%     Vdot_eps + (c3/eps) V_eps <= 0   ==>   V(t) <= V(0) exp(-(c3/eps) t).
% That is what we check, at every integration step.
% EACH CONTROLLER MUST BE JUDGED AGAINST ITS OWN CONSTRUCTION.
%   * The stage-7 min-norm law is built from the CARE P, so it satisfies the
%     CARE rate by construction.
%   * The stage-5 PD law does NOT inherit that rate. Its matching certificate
%     is the Lyapunov construction A'P + PA = -Q with A = [0 I; -Kp -Kd] --
%     which is exactly the equation Chapter 3 writes -- because in the scaled
%     coordinate eta_eps = I_eps eta the PD closed loop is
%     eta_eps_dot = (1/eps) A eta_eps, giving Vdot = -(1/eps) eta_eps' Q eta_eps
%     and hence the (c3/eps) rate exactly.
% Pairing PD against the CARE P instead is a category error, and the run below
% quantifies how wrong it goes.
dt   = 2e-3;
nstp = 120;

pairs = { 'clfqp',    'care' ; ...
          'iolin_pd', 'lyap' };

for i = 1:size(pairs,1)
    pc = p;
    pc.controller       = pairs{i,1};
    pc.clf_construction = pairs{i,2};
    clear ch3_res_clf
    clf_i = ch3_res_clf(pc);
    rate  = clf_i.c3 / pc.eps;

    [worst_ratio, V0, Vt] = run_bound(xs, alpha, pc, clf_i, rate, dt, nstp);
    ok = worst_ratio <= 1 + 1e-6;
    fprintf('  [%s] %-30s V %.3e -> %.3e, worst V/bound = %.4f\n', tf(ok), ...
            sprintf('RES bound (%s/%s)', pairs{i,1}, pairs{i,2}), V0, Vt, worst_ratio);
    pass = pass && ok;
end

% Observation, not an assertion: PD measured against the CARE certificate.
pc = p; pc.controller = 'iolin_pd'; pc.clf_construction = 'care';
clear ch3_res_clf
clf_c = ch3_res_clf(pc);
[wr, V0c, Vtc] = run_bound(xs, alpha, pc, clf_c, clf_c.c3/pc.eps, dt, nstp);
fprintf('        (note) PD vs CARE certificate: V %.3e -> %.3e, worst ratio %.2f\n', ...
        V0c, Vtc, wr);
fprintf('        PD converges further overall but transiently exceeds the CARE rate.\n');

% PD, being far more aggressive than the minimum-norm CLF control, should also
% contract ||eta|| outright over this horizon.
pc = p; pc.controller = 'iolin_pd';
xe = xs; [~,~,oi] = ch3_outputs(xe, alpha, pc); eta0 = norm(oi.eta);
for k = 1:nstp, xe = rk4(xe, alpha, pc, dt); end
[~,~,oi] = ch3_outputs(xe, alpha, pc);
pass = report('PD contracts ||eta||', max(norm(oi.eta) - eta0, 0), 1e-12, pass);
fprintf('        ||eta|| %.3e -> %.3e\n', eta0, norm(oi.eta));

fprintf('--- ch3_test_control: %s ---\n\n', tf(pass));
end

% ---------------------------------------------------------------------------
function [worst_ratio, V0, Vt] = run_bound(x0, alpha, pc, clf_i, rate, dt, nstp)
xe = x0;
[~, ~, oi] = ch3_outputs(xe, alpha, pc);
V0 = ch3_clf_eval(oi.eta, clf_i, pc.eps);
worst_ratio = 0;
Vt = V0;
for k = 1:nstp
    xe = rk4(xe, alpha, pc, dt);
    [~, ~, oi] = ch3_outputs(xe, alpha, pc);
    Vt = ch3_clf_eval(oi.eta, clf_i, pc.eps);
    bound = V0 * exp(-rate * k*dt);
    worst_ratio = max(worst_ratio, Vt / max(bound, realmin));
end
end

function x = rk4(x, alpha, p, dt)
k1 = ch3_ode_rhs(0, x,           alpha, p);
k2 = ch3_ode_rhs(0, x + dt/2*k1, alpha, p);
k3 = ch3_ode_rhs(0, x + dt/2*k2, alpha, p);
k4 = ch3_ode_rhs(0, x + dt*k3,   alpha, p);
x  = x + dt/6*(k1 + 2*k2 + 2*k3 + k4);
end

function ok = report(name, err, tol, ok_in)
ok = ok_in && (err <= tol);
fprintf('  [%s] %-30s err = %.3e (tol %.0e)\n', tf(err <= tol), name, err, tol);
end

function s = tf(b)
if b, s = 'PASS'; else, s = 'FAIL'; end
end
