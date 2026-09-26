function ch4_test_model()
%CH4_TEST_MODEL  The true-vs-nominal model layer and the uncertainty terms.
%
% Every check here is against an INDEPENDENT reference, not against the code's
% own output. The strongest ones exploit a fact that is easy to prove by hand
% and hard to satisfy by accident:
%
%   FOR A UNIFORM MASS/INERTIA SCALE s, THE CONSTRAINED DYNAMICS SPLIT EXACTLY.
%
%   The stance KKT system is
%
%       [ sM  -J' ] [ ddq ]   [ -sV - sG  |  B ]
%       [ J    0  ] [ lam ] = [ -Jdotdq   |  0 ]
%
%   Substituting lam = s*w makes the first row divide through by s and leaves
%   the second untouched, so the DRIFT column reduces to the nominal system
%   exactly:  ddq_drift is INDEPENDENT of s. The input columns carry B, which
%   does not scale, so they pick up exactly one factor:  ddq_in = ddq_in_nom/s.
%
%   Therefore, with no approximation whatsoever,
%
%       Lf^2 y  = Lftil^2 y                  (identical)
%       LgLf y  = Lgtil Lftil y / s
%       Delta2  = (1/s - 1) I                exactly a scalar times identity
%       Delta1  = -(1/s - 1) Lftil^2 y       exactly proportional to the drift
%
% Two consequences worth stating plainly, because they shape the whole chapter:
%
%   * Delta2 being EXACTLY isotropic means the 'scalar' reduction of the max in
%     (4.11) -- the one that turns the min-max into two linear inequalities and
%     keeps (4.12) a genuine QP -- is not an approximation for the perturbation
%     this chapter studies. It is exact.
%   * ||Delta2|| = |1/s - 1| < 1 requires s > 0.5. Below half the nominal mass
%     the worst-case model inside the bound can cancel the control authority
%     entirely, and no worst-case design can help. That is a hard limit of
%     Section 4.1, not a tuning problem.
%
% Checks:
%   1. the nominal path IS Chapter 3, bit for bit
%   2. the KKT split above, for several scales
%   3. Remark 4.1: Delta1 = Delta2 = 0 with no uncertainty
%   4. the analytic Delta1, Delta2 above
%   5. eq (4.3): ydd = mu + Delta1 + Delta2 mu on the TRUE plant
%   6. impact: dq+ scale invariant, impulse linear in s
%   7. torso load: the gravity sign against G(q), and that it DOES move dq+
%   8. a load redrawn every step: the analysed contact force uses the load
%      that step carried, not one load for the whole run
%   9. physical validity: the contact force ch4_step records is the true KKT
%      force under the held torque, and ch4_validity scores lift-off and slip
%  10. the structured model: ch4_link_dynamics reproduces M/V/G at today's
%      link parameters; a structured uncertainty at those parameters is
%      Chapter 3; friction and bias enter exactly through B; prior-set draws
%      are reproducible and genuinely non-uniform (Delta2 not a multiple of I)

fprintf('\n=== ch4_test_model ===\n');
pass = true;

[x0, alpha, p] = ch4_load_gait();

rng(4);
XS = x0 + 0.15*randn(14, 8);          % a spread of test states
XS(:,1) = x0;

%% 1. nominal path is Chapter 3 exactly
e_f = 0; e_g = 0; e_x = 0; e_i = 0;
for k = 1:size(XS,2)
    x = XS(:,k);
    [f3, g3] = ch3_control_affine(x, p);
    [f4, g4] = ch4_control_affine(x, p, []);
    e_f = max(e_f, norm(f3-f4, inf));
    e_g = max(e_g, norm(g3(:)-g4(:), inf));

    [xp3, i3] = ch3_impact(x, p);
    [xp4, i4] = ch4_impact(x, p, []);
    e_x = max(e_x, norm(xp3-xp4, inf));
    e_i = max(e_i, norm(i3-i4, inf));
end
pass = report('nominal f == ch3',      e_f, 0, pass);
pass = report('nominal g == ch3',      e_g, 0, pass);
pass = report('nominal impact == ch3', e_x, 0, pass);
pass = report('nominal impulse == ch3',e_i, 0, pass);

%% 2. the KKT split: ddq_drift invariant, ddq_in ~ 1/s
scales = [0.5 0.7 1.5 3];
e_d = 0; e_in = 0;
for s = scales
    unc = struct('mass_scale', s, 'load_mass', 0);
    for k = 1:size(XS,2)
        [~,~,an] = ch4_control_affine(XS(:,k), p, []);
        [~,~,at] = ch4_control_affine(XS(:,k), p, unc);
        e_d  = max(e_d,  norm(at.ddq_drift - an.ddq_drift, inf));
        e_in = max(e_in, norm(at.ddq_in(:) - an.ddq_in(:)/s, inf));
    end
end
pass = report('ddq_drift scale invariant', e_d,  1e-11, pass);
pass = report('ddq_in scales as 1/s',      e_in, 1e-11, pass);

%% 3. Remark 4.1: no uncertainty, no Delta
e0 = 0;
for k = 1:size(XS,2)
    D = ch4_uncertainty(XS(:,k), alpha, p, []);
    e0 = max(e0, max(D.n1, D.n2));
end
pass = report('Remark 4.1: Delta = 0', e0, 1e-9, pass);

%% 4. the analytic Delta1, Delta2
e_D2 = 0; e_D1 = 0;
for s = scales
    unc = struct('mass_scale', s, 'load_mass', 0);
    d2  = 1/s - 1;
    for k = 1:size(XS,2)
        D       = ch4_uncertainty(XS(:,k), alpha, p, unc);
        Lf2y_n  = ch4_io_lin(XS(:,k), alpha, p, []);
        e_D2 = max(e_D2, norm(D.Delta2 - d2*eye(p.ny), inf));
        e_D1 = max(e_D1, norm(D.Delta1 + d2*Lf2y_n,    inf));
    end
end
pass = report('Delta2 = (1/s-1) I',        e_D2, 1e-10, pass);
pass = report('Delta1 = -(1/s-1) Lf2y',    e_D1, 1e-8,  pass);
fprintf('        ||Delta2|| = |1/s-1| : %s for s = %s\n', ...
        mat2str(abs(1./scales - 1), 4), mat2str(scales));
fprintf('        so worst-case robustness (||Delta2||<1) needs s > 0.5\n');

%% 5. eq (4.3) on the true plant
e_43 = 0;
for s = [0.7 1.5 3]
    unc = struct('mass_scale', s, 'load_mass', 0);
    for k = 1:size(XS,2)
        x = XS(:,k);
        [~, LgLfy_n, u_ff_n] = ch4_io_lin(x, alpha, p, []);
        D  = ch4_uncertainty(x, alpha, p, unc);
        mu = randn(p.ny,1) * 5;

        % apply the NOMINAL pre-control (4.1) to the TRUE plant
        u = u_ff_n + LgLfy_n \ mu;
        [Lf2y_t, LgLfy_t] = ch4_io_lin(x, alpha, p, unc);
        ydd_true = Lf2y_t + LgLfy_t * u;

        ydd_43   = mu + D.Delta1 + D.Delta2 * mu;
        e_43 = max(e_43, norm(ydd_true - ydd_43, inf));
    end
end
pass = report('eq (4.3) ydd = mu+D1+D2 mu', e_43, 1e-8, pass);

%% 6. impact under a uniform scale
e_dq = 0; e_imp = 0;
for s = scales
    unc = struct('mass_scale', s, 'load_mass', 0);
    for k = 1:size(XS,2)
        [xp_n, im_n] = ch4_impact(XS(:,k), p, []);
        [xp_s, im_s] = ch4_impact(XS(:,k), p, unc);
        e_dq  = max(e_dq,  norm(xp_s - xp_n, inf));
        e_imp = max(e_imp, norm(im_s - s*im_n, inf));
    end
end
pass = report('impact dq+ scale invariant', e_dq,  1e-10, pass);
pass = report('impact impulse ~ s',         e_imp, 1e-10, pass);

%% 7. torso load
% The load model claims G gains -mL*g0 in the vertical slot. G itself must
% therefore satisfy G(2) = -m_total*g0 for the unloaded robot, independent of
% q -- which pins the sign convention (U = m g0 z with z = -y). Check that
% against the generated G rather than assuming it.
g_at = zeros(1, size(XS,2));
for k = 1:size(XS,2)
    gv = G(XS(1:p.nq,k));
    g_at(k) = gv(2);
end
pass = report('G(2) constant in q', max(abs(g_at - g_at(1))), 1e-10, pass);
m_implied = -g_at(1) / p.g0;
fprintf('        implied total mass from G(2)/g0 = %.4f kg\n', m_implied);

mL  = 12;
unc = struct('mass_scale', 1, 'load_mass', mL);
[~,~,aL] = ch4_control_affine(x0, p, unc);
[~,~,a0] = ch4_control_affine(x0, p, []);
dM = aL.M - a0.M;
dG = aL.Gv - a0.Gv;
expect_M = zeros(p.nq); expect_M(1,1) = mL; expect_M(2,2) = mL;
expect_G = zeros(p.nq,1); expect_G(2) = -mL*p.g0;
pass = report('load: M += mL on base',  norm(dM(:)-expect_M(:), inf), 1e-10, pass);
pass = report('load: G += -mL g0 on y', norm(dG-expect_G, inf),       1e-10, pass);

% Unlike a uniform scale, a load DOES change the post-impact velocity.
[xp_0] = ch4_impact(x0, p, []);
[xp_L] = ch4_impact(x0, p, unc);
d_load = norm(xp_L - xp_0, inf);
ok = d_load > 1e-6;
fprintf('  [%s] %-30s ||dq+ change|| = %.3e (must be > 0)\n', tf(ok), ...
        'load changes the impact map', d_load);
pass = pass && ok;

%% 8. a load redrawn every step: forces under the load each step carried
% ch4_run_entry analyses a randomized run step by step, so that each step's
% contact force is the one for the mass it carried. Check that against the true
% model directly: at a sample well inside step 2, the reported force must equal
% the KKT force of the robot carrying loads(2) under the torque at that sample,
% and must differ from the force with step 1's load -- which is what one model
% for the whole run would have reported.
pr = p;
pr.controller        = 'clfqp';
pr.uncertainty       = struct('mass_scale', 1, 'load_mass', 0);
pr.load_random_range = [0 70];
e = ch4_run_entry(x0, alpha, pr, struct('n_steps', 2, 'store_traj', true));
ok = e.steps_completed == 2;
if ok
    i2 = find(e.traj.t > e.step_T(1) + 0.25 * e.step_T(2), 1);
    xk = e.traj.x(:, i2);
    uk = e.traj.u(:, i2);
    err_own   = norm(e.traj.lambda(:, i2) - lam_with_load(xk, uk, pr, e.loads(2)), inf);
    gap_other = norm(lam_with_load(xk, uk, pr, e.loads(1)) - ...
                     lam_with_load(xk, uk, pr, e.loads(2)), inf);
    ok = err_own <= 1e-6 && gap_other > 1;
    fprintf(['  [%s] %-30s own load err %.1e; step 1''s load (%.1f vs %.1f kg) ' ...
             'would be off by %.0f N\n'], tf(ok), 'random load: per-step forces', ...
            err_own, e.loads(1), e.loads(2), gap_other);
else
    fprintf('  [FAIL] random-load rollout did not complete 2 steps (%s)\n', e.reason);
end
pass = pass && ok;

%% 9. physical validity: the contact force the simulator records, and its score
% ch4_step records the TRUE contact force at every solver point under the torque
% held over that period. At a point inside step 2 of a 1.5x run it must equal
% the KKT force of the TRUE model there under the held torque. Scoring: the run's
% own score must be bounded by the steps it completed, and the same run with the
% first sample of step 1 forced to lift off, or to slip, must score zero valid
% steps and name the violation.
pv = p;
pv.controller  = 'clfqp';
pv.uncertainty = struct('mass_scale', 1.5, 'load_mass', 0);
sim = ch4_simulate(x0, alpha, pv, 2);
ok = sim.n_ok == 2;
if ok
    st = sim.steps(2);
    j  = find(st.t > 0.4 * st.T, 1);
    ku = find(st.t_u < st.t(j), 1, 'last');
    [~, ~, aux] = ch4_control_affine(st.x(:, j), pv);
    err_lam = norm(st.lambda(:, j) - (aux.lam_drift + aux.lam_in * st.u(:, ku)), inf);
    Vv = ch4_validity(sim, pv);
    lift = sim; lift.steps(1).lambda(:, 1) = [0; -1];
    slip = sim; slip.steps(1).lambda(:, 1) = [0.5; 1];
    Vl = ch4_validity(lift, pv);
    Vs = ch4_validity(slip, pv);
    ok = err_lam <= 1e-8 && numel(st.lambda(1, :)) == numel(st.t) ...
         && Vv.available && Vv.valid_steps >= 0 && Vv.valid_steps <= 2 ...
         && Vl.valid_steps == 0 && strcmp(Vl.first_kind, 'lift-off') ...
         && Vs.valid_steps == 0 && strcmp(Vs.first_kind, 'slip');
    fprintf(['  [%s] %-30s recorded force err %.1e; 1.5x run valid %d/2 ' ...
             '(min Fz %.0f N, max mu %.2f); forced lift-off / slip valid %d / %d\n'], ...
            tf(ok), 'validity: recorded contact', err_lam, Vv.valid_steps, ...
            Vv.Fz_min, Vv.mu_max, Vl.valid_steps, Vs.valid_steps);
    pass = pass && ok;

    % The same run with p.stop_on_invalid: its completed steps must BE its
    % valid steps, and it must say where and how the contact failed.
    ps = pv; ps.stop_on_invalid = true;
    sims = ch4_simulate(x0, alpha, ps, 2);
    if isnan(Vv.first_step)
        ok = sims.n_ok == 2 && ~sims.failed;
    else
        ok = sims.n_ok == Vv.valid_steps && sims.failed ...
             && contains(sims.reason, 'lost contact validity') ...
             && contains(sims.reason, Vv.first_kind);
    end
    fprintf('  [%s] %-30s completed %d = valid %d | %s\n', tf(ok), ...
            'stop_on_invalid: steps = valid', sims.n_ok, Vv.valid_steps, sims.reason);
else
    fprintf('  [FAIL] validity rollout did not complete 2 steps (%s)\n', sim.reason);
end
pass = pass && ok;

% --- torque-box rule ------------------------------------------------------
% The sweeps' box is one actuator rating by default, and a parameter struct
% saved before p.box existed must resolve to the per-case rule it ran with.
pbx = ch4_params();
ok  = strcmp(ch4_box_rule(pbx), 'rating') ...
      && strcmp(ch4_box_rule(rmfield(pbx, 'box')), 'thesis') ...
      && pbx.box.rating_case4 >= pbx.box.rating;
fprintf('  [%s] %-30s default %s (%.0f / %.0f Nm), pre-p.box struct %s\n', ...
        tf(ok), 'box rule: rating by default', ch4_box_rule(pbx), ...
        pbx.box.rating, pbx.box.rating_case4, ch4_box_rule(rmfield(pbx, 'box')));
pass = pass && ok;

%% 10. the structured model: link dynamics, and the true robot built on them
% (a) ch4_link_dynamics at today's link parameters IS M.m / V.m / G.m -- an
%     independent derivation (Newton-Euler on the planar chain) against the
%     generated Lagrangian one.
L0 = ch4_link_params();
e10 = 0;
for k = 1:size(XS, 2)
    q = XS(1:7, k); dq = XS(8:14, k);
    [Ml, Vl, Gl] = ch4_link_dynamics(q, dq, L0, p.g0);
    e10 = max([e10, norm(Ml - M(q), inf), norm(Vl - V([q; dq]), inf), ...
               norm(Gl - G(q), inf)]);
end
pass = report('link dynamics == M/V/G', e10, 1e-10, pass);

% (b) a structured uncertainty at the nominal links, no actuator terms, is the
%     nominal robot
u0 = struct('mass_scale', 1, 'load_mass', 0, 'links', L0);
e10b = 0;
for k = 1:size(XS, 2)
    [f1, g1] = ch4_control_affine(XS(:, k), p, u0);
    [fn, gn] = ch3_control_affine(XS(:, k), p);
    e10b = max([e10b, norm(f1 - fn, inf), norm(g1 - gn, inf)]);
end
pass = report('structured @ nominal == Ch3', e10b, 1e-10, pass);

% (c) friction and bias act through B, like the command: with the links
%     nominal, the drift moves by exactly g * tau_x
ub = u0; ub.tau_bias = [3; -1; 2; 0.5]; ub.b_visc = [0.4; 0.2; 0.3; 0.1];
ub.tau_coulomb = [1; 0.5; 0; 2];
xk = XS(:, 3);
[f1, g1, a1] = ch4_control_affine(xk, p, ub);
[fn, gn]     = ch3_control_affine(xk, p);
tx = ch4_joint_extra(xk(8:14), ch4_structured_part(ub));
e10c = norm(f1 - (fn + gn * tx), inf) + norm(g1 - gn, inf) + norm(a1.tau_x - tx, inf);
pass = report('friction/bias enter through B', e10c, 1e-10, pass);

% (d) a draw from the prior set is reproducible, differs between seeds, and is
%     NOT a uniform scale: Delta2 is no multiple of I
Ua = ch4_uncertainty_set(3, 5); Ub = ch4_uncertainty_set(3, 5); Uc = ch4_uncertainty_set(3, 6);
same   = isequal([Ua.links], [Ub.links]) && isequal([Ua.b_visc], [Ub.b_visc]);
differ = ~isequal([Ua.links], [Uc.links]);
D10    = ch4_uncertainty(x0, alpha, p, Ua(1));
iso    = norm(D10.Delta2 - (trace(D10.Delta2) / p.ny) * eye(p.ny)) / max(D10.n2, eps);
ok     = same && differ && iso > 1e-3 && isfinite(D10.n1) && isfinite(D10.n2);
fprintf(['  [%s] %-30s reproducible %d, seeds differ %d, |Delta2| %.3f with an ' ...
         'anisotropic part %.2f of it, |Delta1| %.1f\n'], tf(ok), 'prior-set draws', ...
        same, differ, D10.n2, iso, D10.n1);
pass = pass && ok;

fprintf('--- ch4_test_model: %s ---\n\n', tf(pass));
end

% ---------------------------------------------------------------------------
function lam = lam_with_load(x, u, p, mL)
%LAM_WITH_LOAD  True stance force of the robot carrying mL, under torque u.
[~, ~, aux] = ch4_control_affine(x, p, struct('mass_scale', 1, 'load_mass', mL));
lam = aux.lam_drift + aux.lam_in * u;
end

% ---------------------------------------------------------------------------
function ok = report(name, err, tol, ok_in)
ok = ok_in && (err <= tol);
fprintf('  [%s] %-30s err = %.3e (tol %.0e)\n', tf(err <= tol), name, err, tol);
end

function s = tf(b)
if b, s = 'PASS'; else, s = 'FAIL'; end
end
