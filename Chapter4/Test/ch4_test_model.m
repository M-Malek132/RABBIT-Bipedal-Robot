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

fprintf('--- ch4_test_model: %s ---\n\n', tf(pass));
end

% ---------------------------------------------------------------------------
function ok = report(name, err, tol, ok_in)
ok = ok_in && (err <= tol);
fprintf('  [%s] %-30s err = %.3e (tol %.0e)\n', tf(err <= tol), name, err, tol);
end

function s = tf(b)
if b, s = 'PASS'; else, s = 'FAIL'; end
end
