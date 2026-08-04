function ch3_test_hzd()
%CH3_TEST_HZD  Verification of the Section 6.3.4 constraint set.
%
% The zero-dynamics construction is the one piece of this chapter that cannot
% be checked by inspection: it reduces a 14-state hybrid system to a scalar
% quadrature, and every step of that reduction is a place to be quietly wrong
% while still producing plausible numbers.  So nothing here is checked against
% itself.  Each claim is checked against the FULL constrained dynamics, or
% against the collocation solution, or against an invariant that only holds if
% the reduction is right.
%
% Checks:
%   1. The closed-form point of Z really is on Z (outputs, phase, contact).
%   2. The velocity direction v_z satisfies ydot = 0 and a fixed stance foot.
%   3. The projection row w annihilates BOTH the actuators and the contact
%      wrench -- the defining property that makes the reduction possible.
%   4. The Coriolis vector is exactly quadratic in qdot, which is what lets
%      b(theta) thetadot^2 be a single coefficient times a single power.
%   5. THE REDUCTION ITSELF: thetaddot from a thetaddot + b thetadot^2 + c = 0
%      against thetaddot read off the full 14-state model under u_ff.
%   6. THE QUADRATURE: (1/2)sigma^2 + V_zero(theta) is CONSTANT along a real
%      rollout of the full robot. This is the one invariant that fails if the
%      integrating factor m is wrong, and it cannot be satisfied by accident.
%   7. delta_zero is a genuine ratio -- Delta is linear in velocity, so
%      thetadot^+ must scale exactly with thetadot^-.
%   8. Against the reference gait: the predicted fixed point matches the rate
%      the collocation actually ends at, and delta_zero^2 matches the Poincare
%      spectral radius computed by 26 independent step simulations.
%   9. Constraint plumbing: length, gating, and that each gate moves only its
%      own row.

fprintf('\n=== ch3_test_hzd ===\n');
p = ch3_params();
pass = true;

rng(7);
% A seed gait, so the suite does not depend on a solved result being present.
alpha = ch3_seed_alpha(p);

th_mid = 0.5*(p.theta_minus + p.theta_plus);
Zp = ch3_zd_point(th_mid, alpha, p);
q  = Zp.q;  v = Zp.v;

%% 1. the point is on Z
[yd_ref] = ch3_yd(alpha, Zp.s, p);
pass = report('q_z on Z: y = 0',        norm(p.H*q - yd_ref, inf),        1e-12, pass);
pass = report('q_z: c_theta q = theta', abs(p.c_theta*q - th_mid),        1e-12, pass);
pass = report('q_z: stance foot at 0',  Zp.res_foot,                      1e-12, pass);

%% 2. the velocity direction
ds_dq = p.c_theta(:).' / (p.theta_plus - p.theta_minus);
[~, dyd] = ch3_yd(alpha, Zp.s, p);
Jy = p.H - dyd * ds_dq;
pass = report('v_z: ydot = 0',          norm(Jy*v, inf),                  1e-10, pass);
pass = report('v_z: c_theta v = 1',     abs(p.c_theta*v - 1),             1e-12, pass);
pass = report('v_z: stance foot fixed', norm(J_st(q)*v, inf),             1e-10, pass);

%% 3. the projection row
% w must see neither the actuators nor the contact wrench, or the scalar
% equation it produces is not the zero dynamics but some other projection with
% u and lambda still in it.
w = Zp.w;
pass = report('w * B = 0',              norm(w*input_matrix(), inf),      1e-10, pass);
pass = report("w * J_st' = 0",          norm(w*J_st(q).', inf),           1e-10, pass);

%% 4. the Coriolis vector is quadratic in qdot
% b(theta) thetadot^2 collects w*M*v_z' and w*V. That is only a single
% coefficient times thetadot^2 if V is homogeneous of degree 2 in qdot.
V1 = V([q; v]);  V2 = V([q; 2*v]);
pass = report('V quadratic in qdot',    norm(V2 - 4*V1, inf),             1e-10, pass);

%% 5. the reduction, against the full 14-state model
% On Z under u_ff the full model must produce exactly the thetaddot the scalar
% ODE predicts. Anything wrong in w, in v_z', or in the sign convention shows
% up here and essentially nowhere else.
e_red = 0;
for thd = [0.5 1.0 2.0 3.5]
    for th = linspace(p.theta_minus, p.theta_plus, 5)
        Zk = ch3_zd_point(th, alpha, p);
        x  = [Zk.q; Zk.v * thd];
        xdot = ch3_col_dynamics(x, alpha, p);
        thdd_full = p.c_theta * xdot(p.nq+1:2*p.nq);
        thdd_red  = -(Zk.b*thd^2 + Zk.c) / Zk.a;
        e_red = max(e_red, abs(thdd_full - thdd_red) / max(1, abs(thdd_full)));
    end
end
pass = report('reduced ODE vs full model', e_red, 1e-5, pass);

%% 6. the quadrature invariant, along a real rollout
% (1/2)sigma^2 + V_zero(theta) is conserved on Z. Checked by integrating the
% ACTUAL robot under u_ff and evaluating both terms from the grid, so it tests
% m and V_zero together against a trajectory neither of them saw.
Zf = ch3_zero_dynamics(alpha, p, 641);
x0 = [Zf.q_z(:,1); Zf.v_z(:,1) * 2.0];
odeopt = odeset('RelTol', 1e-10, 'AbsTol', 1e-12);
[~, XT] = ode45(@(t,xx) ch3_col_dynamics(xx, alpha, p), linspace(0, 0.12, 25), x0, odeopt);
XT = XT.';
th_t  = p.c_theta * XT(1:p.nq, :);
thd_t = p.c_theta * XT(p.nq+1:2*p.nq, :);
inside = th_t >= Zf.theta(1) & th_t <= Zf.theta(end);
m_t = interp1(Zf.theta, Zf.m,      th_t(inside), 'spline');
V_t = interp1(Zf.theta, Zf.V_zero, th_t(inside), 'spline');
Ham = 0.5*(m_t .* thd_t(inside)).^2 + V_t;
pass = report('(1/2)sigma^2 + V_zero const', (max(Ham)-min(Ham))/max(abs(Ham)), 1e-5, pass);

%% 7. delta_zero is a ratio
% Delta solves a linear momentum balance and then permutes indices, so it is
% linear in qdot: doubling thetadot^- must exactly double thetadot^+.
x_end1 = [Zf.q_z(:,end); Zf.v_z(:,end)];
x_end2 = [Zf.q_z(:,end); Zf.v_z(:,end)*2];
xp1 = ch3_impact(x_end1, p);
xp2 = ch3_impact(x_end2, p);
d1 = p.c_theta * xp1(p.nq+1:2*p.nq);
d2 = p.c_theta * xp2(p.nq+1:2*p.nq);
pass = report('Delta linear in qdot',   abs(d2 - 2*d1)/max(1,abs(d1)),    1e-10, pass);

%% 8. against the reference gait, if it is present
ref = fullfile('Results', 'ch3_reference_gait.mat');
if exist(ref, 'file')
    S  = load(ref);
    pr = ch3_upgrade_params(S.p);
    E  = ch3_col_eval(S.z_opt, pr);
    Zr = ch3_zero_dynamics(E.alpha, pr, 321);

    thdN = pr.c_theta * E.X(pr.nq+1:2*pr.nq, E.N);
    zeta_meas = 0.5*(Zr.m(end)*thdN)^2;

    % The gait's own distance from Z bounds how well these can possibly agree:
    % zeta ~ sigma^2, so a relative error in sigma shows up doubled here.
    eta_max = 0;
    for k = 1:E.N
        [yk, ydk] = ch3_outputs(E.X(:,k), E.alpha, pr);
        eta_max = max(eta_max, norm([yk; ydk], inf));
    end
    tol_fp = max(1e-4, 20*eta_max);
    pass = report('fixed point vs collocation', ...
                  abs(zeta_meas - Zr.zeta_star)/abs(zeta_meas), tol_fp, pass);

    % delta_zero^2 is the eigenvalue of the RESTRICTED map. ch3_poincare gets
    % the spectral radius of the FULL step-to-step Jacobian by 26 simulations.
    % For a gait on Z the zero-dynamics mode is the dominant one, so these are
    % two completely independent computations of the same number.
    rho = ch3_poincare(E.X(:,1), E.alpha, pr);
    pass = report('delta^2 vs Poincare rho', abs(Zr.delta2 - rho)/max(rho,eps), 5e-3, pass);
    fprintf('       delta_zero^2 = %.5f   Poincare rho = %.5f\n', Zr.delta2, rho);
else
    fprintf('  [SKIP] reference gait absent; fixed-point cross-check not run\n');
end

%% 9. constraint plumbing
p0 = ch3_params();
z_seed = ch3_col_seed(p0);
c0 = ch3_col_constraints(z_seed, p0);
pass = report('c has 17 rows', abs(numel(c0) - 17), 0, pass);

gates = {'swing_clear', [9 10]; 'liftoff', 11; 'impact', [12 13]; ...
         'hzd', [14 15]; 'phase_mono', 16; 'decoupling', 17};
ok_gate = true;
for i = 1:size(gates,1)
    rows = gates{i,2};
    ok_gate = ok_gate && all(abs(c0(rows) + 1) < 1e-12);   % off => held at -1
    pg = p0;
    pg.limits.enable.(gates{i,1}) = true;
    cg = ch3_col_constraints(z_seed, pg);
    others = setdiff(1:17, rows);
    ok_gate = ok_gate && norm(cg(others) - c0(others), inf) < 1e-12;
end
fprintf('  [%s] %-30s each gate moves only its own rows\n', tf(ok_gate), 'gating isolation');
pass = pass && ok_gate;

fprintf('--- ch3_test_hzd: %s ---\n\n', tf(pass));
if ~pass
    error('ch3_test_hzd:failed', 'ch3_test_hzd had failures.');
end

end

% ---------------------------------------------------------------------------
function alpha = ch3_seed_alpha(p)
%CH3_SEED_ALPHA  The seed gait's alpha, without needing a solved result.
z = ch3_col_seed(p);
[~, ~, alpha] = ch3_col_unpack(z, p);
end

function ok = report(name, err, tol, ok_in)
ok = ok_in && (err <= tol);
fprintf('  [%s] %-30s err = %.3e (tol %.0e)\n', tf(err <= tol), name, err, tol);
end

function s = tf(b)
if b, s = 'PASS'; else, s = 'FAIL'; end
end
