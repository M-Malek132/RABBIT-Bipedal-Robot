function R = ch3_report(z, p, opts)
%CH3_REPORT  Full diagnostic read-out for a collocation result.
%
%   R = ch3_report(z, p)
%   R = ch3_report(z, p, struct('stability', true, 'simulate', 3))
%
% This is the first thing to read after a solve. It answers, in order:
%   * did the transcription converge (residuals)
%   * what gait is it (speed, step length, duration, phase sweep)
%   * is it PHYSICALLY REALIZABLE (all four Table 3.1 quantities, MEASURED and
%     compared against the limits whether or not those limits were enforced)
%   * does it stay on the zero dynamics surface across the step
%   * is it actually STABLE (Poincare rho) -- optional, it costs 26 sims
%   * does it survive being simulated forward under the real controller
%
% MEASURE BEFORE YOU CONSTRAIN.  Every Table 3.1 quantity is reported with its
% measured value AND the limit, with a marker showing whether the limit is
% currently enforced. The thesis limits are ATRIAS numbers (63 kg, 50:1 gears,
% so its 5 Nm is 250 Nm at the joint); RABBIT is ~30 kg direct drive. Copying
% them across without measuring first produces an infeasible problem and a
% solver that fails for reasons that look like bugs.
%
% Inputs
%   z    : collocation decision vector
%   p    : parameter struct
%   opts : optional struct
%            .stability  (false) compute the Poincare spectral radius
%            .simulate   (0)     simulate this many steps under p.controller
%            .verbose    (true)  print
%
% Output
%   R : struct of everything printed, for programmatic use
%
% See also CH3_COL_EVAL, CH3_POINCARE, CH3_FORCES.

p = ch3_upgrade_params(p);   % saved results may predate fields read below

if nargin < 3, opts = struct(); end
if ~isfield(opts,'stability'), opts.stability = false; end
if ~isfield(opts,'simulate'),  opts.simulate  = 0;     end
if ~isfield(opts,'verbose'),   opts.verbose   = true;  end

E = ch3_col_eval(z, p);
[c, ceq] = ch3_col_constraints(z, p);
[X, T, alpha] = ch3_col_unpack(z, p);
N = E.N;

R = struct();
R.alpha = alpha; R.T = T; R.X = X; R.N = N;
R.max_ceq = max(abs(ceq));
R.max_c   = max(c);

% --- residuals by block ---------------------------------------------------
i = 0;
blk = {'node-1 (on Z, phase, contact)', 13; ...
       'Hermite-Simpson defects',       p.nx*(N-1); ...
       'node-N (on guard S)',           2; ...
       'periodicity through Delta',     13; ...
       'NEC1 walking rate',             1};
R.blocks = struct();
for b = 1:size(blk,1)
    n = blk{b,2};
    R.blocks.(matlab.lang.makeValidName(blk{b,1})) = max(abs(ceq(i+1:i+n)));
    i = i + n;
end

% --- gait -----------------------------------------------------------------
R.L_step = E.L_step;
R.speed  = E.L_step / E.T;
R.theta_sweep = [p.c_theta*X(1:p.nq,1), p.c_theta*X(1:p.nq,N)];
R.cost = ch3_col_cost(z, p);
R.int_u2 = E.int_u2;

% --- Table 3.1 quantities, measured --------------------------------------
lam_all = [E.lam, E.lamm];
R.torque_max  = max([abs(E.u(:)); abs(E.um(:))]);
R.Fz_min      = min(lam_all(2,:));
R.Fz_max      = max(lam_all(2,:));
pos = lam_all(2,:) > 1e-6;
if any(pos)
    R.mu_max = max(abs(lam_all(1,pos)) ./ lam_all(2,pos));
else
    R.mu_max = Inf;
end
R.impulse     = norm(E.impulse);
R.clearance   = E.sw_h(max(2, min(N-1, round((N+1)/2))));
R.sw_h_min    = min(E.sw_h(2:N-1));

% --- Section 6.3.4 quantities, measured ----------------------------------
% Same measure-before-you-constrain rule as Table 3.1: every one of these is
% computed whether or not its gate is on, because the only way to pick a
% sensible bound is to see what the gait already does.
h_int = E.sw_h(2:N-1);
if N > 3, h_int = [h_int, E.sw_hm(2:N-2)]; end
R.sw_clear_int = min(h_int);              % NIC3(a)
R.sw_strike    = E.sw_hd(N);              % NIC3(b), negative = descending
R.liftoff      = E.sw_hd_post;            % NEC2
R.impulse_z    = E.impulse(2);            % NEC3
R.impulse_x    = E.impulse(1);
if R.impulse_z > 1e-9
    R.impulse_mu = abs(R.impulse_x) / R.impulse_z;
else
    R.impulse_mu = Inf;
end
R.thetadot_min = min([E.thd, E.thdm]);    % HH6
R.dec_min      = E.dec_min;               % HH2
R.eta_post     = norm(E.eta_post, inf);   % HH4/HH5, implied not imposed

% --- hybrid zero dynamics: NEC4, NEC5, i.e. (5.79)/(5.80) ----------------
% Computed from alpha alone -- no trajectory, no simulation -- so it is an
% INDEPENDENT read on the same gait, and the cross-check printed below (its
% predicted fixed point against the rate the collocation actually ends at) is
% the cheapest evidence that both are right.
R.hzd = [];
if ~isfield(opts,'hzd') || opts.hzd
    try
        R.hzd = ch3_zero_dynamics(alpha, p, p.hzd_grid);
        R.hzd_thetadot_end = p.c_theta * X(p.nq+1:2*p.nq, N);
        R.hzd_zeta_meas    = 0.5 * (R.hzd.m(end) * R.hzd_thetadot_end)^2;
    catch err
        R.hzd_error = err.message;
    end
end

% --- posture: torso pitch and hip height ---------------------------------
% Over nodes AND midpoints, matching what ch3_col_constraints actually
% bounds -- reporting nodes only would under-report both ranges.
qt_all  =  [X(3,:), E.xm(3,:)];
hip_all = -[X(2,:), E.xm(2,:)];     % pz is DOWN-positive; height = -pz
R.qt_lo   = min(qt_all);
R.qt_hi   = max(qt_all);
R.hip_lo  = min(hip_all);
R.hip_hi  = max(hip_all);
R.hip_bob = R.hip_hi - R.hip_lo;
R.hip_dev = max(abs(hip_all - p.limits.hip_h));

% --- zero dynamics drift across the step ---------------------------------
eta_max = 0;
for k = 1:N
    [yk, ydk] = ch3_outputs(X(:,k), alpha, p);
    eta_max = max(eta_max, norm([yk; ydk], inf));
end
R.eta_max = eta_max;

% --- is this an actual trajectory, or a spurious discrete solution? -------
% Checked BEFORE anything else is believed: if the mesh is too coarse, every
% gait number above is fiction.
R.verify = ch3_col_verify(z, p, false);

% --- stance foot drift ----------------------------------------------------
foot0 = P_st(X(1:p.nq,1));
drift = 0;
for k = 1:N
    fk = P_st(X(1:p.nq,k));
    drift = max(drift, norm(fk - foot0, inf));
end
R.stance_drift = drift;

% --- stability ------------------------------------------------------------
R.rho = NaN;
if opts.stability
    [R.rho, ~, pinfo] = ch3_poincare(X(:,1), alpha, p);
    R.poincare = pinfo;
end

% --- forward simulation under the REAL controller -------------------------
R.sim = [];
if opts.simulate > 0
    R.sim = ch3_simulate(X(:,1), alpha, p, opts.simulate);
end

if ~opts.verbose, return; end

%% ------------------------------------------------------------- printing
line = repmat('-', 1, 74);
fprintf('\n%s\n CHAPTER 3 GAIT REPORT   (controller: %s, basis: %s, N = %d)\n%s\n', ...
        line, p.controller, p.basis, N, line);

fprintf(' CONVERGENCE\n');
fprintf('   max |ceq|                    %.3e\n', R.max_ceq);
fprintf('   max c (inequality)           %.3e\n', R.max_c);
fn = fieldnames(R.blocks);
for k = 1:numel(fn)
    fprintf('     %-30s %.3e\n', blk{k,1}, R.blocks.(fn{k}));
end

fprintf('\n GAIT\n');
fprintf('   step length                  %.4f m\n', R.L_step);
fprintf('   step duration                %.4f s\n', R.T);
fprintf('   walking speed                %.4f m/s   (v_des %.4f)\n', R.speed, p.v_des);
fprintf('   theta sweep                  %.4f -> %.4f   (target %.4f -> %.4f)\n', ...
        R.theta_sweep(1), R.theta_sweep(2), p.theta_minus, p.theta_plus);
fprintf('   cost  int||u||^2 / L_step    %.4f      (int||u||^2 = %.4f)\n', R.cost, R.int_u2);

fprintf('\n POSTURE   (design requirements, not Table 3.1)\n');
fprintf('   torso pitch qt               %+.4f .. %+.4f rad  (%+.1f .. %+.1f deg)\n', ...
        R.qt_lo, R.qt_hi, rad2deg(R.qt_lo), rad2deg(R.qt_hi));
fprintf('     box                        [%+.3f %+.3f]  -> leaning %s\n', ...
        p.qt_range(1), p.qt_range(2), ...
        ternary(R.qt_hi < 0, 'BACKWARD', ternary(R.qt_lo > 0, 'FORWARD', 'through vertical')));
fprintf('   hip height                   %.4f .. %.4f m   (bob %.4f m)\n', ...
        R.hip_lo, R.hip_hi, R.hip_bob);
fprintf('     leg extension              %.1f%% .. %.1f%% of the 1.0 m leg\n', ...
        100*R.hip_lo, 100*R.hip_hi);

fprintf('\n PHYSICAL REALIZABILITY  (Table 3.1)   [E] = enforced this solve\n');
prow('peak |torque|',      R.torque_max, p.limits.u_max,       '<=', 'Nm', p.limits.enable.torque);
prow('impact impulse',     R.impulse,    p.limits.impulse_max, '<=', 'Ns', p.limits.enable.impulse);
prow('friction demanded',  R.mu_max,     p.limits.mu_s,        '<=', '-',  p.limits.enable.friction);
prow('min vertical GRF',   R.Fz_min,     p.limits.Fz_min,      '>=', 'N',  p.limits.enable.grf);
prow('mid-step clearance', R.clearance,  p.limits.clearance,   '>=', 'm',  p.limits.enable.clearance);
prow('hip-height dev',     R.hip_dev,    p.limits.hip_h_tol,   '<=', 'm',  p.limits.enable.height);
fprintf('   (max Fz %.1f N, min swing-foot height over the step %.4f m)\n', R.Fz_max, R.sw_h_min);

fprintf('\n SECTION 6.3.4 CONSTRAINT SET   [E] = enforced this solve\n');
prow('NIC1 min normal GRF',   R.Fz_min,       p.limits.Fz_min,        '>=', 'N',   p.limits.enable.grf);
prow('NIC2 friction demanded',R.mu_max,       p.limits.mu_s,          '<=', '-',   p.limits.enable.friction);
prow('NIC3 interior clearance',R.sw_clear_int,p.limits.sw_clear_min,  '>=', 'm',   p.limits.enable.swing_clear);
prow('NIC3 strike rate',      R.sw_strike,   -p.limits.sw_strike_rate,'<=', 'm/s', p.limits.enable.swing_clear);
fprintf('   %s %-22s %10.4f =  %-9.4f %-3s  (residual %.1e)\n', ...
        ternary(~isfield(p,'enforce_nec1') || p.enforce_nec1, '[E]', '[ ]'), ...
        'NEC1 walking rate', R.speed, p.v_des, 'm/s', abs(R.speed - p.v_des));
prow('NEC2 lift-off rate',    R.liftoff,      p.limits.liftoff_rate,  '>=', 'm/s', p.limits.enable.liftoff);
mu_imp = p.limits.mu_s_impact;
if isempty(mu_imp), mu_imp = p.limits.mu_s; end
R.mu_s_impact = mu_imp;
prow('NEC3 impulse Iz',       R.impulse_z,    p.limits.impulse_z_min, '>=', 'Ns',  p.limits.enable.impact);
prow('NEC3 impulse |Ix|/Iz',  R.impulse_mu,   mu_imp,                 '<=', '-',   p.limits.enable.impact);
if ~isempty(R.hzd)
    prow('NEC4 zeta*_2',      R.hzd.zeta_star, R.hzd.V_max/R.hzd.delta2, '>=', '-', p.limits.enable.hzd);
    prow('NEC5 delta_zero^2', R.hzd.delta2,    1,                        '<=', '-', p.limits.enable.hzd);
else
    fprintf('   NEC4/NEC5 not evaluated%s\n', ...
            ternary(isfield(R,'hzd_error'), [': ' R.hzd_error], ''));
end
prow('HH2  sigma_min(LgLfy)', R.dec_min,      p.limits.dec_min,       '>=', '-',     p.limits.enable.decoupling);
prow('HH6  min thetadot',     R.thetadot_min, p.limits.thetadot_min,  '>=', 'rad/s', p.limits.enable.phase_mono);
fprintf('       %-22s %10.3e         %-3s  %s\n', 'HH4/HH5 |eta(Delta x_N)|', R.eta_post, '-', ...
        ternary(R.eta_post < 1e-4, 'ok (implied by periodicity, not imposed)', ...
                                   'CHECK: hybrid invariance is NOT holding'));

if ~isempty(R.hzd)
    Zd = R.hzd;
    fprintf('\n HYBRID ZERO DYNAMICS   (restricted Poincare map, from alpha alone)\n');
    fprintf('   delta_zero                   %+.6f   (sigma^+ / sigma^-)\n', Zd.delta_zero);
    fprintf('   delta_zero^2                 %.6f    -> %s\n', Zd.delta2, ...
            ternary(Zd.stable, 'NEC5 ok, fixed point ATTRACTS', 'NEC5 VIOLATED, fixed point repels'));
    fprintf('   V_zero(theta_f) / V_zero^MAX %+.6f / %+.6f\n', Zd.V_end, Zd.V_max);
    fprintf('   zeta*_2 (fixed point)        %.6f  vs  V^MAX/delta^2 = %.6f  -> %s\n', ...
            Zd.zeta_star, Zd.V_max/Zd.delta2, ternary(Zd.exists, 'NEC4 ok', 'NEC4 VIOLATED, gait stalls mid-step'));
    if isfield(R, 'hzd_thetadot_end')
        rel = abs(R.hzd_zeta_meas - Zd.zeta_star) / max(abs(R.hzd_zeta_meas), eps);
        fprintf('   CROSS-CHECK  thetadot^-       %.6f predicted   %.6f from the collocation\n', ...
                Zd.thetadot_star, R.hzd_thetadot_end);
        fprintf('                zeta^-           %.6f predicted   %.6f from the collocation  (rel %.1e)\n', ...
                Zd.zeta_star, R.hzd_zeta_meas, rel);
        fprintf('     Two independent routes to the same fixed point: a quadrature over\n');
        fprintf('     alpha, and a converged periodicity constraint. Disagreement above\n');
        fprintf('     ~max|eta| means one of them is wrong.\n');
    end
    if opts.stability && isfinite(R.rho)
        fprintf('   vs Poincare rho              %.4f  (full map, %d sims) -- delta^2 is the\n', R.rho, 26);
        fprintf('     ZERO-DYNAMICS eigenvalue, so rho should not fall below it.\n');
    end
end

fprintf('\n TRAJECTORY VALIDITY  (nodes vs a tight rollout -- read this FIRST)\n');
fprintf('   max |X_node - X_true|        %.3e   (tol %.1e)\n', R.verify.max_dev, p.verify_tol);
if R.verify.ok
    fprintf('   -> the nodes lie on a genuine trajectory; the numbers above are real.\n');
else
    fprintf('   -> *** REJECT: SPURIOUS DISCRETE SOLUTION ***\n');
    if ~R.verify.integrated
        fprintf('      The reference rollout could not even complete the step.\n');
    end
    fprintf('      The defects are satisfied but the nodes do not lie on any real\n');
    fprintf('      trajectory, so EVERY gait and Table 3.1 number above is fiction.\n');
    fprintf('      Cause is an under-resolved mesh. Re-solve with larger p.N_nodes\n');
    fprintf('      (currently %d).\n', N);
end

fprintf('\n INVARIANTS\n');
fprintf('   max |eta| over nodes         %.3e   (gait on the zero dynamics surface)\n', R.eta_max);
if R.eta_max > 1e-3
    fprintf('     ^ should be ~0. Large |eta| with a SMALL node-1 residual means the\n');
    fprintf('       mesh is too coarse, not that the invariance is wrong: u_ff makes\n');
    fprintf('       ydd = 0 exactly, so eta can only drift through discretization.\n');
end
fprintf('   stance-foot drift            %.3e m (pinned at every node)\n', R.stance_drift);

if opts.stability
    fprintf('\n ORBITAL STABILITY\n');
    fprintf('   Poincare rho                 %.4f   -> %s\n', R.rho, ...
            ternary(R.rho < 1, 'STABLE (attracting orbit)', 'UNSTABLE (periodic but falls)'));
    fprintf('   fixed-point residual         %.3e\n', R.poincare.fixed_point_residual);
    if R.rho >= 1
        fprintf('   NOTE: periodicity is satisfied but the orbit REPELS. Periodic is not\n');
        fprintf('         stable -- simulate it and the velocities grow step over step.\n');
    end
end

if opts.simulate > 0
    fprintf('\n FORWARD SIMULATION under "%s"  (%d steps requested)\n', p.controller, opts.simulate);
    fprintf('   steps completed              %d\n', R.sim.n_ok);
    fprintf('   outcome                      %s\n', R.sim.reason);
    if R.sim.n_ok > 0
        Ls = arrayfun(@(s) s.L_step, R.sim.steps);
        Ts = arrayfun(@(s) s.T,      R.sim.steps);
        fprintf('   step lengths                 %s\n', num2str(Ls, '%.3f '));
        fprintf('   step durations               %s\n', num2str(Ts, '%.3f '));
    end
end

fprintf('%s\n\n', line);

end

% ---------------------------------------------------------------------------
function prow(name, val, lim, sense, unit, enabled)
if strcmp(sense, '<=')
    ok = val <= lim;
else
    ok = val >= lim;
end
mark = ternary(ok, ' ok ', 'OVER');
if strcmp(sense,'>='), mark = ternary(ok, ' ok ', 'UNDER'); end
tag = ternary(enabled, '[E]', '[ ]');
fprintf('   %s %-22s %10.4f %s %-9.4f %-3s  %s\n', tag, name, val, sense, lim, unit, mark);
end

function o = ternary(cond, a, b)
if cond, o = a; else, o = b; end
end
