function R = ch3_validate_gait(gait_file, opts)
%CH3_VALIDATE_GAIT  Audit stored gaits: their own limits, the declared spec, the true flow.
%
%   R = ch3_validate_gait()               the four reference candidates, compared
%   R = ch3_validate_gait(gait_file)      one gait file (z, z_try or z_opt, and p)
%   R = ch3_validate_gait({f1, f2, ...})  several, with compare.log at the end
%   R = ch3_validate_gait(gait_file, opts)
%
%   opts.n_steps      closed-loop steps per controller (default 10; 0 skips them)
%   opts.controllers  default {'iolin_pd', 'clfqp_con'}
%   opts.control_dt   ZOH period of those runs (default 1e-3, as Ch3 and Ch4)
%   opts.stability    Poincare rho under iolin_pd, 27 steps (default true)
%   opts.n_dense      samples of the true rollout (default 4001)
%   opts.tol          pass tolerance in each quantity's own units (default 1e-6,
%                     fmincon's ConstraintTolerance, as ch3_col_check_limits)
%
% The four candidates, in lineage order: posture_195 (the Ch4/Ch6 reference,
% 195 Nm), speed_ladder/gait_v1200 (NEC1 at the 1.2 m/s design speed),
% torque_march/gait_u162 (162.5 Nm, this family's floor) and
% impulse_march/u162_I1500 (162.5 Nm and 15 Ns, the best so far). The last
% three are git-ignored and live on the machine that solved them; a file this
% machine lacks is skipped with a line in the log.
%
% WHY ANOTHER CHECK. ch3_col_check_limits, and the marches' landing test, ask
% whether a gait meets the limits ITS OWN saved p enables. That is
% self-consistency, and a gait passes it however far its p was relaxed on the
% way. The best gait so far is the case in point: the impulse march reports it
% as meeting every limit, and it does -- the limits of a p with a 162.5 Nm box
% (the spec is 120), NEC1 off (the spec asks 1.2 m/s), a 0.925 +- 0.045 m hip
% band (the spec is 0.90 +- 0.02), and the swing-foot ceiling never on. The
% lineage starts at posture_195, saved 2026-09-10; row 18 arrived 2026-09-12,
% and ch3_upgrade_params fills a gate a saved p never mentions as FALSE, so it
% stayed off for every rung descended from it. posture_195 itself peaks at
% 0.212 m against 0.15.
%
% The rows also see only nodes and Hermite-Simpson midpoints, and midpoints 1
% and N-1 get no lower bound at all: row 9 skips them and says row 2 covers
% them, but row 2 reads the nodes only.
%
% SIX QUESTIONS PER GAIT
%   1. A real trajectory?   ch3_col_verify, then one step of the TRUE flow
%      from node 1 to the guard under u_ff: T, L, impulse and periodicity off
%      the rollout, not off the nodes.
%   2. Its own limits?      Every row of ch3_col_constraints, each measured with
%      its gate on and marked [E] where the solve enforced it.
%   3. The declared spec?   The same rows under ch3_params() -- limit values,
%      gates, qt_range, step_len_min, the box bounds and NEC1 at v_des -- with
%      the gait's own parametrization kept. Physical rows (Table 3.1, NIC, NEC,
%      HH) are summarized apart from design rows (mid-step clearance, ceiling,
%      hip band, torso box, joint-rate rail, speed).
%   4. Between the samples? The same quantities over a dense grid of the true
%      rollout, and at its own guard. '*' marks a limit the collocation meets
%      and the true flow crosses -- between samples, or through the ~1e-4
%      drift ch3_col_verify allows.
%   5. HZD at report grid?  NEC4/NEC5 at p.hzd_grid (161) against the solve's
%      p.hzd_grid_solve (41).
%   6. Closed loop?         n_steps under each controller at a 1 kHz ZOH, scored
%      by ch3_validity (a sample with Fz <= 0 or |Fx|/Fz > mu_s ends the valid
%      run), with the peak torque actually HELD; and the Poincare rho, under
%      iolin_pd because a QP inside ode45 stalls it (ch3_compare_controllers).
%
% Nothing is re-solved and no gait file is written. Output, per gait:
% Results/reruns/validate/<gait>.log and <gait>.mat (struct R); for a list,
% compare.log as well. Every log ends in VALIDATE_GAIT_DONE.
%
% See also CH3_COL_CHECK_LIMITS, CH3_COL_VERIFY, CH3_REPORT, CH3_VALIDITY.

ROOT = fileparts(fileparts(fileparts(fileparts(mfilename('fullpath')))));
OUT  = fullfile(ROOT, 'Results', 'reruns', 'validate');
if ~exist(OUT, 'dir'), mkdir(OUT); end

if nargin < 1 || isempty(gait_file)
    gait_file = {fullfile('Results', 'ch3_gait_posture_195.mat'), ...
                 fullfile('Results', 'reruns', 'speed_ladder', 'gait_v1200.mat'), ...
                 fullfile('Results', 'reruns', 'torque_march', 'gait_u162.mat'), ...
                 fullfile('Results', 'reruns', 'impulse_march', 'u162_I1500.mat')};
end
if nargin < 2, opts = struct(); end
opts = with_defaults(opts);

if ~iscell(gait_file)
    R = validate_one(ROOT, OUT, gait_file, opts);
    return;
end

logf = fullfile(OUT, 'compare.log');
if exist(logf, 'file'), delete(logf); end
R = {};
for i = 1:numel(gait_file)
    if isempty(find_gait(ROOT, gait_file{i}))
        logln(logf, '%s: not on this machine, skipped', gait_file{i});
        continue;
    end
    logln(logf, '%s: validating (%s)', gait_file{i}, datestr(now));
    try
        R{end+1} = validate_one(ROOT, OUT, gait_file{i}, opts); %#ok<AGROW>
    catch err
        % One gait that errors must not cost the others their run.
        logln(logf, '%s: FAILED -- %s', gait_file{i}, err.message);
    end
end
compare(logf, R);
logln(logf, '=== VALIDATE_GAIT_DONE %s', datestr(now));
end

% ---------------------------------------------------------------------------
function R = validate_one(ROOT, OUT, gait_file, opts)
src = find_gait(ROOT, gait_file);
if isempty(src)
    error('ch3_validate_gait:noGait', ...
          ['Gait file not found: %s. The Chapter 3 campaign gaits (Results/reruns/' ...
           'speed_ladder/, torque_march/, impulse_march/) are git-ignored and live ' ...
           'on the machine that solved them; run this there, or copy the file here.'], ...
          gait_file);
end
[~, name] = fileparts(src);
logf = fullfile(OUT, [name '.log']);
if exist(logf, 'file'), delete(logf); end     % appended to below: a rerun must not double it
t_all = tic;
tol   = opts.tol;

S = load(src);
if isfield(S, 'z_try')                        % a march rung
    z = S.z_try;
elseif isfield(S, 'z')
    z = S.z;
elseif isfield(S, 'z_opt')
    z = S.z_opt;
else
    error('ch3_validate_gait:contents', '%s holds no z, z_try or z_opt.', src);
end
if ~isfield(S, 'p'), error('ch3_validate_gait:contents', '%s holds no p.', src); end
p = ch3_upgrade_params(S.p);
d = ch3_params();
logln(logf, '=== validate %s | %s', src, datestr(now));

%% 1. a real trajectory? ----------------------------------------------------
E = ch3_col_eval(z, p);
[X, T, alpha] = ch3_col_unpack(z, p);
N  = E.N;  nq = p.nq;  h = E.h;
c_row = p.c_theta(:).';
V = ch3_col_verify(z, p, false);
[~, ceq_own] = ch3_col_constraints(z, p);

% One step of the true flow under u_ff, to the GUARD rather than to the
% collocation's T, at ch3_col_verify's tolerances. Everything in section 4
% labelled "rollout" is read off this solution.
pf = p;  pf.controller = 'ff';  pf.control_dt = 0;
pf.ode_reltol = 1e-11;  pf.ode_abstol = 1e-13;
st = ch3_step(X(:,1), alpha, pf);
per_true = norm(st.x_next(2:end) - X(2:end,1), inf);
per_node = norm(E.x_next(2:end) - X(2:end,1), inf);
logln(logf, ['1  TRAJECTORY  verify %d (max dev %.2e, tol %.0e) | max|ceq| %.2e | ' ...
             'true flow: guard %s, T %.4f (nodes %.4f), L %.4f (%.4f), ' ...
             'periodicity %.2e (nodes %.2e)'], V.ok, V.max_dev, p.verify_tol, ...
      max(abs(ceq_own)), ternary(st.ok, 'reached', 'NOT REACHED'), st.T, T, ...
      st.L_step, E.L_step, per_true, per_node);

%% 2-3. the constraint rows, own p and spec ------------------------------
[names, gate, phys_row] = row_table();
p_spec = spec_params(p, d);
c_own  = ch3_col_constraints(z, gates_on(p, gate));        % every row measured
c_spec = ch3_col_constraints(z, gates_on(p_spec, gate));
on_own  = cellfun(@(g) gate_on(p, g), gate);
on_spec = cellfun(@(g) gate_on(p_spec, g), gate);

logln(logf, '');
logln(logf, ['2-3  ROWS of ch3_col_constraints, each measured with its gate on ' ...
             '(c <= 0 is met; [E] = enforced; "over" = exceeded with the gate off)']);
logln(logf, '     %-50s %-22s %s', 'row', 'own p', 'spec = ch3_params()');
for r = 1:19
    logln(logf, '  %2d %-50s %s %+11.3e %-4s  %s %+11.3e %-4s', r, names{r}, ...
          etag(on_own(r)), c_own(r), rstr(c_own(r), on_own(r), tol), ...
          etag(on_spec(r)), c_spec(r), rstr(c_spec(r), on_spec(r), tol));
end
bo = bound_violation(z, E, p, N);
bs = bound_violation(z, E, p_spec, N);
logln(logf, ['     box bounds (ch3_col_bounds, <= 0 is inside): own nodes %+.2e, ' ...
             'midpoints %+.2e | spec nodes %+.2e, midpoints %+.2e'], bo(1), bo(2), bs(1), bs(2));

%% 4. the same quantities, at the samples and between them ---------------
t_d = linspace(0, st.T, opts.n_dense);
Xd  = deval(st.sol, t_d);
nd  = numel(t_d);
Ud = zeros(p.nu, nd);  Ld = zeros(2, nd);
Hd = zeros(1, nd);  THd = zeros(1, nd);  DECd = zeros(1, nd);  ETAd = zeros(1, nd);
for k = 1:nd
    [~, LgLfy, u, info] = ch3_io_lin(Xd(:,k), alpha, p);
    Ud(:,k)  = u;
    Ld(:,k)  = info.aux.lam_drift + info.aux.lam_in * u;
    Pk       = P_sw(Xd(1:nq,k));
    Hd(k)    = Pk(2);
    THd(k)   = c_row * Xd(nq+1:end,k);
    DECd(k)  = min(svd(LgLfy));
    ETAd(k)  = norm(info.eta, inf);
end

lam_c = [E.lam, E.lamm];
hip_c = -[X(2,:), E.xm(2,:)];     hip_d = -Xd(2,:);      % pz is down-positive
qt_c  =  [X(3,:), E.xm(3,:)];     qt_d  =  Xd(3,:);
inner = t_d >= h & t_d <= T - h;                          % the span row 9 covers
k_mid = max(2, min(N-1, round((N+1)/2)));
x_mid = deval(st.sol, min((k_mid-1)*h, st.T));
P_mid = P_sw(x_mid(1:nq));
J_end = J_sw(st.x_end(1:nq));    strike_d  = J_end(2,:) * st.x_end(nq+1:end);
J_nxt = J_sw(st.x_next(1:nq));   liftoff_d = J_nxt(2,:) * st.x_next(nq+1:end);
[y_p, yd_p] = ch3_outputs(st.x_next, alpha, p);
mu_imp_o = p.limits.mu_s_impact;  if isempty(mu_imp_o), mu_imp_o = p.limits.mu_s; end
mu_imp_s = d.limits.mu_s_impact;  if isempty(mu_imp_s), mu_imp_s = d.limits.mu_s; end
Zs = ch3_zero_dynamics(alpha, p, p.hzd_grid_solve);
Zr = ch3_zero_dynamics(alpha, p, p.hzd_grid);
eo = p.limits.enable;  es = p_spec.limits.enable;

u_pk = max(abs([E.u(:); E.um(:)]));
mu_c = fric_ratio(lam_c);
Fz_c = min(lam_c(2,:));
h_pk = max([E.sw_h, E.sw_hm]);
v    = E.L_step / E.T;

Q = [];
Q = addq(Q, 'peak |u|  (Table 3.1)', 'Nm', '<=', u_pk, max(abs(Ud(:))), ...
         p.limits.u_max, eo.torque, d.limits.u_max, es.torque, true);
Q = addq(Q, 'impact impulse ||I||  (Table 3.1)', 'Ns', '<=', norm(E.impulse), norm(st.impulse), ...
         p.limits.impulse_max, eo.impulse, d.limits.impulse_max, es.impulse, true);
Q = addq(Q, 'stance friction |Fx|/Fz  (NIC2)', '-', '<=', mu_c, fric_ratio(Ld), ...
         p.limits.mu_s, eo.friction, d.limits.mu_s, es.friction, true);
Q = addq(Q, 'min normal force Fz  (NIC1)', 'N', '>=', Fz_c, min(Ld(2,:)), ...
         p.limits.Fz_min, eo.grf, d.limits.Fz_min, es.grf, true);
Q = addq(Q, 'swing foot height, interior  (NIC3)', 'm', '>=', ...
         min([E.sw_h(2:N-1), E.sw_hm(2:N-2)]), min([Hd(inner), NaN]), ...
         p.limits.sw_clear_min, eo.swing_clear, d.limits.sw_clear_min, es.swing_clear, true);
Q = addq(Q, 'swing foot height, all of (0,T)  (row 2)', 'm', '>=', ...
         min([E.sw_h(2:N-1), E.sw_hm]), min(Hd(2:end-1)), 0, true, 0, true, true);
Q = addq(Q, 'strike velocity  (NIC3)', 'm/s', '<=', E.sw_hd(N), strike_d, ...
         -p.limits.sw_strike_rate, eo.swing_clear, -d.limits.sw_strike_rate, es.swing_clear, true);
Q = addq(Q, 'lift-off velocity  (NEC2)', 'm/s', '>=', E.sw_hd_post, liftoff_d, ...
         p.limits.liftoff_rate, eo.liftoff, d.limits.liftoff_rate, es.liftoff, true);
Q = addq(Q, 'impulse Iz  (NEC3)', 'Ns', '>=', E.impulse(2), st.impulse(2), ...
         p.limits.impulse_z_min, eo.impact, d.limits.impulse_z_min, es.impact, true);
Q = addq(Q, 'impulse |Ix|/Iz  (NEC3)', '-', '<=', abs(E.impulse(1))/E.impulse(2), ...
         abs(st.impulse(1))/st.impulse(2), mu_imp_o, eo.impact, mu_imp_s, es.impact, true);
Q = addq(Q, 'NEC4 residual, grid 41 | 161', '-', '<=', Zs.nec4, Zr.nec4, ...
         0, eo.hzd, 0, es.hzd, true);
Q = addq(Q, 'NEC5 residual, grid 41 | 161', '-', '<=', Zs.nec5, Zr.nec5, ...
         0, eo.hzd, 0, es.hzd, true);
Q = addq(Q, 'min thetadot  (HH6)', 'rad/s', '>=', min([E.thd, E.thdm]), min(THd), ...
         p.limits.thetadot_min, eo.phase_mono, d.limits.thetadot_min, es.phase_mono, true);
Q = addq(Q, 'min sigma(LgLf y)  (HH2)', '-', '>=', E.dec_min, min(DECd), ...
         p.limits.dec_min, eo.decoupling, d.limits.dec_min, es.decoupling, true);
Q = addq(Q, '|eta| after impact  (HH4/HH5, implied)', '-', '<=', norm(E.eta_post, inf), ...
         norm([y_p; yd_p], inf), 1e-4, true, 1e-4, true, true);   % ch3_report's threshold
Q = addq(Q, 'step length  (row 3)', 'm', '>=', E.L_step, st.L_step, ...
         p.step_len_min, true, d.step_len_min, true, true);
Q = addq(Q, 'mid-step clearance  (row 1)', 'm', '>=', E.sw_h(k_mid), P_mid(2), ...
         p.limits.clearance, eo.clearance, d.limits.clearance, es.clearance, false);
Q = addq(Q, 'peak swing-foot height  (row 18)', 'm', '<=', h_pk, max(Hd), ...
         p.limits.clearance_max, eo.clearance_max, d.limits.clearance_max, es.clearance_max, false);
Q = addq(Q, 'hip height, lowest  (row 8)', 'm', '>=', min(hip_c), min(hip_d), ...
         p.limits.hip_h - p.limits.hip_h_tol, eo.height, ...
         d.limits.hip_h - d.limits.hip_h_tol, es.height, false);
Q = addq(Q, 'hip height, highest  (row 8)', 'm', '<=', max(hip_c), max(hip_d), ...
         p.limits.hip_h + p.limits.hip_h_tol, eo.height, ...
         d.limits.hip_h + d.limits.hip_h_tol, es.height, false);
Q = addq(Q, 'torso pitch, lowest  (bound)', 'rad', '>=', min(qt_c), min(qt_d), ...
         p.qt_range(1), true, d.qt_range(1), true, false);
Q = addq(Q, 'torso pitch, highest  (bound)', 'rad', '<=', max(qt_c), max(qt_d), ...
         p.qt_range(2), true, d.qt_range(2), true, false);
Q = addq(Q, 'max |dq|  (bound)', 'rad/s', '<=', ...
         max(max(abs([X(nq+1:end,:), E.xm(nq+1:end,:)]))), max(max(abs(Xd(nq+1:end,:)))), ...
         p.dq_max, true, d.dq_max, true, false);
Q = addq(Q, 'walking speed L/T  (NEC1)', 'm/s', '=', v, st.L_step/st.T, ...
         p.v_des, nec1_on(p), d.v_des, nec1_on(p_spec), false);

logln(logf, '');
logln(logf, ['4  QUANTITIES  colloc = nodes + midpoints; rollout = %d samples of the true ' ...
             'flow and its own guard (for NEC4/NEC5: the 161-point grid); * = met by the ' ...
             'collocation, crossed by the true flow'], nd);
logln(logf, '   %-42s %-5s %11s %11s | %-23s | %s', 'quantity', 'unit', 'colloc', ...
      'rollout', 'own p', 'spec');
for i = 1:numel(Q)
    q = Q(i);
    [ok_o, x_o] = judge(q.col, q.roll, q.lim_o, q.sense, tol);
    [ok_s, x_s] = judge(q.col, q.roll, q.lim_s, q.sense, tol);
    Q(i).ok_o = ok_o;  Q(i).x_o = x_o;  Q(i).ok_s = ok_s;  Q(i).x_s = x_s;
    logln(logf, '   %-42s %-5s %11.5g %11.5g | %s %2s %9.5g %-5s | %s %2s %9.5g %-5s', ...
          q.name, q.unit, q.col, q.roll, ...
          etag(q.on_o), q.sense, q.lim_o, vstr(ok_o, x_o, q.on_o), ...
          etag(q.on_s), q.sense, q.lim_s, vstr(ok_s, x_s, q.on_s));
end
logln(logf, ['   midpoints 1 and N-1, which no row bounds: swing foot at %.3e and %.3e m; ' ...
             'true flow over (0,T): min %.3e m'], E.sw_hm(1), E.sw_hm(N-1), min(Hd(2:end-1)));
logln(logf, '   max|eta| along the true flow %.2e (u_ff holds it at its node-1 value)', max(ETAd));

%% 5. HZD at the report grid ---------------------------------------------
logln(logf, '');
logln(logf, ['5  HZD  delta^2 %.5f (41) / %.5f (161) | zeta* %.4f vs V_max/delta^2 %.4f ' ...
             '(161) | pre-impact thetadot: predicted %.4f, collocation %.4f'], ...
      Zs.delta2, Zr.delta2, Zr.zeta_star, Zr.V_max / Zr.delta2, Zr.thetadot_star, ...
      c_row * X(nq+1:end, N));

%% 6. closed loop ---------------------------------------------------------
logln(logf, '');
CL = struct([]);
cl_parts = {};
if opts.n_steps > 0
    for i = 1:numel(opts.controllers)
        pc = p;
        pc.controller = opts.controllers{i};
        pc.control_dt = opts.control_dt;
        pc.T_max      = 0.5;     % a step not at the guard in 1 s fails (ch3_table_rerun)
        t0  = tic;
        sim = ch3_simulate(X(:,1), alpha, pc, opts.n_steps);
        Vd  = ch3_validity(sim, pc);
        u_held = NaN;
        if sim.n_ok > 0
            Uh = [sim.steps.u];
            u_held = max(abs(Uh(:)));
        end
        first = '';
        if ~isnan(Vd.first_step)
            first = sprintf(', first %s in step %d at t = %.3f s', ...
                            Vd.first_kind, Vd.first_step, Vd.first_t);
        end
        CL(i).controller  = pc.controller;
        CL(i).steps       = sim.n_ok;
        CL(i).reason      = sim.reason;
        CL(i).valid_steps = Vd.valid_steps;
        CL(i).validity    = Vd;
        CL(i).u_held      = u_held;
        logln(logf, ['6  CLOSED LOOP %-9s %2d/%d steps (%s), valid %d%s | min Fz %.1f N, ' ...
                     'max mu %.4f | peak held |u| %.1f Nm (box %.1f) | %.0f s'], ...
              pc.controller, sim.n_ok, opts.n_steps, sim.reason, Vd.valid_steps, first, ...
              Vd.Fz_min, Vd.mu_max, u_held, pc.limits.u_max, toc(t0));
        cl_parts{end+1} = sprintf('%s %d/%d valid', pc.controller, ...
                                  Vd.valid_steps, opts.n_steps); %#ok<AGROW>
    end
end
cl_str = strjoin(cl_parts, ', ');
if isempty(cl_str), cl_str = 'closed loop not run'; end

rho = NaN;
if opts.stability
    pp = p;  pp.controller = 'iolin_pd';  pp.control_dt = 0;  pp.T_max = 0.5;
    t0 = tic;
    [rho, ~, pinfo] = ch3_poincare(X(:,1), alpha, pp);
    logln(logf, ['6  POINCARE  rho %.4f (iolin_pd, continuous; delta^2 %.4f is the ' ...
                 'zero-dynamics eigenvalue) | fixed-point residual %.2e | %.0f s'], ...
          rho, Zr.delta2, pinfo.fixed_point_residual, toc(t0));
end

%% summary -----------------------------------------------------------------
isP = [Q.phys];
onO = [Q.on_o];   okO = [Q.ok_o];   xO = [Q.x_o];
onS = [Q.on_s];   okS = [Q.ok_s];
fail_own  = find(onO & ~okO);
cross_own = find(onO & isfinite(xO));
fail_sP   = find(onS & ~okS & isP);
fail_sD   = find(onS & ~okS & ~isP);
never_on  = find(onS & ~onO);

logln(logf, '');
logln(logf, 'SUMMARY %s', name);
logln(logf, '  real trajectory             %s (verify %.1e, true-flow periodicity %.1e)', ...
      ternary(V.ok && st.ok, 'yes', 'NO'), V.max_dev, per_true);
logln(logf, '  own limits failed           %s', qlist(Q, fail_own, 'own'));
logln(logf, '  crossed by the true flow    %s', qlist(Q, cross_own, 'cross'));
logln(logf, '  spec, physical rows failed  %s', qlist(Q, fail_sP, 'spec'));
logln(logf, '  spec, design rows failed    %s', qlist(Q, fail_sD, 'spec'));
logln(logf, '  on in the spec, off in the solve: %s', qlist(Q, never_on, 'name'));
logln(logf, '  closed loop                 %s; Poincare rho %.4f', cl_str, rho);

R = struct();
R.file   = src;
R.name   = name;
R.verify = V;
R.true_flow = struct('ok', st.ok, 'T', st.T, 'L', st.L_step, 'impulse', st.impulse, ...
                     'periodicity', per_true);
R.rows   = struct('names', {names}, 'c_own', c_own, 'c_spec', c_spec, ...
                  'on_own', on_own, 'on_spec', on_spec, 'phys', phys_row);
R.Q      = Q;
R.hzd    = struct('delta2_41', Zs.delta2, 'delta2_161', Zr.delta2, 'nec4_161', Zr.nec4, ...
                  'nec5_161', Zr.nec5, 'zeta_star', Zr.zeta_star);
R.closed_loop = CL;
R.cl_str = cl_str;
R.rho    = rho;
R.key    = struct('v', v, 'u_pk', u_pk, 'I', norm(E.impulse), 'mu', mu_c, 'Fz', Fz_c, ...
                  'h_pk', h_pk, 'hip', [min(hip_c), max(hip_c)]);
R.fail   = struct('own', fail_own, 'cross', cross_own, 'spec_phys', fail_sP, ...
                  'spec_design', fail_sD, 'never_on', never_on, ...
                  'spec_phys_str', qlist(Q, fail_sP, 'spec'), ...
                  'spec_design_str', qlist(Q, fail_sD, 'spec'));
save(fullfile(OUT, [name '.mat']), 'R');
logln(logf, '=== VALIDATE_GAIT_DONE %s (%.0f s)', datestr(now), toc(t_all));
end

% ---------------------------------------------------------------------------
function compare(logf, R)
%COMPARE  One line per gait against the declared spec, then its failures.
logln(logf, '');
if isempty(R), logln(logf, 'no candidate gait on this machine'); return; end
d = ch3_params();
logln(logf, ['=== COMPARISON against ch3_params: |u| <= %.0f Nm, ||I|| <= %.0f Ns, ' ...
             'mu <= %.2f, Fz >= %.0f N, v = %.2f m/s, swing foot <= %.2f m, hip %.2f +- %.2f m'], ...
      d.limits.u_max, d.limits.impulse_max, d.limits.mu_s, d.limits.Fz_min, d.v_des, ...
      d.limits.clearance_max, d.limits.hip_h, d.limits.hip_h_tol);
logln(logf, '  %-16s %7s %8s %8s %7s %7s %7s %6s %11s %4s %4s %4s  %s', 'gait', 'v m/s', ...
      'verify', 'peak|u|', 'I Ns', 'mu', 'min Fz', 'swing', 'hip m', 'own', 'phys', 'dsgn', ...
      'closed loop; rho');
for i = 1:numel(R)
    r = R{i};  k = r.key;
    logln(logf, '  %-16s %7.4f %8.1e %8.1f %7.2f %7.4f %7.1f %6.3f %5.3f-%5.3f %4d %4d %4d  %s; %.4f', ...
          r.name, k.v, r.verify.max_dev, k.u_pk, k.I, k.mu, k.Fz, k.h_pk, k.hip(1), k.hip(2), ...
          numel(r.fail.own), numel(r.fail.spec_phys), numel(r.fail.spec_design), r.cl_str, r.rho);
end
logln(logf, ['  own / phys / dsgn: limits failed under the gait''s own p, and under the ' ...
             'spec''s physical and design rows']);
for i = 1:numel(R)
    logln(logf, '  %-16s spec physical: %s', R{i}.name, R{i}.fail.spec_phys_str);
    logln(logf, '  %-16s spec design:   %s', R{i}.name, R{i}.fail.spec_design_str);
end
end

% ---------------------------------------------------------------------------
function [names, gate, phys] = row_table()
%ROW_TABLE  The 19 rows of ch3_col_constraints, their gates, and which are
% physical (Table 3.1, contact, hybrid-model hypotheses) rather than design.
names = {'mid-step swing-foot clearance', 'swing foot >= 0 at interior nodes', ...
         'step-length floor', 'peak torque |u| <= u_max  (Table 3.1)', ...
         'friction cone |Fx| <= mu_s Fz  (NIC2)', 'normal force Fz >= Fz_min  (NIC1)', ...
         'impact impulse ||I|| <= impulse_max  (Table 3.1)', 'hip-height band', ...
         'swing foot strictly above ground  (NIC3)', 'transversal strike  (NIC3)', ...
         'post-impact lift-off  (NEC2)', 'impulse compressive Iz >= 0  (NEC3)', ...
         'impulse inside the friction cone  (NEC3)', 'fixed point exists  (NEC4)', ...
         'fixed point stable  (NEC5)', 'theta strictly monotonic  (HH6)', ...
         'decoupling matrix invertible  (HH2)', 'swing-foot height ceiling', ...
         'step-length ceiling'};
gate  = {'clearance', '', '', 'torque', 'friction', 'grf', 'impulse', 'height', ...
         'swing_clear', 'swing_clear', 'liftoff', 'impact', 'impact', 'hzd', 'hzd', ...
         'phase_mono', 'decoupling', 'clearance_max', 'step_len_max'};
phys  = true(1, 19);
phys([1 8 18 19]) = false;
end

function ps = spec_params(p, d)
%SPEC_PARAMS  The gait's own parametrization under the DECLARED requirements.
ps = p;
ps.limits       = d.limits;          % every Table 3.1 / 6.3.4 value AND gate
ps.v_des        = d.v_des;
ps.enforce_nec1 = d.enforce_nec1;
ps.step_len_min = d.step_len_min;
ps.qt_range     = d.qt_range;
ps.T_min  = d.T_min;   ps.T_max = d.T_max;   ps.dq_max = d.dq_max;
ps.hzd_tol = d.hzd_tol;
end

function q = gates_on(q, gate)
for r = 1:numel(gate)
    if ~isempty(gate{r}), q.limits.enable.(gate{r}) = true; end
end
end

function on = gate_on(p, g)
on = isempty(g) || logical(p.limits.enable.(g));
end

function on = nec1_on(p)
on = ~isfield(p, 'enforce_nec1') || logical(p.enforce_nec1);
end

function b = bound_violation(z, E, p, N)
%BOUND_VIOLATION  Worst excursion past ch3_col_bounds: [nodes, midpoints].
% fmincon holds the nodes inside; nothing holds the midpoints.
[lb, ub] = ch3_col_bounds(p, N);
b = [max([lb - z(:); z(:) - ub]), ...
     max(max([E.xm - ub(1:p.nx); lb(1:p.nx) - E.xm]))];
end

function r = fric_ratio(lam)
pos = lam(2,:) > 0;
if any(pos), r = max(abs(lam(1,pos)) ./ lam(2,pos)); else, r = Inf; end
end

function Q = addq(Q, name, unit, sense, col, roll, lim_o, on_o, lim_s, on_s, phys)
q = struct('name', name, 'unit', unit, 'sense', sense, 'col', col, 'roll', roll, ...
           'lim_o', lim_o, 'on_o', logical(on_o), 'lim_s', lim_s, 'on_s', logical(on_s), ...
           'phys', phys);
if isempty(Q), Q = q; else, Q(end+1) = q; end
end

function ok = meets(v, lim, sense, tol)
switch sense
    case '<=', ok = v <= lim + tol;
    case '>=', ok = v >= lim - tol;
    otherwise, ok = abs(v - lim) <= tol;
end
end

function [ok, over] = judge(col, roll, lim, sense, tol)
%JUDGE  Verdict on the collocation value, and how far the rollout crosses the
% limit when the collocation meets it (NaN when it does not cross). The speed
% equality is judged on the collocation only: the rollout's step differs from
% the nodes' by the verify deviation, which would fail a 1e-6 equality anyway.
ok   = meets(col, lim, sense, tol);
over = NaN;
if ok && ~strcmp(sense, '=') && isfinite(roll) && ~meets(roll, lim, sense, tol)
    if strcmp(sense, '<='), over = roll - lim; else, over = lim - roll; end
end
end

function s = vstr(ok, over, on)
if ok, s = 'ok'; elseif on, s = 'FAIL'; else, s = 'over'; end
if isfinite(over), s = [s '*']; end
end

function s = rstr(c, on, tol)
if c <= tol, s = 'ok'; elseif on, s = 'FAIL'; else, s = 'over'; end
end

function s = etag(on)
if on, s = '[E]'; else, s = '[ ]'; end
end

function s = qlist(Q, idx, which)
if isempty(idx), s = 'none'; return; end
parts = cell(1, numel(idx));
for k = 1:numel(idx)
    q = Q(idx(k));
    switch which
        case 'own',   parts{k} = sprintf('%s %.4g %s %.4g', q.name, q.col, q.sense, q.lim_o);
        case 'spec',  parts{k} = sprintf('%s %.4g %s %.4g', q.name, q.col, q.sense, q.lim_s);
        case 'cross', parts{k} = sprintf('%s by %.2e %s', q.name, q.x_o, q.unit);
        otherwise,    parts{k} = q.name;
    end
end
s = strjoin(parts, '; ');
end

function src = find_gait(ROOT, f)
src = fullfile(ROOT, f);
if exist(src, 'file'), return; end
src = f;
if exist(src, 'file'), return; end
src = '';
end

function o = with_defaults(o)
dflt = struct('n_steps', 10, 'control_dt', 1e-3, 'stability', true, ...
              'n_dense', 4001, 'tol', 1e-6);
dflt.controllers = {'iolin_pd', 'clfqp_con'};
fn = fieldnames(dflt);
for i = 1:numel(fn)
    if ~isfield(o, fn{i}) || isempty(o.(fn{i})), o.(fn{i}) = dflt.(fn{i}); end
end
end

function o = ternary(c, a, b)
if c, o = a; else, o = b; end
end

function logln(f, fmt, varargin)
s = sprintf(fmt, varargin{:});
fprintf('%s\n', s);
fid = fopen(f, 'a'); fprintf(fid, '%s\n', s); fclose(fid);
end
