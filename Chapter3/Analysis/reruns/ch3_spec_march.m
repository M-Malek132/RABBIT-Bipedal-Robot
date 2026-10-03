function ch3_spec_march(seed_file, tag, stages, free_theta)
%CH3_SPEC_MARCH  Bring the best physical gait onto ch3_params' design rows, in three stages.
%
%   ch3_spec_march()                    from impulse_march/u162_I1500, stages 1-4
%   ch3_spec_march([], [], 1)           stage 1 only (swing-foot ceiling)
%   ch3_spec_march([], [], 2:3)         resume at stage 2 from the stored rungs
%   ch3_spec_march(seed_file, tag)      from any stored gait (z, z_try or z_opt, and p)
%   ch3_spec_march(seed_file, tag, 2, true)   hip stage with the phase endpoints free
%
% FREE THETA (4th argument). The 2026-10-03 run stalled stage 2 at a hip top of
% 0.9575 m, and the stride identity says why: with theta fixed at
% (-0.15, 0.30), L = 0.4604 h_imp, the stride does not move (HANDOFF, "Stride"),
% so the hip AT IMPACT is pinned at L/0.4604 = 0.949 m and the band's top can
% never go under it. The 0.92 target needs tan(th+) - tan(th-) >= L/0.92, i.e.
% a slightly wider sweep of the stance leg; free_theta makes the pair decision
% variables inside p.theta_bounds. Untested in a solve before this driver.
%
% WHY. Under the 200 Nm box (ch3_params since 2026-10-02) u162_I1500 passes
% every PHYSICAL row of ch3_validate_gait: Table 3.1, NIC, NEC, HH, verified at
% 1.2e-4, Poincare rho 0.72. It fails three DESIGN rows of the spec:
%   row 18  swing-foot ceiling   0.321 m  vs <= 0.15  (gate never on in its lineage)
%   row 8   hip-height band      0.940-0.969 m  vs 0.88-0.92 (its own p: 0.88-0.97)
%   NEC1    walking speed        1.393 m/s  vs = 1.2  (gate off since posture_195)
% Its torque sits on its own 162.5 box; the 37.5 Nm up to 200 is the new room.
%
% THE STAGES, IN THIS ORDER, EACH HOLDING WHAT THE ONES BEFORE IT REACHED.
%   1  ceiling  0.32 -> 0.15 m, 3 cm rungs. Cheapest first: the swing leg is
%      5-17% of the torque-squared cost, the foot goes high only because
%      nothing opposed it (memory: row 1 is inert, row 18 was added for this).
%   2  hip top  0.97 -> 0.92 m, 1 cm rungs, the band's floor held at 0.88. The
%      whole gait sits above the spec band, so the hip has to come down ~5 cm.
%      Row 8's comment says a low hip buys knee torque cheaply, so this may not
%      cost torque at all; with theta fixed the stride shrinks with the hip
%      (L = h_imp (tan th+ - tan th-)), about 0.438 -> 0.43 m.
%   3  speed    1.39 -> 1.20 m/s, 0.05 m/s rungs. NEC1 is switched on AT THE
%      SPEED THE GAIT WALKS (ch3_speed_ladder: with v_des below it, SQP's first
%      subproblem is infeasible and nothing moves), then v_des is marched. Last
%      because it is the risky one: on the torque branch the speed ran UP as
%      limits tightened, and pinning it back collapsed the solve.
%   4  stance friction  0.40 -> 0.32, 0.02 rungs (added 2026-10-03). Not a
%      spec row: a PLANNING margin. Stages 1-3 left the gait exactly on the
%      cone and it slips in step 1 under iolin_pd at 1 kHz (closed-loop demand
%      0.48). Only stance mu_s is tightened; the impact cone stays at 0.40.
%      Judge the result by ch3_validate_gait's closed loop, not by the target.
%
% WHAT ELSE CHANGES FROM THE SEED'S p. The four Table 3.1 values come from
% ch3_params (box 200, impulse 15, mu 0.4, Fz 50) with every physical gate on,
% and the torso box widens to the spec's qt_range (the seed's [0.03 0.45] is a
% leftover of the posture march). Both only loosen the seed's own problem, so
% the seed stays feasible. The mesh stays at the seed's N.
%
% EACH RUNG: a warm-started ch3_col_solve, keeping the lowest-cost iterate that
% already meets every limit if the final iterate misses (ch3_impulse_march,
% "keep the best iterate that already lands"), a polish of a near miss, then
% the landing test max|ceq| <= 1e-6, max c <= 1e-4 and ch3_col_verify. A rung
% that misses is bisected; the step grows back x1.5 after a landing. A stage
% stops when the step falls below its MIN_STEP, and the next stage starts from
% the last landed gait with the stalled stage's limit held where it got to --
% the verdict says which targets were met.
%
% A STALL IS NOT A FLOOR: it is a statement about continuation from this seed
% with this solver. If stage 2 stalls, free theta (p.free_theta, fixed for
% wrench gaits on 2026-09-29) is the lever the stride identity points to.
%
% Output: Results/reruns/spec_march/ -- one .mat per rung, landed or not
% (z_try, p, out, V, landed), spec_march.log, and <tag>_spec.mat (z, p) when
% all three targets are met. A rerun resumes: landed rungs are replayed and a
% stored miss is taken as a miss without re-solving (the solve is
% deterministic); delete that file to retry it.
%
% See also CH3_IMPULSE_MARCH, CH3_SPEED_LADDER, CH3_VALIDATE_GAIT.

if nargin < 1 || isempty(seed_file)
    seed_file = fullfile('Results', 'reruns', 'impulse_march', 'u162_I1500.mat');
end
if nargin < 2 || isempty(tag),    tag    = 'u162_I1500'; end
if nargin < 3 || isempty(stages), stages = 1:4; end
if nargin < 4 || isempty(free_theta), free_theta = false; end

ROOT = fileparts(fileparts(fileparts(fileparts(mfilename('fullpath')))));
SPD  = fullfile(ROOT, 'Results', 'reruns', 'spec_march');
if ~exist(SPD, 'dir'), mkdir(SPD); end
logf = fullfile(SPD, 'spec_march.log');

src = fullfile(ROOT, seed_file);
if ~exist(src, 'file')
    logln(logf, '=== spec march: seed %s not found (Results/reruns/ is git-ignored)', seed_file);
    error('ch3_spec_march:noSeed', 'Seed gait not found: %s', src);
end

S = load(src);
if isfield(S, 'z_try'), z = S.z_try; elseif isfield(S, 'z_opt'), z = S.z_opt; else, z = S.z; end
p = ch3_upgrade_params(S.p);
d = ch3_params();

p.scale_problem = true;                  % off, SQP runs away from a converged gait
% NEC1 stays as the seed saved it: off for u162_I1500 until stage 3 switches it
% on, ON for a seed that already holds the design speed (the 2026-10-03 run's
% s3 gait), so a later hip stage cannot trade the speed away.
if free_theta
    p.free_theta = true;
    z = ch3_col_theta_augment(z, p);     % same gait, theta pair appended to z
end
for g = {'torque', 'impulse', 'friction', 'grf', 'swing_clear', 'liftoff', ...
         'impact', 'hzd', 'phase_mono', 'decoupling'}
    p.limits.enable.(g{1}) = true;
end
p.limits.u_max       = d.limits.u_max;
p.limits.impulse_max = d.limits.impulse_max;
if isempty(p.limits.mu_s_impact)         % the impact cone stays at the spec's
    p.limits.mu_s_impact = d.limits.mu_s; % 0.40 when stage 4 tightens stance
end
p.limits.mu_s        = min(p.limits.mu_s, d.limits.mu_s);   % keep a stage-4 margin
p.limits.Fz_min      = d.limits.Fz_min;
p.qt_range           = d.qt_range;
hip_lo = d.limits.hip_h - d.limits.hip_h_tol;          % 0.88, the band's floor

ST = stage_table(d, hip_lo);

logln(logf, '=== spec march | seed %s | stages %s | box %.1f Nm, impulse %.1f Ns, mu %.2f, Fz %.0f N | %s', ...
      seed_file, mat2str(stages), p.limits.u_max, p.limits.impulse_max, ...
      p.limits.mu_s, p.limits.Fz_min, datestr(now));
report(logf, p, z, 'seed', NaN);

% A stage that is skipped (resume at 2 or 3) still has to leave its limit on
% the p the later stages solve with: replay its landed rungs from disk.
met = false(1, numel(ST));
for k = 1:numel(ST)
    if ismember(k, stages)
        [z, p, met(k)] = run_stage(ST(k), k, z, p, SPD, tag, logf);
    else
        [z, p, met(k)] = replay_stage(ST(k), k, z, p, SPD, tag, logf);
    end
end

names = {ST.name};
if all(met)
    f = fullfile(SPD, sprintf('%s_spec.mat', tag));
    save(f, 'z', 'p');
    logln(logf, ['=== VERDICT: a verified gait meets every physical row at the %.0f Nm ' ...
                 'box AND the design rows (ceiling %.2f m, hip %.2f-%.2f m, %.2f m/s, ' ...
                 'stance friction planned at %.2f): %s'], ...
          p.limits.u_max, ST(1).target, hip_lo, ST(2).target, ST(3).target, ST(4).target, f);
else
    logln(logf, ['=== VERDICT: targets met: %s; not met: %s. The last landed gait holds ' ...
                 'every limit at the value its stage reached (see the rung lines). A ' ...
                 'stall is a continuation result from this seed, not a floor.'], ...
          strjoin(names(met), ', '), strjoin(names(~met), ', '));
end
logln(logf, '=== SPEC_MARCH_DONE %s', datestr(now));
fprintf('SPEC_MARCH_DONE\n');
end

% ---------------------------------------------------------------------------
function ST = stage_table(d, hip_lo)
%STAGE_TABLE  What each stage measures, how it sets its limit, and its steps.
ST(1) = struct('name', 'ceiling', 'unit', 'm', 'target', d.limits.clearance_max, ...
               'step', 0.03, 'min_step', 0.005, 'scale', 1000, ...
               'measure', @(E) max([E.sw_h, E.sw_hm]), ...
               'apply', @set_ceiling, 'equality', false);
ST(2) = struct('name', 'hip', 'unit', 'm', 'target', d.limits.hip_h + d.limits.hip_h_tol, ...
               'step', 0.01, 'min_step', 0.002, 'scale', 1000, ...
               'measure', @(E) max(-[E.X(2,:), E.xm(2,:)]), ...
               'apply', @(p, v) set_hip(p, hip_lo, v), 'equality', false);
ST(3) = struct('name', 'speed', 'unit', 'm/s', 'target', d.v_des, ...
               'step', 0.05, 'min_step', 0.01, 'scale', 1000, ...
               'measure', @(E) E.L_step / E.T, ...
               'apply', @set_speed, 'equality', true);
% 4: a planning margin on STANCE friction. The 2026-10-03 validation: s1 and s3
% sit exactly on mu 0.40 and slip in step 1 under iolin_pd at 1 kHz, at a
% closed-loop demand of 0.48 -- tracking adds ~0.08 late in stance. 0.32 plans
% for that; the impact cone (row 13) stays at the spec's 0.40.
ST(4) = struct('name', 'friction', 'unit', '-', 'target', 0.32, ...
               'step', 0.02, 'min_step', 0.005, 'scale', 1000, ...
               'measure', @stance_mu, ...
               'apply', @(p, v) set_friction(p, v, d.limits.mu_s), 'equality', false);
end

function m = stance_mu(E)
lam = [E.lam, E.lamm];
m   = max(abs(lam(1,:)) ./ lam(2,:));
end

function p = set_friction(p, v, mu_impact)
p.limits.enable.friction = true;
p.limits.mu_s            = v;
if isempty(p.limits.mu_s_impact), p.limits.mu_s_impact = mu_impact; end
end

function p = set_ceiling(p, v)
p.limits.enable.clearance_max = true;
p.limits.clearance_max        = v;
end

function p = set_hip(p, lo, hi)
p.limits.enable.height = true;
p.limits.hip_h         = (lo + hi) / 2;
p.limits.hip_h_tol     = (hi - lo) / 2;
end

function p = set_speed(p, v)
p.enforce_nec1 = true;
p.v_des        = v;
end

% ---------------------------------------------------------------------------
function [z, p, met] = run_stage(s, k, z, p, SPD, tag, logf)
%RUN_STAGE  March one limit from the gait's own value to the stage target.
cur = s.measure(ch3_col_eval(z, p));
logln(logf, '--- stage %d: %s | gait %.4f %s -> target %.4f, step %.3f (min %.3f)', ...
      k, s.name, cur, s.unit, s.target, s.step, s.min_step);
sgn = sign(s.target - cur);

% GATE ON AT THE GAIT'S OWN VALUE before any step: the gait meets it exactly,
% so no solve is needed, and a stage whose first rung misses still leaves its
% limit on for the stages after it. For NEC1 this is also the switch-on at the
% achieved speed that ch3_speed_ladder found necessary.
p = s.apply(p, cur);
logln(logf, '  %s gate on at the gait''s own value %.4f %s', s.name, cur, s.unit);
if (sgn >= 0 && ~s.equality) || abs(s.target - cur) < 1e-4
    p   = s.apply(p, s.target);
    met = true;
    logln(logf, '  already within the target; %s limit set to %.4f %s', s.name, s.target, s.unit);
    return;
end

step = s.step;
while sgn * (s.target - cur) > 1e-9
    v_try = cur + sgn * min(step, sgn * (s.target - cur));
    f     = rung_file(SPD, tag, k, s, v_try);
    if exist(f, 'file')
        L = load(f);
        if L.landed
            z = L.z_try; p = L.p; cur = v_try;
            report(logf, p, z, sprintf('%s %.4f (stored)', s.name, v_try), NaN);
            continue;
        end
        logln(logf, '  stored rung %s %.4f is a miss; not re-solving it (delete %s to retry)', ...
              s.name, v_try, f);
        landed = false;
    else
        [z_try, p_try, out, V, landed] = solve_rung(s.apply(p, v_try), z, ...
                                                    sprintf('%s %.4f', s.name, v_try), logf);
        save_rung(f, z_try, p_try, out, V, landed);
    end
    if landed
        z = z_try; p = p_try; cur = v_try;
        step = min(step * 1.5, s.step);          % grow again once it is moving
        % SKIP AHEAD when the gait overshot its own limit. The 2026-10-03 run's
        % first friction rung (limit 0.38) landed at a stance mu of 0.240, and
        % the march still spent 2 h solving 0.36 and 0.34, both inactive.
        m = s.measure(ch3_col_eval(z, p));
        if ~s.equality && sgn * (m - cur) > 0
            logln(logf, '  the gait is already at %.4f %s, past its %.4f limit; marching on from there', ...
                  m, s.unit, cur);
            cur = m;
        end
    else
        step = step / 2;
        if step < s.min_step
            logln(logf, '  stage %d stalls at %s %.4f %s: the step fell below %.3f', ...
                  k, s.name, cur, s.unit, s.min_step);
            break;
        end
        logln(logf, '  bisecting: next step %.4f %s', step, s.unit);
    end
end
met = sgn * (s.target - cur) <= 1e-9;
if met && ~s.equality
    p = s.apply(p, s.target);                % the gait is at or past it
end
logln(logf, '--- stage %d %s: %s at %.4f %s', k, s.name, ...
      ternary(met, 'TARGET MET', 'stalled'), cur, s.unit);
end

% ---------------------------------------------------------------------------
function [z, p, met] = replay_stage(s, k, z, p, SPD, tag, logf)
%REPLAY_STAGE  A stage not asked for: take its last landed rung off disk.
files = dir(fullfile(SPD, sprintf('%s_s%d_%s_*.mat', tag, k, s.name)));
best  = []; cur = s.measure(ch3_col_eval(z, p)); dir_ = sign(s.target - cur);
for i = 1:numel(files)
    L = load(fullfile(files(i).folder, files(i).name));
    if ~L.landed, continue; end
    v = s.measure(ch3_col_eval(L.z_try, L.p));
    if isempty(best) || dir_ * (v - best.v) > 0
        best = struct('z', L.z_try, 'p', L.p, 'v', v);
    end
end
if isempty(best)
    met = abs(s.target - cur) <= 1e-4 || (~s.equality && dir_ >= 0);
    if met, p = s.apply(p, s.target); else, p = s.apply(p, cur); end
    logln(logf, '--- stage %d %s: skipped, no landed rung on disk; gait at %.4f %s', ...
          k, s.name, cur, s.unit);
    return;
end
z = best.z; p = best.p;
met = abs(s.target - best.v) <= 1e-3 * max(1, abs(s.target));
report(logf, p, z, sprintf('%s (replayed)', s.name), NaN);
end

% ---------------------------------------------------------------------------
function [z_try, p, out, V, landed] = solve_rung(p, z, label, logf)
%SOLVE_RUNG  One rung, with the harvest and polish of ch3_impulse_march.
tracker = containers.Map();
tracker('best')  = struct('z', [], 'J', inf, 'iter', NaN, 'max_c', NaN, 'max_ceq', NaN);
tracker('since') = 0;
opts = struct('MaxFunctionEvaluations', 4e5, 'MaxIterations', 400, ...
              'OutputFcn', @(zz, ov, state) track_landed(zz, ov, state, p, tracker));
t0 = tic;
[z_try, out] = ch3_col_solve(p, z, opts);
secs = toc(t0);
logln(logf, '%-18s exitflag %d | max|c| %.2e | max|ceq| %.2e | %.0f s', ...
      label, out.exitflag, out.max_c, out.max_ceq, secs);

best = tracker('best');
out.harvested = NaN;
if ~(out.max_ceq <= 1e-6 && out.max_c <= 1e-4) && ~isempty(best.z)
    logln(logf, ['  the final iterate misses; iterate %d already met every limit ' ...
                 '(max|c| %.2e, max|ceq| %.2e, J %.2f) and is used instead'], ...
          best.iter, best.max_c, best.max_ceq, best.J);
    z_try = best.z;
    out.max_c = best.max_c;  out.max_ceq = best.max_ceq;  out.harvested = best.iter;
end
V = report(logf, p, z_try, label, secs);

if isnan(out.harvested) && out.max_c > 1e-6 && out.max_c < 1e-4 && V.ok
    logln(logf, '  polishing %s (max|c| %.2e)', label, out.max_c);
    t1 = tic;
    [z_p, out_p] = ch3_col_solve(p, z_try, struct('MaxFunctionEvaluations', 4e5, ...
                                                  'MaxIterations', 400));
    V_p = report(logf, p, z_p, [label ' polished'], toc(t1));
    if V_p.ok && out_p.max_c < out.max_c
        out_p.harvested = NaN;
        z_try = z_p; out = out_p; V = V_p;
    end
end
landed = V.ok && out.max_ceq <= 1e-6 && out.max_c <= 1e-4;
end

% ---------------------------------------------------------------------------
function stop = track_landed(z, ov, state, p, tracker)
%TRACK_LANDED  fmincon OutputFcn: remember the lowest-cost iterate that meets
% the landing test in true units; stop once PATIENCE iterations bring no better.
PATIENCE = 30;
stop = false;
if ~strcmp(state, 'iter'), return; end
best = tracker('best');
if ov.fval < best.J
    [c, ceq] = ch3_col_constraints(z, p);
    if max(c) <= 1e-4 && max(abs(ceq)) <= 1e-6
        tracker('best')  = struct('z', z, 'J', ov.fval, 'iter', ov.iteration, ...
                                  'max_c', max(c), 'max_ceq', max(abs(ceq)));
        tracker('since') = ov.iteration; %#ok<NASGU> containers.Map is a handle
        return;
    end
end
if ~isempty(best.z) && ov.iteration - tracker('since') >= PATIENCE
    stop = true;
end
end

% ---------------------------------------------------------------------------
function save_rung(f, z_try, p, out, V, landed)
save(f, 'z_try', 'p', 'out', 'V', 'landed');
end

function f = rung_file(SPD, tag, k, s, v)
f = fullfile(SPD, sprintf('%s_s%d_%s_%04.0f.mat', tag, k, s.name, s.scale * v));
end

function V = report(logf, p, z, label, secs)
E   = ch3_col_eval(z, p);
c   = ch3_col_constraints(z, p);
chk = ch3_col_check_limits(z, p);
V   = ch3_col_verify(z, p, false);
lam = [E.lam, E.lamm];
pk  = max([abs(E.u(:)); abs(E.um(:))]);
h   = -[E.X(2,:), E.xm(2,:)];
[cw, iw] = max(c);
[~, ~, ~, th] = ch3_col_unpack(z, p);
label = sprintf('%s [th %+.3f %+.3f]', label, th(1), th(2));
logln(logf, ['%-22s v %.4f | T %.4f | L %.4f | peak|u| %6.1f of %5.1f | impulse %5.2f | ' ...
             'mu %.4f | minFz %6.1f | foot apex %.3f | hip %.3f-%.3f | limits ok %d ' ...
             '(max c %.2e, worst row %d) | verify %d (%.1e) | %.0f s'], ...
      label, E.L_step/E.T, E.T, E.L_step, pk, p.limits.u_max, norm(E.impulse), ...
      max(abs(lam(1,:))./lam(2,:)), min(lam(2,:)), max([E.sw_h, E.sw_hm]), ...
      min(h), max(h), chk.ok, cw, iw, V.ok, V.max_dev, secs);
end

function s = ternary(b, x, y)
if b, s = x; else, s = y; end
end

function logln(f, fmt, varargin)
s = sprintf(fmt, varargin{:});
fprintf('%s\n', s);
fid = fopen(f, 'a'); fprintf(fid, '%s\n', s); fclose(fid);
end
