function ch3_speed_march(iters, max_stage)
%CH3_SPEED_MARCH  Sweep walking speed 0.35 -> 1.2 m/s holding lean + tall.
%
%   ch3_speed_march              full sweep, per-rung default budgets
%   ch3_speed_march(2, 1)        smoke test: 2 iterations, one rung
%
% Produces a gait FAMILY spanning 0.35 to 1.2 m/s in which every member keeps
% the posture won by ch3_lean_tall_march: torso pitched FORWARD (qt box
% [0.08 0.25] rad) and hip held TALL in the band [0.915 0.950] m.
%
% WHERE THE SWEEP STARTS, AND WHY IT IS NOT 0.35.  The seed is stage 3 of
% Results/ch3_talllean_state.mat -- an N = 61 gait verified at 2.9e-05 that
% already walks at 0.9781 m/s with qt in [+0.080 +0.142] and the hip in
% [0.9280 0.9500]. It sits in the MIDDLE of the requested range, not at its
% bottom, so the sweep runs in two phases OUT from it rather than in one pass
% along it:
%
%   DOWN  0.978 -> 0.35   fills the slow half
%   UP    0.978 -> 1.20   fills the fast half (restarts from the same seed)
%
% Marching down to 0.35 and then back up through territory already solved
% would double the cost of the slow half to buy nothing: continuation cares
% that each solve starts near-feasible, not that the targets are monotone.
% The two phases are stitched into one speed-ordered library at the end.
%
% ENFORCE_NEC1 MUST BE ON HERE, AND THE SEED HAS IT OFF.  Speed only becomes a
% constraint when NEC1 (L_step/T = v_des) is gated on; ch3_continuation sets
% p.v_des but never touches the gate. Handed the seed's params verbatim it
% would run every rung as the same free-speed problem, report the same 0.978
% m/s each time, and look like a march that simply refused to move. The gate
% is switched on below, at the achieved speed, which is exactly the handoff
% ch3_col_constraints documents: solve cold with NEC1 off, then turn it on and
% march v_des.
%
% SCALEPROBLEM MUST BE ON, AND ch3_col_solve TURNS IT OFF.  This is the single
% thing that makes the sweep work at all, so do not quietly drop it. Measured on
% this seed, with scaling OFF:
%
%   re-solving the gait against its OWN UNCHANGED params diverges. Not a speed
%   march, not NEC1, nothing altered: fmincon's first step has norm 15.3,
%   feasibility goes 1.9e-08 -> 3.1e-01, the hip climbs out of its own band to
%   0.973 m and J reaches 8311 in 15 iterations.
%
% The mechanism is scaling. SQP starts from B = I, and ||grad J|| ~ 6.8e3 here,
% so the first QP step wants a norm in the thousands and the line search only
% cuts it to ~15 -- far outside the region where the linearisation holds. Once
% there, gradients read ~1e10 and the solve never returns. The same seed with
% ScaleProblem = true takes first steps of norm ~0.3, drops feasibility
% monotonically, and lands a 0.02 m/s rung dead on target (J 2232 -> 2227,
% verify 3.7e-04) in 15 iterations and 177 s.
%
% ch3_col_solve sets ScaleProblem = false on purpose, but its stated reason is
% diagnostic readability -- keeping the reported Feasibility comparable to
% ConstraintTolerance rather than to a rescaled surrogate. That is a real cost:
% the Feasibility and first-order optimality columns in the log below are
% scaled surrogates and should not be read against p.verify_tol. Read
% ch3_col_verify's max_dev instead, which is computed from the trajectory and
% is unaffected. Convergence beats readability here.
%
% Two other cures were tried on this seed and did NOT work, so they are not
% worth re-running: interior-point (J -> 15774, verify 2.6) and a remesh to
% N = 81 (verify 5.3e-01, and the remesh alone already broke verification at
% 2.9e-03). The N=81 recipe in the project notes fixes UNDER-RESOLUTION, which
% is a different failure from this one.
%
% A MISSED RUNG BISECTS RATHER THAN STOPPING.  ch3_continuation halts the whole
% march the moment a speed fails to verify, on the sound argument that
% continuing from a fictional gait propagates it. That is right for its own
% one-shot use and too brittle for an 18-rung sweep, where one awkward speed
% should cost a smaller step rather than the run. Each rung here retries at the
% midpoint between the last ACHIEVED speed and the target, up to MAX_BISECT
% times; only if a halved step still misses does the phase stop, and the other
% phase still runs.
%
% EXPECT THE FAST END TO BE THE HARD ONE. Tall and fast pull against each
% other: a longer step at fixed leg length wants a LOWER hip, and the band is
% already parked at 0.94-0.95 m against a 1.0 m leg (94-95% extension, right at
% the conditioning limit ch3_params warns about). The height band is therefore
% held FIXED at the one ch3_lean_tall_march actually verified -- its stage 4
% attempt to raise the ceiling to 0.955 failed twice (dev 7.9e+01, then
% 1.4e+01) -- so this driver moves ONE knob, the speed, and leaves posture
% alone. Raising height and speed together is what that failed stage already
% showed does not work.
%
% Resumable: state is saved after every rung to
% Results/ch3_speed_march_state.mat and progress appended to
% Results/ch3_speed_march.log, so a crash costs one rung rather than the run.
% Re-running picks up at the first unfinished rung.
%
% Writes Results/ch3_gait_library_lean_tall.mat (the whole family),
% Results/ch3_gait_lean_tall_v120.mat (the fast end, with a full ch3_report)
% and Results/ch3_walk_lean_tall_v120.gif.
%
% See also CH3_CONTINUATION, CH3_LEAN_TALL_MARCH, CH3_COL_VERIFY.

if nargin < 1, iters = []; end

MAX_BISECT = 3;
V_LO = 0.35;
V_HI = 1.20;
DV   = 0.05;        % nominal rung size [m/s]

% WHERE A RUNG ACTUALLY CONVERGES. Measured on the first 0.05 m/s rung run at a
% 120-iteration cap: feasibility reaches 1.7e-06 by iteration 20 and sits at
% ~1e-07 from there on, while J drifts down only 2% (2303 -> 2256) over the
% remaining 100 iterations. The constraint work -- the part a continuation
% needs -- is done early; the rest is polish, and it tripled the rung's cost
% (1254 s). 60 keeps 3x margin over where feasibility settles while cutting the
% sweep from ~6 h to ~3 h. It also stays well clear of letting a diverging rung
% run long, which is what the unscaled failure did for 400 iterations / 3 h.
ITERS_DEFAULT = 60;

% J_GUARD stops a rung the moment the cost runs away, which is what a diverging
% solve does here long before it exhausts its iterations. A healthy rung barely
% moves J (2232 -> 2227 measured); a diverging one passed 6e4 by iteration 40.
J_GUARD = 10;       % multiples of the rung's starting cost

root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
resd = fullfile(root, 'Results');
if ~exist(resd, 'dir'), mkdir(resd); end

STATE = fullfile(resd, 'ch3_speed_march_state.mat');
LOG   = fullfile(resd, 'ch3_speed_march.log');
SEED  = fullfile(resd, 'ch3_talllean_state.mat');
LIB   = fullfile(resd, 'ch3_gait_library_lean_tall.mat');
FINAL = fullfile(resd, 'ch3_gait_lean_tall_v120.mat');
GIF   = fullfile(resd, 'ch3_walk_lean_tall_v120.gif');

% --- start or resume ------------------------------------------------------
if exist(STATE, 'file')
    S = load(STATE);
    k0 = S.k_done + 1;  z = S.z;  p = S.p;  p0 = S.p0;  z0 = S.z0;
    lib = S.lib;  stages = S.stages;
    ch3_logln(LOG, sprintf('RESUME at rung %d/%d', k0, numel(stages)));
else
    if ~exist(SEED, 'file')
        error('ch3_speed_march:seed', ...
              ['Cannot find the lean+tall seed "%s". Run ch3_lean_tall_march ' ...
               'first.'], SEED);
    end
    S  = load(SEED);
    z0 = S.z;
    p0 = ch3_upgrade_params(S.p);

    % The posture this sweep is required to hold. Read off the seed rather
    % than hard-coded so a re-run of the posture march upstream carries
    % through, but asserted below so a seed that quietly lost the lean fails
    % loudly here instead of producing a sweep that is not lean+tall at all.
    p0.enforce_nec1         = true;      % <-- without this the march is a no-op
    p0.limits.enable.height = true;

    if p0.qt_range(1) < 0
        error('ch3_speed_march:notlean', ...
              ['Seed pitch box is [%+.3f %+.3f]; a box straddling zero comes ' ...
               'back leaning BACKWARD (ch3_params SIGN note). Expected an ' ...
               'entirely positive box such as [0.08 0.25].'], ...
              p0.qt_range(1), p0.qt_range(2));
    end

    E0 = ch3_col_eval(z0, p0);
    V0 = ch3_col_verify(z0, p0, false);
    v0 = E0.L_step / E0.T;

    if ~V0.ok
        error('ch3_speed_march:seedunverified', ...
              ['Seed gait does not verify (%.3e > %.1e). Marching from a gait ' ...
               'that is not a real trajectory only propagates it.'], ...
              V0.max_dev, p0.verify_tol);
    end

    stages = build_stages(v0, V_LO, V_HI, DV);

    ch3_logln(LOG, sprintf('=== SPEED sweep %.2f -> %.2f m/s, lean+tall held, %d rungs ===', ...
                           V_LO, V_HI, numel(stages)));
    ch3_logln(LOG, sprintf(['    seed: N=%d v=%.4f qt=[%+.4f %+.4f] hip=[%.4f %.4f] ' ...
                            'verify %.3e (ok=%d)'], size(E0.X,2), v0, ...
                           min(E0.X(3,:)), max(E0.X(3,:)), ...
                           min(-E0.X(2,:)), max(-E0.X(2,:)), V0.max_dev, V0.ok));
    ch3_logln(LOG, sprintf('    holding qt box [%+.3f %+.3f] rad, hip band [%.3f %.3f] m', ...
                           p0.qt_range(1), p0.qt_range(2), ...
                           p0.limits.hip_h - p0.limits.hip_h_tol, ...
                           p0.limits.hip_h + p0.limits.hip_h_tol));

    p = p0;  z = z0;
    lib = {snap(z0, p0, v0, 'seed')};
    k0 = 1;
end

if nargin >= 2 && ~isempty(max_stage)
    stages = stages(1:min(max_stage, numel(stages)));
end

% --- the sweep ------------------------------------------------------------
skip_phase = '';
for k = k0:numel(stages)
    st = stages(k);
    t0 = tic;

    if strcmp(st.phase, skip_phase)
        ch3_logln(LOG, sprintf('--- rung %d/%d [%s] v=%.4f  SKIPPED (phase stopped)', ...
                               k, numel(stages), st.phase, st.v));
        continue;
    end

    % A phase restarts from the SEED, not from wherever the previous phase
    % ended: 'up' begins at the same 0.978 m/s gait 'down' did, and inheriting
    % the 0.35 m/s gait instead would make its first rung a 0.85 m/s jump.
    if st.reseed
        z = z0;  p = p0;
        ch3_logln(LOG, sprintf('=== phase %s: restarting from the seed gait ===', st.phase));
    end

    % Set the budget every rung rather than inheriting it: a resumed run picks
    % p up from the state file, which carries whatever the SEED was solved with
    % (400 here) -- far more than a scaled rung needs.
    if ~isempty(iters), p.max_iter = iters; else, p.max_iter = ITERS_DEFAULT; end
    pk = p;

    % Scaling is what makes these solves converge (see the header). The guard
    % is cheap insurance on top: it ends a runaway in minutes instead of at the
    % iteration cap. Note this OutputFcn REPLACES the checkpoint one
    % ch3_col_solve installs -- acceptable because the march already saves
    % state after every rung, so a crash costs one rung either way.
    J_start = ch3_col_cost(z, pk);
    opts = struct('ScaleProblem', true, ...
                  'OutputFcn', @(zz, ov, state) cost_guard(ov, state, J_start, J_GUARD));

    ch3_logln(LOG, sprintf('--- rung %d/%d [%s] target v = %.4f m/s', ...
                           k, numel(stages), st.phase, st.v));

    v_entry = measure(z, pk);
    v_from  = v_entry;
    v_tgt   = st.v;
    ok      = false;
    nbis    = 0;

    while true
        [z_try, h] = ch3_continuation(pk, v_tgt, pk.max_iter, z, LOG, opts);

        if ~isempty(h) && h(end).verify_ok
            z      = z_try;
            v_from = measure(z, pk);
            if abs(v_tgt - st.v) < 1e-6
                ok = true;
                break;                      % reached the rung's own target
            end
            % A BISECTED STEP LANDED, SO RE-AIM AT THE FULL TARGET FROM HERE.
            % Banking the shortfall instead would silently hand the next rung a
            % step larger than DV -- the sweep would quietly coarsen exactly
            % where it had just proved it needed to be finer.
            ch3_logln(LOG, sprintf('    intermediate v=%.4f verified; re-aiming at %.4f', ...
                                   v_from, st.v));
            v_tgt = st.v;
            continue;
        end

        if isempty(h), dev = NaN; else, dev = h(end).verify_dev; end
        nbis = nbis + 1;
        if nbis > MAX_BISECT
            ch3_logln(LOG, sprintf(['    rung missed at v=%.4f (dev %.3e) after ' ...
                                    '%d bisections'], v_tgt, dev, nbis - 1));
            break;
        end

        % Halve the step from the last ACHIEVED speed, and restart from the
        % PRE-rung gait -- never from the failed one, which is by definition
        % not a real trajectory.
        v_tgt = 0.5 * (v_from + v_tgt);
        ch3_logln(LOG, sprintf(['    missed (dev %.3e); bisecting to v=%.4f ' ...
                                '(attempt %d/%d)'], dev, v_tgt, nbis + 1, MAX_BISECT + 1));
    end

    % PARTIAL PROGRESS STILL COUNTS. A rung that got most of the way to its
    % target left the gait verified at a NEW speed nearer the goal, and the
    % next rung pulls from there; only a rung that moved essentially nowhere
    % means the sweep is genuinely stuck and the phase should stop.
    if ~ok
        moved = abs(v_from - v_entry);
        want  = abs(st.v   - v_entry);
        if want > 0 && moved / want >= 0.2
            ch3_logln(LOG, sprintf(['    PARTIAL: reached %.4f of %.4f (%.0f%% of the ' ...
                                    'rung); continuing from there'], ...
                                   v_from, st.v, 100 * moved / want));
        else
            ch3_logln(LOG, sprintf('STOP phase %s at rung %d: no verified step toward %.4f', ...
                                   st.phase, k, st.v));
            skip_phase = st.phase;
            k_done = k; %#ok<NASGU>
            save(STATE, 'k_done', 'z', 'p', 'z0', 'p0', 'lib', 'stages');
            continue;
        end
    end

    % Record the speed the gait ACTUALLY walks at, not the one that was asked
    % for: after a partial rung those differ, and a p whose v_des disagrees
    % with its own z is what makes a stored library member unreproducible.
    p       = pk;
    p.v_des = v_from;

    E   = ch3_col_eval(z, p);
    V   = ch3_col_verify(z, p, false);
    qt  = E.X(3,:);
    hip = -E.X(2,:);                       % pz is DOWN-positive
    pku = max(max(abs(E.u(:))), max(abs(E.um(:))));

    ch3_logln(LOG, sprintf(['    N=%d  T=%.4f  L=%.4f  v=%.4f m/s (target %.4f)\n' ...
                            '    qt  [%+.4f %+.4f] rad (%+.1f .. %+.1f deg)\n' ...
                            '    hip [ %.4f  %.4f] m  (bob %.4f)\n' ...
                            '    peak|u|=%.1f Nm  ||impulse||=%.3f Ns\n' ...
                            '    verify %.3e (ok=%d)   %.0f s'], ...
                           size(E.X,2), E.T, E.L_step, E.L_step/E.T, st.v, ...
                           min(qt), max(qt), rad2deg(min(qt)), rad2deg(max(qt)), ...
                           min(hip), max(hip), max(hip)-min(hip), ...
                           pku, norm(E.impulse), V.max_dev, V.ok, toc(t0)));

    lib{end+1} = snap(z, p, E.L_step/E.T, st.phase); %#ok<AGROW>

    k_done = k; %#ok<NASGU>
    save(STATE, 'k_done', 'z', 'p', 'z0', 'p0', 'lib', 'stages');
end

% --- assemble ------------------------------------------------------------
speeds = cellfun(@(s) s.speed, lib);
[speeds, ord] = sort(speeds);
lib = lib(ord);
library = lib; %#ok<NASGU>
save(LIB, 'library');
ch3_logln(LOG, sprintf('library: %d gaits, %.4f .. %.4f m/s -> %s', ...
                       numel(lib), min(speeds), max(speeds), LIB));

% The fastest VERIFIED member is the deliverable: the 1.2 m/s gait when the up
% phase completed, and the fastest rung that did otherwise -- say which, rather
% than silently shipping a slower gait under the v120 name.
[v_best, i_best] = max(speeds);
z = lib{i_best}.z;  p = lib{i_best}.p;
if abs(v_best - V_HI) > 1e-3
    ch3_logln(LOG, sprintf(['NOTE: sweep topped out at %.4f m/s, short of the ' ...
                            '%.2f m/s target; saving that gait as the fast end.'], ...
                           v_best, V_HI));
end

R = ch3_report(z, p, struct('stability', true, 'simulate', 5)); %#ok<NASGU>
save(FINAL, 'z', 'p', 'R');
[X, ~, alpha] = ch3_col_unpack(z, p);
try
    ch3_animate(X(:,1), alpha, p, 4, GIF);
catch ME
    ch3_logln(LOG, sprintf('animate failed: %s', ME.message));
end
ch3_logln(LOG, sprintf('FINAL v=%.4f m/s  rho=%.4f  -> %s', v_best, R.rho, FINAL));
ch3_logln(LOG, 'MARKER_ALLDONE');

end

% ---------------------------------------------------------------- helpers
function stages = build_stages(v0, v_lo, v_hi, dv)
% Two phases out from the seed speed. Rungs step by dv and the LAST of each
% phase absorbs the remainder, so no rung is ever larger than dv.
down = v0 - dv : -dv : v_lo;
if isempty(down) || abs(down(end) - v_lo) > 1e-9, down = [down, v_lo]; end
up = v0 + dv : dv : v_hi;
if isempty(up) || abs(up(end) - v_hi) > 1e-9, up = [up, v_hi]; end

stages = struct('v', {}, 'phase', {}, 'reseed', {});
for i = 1:numel(down)
    stages(end+1) = struct('v', down(i), 'phase', 'down', 'reseed', i == 1); %#ok<AGROW>
end
for i = 1:numel(up)
    stages(end+1) = struct('v', up(i), 'phase', 'up', 'reseed', i == 1); %#ok<AGROW>
end
end

function v = measure(z, p)
E = ch3_col_eval(z, p);
v = E.L_step / E.T;
end

function stop = cost_guard(ov, state, J_start, factor)
% Stop a rung whose cost has run away. A diverging solve here climbs J by
% orders of magnitude while feasibility crawls, and fmincon reports the same
% exitflag 0 for that as for "still converging, out of iterations" -- so
% without this the only way to tell them apart is to wait for the cap.
stop = false;
if ~strcmp(state, 'iter'), return; end
if J_start > 0 && ov.fval > factor * J_start
    stop = true;
end
end

function s = snap(z, p, v, phase)
E = ch3_col_eval(z, p);
V = ch3_col_verify(z, p, false);
s = struct('z', z, 'p', p, 'speed', v, 'phase', phase, ...
           'T', E.T, 'L_step', E.L_step, 'N', size(E.X, 2), ...
           'qt_lo', min(E.X(3,:)), 'qt_hi', max(E.X(3,:)), ...
           'hip_lo', min(-E.X(2,:)), 'hip_hi', max(-E.X(2,:)), ...
           'peak_u', max(max(abs(E.u(:))), max(abs(E.um(:)))), ...
           'impulse', norm(E.impulse), ...
           'verify_dev', V.max_dev, 'verify_ok', V.ok);
end
