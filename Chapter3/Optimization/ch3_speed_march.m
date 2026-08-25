function ch3_speed_march(iters, max_stage, iters_posture)
%CH3_SPEED_MARCH  Forward-lean, tall, fully limited, swept 0.35 -> 1.20 m/s.
%
%   ch3_speed_march              full campaign, 40 iters/speed rung, 80/posture
%   ch3_speed_march(60)          bigger speed budget
%   ch3_speed_march(40, [], 120) bigger POSTURE budget, speed budget unchanged
%   ch3_speed_march(2, 1)        smoke test: 2 iterations, one rung
%
% Builds a FAMILY of gaits across the speed range that hold the forward lean
% and the raised hip AND stay inside all four Table 3.1 limits -- GRF,
% friction, torque and impulse.
%
% WHY THE SEED IS ch3_gait_full_constrained AND NOT ch3_gait_forward_lean_tall.
% The obvious seed is the one already named for this posture. It cannot be
% used: it has min Fz = -192.6 N. A negative normal force is the ground PULLING
% the foot down, so that gait only exists with the stance foot glued to the
% floor, and its |Fx|/Fz reads 1.4e10 because the ratio is dividing through
% zero as Fz changes sign. Enabling friction there is exactly the
% division-by-near-zero ch3_col_solve warns about, and every torque and impulse
% number measured on it describes a trajectory that is not realizable anyway.
%
% ch3_gait_full_constrained is the opposite: all twelve gates on, all four
% limits met exactly (Fz 50.0 N, |Fx|/Fz 0.400, peak |u| 120.0 Nm,
% ||I|| 15.00 Ns), verified, N = 61, walking at 0.9747 m/s with the hip at
% 0.9125-0.9250 m. Everything the request asks for except the lean: it sits at
% qt in [-0.100 -0.019], leaning slightly BACKWARD. So the campaign starts from
% a valid gait and moves the one property that is wrong, rather than starting
% from the right posture and trying to repair the physics.
%
% THREE PHASES, ONE REQUIREMENT AT A TIME.
%
%   Phase 1 (lean).  The optimizer RIDES THE LOWER BOUND of qt_range -- that is
%   what the achieved -0.100 against a box of [-0.100 +0.450] means -- so the
%   lean is moved by walking that floor forward, -0.05, 0.00, +0.04, +0.08,
%   with the ceiling closing in behind it. ch3_posture_march does the rungs; it
%   reposes the torso into each new box first, because letting fmincon clamp
%   the old solution into it shifts theta at every node and detonates the
%   defects.
%
%   Phase 2 (height).  The gait rides the TOP of its band ([0.845 0.925],
%   achieved 0.9125-0.9250), so "tall" is currently a preference, not a
%   constraint -- nothing stops the hip sagging to 0.845 once the torso pitches
%   forward. One rung raises the floor to 0.885 and locks the height in before
%   speed starts moving.
%
%   Phase 3 (speed).  Only now does NEC1 come on and the sweep run, from
%   whatever speed phases 1-2 settled at.
%
% NO LIMIT IS EVER TOUCHED.  All four are on in the seed and stay on, at their
% final values, for the whole campaign. There is no limit ladder here and no
% ordering question: the GRF-before-friction rule matters when a gate is being
% SWITCHED ON against a violated constraint, and nothing is switched on.
%
% THE SEED HAS NO SLACK, WHICH IS THE RISK.  It meets all four limits EXACTLY,
% so every rung of phase 1 asks for a posture change with nothing to trade. If
% the lean stalls, that is the honest answer -- the forward lean and the full
% Table 3.1 set are not simultaneously reachable at this speed -- and the log
% says which rung and which constraint.
%
% EXPECT THE FAST END TO BE THE HARD ONE.  Impact impulse grows with speed and
% ||I|| is already at its 15 Ns cap at 0.97 m/s, so the up branch is likely to
% stall well before 1.20. That stall locates the fastest gait reachable at this
% posture inside Table 3.1, which is a result rather than a failure, and it
% cannot contaminate the down branch -- that restarts from the phase-2 gait.
%
% CHECKPOINTS EVERY RUNG.  This is ~22 solves of several minutes each and this
% MATLAB install aborts from its add-on registry often enough that an
% unattended multi-hour run will probably not survive intact (see ch3_col_resume
% for the stack). State is written after every rung, so re-running resumes.
%
% Inputs
%   iters     : fmincon iterations per rung (default 40)
%   max_stage : stop after this many rungs (smoke testing)
%
% Outputs (written to Results/)
%   ch3_speed_march_state.mat  k_done, z, p, hist, stages -- the resume point
%   ch3_speed_family.mat       fam: one entry per converged speed rung,
%                              ascending in speed, each with .v .z .p and its
%                              measurements
%   ch3_speed_march.log        per-rung progress, timestamped and flushed
%
% See also CH3_POSTURE_MARCH, CH3_CONTINUATION, CH3_REALIZABILITY_MARCH,
%          CH3_COL_SOLVE, CH3_COL_BUDGET, CH3_LOGLN.

if nargin < 1 || isempty(iters), iters = 40; end
if nargin < 2, max_stage = []; end
% The lean rungs are the hard phase -- the seed meets all four limits exactly,
% so each one asks for a posture change with nothing to trade. The first
% attempt at this campaign died on lean rung 1 with exitflag -2 at
% max|ceq| = 1.0e-02, which is SQP stalling a hair outside feasibility rather
% than a genuinely infeasible problem. Smaller steps (below) and a bigger
% budget (here) are the documented cure.
if nargin < 3 || isempty(iters_posture), iters_posture = 2*iters; end

root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
resd = fullfile(root, 'Results');
if ~exist(resd, 'dir'), mkdir(resd); end

STATE = fullfile(resd, 'ch3_speed_march_state.mat');
LOG   = fullfile(resd, 'ch3_speed_march.log');
SEED  = fullfile(resd, 'ch3_gait_full_constrained.mat');
FINAL = fullfile(resd, 'ch3_speed_family.mat');

V_GRID = 0.35:0.05:1.20;      % the requested range, absolute

% --- start or resume ------------------------------------------------------
if exist(STATE, 'file')
    S = load(STATE);
    k0 = S.k_done + 1;  z = S.z;  p = S.p;  hist = S.hist;
    stages = S.stages;  z_split = S.z_split;
    if isfield(S,'dead'), dead = S.dead; else, dead = ''; end
    ch3_logln(LOG, sprintf('RESUME at rung %d/%d', k0, numel(stages)));
else
    if ~exist(SEED, 'file')
        error('ch3_speed_march:seed', 'Cannot find the seed gait "%s".', SEED);
    end
    S = load(SEED);
    z = S.z;
    if isfield(S, 'p'), p = S.p; else, p = S.pf; end
    p = ch3_upgrade_params(p);
    p.enforce_nec1 = false;         % phases 1-2 let the speed float
    k0 = 1;  hist = {};  dead = '';  z_split = [];
    stages = posture_ladder();      % phase 3 is appended once posture lands

    E0 = ch3_col_eval(z, p);
    V0 = ch3_col_verify(z, p, false);
    lam0 = [E0.lam, E0.lamm];
    ch3_logln(LOG, sprintf(['=== SPEED march: %d iters/speed rung, %d/posture rung, ' ...
                            'all four limits ON ==='], iters, iters_posture));
    ch3_logln(LOG, sprintf(['    seed: N=%d v=%.4f T=%.4f qt=[%+.4f %+.4f] ' ...
                            'hip=[%.4f %.4f] verify %.3e (ok=%d)'], ...
                           E0.N, E0.L_step/E0.T, E0.T, ...
                           min(E0.X(3,:)), max(E0.X(3,:)), ...
                           min(-E0.X(2,:)), max(-E0.X(2,:)), V0.max_dev, V0.ok));
    ch3_logln(LOG, sprintf(['    seed limits: Fz>=%.1f (cap %.1f)  |Fx|/Fz=%.4f (cap %.2f)  ' ...
                            'peak|u|=%.1f (cap %.0f)  ||I||=%.2f (cap %.1f)'], ...
                           min(lam0(2,:)), p.limits.Fz_min, ...
                           max(abs(lam0(1,:))./max(lam0(2,:),1e-9)), p.limits.mu_s, ...
                           max(max(abs(E0.u(:))),max(abs(E0.um(:)))), p.limits.u_max, ...
                           norm(E0.impulse), p.limits.impulse_max));
    ch3_logln(LOG, sprintf('    target: qt box -> [0.080 0.250], hip floor -> 0.885'));
end

% --- the campaign ---------------------------------------------------------
% A while loop, not for: phase 3 is APPENDED to stages partway through, and a
% for loop fixes its bound at entry and would never see the new rungs.
k = k0;
while k <= numel(stages)
    if ~isempty(max_stage) && k > max_stage, break; end
    st = stages(k);

    if ~isempty(st.branch) && strcmp(st.branch, dead)
        k = k + 1;  continue;
    end

    t0 = tic;
    ch3_logln(LOG, sprintf('--- rung %d/%d [%s] %s', k, numel(stages), st.kind, st.desc));

    ok = true;
    % Assigned by every branch below, but pre-seeded so the measurement block
    % below can report a rung that FAILED before producing solver output.
    out = struct('fval', NaN, 'exitflag', NaN, 'max_ceq', NaN, 'max_c', NaN);
    dev = NaN;
    switch st.kind
        case 'posture'
            % Delegate to ch3_posture_march: it reposes the torso into the new
            % box before solving, which is what keeps the phase clock and both
            % feet intact across a pitch change.
            p = ch3_col_budget(p, iters_posture, z);
            [z_new, h] = ch3_posture_march(p, st.opt, iters_posture, z, LOG);
            if isempty(h) || ~h(end).verify_ok
                ok = false;
                if ~isempty(h)
                    dev = h(end).verify_dev;
                    out = struct('fval', h(end).fval, 'exitflag', h(end).exitflag, ...
                                 'max_ceq', h(end).max_ceq, 'max_c', h(end).max_c);
                end
            else
                p.qt_range             = h(end).qt_range;
                p.limits.hip_h         = h(end).hip_h;
                p.limits.hip_h_tol     = h(end).hip_h_tol;
                p.limits.enable.height = h(end).height;
                dev = h(end).verify_dev;
                out = struct('fval', h(end).fval, 'exitflag', h(end).exitflag, ...
                             'max_ceq', h(end).max_ceq, 'max_c', h(end).max_c);
            end

        case 'speed'
            if st.restart
                z = z_split;
                ch3_logln(LOG, sprintf('--- restarting from the phase-2 gait for the %s branch', st.branch));
            end
            p.enforce_nec1 = true;
            p.v_des = st.target;
            p = ch3_col_budget(p, iters, z);

            [z_new, out] = ch3_col_solve(p, z);
            V = ch3_col_verify(z_new, p, false);

            % ONE RETRY on a bigger budget. At N = 61 the mesh is not the
            % suspect -- a miss is exhausted iterations, the same reasoning
            % ch3_realizability_march uses. Restart from the PRE-rung gait.
            if ~V.ok
                ch3_logln(LOG, sprintf(['    rung missed (verify %.3e); retrying ' ...
                                        'with %d iterations'], V.max_dev, 2*iters));
                p2 = ch3_col_budget(p, 2*iters, z);
                [z_new, out] = ch3_col_solve(p2, z);
                V = ch3_col_verify(z_new, p2, false);
                if V.ok, p = p2; end
            end
            ok  = V.ok;
            dev = V.max_dev;
    end

    % --- measure whatever the rung produced --------------------------------
    E   = ch3_col_eval(z_new, p);
    qt  = E.X(3,:);
    hip = -E.X(2,:);
    lam = [E.lam, E.lamm];
    peak_u = max(max(abs(E.u(:))), max(abs(E.um(:))));

    ch3_logln(LOG, sprintf([ ...
        '    v=%.4f  T=%.4f  L=%.4f  J=%.2f  exitflag=%d\n' ...
        '    qt [%+.4f %+.4f] rad   hip [%.4f %.4f] m (bob %.4f)\n' ...
        '    Fz>=%.1f (%.1f)  |Fx|/Fz=%.4f (%.2f)  |u|=%.1f (%.0f)  ||I||=%.2f (%.1f)\n' ...
        '    max|ceq|=%.2e  max c=%.2e  verify %.3e (ok=%d)   %.0f s'], ...
        E.L_step/E.T, E.T, E.L_step, out.fval, out.exitflag, ...
        min(qt), max(qt), min(hip), max(hip), max(hip)-min(hip), ...
        min(lam(2,:)), p.limits.Fz_min, ...
        max(abs(lam(1,:))./max(lam(2,:),1e-9)), p.limits.mu_s, ...
        peak_u, p.limits.u_max, norm(E.impulse), p.limits.impulse_max, ...
        out.max_ceq, out.max_c, dev, ok, toc(t0)));

    if ~ok
        if strcmp(st.kind, 'speed')
            dead = st.branch;
            ch3_logln(LOG, sprintf(['    STOP on the %s branch at v=%.3f: did not ' ...
                                    'verify (%.3e). Skipping the rest of it.'], ...
                                   st.branch, st.target, dev));
            k_done = k; %#ok<NASGU>
            save(STATE, 'k_done', 'z', 'p', 'hist', 'stages', 'z_split', 'dead');
            k = k + 1;  continue;
        end
        ch3_logln(LOG, sprintf('    STOP in %s: did not verify (%.3e)', st.kind, dev));
        ch3_logln(LOG, 'MARKER_STOPPED');
        k_done = k - 1; %#ok<NASGU>
        save(STATE, 'k_done', 'z', 'p', 'hist', 'stages', 'z_split', 'dead');
        return;
    end

    if strcmp(st.kind, 'speed')
        hist{end+1} = rec(E, p, out, dev, st.target, st.branch, qt, hip, lam, peak_u, z_new); %#ok<AGROW>
    end

    z = z_new;
    k_done = k; %#ok<NASGU>

    % Posture phases just finished: split the speed grid where they landed and
    % append phase 3. z_split is the warm start BOTH branches start from.
    if k == numel(stages) && isempty(z_split)
        v0 = E.L_step / E.T;
        z_split = z;
        hist{end+1} = rec(E, p, out, dev, v0, 'posture', qt, hip, lam, peak_u, z); %#ok<AGROW>
        stages = [stages, speed_ladder(V_GRID, v0)]; %#ok<AGROW>
        ch3_logln(LOG, sprintf(['=== posture done at v=%.4f, qt [%+.4f %+.4f], ' ...
                                'hip [%.4f %.4f]. Phase 3: %d speed rungs ==='], ...
                               v0, min(qt), max(qt), min(hip), max(hip), ...
                               numel(stages)-k));
    end

    save(STATE, 'k_done', 'z', 'p', 'hist', 'stages', 'z_split', 'dead');
    k = k + 1;
end

% --- collect the family, ascending in speed -------------------------------
if isempty(hist)
    ch3_logln(LOG, 'no rung converged; nothing to write');
    ch3_logln(LOG, 'MARKER_ALLDONE');
    return;
end

fam = [hist{:}];
[~, ord] = sort([fam.v]);
fam = fam(ord);
save(FINAL, 'fam', '-v7.3');

ch3_logln(LOG, sprintf('FAMILY: %d gaits, v = %.4f .. %.4f m/s -> %s', ...
                       numel(fam), min([fam.v]), max([fam.v]), FINAL));
ch3_logln(LOG, 'MARKER_ALLDONE');

end

% ---------------------------------------------------------------- helpers
function r = rec(E, p, out, dev, target, branch, qt, hip, lam, peak_u, z)
r = struct('v', E.L_step/E.T, 'v_target', target, 'branch', branch, ...
           'z', z, 'p', p, 'N', E.N, 'T', E.T, 'L_step', E.L_step, ...
           'fval', out.fval, 'exitflag', out.exitflag, ...
           'max_ceq', out.max_ceq, 'max_c', out.max_c, ...
           'qt_lo', min(qt), 'qt_hi', max(qt), ...
           'hip_lo', min(hip), 'hip_hi', max(hip), ...
           'Fz_min', min(lam(2,:)), ...
           'fric', max(abs(lam(1,:))./max(lam(2,:),1e-9)), ...
           'peak_u', peak_u, 'impulse', norm(E.impulse), 'verify_dev', dev);
end

function stages = posture_ladder()
% Phase 1 walks the qt box FLOOR forward from -0.100 to +0.080 -- the floor is
% what the optimizer rides -- closing the ceiling in behind it. Phase 2 then
% raises the hip band floor so the height cannot sag once the torso is forward.
stages = [ ...
    mkp('lean: qt box [-0.080 0.450]', struct('qt_range', [-0.080 0.450])), ...
    mkp('lean: qt box [-0.060 0.450]', struct('qt_range', [-0.060 0.450])), ...
    mkp('lean: qt box [-0.040 0.450]', struct('qt_range', [-0.040 0.450])), ...
    mkp('lean: qt box [-0.020 0.420]', struct('qt_range', [-0.020 0.420])), ...
    mkp('lean: qt box [ 0.000 0.400]', struct('qt_range', [ 0.000 0.400])), ...
    mkp('lean: qt box [ 0.020 0.370]', struct('qt_range', [ 0.020 0.370])), ...
    mkp('lean: qt box [ 0.040 0.320]', struct('qt_range', [ 0.040 0.320])), ...
    mkp('lean: qt box [ 0.060 0.280]', struct('qt_range', [ 0.060 0.280])), ...
    mkp('lean: qt box [ 0.080 0.250]  TARGET', struct('qt_range', [0.080 0.250])), ...
    mkp('height: hip band [0.885 0.925]', ...
        struct('hip_h', 0.9050, 'hip_h_tol', 0.0200, 'height', true)) ];
end

function stages = speed_ladder(grid, v0)
% Split the absolute speed grid at wherever the posture phases landed:
% everything below marched downward, everything above upward, each branch
% restarting from the posture gait.
down = sort(grid(grid < v0 - 1e-9), 'descend');
up   = sort(grid(grid > v0 + 1e-9), 'ascend');

stages = struct('kind', {}, 'desc', {}, 'opt', {}, 'target', {}, ...
                'branch', {}, 'restart', {});
for i = 1:numel(down)
    stages(end+1) = mks(down(i), 'down', i == 1); %#ok<AGROW>
end
for i = 1:numel(up)
    stages(end+1) = mks(up(i), 'up', i == 1); %#ok<AGROW>
end
end

function s = mkp(desc, opt)
s = struct('kind', 'posture', 'desc', desc, 'opt', opt, 'target', [], ...
           'branch', '', 'restart', false);
end

function s = mks(target, branch, restart)
s = struct('kind', 'speed', ...
           'desc', sprintf('v_des = %.3f m/s [%s]', target, branch), ...
           'opt', struct(), 'target', target, 'branch', branch, 'restart', restart);
end
