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
SEED  = fullfile(resd, 'ch3_gait_lean_tall_constrained.mat');
FINAL = fullfile(resd, 'ch3_speed_family.mat');

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
    p.enforce_nec1 = false;         % the lean leg lets the speed float

    % HZD OFF AS A CONSTRAINT, MEASURED AS A DIAGNOSTIC INSTEAD.
    % NEC4/NEC5 come from ch3_zero_dynamics, which ch3_col_constraints itself
    % documents as the expensive row: one evaluation is p.hzd_grid_solve (41)
    % points of roughly four 7x7 solves, and fmincon finite-differences it once
    % per decision variable. Measured here it costs 2 min/iteration against
    % 18 s with the gate off -- a 7x tax that puts this campaign past 30 hours.
    % That cost is why the gate ships defaulting to false.
    %
    % So it is switched off for the marching and evaluated ONCE PER RUNG below,
    % which is free by comparison and is the same "report it rather than impose
    % it" treatment ch3_col_eval gives hybrid invariance. If a rung ever comes
    % back with nec4/nec5 violated the log says so immediately, rather than the
    % campaign discovering it at the end.
    p.limits.enable.hzd = false;

    % RAISE THE HIP CEILING BEFORE MARCHING ANYTHING.
    % The seed rides its band ceiling exactly: band [0.8450 0.9250], gait hip
    % max 0.9250. Together with the four limit caps -- Fz 50.0, |Fx|/Fz 0.4000,
    % |u| 120.0, ||I|| 15.00, every one of them exactly active -- that is FIVE
    % active constraints, which pins the gait at a vertex of the feasible set
    % with no direction left to move in. It showed up as a 0.05 m/s speed rung
    % collapsing: objective and feasibility frozen at 3.286e-02 while the step
    % length fell to 3.1e-07.
    %
    % Four of those five are physics or Table 3.1 and stay. The hip ceiling is
    % the one that is pure design choice, and raising it is MORE aligned with
    % wanting a tall gait, not less -- 0.955 is the ceiling ch3_lean_tall_march
    % already targets. The floor is left where it is; the final height rung
    % raises that instead, once the marching is done.
    p.limits.hip_h     = 0.9000;
    p.limits.hip_h_tol = 0.0550;      % band [0.8450 0.9550]
    k0 = 1;  hist = {};  dead = '';  z_split = [];

    E0 = ch3_col_eval(z, p);
    V0 = ch3_col_verify(z, p, false);
    stages = ladder(E0.L_step / E0.T);   % ladder is anchored to the seed speed
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
    ch3_logln(LOG, sprintf(['    hip ceiling raised 0.9250 -> 0.9550 (it was a 5th ' ...
                            'active constraint, pinning the gait)']));
    ch3_logln(LOG, sprintf(['    plan: lean is DONE (qt floor +0.080); march speed UP ' ...
                            'from %.4f to 1.20 m/s until it breaks'], E0.L_step/E0.T));
end

% --- the campaign ---------------------------------------------------------
% A while loop, not for: phase 3 is APPENDED to stages partway through, and a
% for loop fixes its bound at entry and would never see the new rungs.
k = k0;
stopped = false;
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
            p.enforce_nec1 = true;
            p.v_des = st.target;
            p = ch3_col_budget(p, iters, z);

            % ch3_col_solve ALREADY ran ch3_col_verify and put it in out.verify.
            % Re-running it here cost a second ode45 at RelTol 1e-11 -- about
            % 13 minutes a rung, comparable to the solve itself -- for an answer
            % already in hand.
            [z_new, out] = ch3_col_solve(p, z);
            V = out.verify;

            % ONE RETRY on a bigger budget, on EITHER failure mode. At N = 61
            % the mesh is not the suspect -- a miss is exhausted iterations,
            % the same reasoning ch3_realizability_march uses. Restart from the
            % PRE-rung gait.
            %
            % "Verified but did not move" is its own failure and needs the
            % retry just as much: fmincon can return exitflag -2 after two
            % iterations having changed nothing, and the unchanged gait then
            % verifies perfectly. Gating the retry on verification alone let
            % that case through untried.
            if ~V.ok || ~target_met(st, ch3_col_eval(z_new, p), p)
                ch3_logln(LOG, sprintf(['    rung missed (verify %.3e); retrying ' ...
                                        'with %d iterations'], V.max_dev, 2*iters));
                p2 = ch3_col_budget(p, 2*iters, z);
                [z_new, out] = ch3_col_solve(p2, z);
                V = out.verify;
                if V.ok, p = p2; end
            end
            ok  = V.ok;
            dev = V.max_dev;
    end

    % DID THE RUNG ACTUALLY DO ANYTHING?  Verification asks whether the nodes
    % lie on a real trajectory. It does NOT ask whether this rung reached the
    % target it was set. When fmincon bails immediately -- exitflag -2, z
    % returned essentially unchanged -- the UNCHANGED gait still verifies,
    % because it is the previous rung's already-verified result. The rung then
    % records as a success having achieved nothing.
    %
    % That is not hypothetical: five up-leg rungs asking for 0.65-0.85 m/s all
    % "succeeded" in 15 seconds each while sitting at 1.0616 m/s, and were
    % written into the family labelled with speeds they never reached. The
    % family is the deliverable, so a mislabelled member is worse than a
    % missing one.
    Etmp = ch3_col_eval(z_new, p);
    if ok && ~target_met(st, Etmp, p)
        ok = false;
        ch3_logln(LOG, '    rung did NOT reach its target (solver made no useful progress)');
    end

    % --- measure whatever the rung produced --------------------------------
    E   = ch3_col_eval(z_new, p);
    qt  = E.X(3,:);
    hip = -E.X(2,:);
    lam = [E.lam, E.lamm];
    peak_u = max(max(abs(E.u(:))), max(abs(E.um(:))));

    % NEC4/NEC5 once per rung -- cheap here, ruinous inside a gradient.
    nec = [NaN NaN];
    try
        Zd = ch3_zero_dynamics(E.alpha, p, p.hzd_grid_solve);
        nec = [Zd.nec4, Zd.nec5];
    catch
    end

    ch3_logln(LOG, sprintf([ ...
        '    v=%.4f  T=%.4f  L=%.4f  J=%.2f  exitflag=%d\n' ...
        '    qt [%+.4f %+.4f] rad   hip [%.4f %.4f] m (bob %.4f)\n' ...
        '    Fz>=%.1f (%.1f)  |Fx|/Fz=%.4f (%.2f)  |u|=%.1f (%.0f)  ||I||=%.2f (%.1f)\n' ...
        '    max|ceq|=%.2e  max c=%.2e  verify %.3e (ok=%d)\n' ...
        '    NEC4=%+.4f  NEC5=%+.4f  (both must stay <= 0)   %.0f s'], ...
        E.L_step/E.T, E.T, E.L_step, out.fval, out.exitflag, ...
        min(qt), max(qt), min(hip), max(hip), max(hip)-min(hip), ...
        min(lam(2,:)), p.limits.Fz_min, ...
        max(abs(lam(1,:))./max(lam(2,:),1e-9)), p.limits.mu_s, ...
        peak_u, p.limits.u_max, norm(E.impulse), p.limits.impulse_max, ...
        out.max_ceq, out.max_c, dev, ok, nec(1), nec(2), toc(t0)));

    if ~ok
        % THE CHAIN HALTS HERE, whichever leg it was. Every leg warm-starts
        % from the one before -- the lean rungs expect the 0.35 m/s gait, the
        % up rungs expect the leaned one -- so skipping a failed rung would run
        % the next leg from a gait it was not designed for, and a failed rung's
        % z is not a real trajectory in the first place.
        %
        % A stall on the UP leg is the expected outcome rather than a fault: it
        % locates the fastest gait reachable at this posture inside Table 3.1.
        % Either way the family collected so far is written out below.
        if strcmp(st.branch, 'up')
            ch3_logln(LOG, sprintf(['    CEILING FOUND: v = %.3f m/s does not ' ...
                                    'verify (%.3e). The gait below it is the ' ...
                                    'fastest reachable at this posture.'], ...
                                   st.target, dev));
        else
            ch3_logln(LOG, sprintf('    STOP in %s (%s): did not verify (%.3e)', ...
                                   st.kind, st.desc, dev));
        end
        k_done = k - 1; %#ok<NASGU>
        save(STATE, 'k_done', 'z', 'p', 'hist', 'stages', 'z_split', 'dead');
        stopped = true;
        break;
    end

    if strcmp(st.kind, 'speed')
        hist{end+1} = rec(E, p, out, dev, st.target, st.branch, qt, hip, lam, peak_u, z_new); %#ok<AGROW>
    end

    z = z_new;
    k_done = k; %#ok<NASGU>

    save(STATE, 'k_done', 'z', 'p', 'hist', 'stages', 'z_split', 'dead');
    k = k + 1;
end

% --- collect the family, ascending in speed -------------------------------
if isempty(hist)
    ch3_logln(LOG, 'no rung converged; nothing to write');
    ch3_logln(LOG, marker(stopped));
    return;
end

fam = [hist{:}];
[~, ord] = sort([fam.v]);
fam = fam(ord);
save(FINAL, 'fam', '-v7.3');

ch3_logln(LOG, sprintf('FAMILY: %d gaits, v = %.4f .. %.4f m/s -> %s', ...
                       numel(fam), min([fam.v]), max([fam.v]), FINAL));
ch3_logln(LOG, sprintf('    lean at the fast end: qt [%+.4f %+.4f] rad', ...
                       fam(end).qt_lo, fam(end).qt_hi));
ch3_logln(LOG, marker(stopped));

end

% ---------------------------------------------------------------- helpers
function tf = target_met(st, E, p)
%TARGET_MET  Did this rung achieve what it was set, not merely stay valid?
tol = 2e-3;
switch st.kind
    case 'speed'
        tf = abs(E.L_step/E.T - st.target) <= tol;
    case 'posture'
        tf = true;
        if isfield(st.opt,'qt_range') && ~isempty(st.opt.qt_range)
            tf = tf && min(E.X(3,:)) >= st.opt.qt_range(1) - tol;
        end
        if isfield(st.opt,'hip_h') && ~isempty(st.opt.hip_h)
            tf = tf && min(-E.X(2,:)) >= (p.limits.hip_h - p.limits.hip_h_tol) - tol;
        end
    otherwise
        tf = true;
end
end

function m = marker(stopped)
% The wrapper greps for these: ALLDONE ends the run, STOPPED ends it too but
% says the campaign halted itself rather than finishing the ladder.
if stopped, m = 'MARKER_STOPPED'; else, m = 'MARKER_ALLDONE'; end
end

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

function stages = ladder(v0)
% ONE CONTINUOUS CHAIN, ordered so each leg is solved where it is cheapest.
%
% The first attempt marched the lean at ~1.0 m/s and hit a wall: rung 1
% converged (exitflag 1), rung 2 only just made it on iterations (exitflag 0,
% 1806 s), rung 3 returned exitflag -2 after 12 of its 100 iterations. That
% gradient is a feasibility boundary, not a budget shortfall -- the gait was
% pinned at all four limit caps AND riding the top of the hip band, with
% nothing left to trade for torso pitch.
%
% STEP SIZE IS 0.05 m/s, NOT 0.10. A 0.10 step was tried and stalled: the
% solve ended at max|ceq| = 8.6e-02 with the retry making no progress --
% feasibility stuck at 9.05e-02, step length 2.3e-03, and first-order
% optimality GROWING. With all four limits pinned at their caps the gait has no
% slack, so even a speed change is expensive and the step has to be small.
%
% Speed is the slack. Torque and impact impulse both fall with speed, so the
% same lean that is infeasible at 1.0 m/s has room at 0.35. Hence: walk the
% speed DOWN first at the partial lean, do the lean where it is cheap, then
% walk the speed back UP at the full lean until it breaks. Where it breaks is
% the answer to "how fast can this robot walk leaning forward, inside
% Table 3.1" -- a measurement, not a failure.
%
% Every leg continues from the one before it, so there is no restart and no
% branch: it is a single warm-start chain from 1.00 down to 0.35, through the
% lean, and back up.
stages = struct('kind', {}, 'desc', {}, 'opt', {}, 'target', {}, ...
                'branch', {}, 'restart', {});

% THE LEAN IS DONE. Results/ch3_gait_lean_tall_constrained.mat holds it:
% qt [+0.0800 +0.1561] rad, hip [0.9282 0.9550] m, v = 1.0616 m/s, with
% Fz 50.0, |Fx|/Fz 0.4000, |u| 85.6 and ||I|| 15.00 all inside their caps and
% NEC4/NEC5 at -1.19/-0.49. Only the speed ceiling is left to find.
%
% THE GRID IS RELATIVE TO THE SEED, NOT ABSOLUTE. The previous ladder hardcoded
% an up leg starting at 0.65 m/s on the assumption the lean would finish near
% 0.600. It finished at 1.0616 -- speed floats during posture rungs, and the
% forward lean made the gait want to go faster -- so the first up rung demanded
% a 0.41 m/s DECREASE and fmincon gave up in 15 seconds. Anchoring the ladder
% to the seed's measured speed removes that whole class of failure.
% 0.02 STEPS, NOT 0.05. The first attempt jumped 1.0616 -> 1.100 and fmincon
% called the QP infeasible in 45 seconds. A 0.038 m/s step failing that fast is
% suggestive of a wall but not proof of one, and near a ceiling the honest
% thing is to walk up to it in small steps and let it stop the march itself.
for v = (ceil(v0/0.02)*0.02) : 0.02 : 1.20
    stages(end+1) = mks(v, 'up'); %#ok<AGROW>
end
end

function s = mks(target, branch)
s = struct('kind', 'speed', ...
           'desc', sprintf('v_des = %.3f m/s [%s]', target, branch), ...
           'opt', struct(), 'target', target, 'branch', branch, 'restart', false);
end
