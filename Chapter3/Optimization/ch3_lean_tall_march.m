function ch3_lean_tall_march(route, iters, max_stage)
%CH3_LEAN_TALL_MARCH  March to the forward-lean gait with the hip at ~0.94 m.
%
%   ch3_lean_tall_march                  warm route, per-stage iteration budgets
%   ch3_lean_tall_march('cold')          cold route, all the way from a seed
%   ch3_lean_tall_march('warm', 2, 1)    smoke test: 2 iterations, one rung
%
% Two routes to the same gait. They share every stage mechanic below and differ
% only in where they start and how long the ladder is:
%
%   'warm' (default)  Starts from Results/ch3_gait_forward_lean_tall_nec3.mat,
%                     an N = 61 gait ALREADY at the target pitch (qt +0.080 ..
%                     +0.163 rad) with the hip at 0.909-0.925 m, and walks only
%                     the height band up: ceiling 0.925 -> 0.955, four rungs.
%
%                     THIS SEED CARRIES NEC3 ENFORCED (enable.impact, with
%                     mu_s_impact = 0.4) and sits at |Ix|/Iz = 0.3667, inside
%                     the impulse friction cone by 0.033. The gate is stored
%                     true, so ch3_upgrade_params keeps it on and every rung
%                     below must HOLD the cone rather than merely inherit it.
%                     That is the point of seeding from this file: its
%                     predecessor, ch3_gait_forward_lean_tall.mat, sat at
%                     0.4628 -- OUTSIDE a 0.4 cone -- so a warm run that ever
%                     enabled impact began infeasible.
%
%                     THE PRICE IS TORQUE. Marching NEC3 in cost peak |u|
%                     107.4 -> 140.4 Nm (+31%), which is ABOVE the 120 Nm
%                     ch3_realizability_march marches down to. Run that march
%                     AFTER this one, and expect its torque rungs to start
%                     further out than its own header's 191.4 -> 120 Nm note
%                     implies. The older seed is still on disk if a campaign
%                     wants the cheaper-torque gait and does not need NEC3.
%
%   'cold'            Starts from ch3_col_seed at N = 21 with the pitch box
%                     wide open, and walks the whole way: refine, seven pitch
%                     rungs, seven height rungs. Kept because it is the route
%                     that does not presuppose a converged gait on disk.
%
% PREFER THE WARM ROUTE. Measured on this project: one pitch rung starting from
% the cold gait's qt = -0.77 cost 3.1 hours across three solves (N = 41, refine
% to 61, retry) and still missed the verify tolerance by 1.5x -- not because the
% mesh was coarse (N = 61 verified the same shape at 5.3e-05) nor because the
% step was large (repose +0.129 rad), but because it hit MaxIterations with J
% still falling. Fourteen such rungs is 25-35 hours. The warm route's rungs
% instead start from a CONVERGED, VERIFIED gait and move one bound a little,
% which is the regime SQP is good at, and is exactly what the two successful
% rungs in that gait's own history did (they verified at 2.2e-05 and 1.3e-05).
%
% THE HEIGHT BAND'S CEILING IS THE KNOB, NOT ITS FLOOR. Measured on both of the
% earlier successful height stages, the gait climbs and parks on the band's
% UPPER bound (0.925 in each), so raising the ceiling is what makes it taller;
% the floor follows to keep the band tight. Move the ceiling ~0.010 m per rung:
% the one attempt that jumped it 0.020 m in a single move (centre 0.915, band
% [0.885 0.945]) came back at max|ceq| 0.12 and verify 3.0, and a jump straight
% to [0.905 0.995] was infeasible. Target ~0.94 m of hip height, 94% of the
% 1.0 m full leg extension -- close to the ~0.95 conditioning limit ch3_params
% warns about, so expect J to rise and the last rungs to be the hard ones.
%
% The cold route's pitch ladder is the one that verified in
% Results/ch3_posture_march_hist:
%   [-0.30 0.45] -> [-0.10 0.45] -> [0 0.40] -> [0.08 0.25]
%
% Resumable: state is saved after every stage and progress appended to a log,
% both named per route (see STATE/LOG below), so a MATLAB crash costs one stage
% rather than the run. Re-running picks up at the first unfinished stage.
%
% Both routes finish by writing Results/ch3_gait_lean_tall.mat and
% Results/ch3_walk_lean_tall.gif.

if nargin < 1 || isempty(route),  route = 'warm'; end
if nargin < 2 || isempty(iters),  iters = [];     end
route = validatestring(route, {'warm', 'cold'});

% Resolve Results/ from THIS FILE's location, not from the working directory.
% A relative 'Results/...' silently becomes an invalid file identifier the
% moment the session is cd'd anywhere else, and the first fprintf then dies
% with "Invalid file identifier" -- which reads like a logging bug rather than
% a path one. This file lives in <root>/Chapter3/Optimization.
root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
resd = fullfile(root, 'Results');
if ~exist(resd, 'dir'), mkdir(resd); end

% The per-route state/log names are the ones the two predecessor drivers used
% (ch3_leantall_drive and ch3_tall_from_lean). They are kept verbatim so a
% campaign already in flight when this consolidation landed still resumes.
switch route
    case 'cold'
        STATE = fullfile(resd, 'ch3_leantall_state.mat');
        LOG   = fullfile(resd, 'ch3_leantall.log');
    case 'warm'
        STATE = fullfile(resd, 'ch3_talllean_state.mat');
        LOG   = fullfile(resd, 'ch3_talllean.log');
end
% The NEC3-enforced gait, not its predecessor ch3_gait_forward_lean_tall.mat.
% See the warm-route note in the header for what that buys and what it costs.
SEED = fullfile(resd, 'ch3_gait_forward_lean_tall_nec3.mat');

stages = ladder(route);

if nargin >= 3 && ~isempty(max_stage)
    stages = stages(1:min(max_stage, numel(stages)));
end

% --- start or resume ------------------------------------------------------
if exist(STATE, 'file')
    S = load(STATE);
    k0 = S.k_done + 1;  z = S.z;  p = S.p;  hist = S.hist;
    p = ch3_upgrade_params(p);
    ch3_logln(LOG, sprintf('RESUME at stage %d/%d', k0, numel(stages)));

    % A RESUME CARRIES THE STATE'S GATES, NOT THE SEED'S -- the seed is not
    % read at all on this path. So changing SEED (or the gates inside it) has
    % NO EFFECT until this state file is gone, and the run would quietly
    % enforce less than the seed promises while looking like a normal resume.
    % Name both sets and say so, rather than leaving it to be discovered.
    gs = fieldnames(p.limits.enable)';
    gs = gs(cellfun(@(g) logical(p.limits.enable.(g)), gs));
    ch3_logln(LOG, sprintf('    state gates: %s', strjoin(gs, ' ')));

    if exist(SEED, 'file')
        Q = load(SEED);
        if     isfield(Q, 'p'),  qp = ch3_upgrade_params(Q.p);
        elseif isfield(Q, 'pf'), qp = ch3_upgrade_params(Q.pf);
        else,                    qp = [];
        end
        if ~isempty(qp)
            gq = fieldnames(qp.limits.enable)';
            gq = gq(cellfun(@(g) logical(qp.limits.enable.(g)), gq));
            missing = setdiff(gq, gs);
            if ~isempty(missing)
                ch3_logln(LOG, sprintf([ ...
                    '    NOTE: the seed enforces %s, which this state does NOT.\n' ...
                    '    The seed is not read when resuming, so those gates stay ' ...
                    'off for the\n    rest of the march. Delete %s to restart ' ...
                    'from the seed, or set them\n    on the state''s p if you ' ...
                    'want them held from here.'], ...
                    strjoin(missing, ' '), STATE));
            end
        end
    end

    Cr = ch3_col_check_limits(z, p);
    ch3_logln(LOG, sprintf('    state gait: limits max c=%.2e (ok=%d)', Cr.max_c, Cr.ok));
    if ~Cr.ok
        ch3_logln(LOG, sprintf('    %s', Cr.report));
        error('ch3_lean_tall_march:stateInfeasible', ...
              'Resume state "%s" violates a limit its own params enable:\n%s', ...
              STATE, Cr.report);
    end
elseif strcmp(route, 'cold')
    k0 = 1;  z = [];  hist = {};
    p = ch3_params('N_nodes', 21, 'enforce_nec1', false, 'qt_range', [-1 1]);

    % THE COLD STAGE RUNS BARE, which is what its own stage description says
    % ("qt box wide, height off, NEC1 off") and what the staged workflow in
    % ch3_params prescribes: solve once with the gates down, then enable them
    % one at a time from a converged gait. ch3_params has shipped every gate
    % defaulting to TRUE since e7e101b, so without this the cold solve starts
    % from ch3_col_seed carrying the full NIC/NEC set -- including NEC3, which
    % ch3_impact_march exists precisely because it cannot be switched on in one
    % move (a fresh seed sits outside the impulse friction cone by a factor of
    % six). Only clearance stays on; it is the one gate this pipeline has
    % always solved cold with.
    %
    % The height rungs below switch enable.height back on themselves, and
    % realizability is a separate march (ch3_realizability_march) that runs
    % after this one -- so nothing here silently loses a constraint it needs.
    gates = fieldnames(p.limits.enable);
    for i = 1:numel(gates)
        p.limits.enable.(gates{i}) = false;
    end
    p.limits.enable.clearance = true;

    ch3_logln(LOG, sprintf('=== LEAN+TALL march (cold), %d stages ===', numel(stages)));
else
    if ~exist(SEED, 'file')
        error('ch3_lean_tall_march:seed', ...
              ['Cannot find the warm-route seed gait "%s". Run the cold ' ...
               'route instead: ch3_lean_tall_march(''cold'').'], SEED);
    end
    S = load(SEED);
    z = S.z;

    % ACCEPT EITHER NAME. Gaits written by the solve drivers store their
    % parameters as 'p'; the older hand-built files in Results/ store them as
    % 'pf'. Hard-coding one made the seed filename and the variable name a
    % matched pair, so re-pointing SEED at a differently-written file failed
    % with "Unrecognized field name" rather than anything about seeds.
    if     isfield(S, 'p'),  p = ch3_upgrade_params(S.p);
    elseif isfield(S, 'pf'), p = ch3_upgrade_params(S.pf);
    else
        error('ch3_lean_tall_march:seedVars', ...
              ['Seed "%s" contains neither p nor pf, so there are no ' ...
               'parameters to warm-start from (has: %s).'], ...
              SEED, strjoin(fieldnames(S)', ', '));
    end

    p.qt_range             = [0.08 0.25];
    p.limits.enable.height = true;
    k0 = 1;  hist = {};
    ch3_logln(LOG, '=== LEAN+TALL march (warm): hip band 0.925 -> 0.955 ceiling ===');

    E0 = ch3_col_eval(z, p);
    V0 = ch3_col_verify(z, p, false);
    C0 = ch3_col_check_limits(z, p);

    % SAY WHICH GATES THIS RUN IS HOLDING. The seed's own enable struct decides
    % that -- enforcing NEC3 through the height rungs is inherited from the
    % file, not set here -- so a log that does not name them cannot be read
    % back later to tell which constraints a result actually respects.
    gon = fieldnames(p.limits.enable)';
    gon = gon(cellfun(@(g) logical(p.limits.enable.(g)), gon));

    ch3_logln(LOG, sprintf(['    seed: N=%d v=%.4f qt=[%+.4f %+.4f] hip=[%.4f %.4f] ' ...
                        'verify %.3e (ok=%d)\n' ...
                        '    gates: %s\n' ...
                        '    |Ix|/Iz=%.4f  peak|u|=%.1f Nm  limits max c=%.2e (ok=%d)'], ...
                       size(E0.X,2), E0.L_step/E0.T, ...
                       min(E0.X(3,:)), max(E0.X(3,:)), ...
                       min(-E0.X(2,:)), max(-E0.X(2,:)), V0.max_dev, V0.ok, ...
                       strjoin(gon, ' '), ...
                       abs(E0.impulse(1))/E0.impulse(2), ...
                       max(max(abs(E0.u(:))), max(abs(E0.um(:)))), ...
                       C0.max_c, C0.ok));

    % A seed that already violates something it declares is not a warm start,
    % it is the bug this guard exists for -- fail here rather than after the
    % first rung has spent an hour inheriting it.
    if ~C0.ok
        ch3_logln(LOG, sprintf('    %s', C0.report));
        error('ch3_lean_tall_march:seedInfeasible', ...
              ['Warm-route seed "%s" violates a limit its own params ' ...
               'enable:\n%s'], SEED, C0.report);
    end
end

% --- the march ------------------------------------------------------------
for k = k0:numel(stages)
    st = stages(k);
    t0 = tic;
    if isfield(st.opt, 'iters') && ~isempty(st.opt.iters), p.max_iter = st.opt.iters; end
    if ~isempty(iters), p.max_iter = iters; end
    ch3_logln(LOG, sprintf('--- stage %d/%d [%s] %s', k, numel(stages), st.kind, st.desc));

    switch st.kind
        case 'cold'
            p = ch3_col_budget(p);
            z = ch3_col_solve(p, ch3_col_seed(p));

        case {'pitch', 'height'}
            p = ch3_col_budget(p);
            [z2, h] = ch3_posture_march(p, st.opt, p.max_iter, z, LOG);

            % ONE RETRY, AND ITS CURE DEPENDS ON WHY THE RUNG MISSED.
            % ch3_posture_march's own advice when a rung fails to verify is
            % "refine the mesh or take smaller posture steps"; the rungs here
            % are already small, so that leaves the mesh. But once the mesh is
            % at 61 a miss is no longer coarseness -- the cold march's rungs
            % failed there with exitflag 0, out of iterations rather than out
            % of progress -- so the cure becomes a bigger iteration budget.
            % Below 61 refine, at 61 buy more iterations. Either way restart
            % from the PRE-stage gait, never from the failed one, which is by
            % definition not a real trajectory.
            if isempty(h) || ~h(end).verify_ok
                if isempty(h), dev = NaN; else, dev = h(end).verify_dev; end
                if p.N_nodes < 61
                    ch3_logln(LOG, sprintf(['    rung did not verify at N=%d ' ...
                                        '(dev %.3e); refining to 61 and retrying'], ...
                                       p.N_nodes, dev));
                    [z, p] = ch3_col_remesh(z, p, 61);
                    p = ch3_col_budget(p);
                    z = ch3_col_solve(p, z);
                else
                    ch3_logln(LOG, sprintf(['    rung missed (dev %.3e); retrying ' ...
                                        'with %d iterations'], dev, 2*p.max_iter));
                    p = ch3_col_budget(p, 2 * p.max_iter);
                end
                [z2, h] = ch3_posture_march(p, st.opt, p.max_iter, z, LOG);
            end

            if isempty(h) || ~h(end).verify_ok
                if isempty(h), dev = NaN; else, dev = h(end).verify_dev; end
                ch3_logln(LOG, sprintf('STOP at stage %d: did not verify (dev %.3e)', k, dev));
                ch3_logln(LOG, 'MARKER_STOPPED');
                return;
            end

            z = z2;
            p.qt_range             = h(end).qt_range;
            p.limits.hip_h         = h(end).hip_h;
            p.limits.hip_h_tol     = h(end).hip_h_tol;
            p.limits.enable.height = h(end).height;

        case 'remesh'
            % Climb the mesh until the solution verifies as a real trajectory.
            % alpha and T carry across untouched, so each rung warm-starts from
            % the last rather than restarting the structural work.
            for Nn = st.opt.N
                [z, p] = ch3_col_remesh(z, p, Nn);
                p = ch3_col_budget(p);
                z = ch3_col_solve(p, z);
                Vr = ch3_col_verify(z, p, false);
                ch3_logln(LOG, sprintf('    N=%d: verify %.3e (ok=%d)', ...
                                   Nn, Vr.max_dev, Vr.ok));
                if Vr.ok, break; end
            end

        case 'final'
            % A DELIVERABLE GAIT IS NEVER WRITTEN VIOLATING ITS OWN LIMITS.
            % Not a warning: this file is what everything downstream loads and
            % warm-starts from, and a gait whose params claim a constraint it
            % misses propagates that claim to every consumer.
            FINAL = fullfile(resd, 'ch3_gait_lean_tall.mat');
            ch3_assert_limits(z, p, FINAL, LOG);
            R = ch3_report(z, p, struct('stability', true, 'simulate', 5));
            save(FINAL, 'z', 'p', 'R');
            [X, ~, alpha] = ch3_col_unpack(z, p);
            try
                ch3_animate(X(:,1), alpha, p, 4, ...
                            fullfile(resd, 'ch3_walk_lean_tall.gif'));
            catch ME
                ch3_logln(LOG, sprintf('animate failed: %s', ME.message));
            end
            ch3_logln(LOG, sprintf('FINAL rho=%.4f  speed=%.4f m/s', R.rho, R.speed));
    end

    % --- measure whatever this stage produced -----------------------------
    E   = ch3_col_eval(z, p);
    V   = ch3_col_verify(z, p, false);
    chk = ch3_col_check_limits(z, p);
    qt  = E.X(3,:);
    hip = -E.X(2,:);                       % pz is DOWN-positive
    ch3_logln(LOG, sprintf(['    N=%d  T=%.4f  L=%.4f  v=%.4f m/s\n' ...
                        '    qt  [%+.4f %+.4f] rad (%+.1f .. %+.1f deg)\n' ...
                        '    hip [ %.4f  %.4f] m  (bob %.4f)\n' ...
                        '    verify %.3e (ok=%d)  limits max c %.2e (ok=%d)   %.0f s'], ...
                       size(E.X,2), E.T, E.L_step, E.L_step/E.T, ...
                       min(qt), max(qt), rad2deg(min(qt)), rad2deg(max(qt)), ...
                       min(hip), max(hip), max(hip)-min(hip), ...
                       V.max_dev, V.ok, chk.max_c, chk.ok, toc(t0)));

    hist{end+1} = struct('k', k, 'kind', st.kind, 'desc', st.desc, ...
                         'z', z, 'N', size(E.X,2), 'T', E.T, ...
                         'speed', E.L_step/E.T, 'qt_lo', min(qt), 'qt_hi', max(qt), ...
                         'hip_lo', min(hip), 'hip_hi', max(hip), ...
                         'verify_dev', V.max_dev, 'verify_ok', V.ok, ...
                         'limits_ok', chk.ok, 'limits_max_c', chk.max_c); %#ok<AGROW>

    % chk TRAVELS WITH THE STATE FILE.  State is a resume point, so it is
    % written even when the stage went wrong -- but it is also, in practice, a
    % warm-start seed for other marches (ch3_talllean_state.mat is one), and a
    % seed that cannot say whether it was feasible is how the eight-file
    % incident in ch3_upgrade_params/merge_gates stayed invisible.
    k_done = k; %#ok<NASGU>
    save(STATE, 'k_done', 'z', 'p', 'hist', 'chk');

    if ~chk.ok
        ch3_logln(LOG, sprintf('    %s', chk.report));
        ch3_logln(LOG, sprintf(['STOP at stage %d: the gait violates a limit its ' ...
                            'own params enable. Not a mesh problem -- verify ' ...
                            'passed at %.3e.'], k, V.max_dev));
        ch3_logln(LOG, 'MARKER_STOPPED');
        return;
    end

    if ~V.ok && ~strcmp(st.kind, 'final')
        if st.soft
            ch3_logln(LOG, sprintf(['    coarse (%.3e > %.1e) -- expected here; ' ...
                                'the next stage refines'], V.max_dev, p.verify_tol));
        else
            ch3_logln(LOG, sprintf('STOP at stage %d: verify failed (%.3e > %.1e)', ...
                               k, V.max_dev, p.verify_tol));
            ch3_logln(LOG, 'MARKER_STOPPED');
            return;
        end
    end
end

ch3_logln(LOG, 'MARKER_ALLDONE');

end

% ---------------------------------------------------------------- helpers
function stages = ladder(route)
% NOTE: keep comments OUT of the bracketed stage lists below. Inside [ ], a
% comment following a comma acts as a ROW separator and silently turns the
% 1 x N stage array into something else.
switch route
    case 'cold'
        % REFINEMENT COMES SECOND, NOT SIXTH. A coarse cold solve failing
        % ch3_col_verify is the EXPECTED outcome, not a reason to stop: the
        % documented cure is ch3_col_remesh warm-started from it. Measured
        % here, the N = 21 cold solve landed 7.7e-03 from a true rollout
        % against a 1e-03 tolerance. So the cold stage is marked soft (it may
        % fail verification) and the very next stage refines; every stage after
        % that must verify, because from there on a failure means the POSTURE
        % step was too big, which refinement does not fix.
        stages = [ ...
            mk('cold',   'cold solve, qt box wide, height off, NEC1 off', ...
               struct('iters', 500), true), ...
            mk('remesh', 'remesh 21 -> 41 (-> 61 if still coarse)', ...
               struct('N', [41 61], 'iters', 400)), ...
            mk('pitch',  'pitch box [-0.62 0.45]', struct('qt_range', [-0.62 0.45], 'iters', 300)), ...
            mk('pitch',  'pitch box [-0.52 0.45]', struct('qt_range', [-0.52 0.45], 'iters', 300)), ...
            mk('pitch',  'pitch box [-0.42 0.45]', struct('qt_range', [-0.42 0.45], 'iters', 300)), ...
            mk('pitch',  'pitch box [-0.30 0.45]', struct('qt_range', [-0.30 0.45], 'iters', 300)), ...
            mk('pitch',  'pitch box [-0.10 0.45]', struct('qt_range', [-0.10 0.45], 'iters', 300)), ...
            mk('pitch',  'pitch box [ 0.00 0.40]', struct('qt_range', [ 0.00 0.40], 'iters', 300)), ...
            mk('pitch',  'pitch box [ 0.08 0.25]  FORWARD LEAN', ...
               struct('qt_range', [0.08 0.25], 'iters', 300)), ...
            mk('height', 'hip band [0.835 0.925]', ...
               struct('hip_h', 0.8800, 'hip_h_tol', 0.0450, 'height', true, 'iters', 300)), ...
            mk('height', 'hip band [0.868 0.925]', ...
               struct('hip_h', 0.8965, 'hip_h_tol', 0.0285, 'height', true, 'iters', 300)), ...
            mk('height', 'hip band [0.888 0.925]', ...
               struct('hip_h', 0.9065, 'hip_h_tol', 0.0185, 'height', true, 'iters', 300)), ...
            mk('height', 'hip band [0.900 0.935]', ...
               struct('hip_h', 0.9175, 'hip_h_tol', 0.0175, 'height', true, 'iters', 300)), ...
            mk('height', 'hip band [0.912 0.945]', ...
               struct('hip_h', 0.9285, 'hip_h_tol', 0.0165, 'height', true, 'iters', 300)), ...
            mk('height', 'hip band [0.921 0.950]', ...
               struct('hip_h', 0.9355, 'hip_h_tol', 0.0145, 'height', true, 'iters', 300)), ...
            mk('height', 'hip band [0.930 0.955]  TARGET ~0.94 m', ...
               struct('hip_h', 0.9425, 'hip_h_tol', 0.0125, 'height', true, 'iters', 300)), ...
            mk('final',  'report + save', struct())];

    case 'warm'
        % Bigger per-rung budgets than the cold ladder: these rungs start from
        % a converged N = 61 gait, so iterations are the binding resource
        % rather than mesh quality.
        stages = [ ...
            mk('height', 'hip band [0.888 0.935]', ...
               struct('hip_h', 0.9115, 'hip_h_tol', 0.0235, 'height', true, 'iters', 400)), ...
            mk('height', 'hip band [0.900 0.945]', ...
               struct('hip_h', 0.9225, 'hip_h_tol', 0.0225, 'height', true, 'iters', 400)), ...
            mk('height', 'hip band [0.915 0.950]', ...
               struct('hip_h', 0.9325, 'hip_h_tol', 0.0175, 'height', true, 'iters', 400)), ...
            mk('height', 'hip band [0.930 0.955]  TARGET ~0.94 m', ...
               struct('hip_h', 0.9425, 'hip_h_tol', 0.0125, 'height', true, 'iters', 400)), ...
            mk('final',  'report + save', struct())];
end
end

function s = mk(kind, desc, opt, soft)
% soft = this stage is allowed to fail ch3_col_verify without stopping the
% march (only the cold solve, whose failure refinement is meant to fix).
if nargin < 4 || isempty(soft), soft = false; end
s = struct('kind', kind, 'desc', desc, 'opt', opt, 'soft', soft);
end
