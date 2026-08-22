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
%   'warm' (default)  Starts from Results/ch3_gait_forward_lean_tall.mat, an
%                     N = 61 gait ALREADY at the target pitch (qt +0.080 ..
%                     +0.166 rad) with the hip at 0.905-0.925 m, and walks only
%                     the height band up: ceiling 0.925 -> 0.955, four rungs.
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
SEED = fullfile(resd, 'ch3_gait_forward_lean_tall.mat');

stages = ladder(route);

if nargin >= 3 && ~isempty(max_stage)
    stages = stages(1:min(max_stage, numel(stages)));
end

% --- start or resume ------------------------------------------------------
if exist(STATE, 'file')
    S = load(STATE);
    k0 = S.k_done + 1;  z = S.z;  p = S.p;  hist = S.hist;
    ch3_logln(LOG, sprintf('RESUME at stage %d/%d', k0, numel(stages)));
elseif strcmp(route, 'cold')
    k0 = 1;  z = [];  hist = {};
    p = ch3_params('N_nodes', 21, 'enforce_nec1', false, 'qt_range', [-1 1]);
    ch3_logln(LOG, sprintf('=== LEAN+TALL march (cold), %d stages ===', numel(stages)));
else
    if ~exist(SEED, 'file')
        error('ch3_lean_tall_march:seed', ...
              ['Cannot find the warm-route seed gait "%s". Run the cold ' ...
               'route instead: ch3_lean_tall_march(''cold'').'], SEED);
    end
    S = load(SEED);
    z = S.z;
    p = ch3_upgrade_params(S.pf);       % this file stores its params as pf
    p.qt_range             = [0.08 0.25];
    p.limits.enable.height = true;
    k0 = 1;  hist = {};
    ch3_logln(LOG, '=== LEAN+TALL march (warm): hip band 0.925 -> 0.955 ceiling ===');

    E0 = ch3_col_eval(z, p);
    V0 = ch3_col_verify(z, p, false);
    ch3_logln(LOG, sprintf(['    seed: N=%d v=%.4f qt=[%+.4f %+.4f] hip=[%.4f %.4f] ' ...
                        'verify %.3e (ok=%d)'], size(E0.X,2), E0.L_step/E0.T, ...
                       min(E0.X(3,:)), max(E0.X(3,:)), ...
                       min(-E0.X(2,:)), max(-E0.X(2,:)), V0.max_dev, V0.ok));
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
            R = ch3_report(z, p, struct('stability', true, 'simulate', 5));
            save(fullfile(resd, 'ch3_gait_lean_tall.mat'), 'z', 'p', 'R');
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
    qt  = E.X(3,:);
    hip = -E.X(2,:);                       % pz is DOWN-positive
    ch3_logln(LOG, sprintf(['    N=%d  T=%.4f  L=%.4f  v=%.4f m/s\n' ...
                        '    qt  [%+.4f %+.4f] rad (%+.1f .. %+.1f deg)\n' ...
                        '    hip [ %.4f  %.4f] m  (bob %.4f)\n' ...
                        '    verify %.3e (ok=%d)   %.0f s'], ...
                       size(E.X,2), E.T, E.L_step, E.L_step/E.T, ...
                       min(qt), max(qt), rad2deg(min(qt)), rad2deg(max(qt)), ...
                       min(hip), max(hip), max(hip)-min(hip), ...
                       V.max_dev, V.ok, toc(t0)));

    hist{end+1} = struct('k', k, 'kind', st.kind, 'desc', st.desc, ...
                         'z', z, 'N', size(E.X,2), 'T', E.T, ...
                         'speed', E.L_step/E.T, 'qt_lo', min(qt), 'qt_hi', max(qt), ...
                         'hip_lo', min(hip), 'hip_hi', max(hip), ...
                         'verify_dev', V.max_dev, 'verify_ok', V.ok); %#ok<AGROW>

    k_done = k; %#ok<NASGU>
    save(STATE, 'k_done', 'z', 'p', 'hist');

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
