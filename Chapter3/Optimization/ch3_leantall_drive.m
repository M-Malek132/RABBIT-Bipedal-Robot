function ch3_leantall_drive(iters, max_stage)
%CH3_LEANTALL_DRIVE  Cold solve -> forward lean -> hip band with floor 0.888.
%
%   ch3_leantall_drive              full march, p.max_iter per stage
%   ch3_leantall_drive(2, 2)        smoke test: 2 iterations, stop after stage 2
%
% Resumable: state is saved after every stage to Results/ch3_leantall_state.mat
% and progress is appended to Results/ch3_leantall.log. Re-running picks up at
% the first unfinished stage, so a MATLAB crash costs one stage, not the run.
%
% The pitch ladder is the one that verified in Results/ch3_posture_march_hist:
%   [-0.30 0.45] -> [-0.10 0.45] -> [0 0.40] -> [0.08 0.25]
% The height ladder mirrors ch3_gait_forward_lean_tall's, with the band FLOOR
% raised to 0.888 at the last rung (its ceiling 0.925 is where both successful
% stages settled, so it is held fixed and only the floor moves).

if nargin < 1 || isempty(iters), iters = []; end

% Resolve Results/ from THIS FILE's location, not from the working directory.
% A relative 'Results/...' silently becomes an invalid file identifier the
% moment the session is cd'd anywhere else, and the first fprintf then dies
% with "Invalid file identifier" -- which reads like a logging bug rather than
% a path one. This file lives in <root>/Chapter3/Optimization.
root  = fileparts(fileparts(fileparts(mfilename('fullpath'))));
resd  = fullfile(root, 'Results');
if ~exist(resd, 'dir'), mkdir(resd); end
STATE = fullfile(resd, 'ch3_leantall_state.mat');
LOG   = fullfile(resd, 'ch3_leantall.log');

% REFINEMENT COMES SECOND, NOT SIXTH. A coarse cold solve failing
% ch3_col_verify is the EXPECTED outcome, not a reason to stop: the documented
% cure is ch3_col_remesh warm-started from it. Measured here, the N = 21 cold
% solve landed 7.7e-03 from a true rollout against a 1e-03 tolerance. So the
% cold stage is marked soft (it may fail verification) and the very next stage
% refines; every stage after that must verify, because from there on a failure
% means the POSTURE step was too big, which refinement does not fix.
%
% THE HEIGHT BAND'S CEILING IS THE KNOB, NOT ITS FLOOR. Measured on both of the
% earlier successful height stages, the gait climbs and parks on the band's
% UPPER bound (0.925 in each), so raising the ceiling is what makes it taller.
% Move it ~0.010 m per rung: the one attempt that jumped it 0.020 m in a single
% move (centre 0.915, band [0.885 0.945]) came back at max|ceq| 0.12 and verify
% 3.0, and a jump straight to [0.905 0.995] was infeasible. Target ~0.94 m of
% hip height, 94% of the 1.0 m full leg extension -- close to the ~0.95
% conditioning limit ch3_params warns about, so expect J to rise and some rungs
% to need a second pass.
%
% NOTE: keep comments OUT of the bracketed stage list below. Inside [ ], a
% comment following a comma acts as a ROW separator and silently turns the
% 1 x N stage array into something else.
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

if exist(STATE, 'file')
    S = load(STATE);
    k0 = S.k_done + 1;  z = S.z;  p = S.p;  hist = S.hist;
    logln(LOG, sprintf('RESUME at stage %d/%d', k0, numel(stages)));
else
    k0 = 1;
    p = ch3_params('N_nodes', 21, 'enforce_nec1', false, 'qt_range', [-1 1]);
    z = [];  hist = {};
    logln(LOG, sprintf('=== LEAN+TALL march, %d stages, cold from seed ===', numel(stages)));
end

if nargin >= 2 && ~isempty(max_stage)
    stages = stages(1:min(max_stage, numel(stages)));
end

for k = k0:numel(stages)
    st = stages(k);
    t0 = tic;
    if isfield(st.opt, 'iters') && ~isempty(st.opt.iters), p.max_iter = st.opt.iters; end
    if ~isempty(iters), p.max_iter = iters; end
    logln(LOG, sprintf('--- stage %d/%d [%s] %s', k, numel(stages), st.kind, st.desc));

    switch st.kind
        case 'cold'
            p = budget(p);
            z = ch3_col_solve(p, ch3_col_seed(p));

        case {'pitch', 'height'}
            p = budget(p);
            [z2, h] = ch3_posture_march(p, st.opt, p.max_iter, z);

            % ch3_posture_march's own advice when a rung fails to verify is
            % "refine the mesh or take smaller posture steps". The rungs here
            % are already small, so refine and retry the SAME target once --
            % from the pre-stage gait, never from the failed one, which is by
            % definition not a real trajectory.
            if (isempty(h) || ~h(end).verify_ok) && p.N_nodes < 61
                if isempty(h), dev = NaN; else, dev = h(end).verify_dev; end
                logln(LOG, sprintf(['    rung did not verify at N=%d ' ...
                                    '(dev %.3e); refining to 61 and retrying'], ...
                                   p.N_nodes, dev));
                [z, p] = ch3_col_remesh(z, p, 61);
                p = budget(p);
                z = ch3_col_solve(p, z);
                [z2, h] = ch3_posture_march(p, st.opt, p.max_iter, z);
            end

            if isempty(h) || ~h(end).verify_ok
                logln(LOG, sprintf('STOP at stage %d: stage did not verify', k));
                logln(LOG, 'MARKER_STOPPED');
                return;
            end
            z = z2;
            p.qt_range           = h(end).qt_range;
            p.limits.hip_h       = h(end).hip_h;
            p.limits.hip_h_tol   = h(end).hip_h_tol;
            p.limits.enable.height = h(end).height;

        case 'remesh'
            % Climb the mesh until the solution verifies as a real trajectory.
            % alpha and T carry across untouched, so each rung warm-starts from
            % the last rather than restarting the structural work.
            for Nn = st.opt.N
                [z, p] = ch3_col_remesh(z, p, Nn);
                p = budget(p);
                z = ch3_col_solve(p, z);
                Vr = ch3_col_verify(z, p, false);
                logln(LOG, sprintf('    N=%d: verify %.3e (ok=%d)', ...
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
                logln(LOG, sprintf('animate failed: %s', ME.message));
            end
            logln(LOG, sprintf('FINAL rho=%.4f  speed=%.4f m/s', R.rho, R.speed));
    end

    % --- measure whatever this stage produced -----------------------------
    E   = ch3_col_eval(z, p);
    V   = ch3_col_verify(z, p, false);
    qt  = E.X(3,:);
    hip = -E.X(2,:);                       % pz is DOWN-positive
    line = sprintf(['    N=%d  T=%.4f  L=%.4f  v=%.4f m/s\n' ...
                    '    qt  [%+.4f %+.4f] rad (%+.1f .. %+.1f deg)\n' ...
                    '    hip [ %.4f  %.4f] m\n' ...
                    '    verify %.3e (ok=%d)   %.0f s'], ...
                   size(E.X,2), E.T, E.L_step, E.L_step/E.T, ...
                   min(qt), max(qt), rad2deg(min(qt)), rad2deg(max(qt)), ...
                   min(hip), max(hip), V.max_dev, V.ok, toc(t0));
    logln(LOG, line);

    hist{end+1} = struct('k', k, 'kind', st.kind, 'desc', st.desc, ...
                         'z', z, 'N', size(E.X,2), 'T', E.T, ...
                         'speed', E.L_step/E.T, 'qt_lo', min(qt), 'qt_hi', max(qt), ...
                         'hip_lo', min(hip), 'hip_hi', max(hip), ...
                         'verify_dev', V.max_dev, 'verify_ok', V.ok); %#ok<AGROW>

    k_done = k; %#ok<NASGU>
    save(STATE, 'k_done', 'z', 'p', 'hist');

    if ~V.ok && ~strcmp(st.kind, 'final')
        if st.soft
            logln(LOG, sprintf(['    coarse (%.3e > %.1e) -- expected here; ' ...
                                'the next stage refines'], V.max_dev, p.verify_tol));
        else
            logln(LOG, sprintf('STOP at stage %d: verify failed (%.3e > %.1e)', ...
                               k, V.max_dev, p.verify_tol));
            logln(LOG, 'MARKER_STOPPED');
            return;
        end
    end
end

logln(LOG, 'MARKER_ALLDONE');

end

% ---------------------------------------------------------------- helpers
function s = mk(kind, desc, opt, soft)
% soft = this stage is allowed to fail ch3_col_verify without stopping the
% march (only the cold solve, whose failure refinement is meant to fix).
if nargin < 4 || isempty(soft), soft = false; end
s = struct('kind', kind, 'desc', desc, 'opt', opt, 'soft', soft);
end

function p = budget(p)
% Make MaxIterations the binding cap rather than MaxFunctionEvaluations.
% Central differences cost ~2*n_vars evaluations per gradient, so the default
% 3e5 silently caps a 599-variable solve at ~250 iterations.
n_vars = 14 * p.N_nodes + 1 + p.ny * p.n_ctrl;
p.max_fun_evals = ceil(2.2 * n_vars * p.max_iter);
end

function logln(f, msg)
fid = fopen(f, 'a');
if fid < 0
    error('ch3_leantall_drive:log', ...
          'Cannot open the log file "%s" for append.', f);
end
fprintf(fid, '%s  %s\n', datestr(now, 'HH:MM:SS'), msg); %#ok<TNOW1,DATST>
fclose(fid);          % close every time -- that is what makes it readable live
fprintf('%s\n', msg);
end
