function ch3_tall_from_lean(iters, max_stage)
%CH3_TALL_FROM_LEAN  Raise the forward-lean gait's hip to ~0.94 m.
%
%   ch3_tall_from_lean              full march, per-stage iteration budgets
%   ch3_tall_from_lean(2, 1)        smoke test: 2 iterations, one rung
%
% Starts from Results/ch3_gait_forward_lean_tall.mat -- an N = 61 gait already
% at the target pitch (qt +0.080 .. +0.166 rad) with the hip at 0.905-0.925 m --
% and walks ONLY the height band upward, ceiling 0.925 -> 0.955.
%
% WHY NOT MARCH POSTURE FROM THE COLD SOLVE. Measured on this project: one
% pitch rung starting from the cold gait's qt = -0.77 cost 3.1 hours across
% three solves (N = 41, refine to 61, retry) and still missed the verify
% tolerance by 1.5x -- not because the mesh was coarse (N = 61 verified the
% same shape at 5.3e-05) nor because the step was large (repose +0.129 rad),
% but because it hit MaxIterations with J still falling. Fourteen such rungs is
% 25-35 hours. The rungs here instead start from a CONVERGED, VERIFIED gait and
% move one bound a little, which is the regime SQP is good at, and is exactly
% what the two successful rungs in this gait's own history did (they verified
% at 2.2e-05 and 1.3e-05).
%
% THE CEILING IS THE KNOB. On both of those earlier rungs the gait climbed and
% parked on the band's UPPER bound, so raising the ceiling is what makes it
% taller; the floor follows to keep the band tight. ~0.94 m is 94% of the 1.0 m
% full leg extension, close to the ~0.95 limit ch3_params warns about, so expect
% J to rise and the last rungs to be the hard ones.
%
% Resumable: state after every rung in Results/ch3_talllean_state.mat, progress
% appended to Results/ch3_talllean.log.

if nargin < 1 || isempty(iters), iters = []; end

root  = fileparts(fileparts(fileparts(mfilename('fullpath'))));
resd  = fullfile(root, 'Results');
if ~exist(resd, 'dir'), mkdir(resd); end
STATE = fullfile(resd, 'ch3_talllean_state.mat');
LOG   = fullfile(resd, 'ch3_talllean.log');
SEED  = fullfile(resd, 'ch3_gait_forward_lean_tall.mat');

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

if exist(STATE, 'file')
    S = load(STATE);
    k0 = S.k_done + 1;  z = S.z;  p = S.p;  hist = S.hist;
    logln(LOG, sprintf('RESUME at rung %d/%d', k0, numel(stages)));
else
    if ~exist(SEED, 'file')
        error('ch3_tall_from_lean:seed', 'Cannot find the seed gait "%s".', SEED);
    end
    S = load(SEED);
    z = S.z;
    p = ch3_upgrade_params(S.pf);       % this file stores its params as pf
    p.qt_range             = [0.08 0.25];
    p.limits.enable.height = true;
    k0 = 1;  hist = {};
    logln(LOG, '=== TALL-from-LEAN march: hip band 0.925 -> 0.955 ceiling ===');

    E0 = ch3_col_eval(z, p);
    V0 = ch3_col_verify(z, p, false);
    logln(LOG, sprintf(['    seed: N=%d v=%.4f qt=[%+.4f %+.4f] hip=[%.4f %.4f] ' ...
                        'verify %.3e (ok=%d)'], size(E0.X,2), E0.L_step/E0.T, ...
                       min(E0.X(3,:)), max(E0.X(3,:)), ...
                       min(-E0.X(2,:)), max(-E0.X(2,:)), V0.max_dev, V0.ok));
end

if nargin >= 2 && ~isempty(max_stage)
    stages = stages(1:min(max_stage, numel(stages)));
end

for k = k0:numel(stages)
    st = stages(k);
    t0 = tic;
    if isfield(st.opt, 'iters') && ~isempty(st.opt.iters), p.max_iter = st.opt.iters; end
    if ~isempty(iters), p.max_iter = iters; end
    logln(LOG, sprintf('--- rung %d/%d [%s] %s', k, numel(stages), st.kind, st.desc));

    switch st.kind
        case 'height'
            p = budget(p);
            [z2, h] = ch3_posture_march(p, st.opt, p.max_iter, z);

            % The cold march's rungs failed with exitflag 0 -- out of
            % iterations, not out of progress -- so buy more before giving up,
            % restarting from the PRE-rung gait rather than the failed one.
            if isempty(h) || ~h(end).verify_ok
                if isempty(h), dev = NaN; else, dev = h(end).verify_dev; end
                logln(LOG, sprintf(['    rung missed (dev %.3e); retrying with ' ...
                                    '%d iterations'], dev, 2*p.max_iter));
                p.max_iter = 2 * p.max_iter;
                p = budget(p);
                [z2, h] = ch3_posture_march(p, st.opt, p.max_iter, z);
            end

            if isempty(h) || ~h(end).verify_ok
                if isempty(h), dev = NaN; else, dev = h(end).verify_dev; end
                logln(LOG, sprintf('STOP at rung %d: did not verify (dev %.3e)', k, dev));
                logln(LOG, 'MARKER_STOPPED');
                return;
            end

            z = z2;
            p.limits.hip_h         = h(end).hip_h;
            p.limits.hip_h_tol     = h(end).hip_h_tol;
            p.limits.enable.height = h(end).height;

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
            logln(LOG, sprintf('FINAL rho=%.4f speed=%.4f m/s', R.rho, R.speed));
    end

    E   = ch3_col_eval(z, p);
    V   = ch3_col_verify(z, p, false);
    qt  = E.X(3,:);
    hip = -E.X(2,:);                       % pz is DOWN-positive
    logln(LOG, sprintf(['    N=%d  T=%.4f  L=%.4f  v=%.4f m/s\n' ...
                        '    qt  [%+.4f %+.4f] rad (%+.1f .. %+.1f deg)\n' ...
                        '    hip [ %.4f  %.4f] m  (bob %.4f)\n' ...
                        '    verify %.3e (ok=%d)   %.0f s'], ...
                       size(E.X,2), E.T, E.L_step, E.L_step/E.T, ...
                       min(qt), max(qt), rad2deg(min(qt)), rad2deg(max(qt)), ...
                       min(hip), max(hip), max(hip)-min(hip), ...
                       V.max_dev, V.ok, toc(t0)));

    hist{end+1} = struct('k', k, 'desc', st.desc, 'z', z, ...
                         'N', size(E.X,2), 'speed', E.L_step/E.T, ...
                         'qt_lo', min(qt), 'qt_hi', max(qt), ...
                         'hip_lo', min(hip), 'hip_hi', max(hip), ...
                         'verify_dev', V.max_dev, 'verify_ok', V.ok); %#ok<AGROW>

    k_done = k; %#ok<NASGU>
    save(STATE, 'k_done', 'z', 'p', 'hist');
end

logln(LOG, 'MARKER_ALLDONE');

end

% ---------------------------------------------------------------- helpers
function s = mk(kind, desc, opt)
s = struct('kind', kind, 'desc', desc, 'opt', opt);
end

function p = budget(p)
% Make MaxIterations the binding cap rather than MaxFunctionEvaluations.
% Central differences cost ~2*n_vars evaluations per gradient, so the default
% 3e5 silently caps an 879-variable solve at ~170 iterations.
n_vars = 14 * p.N_nodes + 1 + p.ny * p.n_ctrl;
p.max_fun_evals = ceil(2.2 * n_vars * p.max_iter);
end

function logln(f, msg)
fid = fopen(f, 'a');
if fid < 0
    error('ch3_tall_from_lean:log', 'Cannot open the log file "%s".', f);
end
fprintf(fid, '%s  %s\n', datestr(now, 'HH:MM:SS'), msg); %#ok<TNOW1,DATST>
fclose(fid);          % close every time -- that is what makes it readable live
fprintf('%s\n', msg);
end
