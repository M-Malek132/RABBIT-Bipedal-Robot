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

STATE = fullfile('Results', 'ch3_leantall_state.mat');
LOG   = fullfile('Results', 'ch3_leantall.log');

stages = [ ...
    mk('cold',   'cold solve, qt box wide, height off, NEC1 off', struct()), ...
    mk('pitch',  'pitch box [-0.30 0.45]', struct('qt_range', [-0.30 0.45])), ...
    mk('pitch',  'pitch box [-0.10 0.45]', struct('qt_range', [-0.10 0.45])), ...
    mk('pitch',  'pitch box [ 0.00 0.40]', struct('qt_range', [ 0.00 0.40])), ...
    mk('pitch',  'pitch box [ 0.08 0.25]  FORWARD LEAN', struct('qt_range', [0.08 0.25])), ...
    mk('remesh', 'remesh 21 -> 41', struct('N', 41)), ...
    mk('height', 'hip band [0.835 0.925]', struct('hip_h', 0.8800, 'hip_h_tol', 0.0450, 'height', true)), ...
    mk('height', 'hip band [0.868 0.925]', struct('hip_h', 0.8965, 'hip_h_tol', 0.0285, 'height', true)), ...
    mk('height', 'hip band [0.888 0.925]  TARGET', struct('hip_h', 0.9065, 'hip_h_tol', 0.0185, 'height', true)), ...
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
    if ~isempty(iters), p.max_iter = iters; end
    logln(LOG, sprintf('--- stage %d/%d [%s] %s', k, numel(stages), st.kind, st.desc));

    switch st.kind
        case 'cold'
            p = budget(p);
            z = ch3_col_solve(p, ch3_col_seed(p));

        case {'pitch', 'height'}
            p = budget(p);
            [z2, h] = ch3_posture_march(p, st.opt, p.max_iter, z);
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
            [z, p] = ch3_col_remesh(z, p, st.opt.N);
            p = budget(p);
            z = ch3_col_solve(p, z);

        case 'final'
            R = ch3_report(z, p, struct('stability', true, 'simulate', 5));
            save(fullfile('Results','ch3_gait_lean_tall_888.mat'), 'z', 'p', 'R');
            [X, ~, alpha] = ch3_col_unpack(z, p);
            try
                ch3_animate(X(:,1), alpha, p, 4, fullfile('Results','ch3_walk_lean_tall_888.gif'));
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
        logln(LOG, sprintf('STOP at stage %d: verify failed (%.3e > %.1e)', ...
                           k, V.max_dev, p.verify_tol));
        logln(LOG, 'MARKER_STOPPED');
        return;
    end
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
% 3e5 silently caps a 599-variable solve at ~250 iterations.
n_vars = 14 * p.N_nodes + 1 + p.ny * p.n_ctrl;
p.max_fun_evals = ceil(2.2 * n_vars * p.max_iter);
end

function logln(f, msg)
fid = fopen(f, 'a');
fprintf(fid, '%s  %s\n', datestr(now, 'HH:MM:SS'), msg); %#ok<TNOW1,DATST>
fclose(fid);          % close every time -- that is what makes it readable live
fprintf('%s\n', msg);
end
