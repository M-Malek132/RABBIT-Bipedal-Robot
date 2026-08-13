function ch3_realizability_march()
%CH3_REALIZABILITY_MARCH  Finish Table 3.1: enable torque, then impulse.
%
%   ch3_realizability_march()
%
% Starting point is Results/ch3_gait_fix.mat -- an N = 61 gait already
% verified with friction, GRF, clearance, hip-height, NEC3 impulse validity
% and the HZD/HH gates all ENFORCED (see ch3_report on that file). Only two
% Table 3.1 limits are still off there:
%
%   peak |torque|     191.4 Nm measured  vs  u_max        = 120 Nm
%   impact impulse     18.8 Ns measured  vs  impulse_max  =  15 Ns
%
% Both are plain box/linear inequalities on quantities evaluated at nodes and
% midpoints (torque) or at the single impact endpoint (impulse) -- neither has
% the division-by-near-zero pathology that makes friction-before-GRF unsafe, so
% there is no ordering hazard here. They ARE marched rather than switched on at
% their final value in one move anyway, for the same reason every other limit
% in this pipeline is: torque needs to come down 37%, which is a big step for a
% warm start to survive in one SQP solve, and a failed jump wastes the entire
% solve rather than just the last rung of it.
%
% Torque goes first, impulse second -- torque is the larger relative move
% (37% vs 20%) and enabling it can itself move the impulse (both depend on the
% same orbit), so there is no point spending the smaller, harder-to-move
% impulse budget before the larger one has settled.
%
% Resumable exactly like ch3_lean_tall_march: state is saved after every stage
% to Results/ch3_realizability_state.mat and progress appended to
% Results/ch3_realizability.log, so a MATLAB crash (this install's -batch has a
% ~40% chance of dying in libcurl before user code even runs -- unrelated to
% this code) costs one stage, not the run. Re-running picks up where it left
% off.
%
% Finishes by writing Results/ch3_gait_full_constrained.mat (z, p, R) and
% Results/ch3_walk_full_constrained.gif.
%
% See also CH3_IMPACT_MARCH, CH3_POSTURE_MARCH, CH3_LEAN_TALL_MARCH.

root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
resd = fullfile(root, 'Results');
STATE = fullfile(resd, 'ch3_realizability_state.mat');
LOG   = fullfile(resd, 'ch3_realizability.log');
SEED  = fullfile(resd, 'ch3_gait_fix.mat');
FINAL = fullfile(resd, 'ch3_gait_full_constrained.mat');

stages = ladder();

if exist(STATE, 'file')
    S = load(STATE);
    k0 = S.k_done + 1;  z = S.z;  p = S.p;  hist = S.hist;
    logln(LOG, sprintf('RESUME at stage %d/%d', k0, numel(stages)));
else
    if ~exist(SEED, 'file')
        error('ch3_realizability_march:seed', 'Cannot find seed gait "%s".', SEED);
    end
    S = load(SEED);
    z = S.z;
    p = ch3_upgrade_params(S.pf);
    k0 = 1;  hist = {};
    logln(LOG, sprintf('=== REALIZABILITY march: %d stages ===', numel(stages)));

    E0 = ch3_col_eval(z, p);
    V0 = ch3_col_verify(z, p, false);
    peak_u0 = max(max(abs(E0.u(:))), max(abs(E0.um(:))));
    logln(LOG, sprintf(['    seed: N=%d peak|u|=%.2f Nm  ||impulse||=%.3f Ns  ' ...
                        'verify %.3e (ok=%d)'], size(E0.X,2), peak_u0, ...
                       norm(E0.impulse), V0.max_dev, V0.ok));
end

for k = k0:numel(stages)
    st = stages(k);
    t0 = tic;
    p = budget(p, st.iters);
    logln(LOG, sprintf('--- stage %d/%d [%s] %s', k, numel(stages), st.kind, st.desc));

    switch st.kind
        case 'torque'
            p.limits.enable.torque = true;
            p.limits.u_max = st.target;
        case 'impulse'
            p.limits.enable.impulse = true;
            p.limits.impulse_max = st.target;
        case 'final'
            R = ch3_report(z, p, struct('stability', true, 'simulate', 5));
            save(FINAL, 'z', 'p', 'R');
            [X, ~, alpha] = ch3_col_unpack(z, p);
            try
                ch3_animate(X(:,1), alpha, p, 4, ...
                            fullfile(resd, 'ch3_walk_full_constrained.gif'));
            catch ME
                logln(LOG, sprintf('animate failed: %s', ME.message));
            end
            logln(LOG, sprintf('FINAL rho=%.4f  speed=%.4f m/s  file=%s', ...
                               R.rho, R.speed, FINAL));
            k_done = k; %#ok<NASGU>
            save(STATE, 'k_done', 'z', 'p', 'hist');
            logln(LOG, 'MARKER_ALLDONE');
            return;
    end

    [z_new, out] = ch3_col_solve(p, z);
    V = ch3_col_verify(z_new, p, false);

    % ONE RETRY: bigger iteration budget first (this gait is already on a
    % fine-enough mesh, N=61, so a miss here is exhausted iterations rather
    % than mesh coarseness -- see ch3_lean_tall_march for the same reasoning).
    % If it still misses, refine the mesh next.
    if ~V.ok
        logln(LOG, sprintf(['    stage missed (verify %.3e); retrying with ' ...
                            '%d iterations'], V.max_dev, 2*p.max_iter));
        p2 = budget(p, 2*st.iters);
        [z_new, out] = ch3_col_solve(p2, z);
        V = ch3_col_verify(z_new, p2, false);
        if V.ok, p = p2; end
    end

    if ~V.ok
        logln(LOG, sprintf(['    stage missed again (verify %.3e); refining ' ...
                            'mesh 61 -> 81 and retrying'], V.max_dev));
        [z_r, p_r] = ch3_col_remesh(z, p, 81);
        p_r = budget(p_r, st.iters);
        [z_new, out] = ch3_col_solve(p_r, z_r);
        V = ch3_col_verify(z_new, p_r, false);
        if V.ok, p = p_r; end
    end

    E = ch3_col_eval(z_new, p);
    peak_u = max(max(abs(E.u(:))), max(abs(E.um(:))));
    logln(LOG, sprintf(['    N=%d  J=%.2f  exitflag=%d  max|ceq|=%.2e  max c=%.2e\n' ...
                        '    peak|u|=%.3f Nm  ||impulse||=%.3f Ns  speed=%.4f m/s\n' ...
                        '    verify %.3e (ok=%d)   %.0f s'], ...
                       size(E.X,2), out.fval, out.exitflag, out.max_ceq, out.max_c, ...
                       peak_u, norm(E.impulse), E.L_step/E.T, V.max_dev, V.ok, toc(t0)));

    hist{end+1} = struct('k', k, 'kind', st.kind, 'target', st.target, ...
                         'z', z_new, 'N', size(E.X,2), 'peak_u', peak_u, ...
                         'impulse_norm', norm(E.impulse), 'speed', E.L_step/E.T, ...
                         'verify_dev', V.max_dev, 'verify_ok', V.ok); %#ok<AGROW>

    if ~V.ok
        logln(LOG, sprintf('STOP at stage %d: did not verify (%.3e)', k, V.max_dev));
        logln(LOG, 'MARKER_STOPPED');
        k_done = k - 1; %#ok<NASGU>
        save(STATE, 'k_done', 'z', 'p', 'hist');
        return;
    end

    z = z_new;
    k_done = k; %#ok<NASGU>
    save(STATE, 'k_done', 'z', 'p', 'hist');
end

end

% ---------------------------------------------------------------- helpers
function stages = ladder()
% Torque: 191.4 -> 120 Nm, four rungs. Impulse: 18.8 -> 15 Ns, three rungs.
% Step sizes follow the same "move the bound a little, warm-start" logic as
% every other march in this package.
stages = [ ...
    mk('torque',  'u_max = 170 Nm', 170, 300), ...
    mk('torque',  'u_max = 150 Nm', 150, 300), ...
    mk('torque',  'u_max = 135 Nm', 135, 300), ...
    mk('torque',  'u_max = 120 Nm  TARGET', 120, 300), ...
    mk('impulse', 'impulse_max = 17.5 Ns', 17.5, 300), ...
    mk('impulse', 'impulse_max = 16.0 Ns', 16.0, 300), ...
    mk('impulse', 'impulse_max = 15.0 Ns  TARGET', 15.0, 300), ...
    mk('final',   'report + save', [], [])];
end

function s = mk(kind, desc, target, iters)
s = struct('kind', kind, 'desc', desc, 'target', target, 'iters', iters);
end

function p = budget(p, iters)
if ~isempty(iters), p.max_iter = iters; end
n_vars = 14 * p.N_nodes + 1 + p.ny * p.n_ctrl;
p.max_fun_evals = ceil(2.2 * n_vars * p.max_iter);
end

function logln(f, msg)
fid = fopen(f, 'a');
if fid < 0
    error('ch3_realizability_march:log', 'Cannot open log file "%s".', f);
end
fprintf(fid, '%s  %s\n', datestr(now, 'HH:MM:SS'), msg); %#ok<TNOW1,DATST>
fclose(fid);
fprintf('%s\n', msg);
end
