function lib = ch6_lib_build(p, seed_file, out_file)
%CH6_LIB_BUILD  Solve the gait library of Section 6.4.  (6.24)
%
%   lib = ch6_lib_build(p)
%   lib = ch6_lib_build(p, seed_file, out_file)
%
% Solves one gait per entry of p.lib.targets and writes them to a .mat that
% ch6_lib_load reads. Each solve is ch6_lib_solve -- a Chapter-3 collocation
% solve with the step length pinned.
%
% ============================================================ MARCH, DO NOT SWEEP
% The solves are ordered OUTWARD FROM THE SEED and each one warm-starts from
% its neighbour, not from the seed:
%
%       ... <- L_(i-1) <- L_i = seed -> L_(i+1) -> ...
%
% This is the same continuation argument ch3_continuation makes for speed, and
% it is not optional here either. A pinned step length is an equality on a
% quantity that depends on the whole orbit; starting a 0.20 m solve from a
% 0.35 m gait makes that equality the dominant residual from iteration one and
% SQP spends its budget fighting it. Marching keeps each individual solve nearly
% feasible at its start, which is the regime SQP is good at.
%
% ------------------------------------------------------ a gait that fails is DROPPED
% A solve whose mesh check fails is not a trajectory, and interpolating between
% a real gait and a fictional one produces a fictional gait everywhere in
% between. Such an entry is left out of the library and the march STOPS in that
% direction rather than stepping past it -- continuing from a solution that is
% not a real trajectory only propagates the problem, which is the rule
% ch3_continuation already follows.
%
% The consequence is that the library's ACHIEVED grid can be shorter than
% p.lib.targets. That grid is stored in the file and read back by
% ch6_lib_alpha, never assumed, so a short library gives more extrapolation
% rather than a wrong answer -- and Remark 6.6 is precisely about what that
% costs.
%
% ------------------------------------------------------------------- the seed
% seed_file must contain a converged Chapter-3 result (`z_opt` and `p`). The
% repo's reference gait, Results/ch3_reference_gait.mat, is the default: it is
% periodic, mesh-verified and stable, and its own step length (0.353 m) sits
% inside the default target grid so the march runs both ways.
%
% Inputs
%   p         : parameter struct (uses p.lib.targets, p.lib.iters)
%   seed_file : .mat with z_opt (default Results/ch3_reference_gait.mat)
%   out_file  : where to write (default Results/ch6_gait_library.mat)
%
% Output
%   lib : struct
%     .L        1 x m achieved step lengths, ascending
%     .alpha    ny x n_ctrl x m
%     .T .speed .fval .exitflag .verify_dev   1 x m
%     .targets  the requested grid
%     .p        the parameter struct used
%     .dropped  targets that failed and why
%
% See also CH6_LIB_SOLVE, CH6_LIB_ALPHA, CH6_LIB_LOAD, CH3_CONTINUATION.

root = fileparts(fileparts(mfilename('fullpath')));   % .../Chapter6
root = fileparts(root);                               % repo root

if nargin < 2 || isempty(seed_file)
    seed_file = fullfile(root, 'Results', 'ch3_reference_gait.mat');
end
if nargin < 3 || isempty(out_file)
    out_file = fullfile(root, 'Results', 'ch6_gait_library.mat');
end

S = load(seed_file);
if ~isfield(S, 'z_opt')
    error('ch6_lib_build:seed', ...
          'Seed file "%s" has no z_opt (fields: %s).', ...
          seed_file, strjoin(fieldnames(S).', ', '));
end
z_seed = S.z_opt;

% The seed's own parameters define the MESH, and the mesh is baked into the
% length of z. Take p's controller/limit settings but the seed's transcription,
% or the warm start is a vector of the wrong size.
p_solve = p;
if isfield(S, 'p')
    % Copy only the fields the seed actually carries. A result saved before a
    % parameter existed is still a perfectly good warm start, and the alternative
    % -- assuming the field list -- makes every older .mat in Results/ unusable
    % as a seed for a reason that has nothing to do with the gait in it.
    for f = {'N_nodes', 'bez_deg', 'n_ctrl', 'qt_range'}
        if isfield(S.p, f{1}), p_solve.(f{1}) = S.p.(f{1}); end
    end
end
p_solve = ch3_upgrade_params(p_solve);
p_solve.controller = 'ff';        % the gait is DESIGNED on Z; see ch3_main

E0 = ch3_col_eval(z_seed, p_solve);
L0 = E0.L_step;
fprintf('[ch6_lib_build] seed %s: L_step = %.4f m, T = %.3f s\n', ...
        seed_file, L0, E0.T);

% MESH. A longer step needs a finer mesh than the seed's, and the failure mode
% is silent: MEASURED, the L* = 0.40 m solve hit its target to 3e-15 on the
% seed's N = 21 mesh and still missed a true trajectory by 3.19e-02 (tol 1e-03).
% Small defects do not imply a real trajectory -- ch3_col_verify's whole reason
% for existing. Remesh once, up front, so every library gait is solved where it
% can be verified rather than discovering it one gait at a time.
if isfield(p.lib, 'N_nodes') && ~isempty(p.lib.N_nodes) && ...
        p.lib.N_nodes ~= p_solve.N_nodes
    fprintf('[ch6_lib_build] remeshing the seed %d -> %d nodes\n', ...
            p_solve.N_nodes, p.lib.N_nodes);
    [z_seed, p_solve] = ch3_col_remesh(z_seed, p_solve, p.lib.N_nodes);
end

% EVERY MARCHED GAIT IS KEPT. p.lib.targets is the march itself, not a grid to
% be reached through hidden intermediate solves. Two reasons: a denser library
% is strictly better for (6.22), since interpolation error falls with the
% spacing; and a solve that is only a stepping stone to another one still costs
% five to twenty minutes, so discarding it would be paying full price for
% nothing. Choose targets close enough together that each solve starts nearly
% feasible -- that is what the march needs -- and take the density for free.
targets = sort(p.lib.targets(:).');
m       = numel(targets);

% Split the grid at the seed and march outward in both directions.
i_up   = find(targets >= L0);
i_down = fliplr(find(targets < L0));

alpha_all = nan(p.ny, p_solve.n_ctrl, m);
L_all     = nan(1, m);
T_all     = nan(1, m);
v_all     = nan(1, m);
J_all     = nan(1, m);
ef_all    = nan(1, m);
dev_all   = nan(1, m);
ok        = false(1, m);
dropped   = {};

for dir = 1:2
    if dir == 1, idx = i_up; else, idx = i_down; end
    z = z_seed;
    for j = idx
        fprintf('\n===== library gait %d/%d : L* = %.3f m =====\n', ...
                find(targets == targets(j), 1), m, targets(j));
        [z_new, o] = ch6_lib_solve(p_solve, targets(j), z, p.lib.iters);

        if ~o.verify.ok
            dropped{end+1} = sprintf('L* = %.3f (mesh dev %.2e)', ...
                                     targets(j), o.verify.max_dev); %#ok<AGROW>
            fprintf('[ch6_lib_build] dropping L* = %.3f and stopping this direction.\n', ...
                    targets(j));
            break;
        end

        alpha_all(:,:,j) = o.alpha;
        L_all(j)   = o.L_step;
        T_all(j)   = o.T;
        v_all(j)   = o.speed;
        J_all(j)   = o.fval;
        ef_all(j)  = o.exitflag;
        dev_all(j) = o.verify.max_dev;
        ok(j)      = true;

        z = z_new;
    end
end

if ~any(ok)
    error('ch6_lib_build:empty', ...
          'No library gait converged. Check the seed and p.lib.targets.');
end

% Sort by ACHIEVED length, not by target: a solve that lands 3 mm off its target
% must not put the grid out of order, because ch6_lib_alpha's bracketing search
% assumes ascending L and would silently interpolate the wrong pair.
[Ls, ord] = sort(L_all(ok));
sel  = find(ok);
sel  = sel(ord);

lib = struct();
lib.L          = Ls;
lib.alpha      = alpha_all(:,:,sel);
lib.T          = T_all(sel);
lib.speed      = v_all(sel);
lib.fval       = J_all(sel);
lib.exitflag   = ef_all(sel);
lib.verify_dev = dev_all(sel);
lib.targets    = targets;
lib.seed_file  = seed_file;
lib.p          = p_solve;
lib.dropped    = dropped;
lib.built      = datestr(now, 'yyyy-mm-dd HH:MM:SS'); %#ok<TNOW1,DATST>

save(out_file, '-struct', 'lib');
fprintf('\n[ch6_lib_build] %d/%d gaits, L in [%.3f, %.3f] m -> %s\n', ...
        numel(Ls), m, Ls(1), Ls(end), out_file);
if ~isempty(dropped)
    fprintf('[ch6_lib_build] dropped: %s\n', strjoin(dropped, ', '));
end

end
