function [z_opt, fval, exitflag, p] = rabbit_hzd_collocation(N, warm_start_file, max_iter, p_over)
%RABBIT_HZD_COLLOCATION  HZD gait design by Hermite-Simpson DIRECT COLLOCATION.
%
%   [z_opt, fval, exitflag, p] = rabbit_hzd_collocation(N, warm_start_file, max_iter)
%
% The direct-collocation counterpart of rabbit_hzd_trajectory_optimization.m
% (which is single shooting). This matches the thesis method: the dynamics are
% enforced as Hermite-Simpson defect CONSTRAINTS (no ode45 in the loop) and the
% decision variables are node states/torques/contact-forces, the B-spline
% (Bezier) coefficients, and the step time -- NOT a single integrated initial
% state. The virtual constraints y = q(4:7) - yd(s, coeffs) = 0 are enforced at
% every node, so the gait lives on the hybrid zero-dynamics manifold and
% `coeffs` parametrizes it. See Collocation/DOCS.md.
%
% Args (all optional):
%   N               node count (default 20)
%   warm_start_file Results/*.mat with z_opt,p to seed from a previous PD gait
%                   (default: newest hzd_result_*.mat). The PD gait is
%                   simulated and resampled to N nodes (col_seed_from_pd).
%   max_iter        override fmincon MaxIterations (default from setup)
%   p_over          struct of p fields to override (e.g. enabling physical
%                   constraints): struct('enforce_grf',true,'enforce_impulse',true)
%
% Physical Table-3.1 limits are inherited from hzd_problem_setup toggles
% (p.enforce_torque as a torque box on U; friction/GRF/impulse as inequalities);
% flip them here via p_over rather than editing hzd_problem_setup.
%
% Warm start: if warm_start_file is a COLLOCATION result (col_result_*.mat, has
% an 'N' field) its z_opt is used directly (N must match); if it is a PD /
% shooting result (hzd_result_*.mat) the step is simulated and resampled to N
% nodes. This lets you add physical constraints ONE PHASE AT A TIME, warm-
% starting each solve from the previous converged collocation gait.

    here = fileparts(mfilename('fullpath'));
    root = fileparts(fileparts(here));
    addpath(genpath(root));

    if nargin < 1 || isempty(N), N = 20; end
    if nargin < 4, p_over = []; end

    [p, lb, ub, options] = col_problem_setup(N, p_over);
    if nargin >= 3 && ~isempty(max_iter)
        options = optimoptions(options, 'MaxIterations', max_iter);
    end

    % ---- resolve warm-start file ----------------------------------------
    if nargin < 2 || isempty(warm_start_file)
        d = dir(fullfile(root,'Results','hzd_result_*.mat'));
        [~,i] = max([d.datenum]);
        warm_start_file = fullfile(d(i).folder, d(i).name);
    elseif ~isfile(warm_start_file)
        warm_start_file = fullfile(root,'Results', warm_start_file);
    end
    S = load(warm_start_file);

    % ---- build the initial decision vector z0 ---------------------------
    if isfield(S, 'N')      % a collocation result: reuse its vector directly
        assert(numel(S.z_opt) == numel(lb), ...
            'collocation warm-start N mismatch: saved N=%d, requested N=%d.', S.N, N);
        z0 = S.z_opt;
        fprintf('Warm-starting collocation (N=%d) from COLLOCATION gait:\n  %s\n', ...
                N, warm_start_file);
    else                    % a PD/shooting result: seed by simulation
        [coeffs0, x_start0] = unpack_z(S.z_opt, p);
        z0 = col_seed_from_pd(coeffs0, x_start0, p);
        fprintf('Seeding collocation (N=%d) from PD gait:\n  %s\n', N, warm_start_file);
    end

    % ---- solve -----------------------------------------------------------
    [z_opt, fval, exitflag] = fmincon(@(z) col_cost(z, p), z0, ...
        [], [], [], [], lb, ub, @(z) col_constraints(z, p), options);

    % ---- report ----------------------------------------------------------
    [c, ceq] = col_constraints(z_opt, p);
    fprintf('\nexitflag=%d  cost=%.4g  max|ceq|=%.3e  max c=%.3e\n', ...
            exitflag, fval, max(abs(ceq)), max(c));

    % ---- save (mirrors the shooting driver's convention) -----------------
    results_dir = fullfile(root, 'Results');
    ts = datestr(now, 'yyyy-mm-dd_HH-MM-SS');
    result_file = fullfile(results_dir, ['col_result_', ts, '.mat']);
    save(result_file, 'z_opt', 'p', 'fval', 'exitflag', 'N');
    fprintf('Saved: %s\n', result_file);
end
