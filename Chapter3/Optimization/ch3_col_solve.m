function [z_opt, out] = ch3_col_solve(p, z0, opts_override)
%CH3_COL_SOLVE  Stage 3: solve for the Bezier coefficients alpha.
%
%   [z_opt, out] = ch3_col_solve(p)
%   [z_opt, out] = ch3_col_solve(p, z0)
%   [z_opt, out] = ch3_col_solve(p, z0, opts_override)
%
% Runs fmincon on the direct-collocation transcription: minimize
% torque-squared per distance subject to the Hermite-Simpson defects,
% periodicity through Delta, and whichever of the Table 3.1 limits are
% enabled.  The output of this stage is one periodic orbit and the alpha that
% defines it.
%
% SOLVER CHOICES, AND WHY:
%   * SQP.  The transcription is mostly equality constrained with a smooth
%     objective, which is where SQP is strongest; interior-point spends its
%     effort on a barrier that has little to do here.
%   * Central finite differences.  The constraint values pass through a KKT
%     solve and a matrix inversion, so forward differences at the default step
%     lose too many digits near a well-conditioned solution.
%   * ScaleProblem off, so the Feasibility number fmincon reports is directly
%     comparable to ConstraintTolerance rather than to a rescaled surrogate.
%
% WARM STARTING.  Pass z0 from a previous solve to enable a Table 3.1 limit
% incrementally.  The recommended order is GRF first (get Fz > 0), then
% friction, then torque, then impulse: friction is |Fx|/Fz, so enabling it
% while Fz still crosses zero sets the optimizer chasing a division by zero.
%
% Inputs
%   p             : parameter struct
%   z0            : optional warm start (default: ch3_col_seed)
%   opts_override : optional struct of optimoptions name/value overrides
%
% Outputs
%   z_opt : optimized decision vector
%   out   : struct .fval .exitflag .output .alpha .X .T .z0 .seed_info
%           .wall_time .max_ceq .max_c
%
% See also CH3_COL_SEED, CH3_COL_CONSTRAINTS, CH3_REPORT.

p = ch3_upgrade_params(p);   % a warm start may carry an older parameter struct

seed_info = [];
if nargin < 2 || isempty(z0)
    [z0, seed_info] = ch3_col_seed(p);
end

N = (numel(z0) - 1 - p.ny*p.n_ctrl) / p.nx;
[lb, ub] = ch3_col_bounds(p, N);

options = optimoptions('fmincon', ...
    'Algorithm',               'sqp', ...
    'Display',                 'iter', ...
    'MaxIterations',           p.max_iter, ...
    'MaxFunctionEvaluations',  p.max_fun_evals, ...
    'OptimalityTolerance',     1e-6, ...
    'ConstraintTolerance',     1e-6, ...
    'StepTolerance',           1e-10, ...
    'FiniteDifferenceType',    'central', ...
    'ScaleProblem',            false);

% CHECKPOINTING.  A solve here runs for minutes and the only artifact is the
% final z, so anything that kills the process -- including causes that have
% nothing to do with this code -- costs the entire run. (Observed: MATLAB's
% add-on registry crashed in a network call and took down a solve at iteration
% 76.) Writing z every few iterations makes a crash cost minutes, not a run,
% and the checkpoint doubles as a warm start.
if ~isempty(p.checkpoint_file)
    options = optimoptions(options, 'OutputFcn', ...
        @(z, ov, state) checkpoint_fcn(z, ov, state, p));
end

if nargin >= 3 && ~isempty(opts_override)
    fn = fieldnames(opts_override);
    for i = 1:numel(fn)
        options = optimoptions(options, fn{i}, opts_override.(fn{i}));
    end
end

t_start = tic;
[z_opt, fval, exitflag, output] = fmincon( ...
    @(z) ch3_col_cost(z, p), z0, [], [], [], [], lb, ub, ...
    @(z) ch3_col_constraints(z, p), options);
wall_time = toc(t_start);

[c, ceq] = ch3_col_constraints(z_opt, p);
[X, T, alpha] = ch3_col_unpack(z_opt, p);

out = struct('fval', fval, 'exitflag', exitflag, 'output', output, ...
             'alpha', alpha, 'X', X, 'T', T, ...
             'z0', z0, 'seed_info', seed_info, ...
             'wall_time', wall_time, ...
             'max_ceq', max(abs(ceq)), 'max_c', max(c));

fprintf('\n[ch3_col_solve] exitflag %d, J = %.4f, max|ceq| = %.3e, max c = %.3e, %.1f s\n', ...
        exitflag, fval, out.max_ceq, out.max_c, wall_time);

% A converged solve with tiny residuals can still be a spurious discrete
% solution, so check it here rather than leaving it to whoever reports later.
% Reporting "converged" without this is how a fictional gait gets believed.
out.verify = ch3_col_verify(z_opt, p, false);
if out.verify.ok
    fprintf('[ch3_col_solve] mesh check OK: nodes within %.2e of a true trajectory\n', ...
            out.verify.max_dev);
else
    fprintf(['[ch3_col_solve] *** WARNING: SPURIOUS DISCRETE SOLUTION ***\n' ...
             '   nodes deviate %.3e from a true rollout (tol %.1e).\n' ...
             '   Small defects do NOT imply a real trajectory -- the mesh is too\n' ...
             '   coarse. Increase p.N_nodes (currently %d) and re-solve; do not\n' ...
             '   trust any gait number from this result.\n'], ...
            out.verify.max_dev, p.verify_tol, N);
end

end

% ---------------------------------------------------------------------------
function stop = checkpoint_fcn(z, ov, state, p)
stop = false;
if ~strcmp(state, 'iter'), return; end
if mod(ov.iteration, p.checkpoint_every) ~= 0, return; end
ck = struct('z', z, 'iteration', ov.iteration, 'fval', ov.fval, ...
            'feasibility', ov.constrviolation, 'p', p); %#ok<NASGU>
save(p.checkpoint_file, '-struct', 'ck');
end
