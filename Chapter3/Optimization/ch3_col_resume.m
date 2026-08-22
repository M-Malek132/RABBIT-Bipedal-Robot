function [z, out] = ch3_col_resume(ckpt_file, n_iters, p_override)
%CH3_COL_RESUME  Run a bounded chunk of a solve, resuming from a checkpoint.
%
%   [z, out] = ch3_col_resume(ckpt_file, n_iters)
%   [z, out] = ch3_col_resume(ckpt_file, n_iters, p_override)
%
% Runs at most n_iters fmincon iterations starting from the checkpoint (or
% from a fresh seed if the checkpoint does not exist yet), then writes the
% checkpoint back and returns.  Call it repeatedly -- each call in its own
% MATLAB process -- to complete a long solve in bounded chunks.
%
% WHY THIS EXISTS.  This MATLAB installation aborts mid-run from a heap
% corruption inside its own Add-On registry:
%
%     abort() <- malloc_zone_error <- libcurl <- libmwaddons_registry_core
%     Current Thread: 'OrderedWorkQ'
%
% That is a background thread doing network fetches for add-on metadata. It is
% unrelated to this code, it takes the whole process with it, and -nojvm does
% not prevent it because the registry is native rather than Java. Observed
% three times, after 5 to 10 minutes of solving.
%
% Rather than hope a long solve survives, chunk it: each process does a few
% iterations and persists its state, so a crash costs one chunk and the next
% process picks up exactly where it left off. Progress becomes monotone
% regardless of the platform's stability.
%
% WHAT CHUNKING COSTS.  Only z is carried across a chunk boundary, not
% fmincon's internal state -- the BFGS Hessian approximation and the merit
% function's penalty parameters are rebuilt from scratch each time. SQP is
% therefore slower per iteration count than one long run would be, and very
% small chunks are actively wasteful. Prefer the largest chunk that reliably
% finishes: on this machine crashes were seen after 5-10 minutes, so chunks of
% roughly 2-4 minutes are a reasonable compromise.
%
% Inputs
%   ckpt_file  : path to the checkpoint .mat (created if absent)
%   n_iters    : fmincon iterations to run in this chunk
%   p_override : optional struct of parameter fields to set (applied only when
%                seeding fresh; a resumed checkpoint keeps its own p)
%
% Outputs
%   z   : decision vector after this chunk
%   out : ch3_col_solve output for this chunk
%
% See also CH3_COL_SOLVE, CH3_COL_VERIFY, CH3_COL_BUDGET.

if exist(ckpt_file, 'file')
    S = load(ckpt_file);
    p = ch3_upgrade_params(S.p);
    z0 = S.z;
    done = 0;
    if isfield(S, 'iters_done'), done = S.iters_done; end
    fprintf('[resume] checkpoint found: %d iterations done so far\n', done);
else
    p = ch3_params();
    if nargin >= 3 && ~isempty(p_override)
        fn = fieldnames(p_override);
        for i = 1:numel(fn), p.(fn{i}) = p_override.(fn{i}); end
    end
    z0 = ch3_col_seed(p);
    done = 0;
    fprintf('[resume] no checkpoint; starting from a fresh seed\n');
end

p.checkpoint_file = ckpt_file;

% Size the evaluation cap to this chunk's iteration count, taking the mesh from
% z0 rather than p.N_nodes -- a checkpoint may have been written after a
% remesh. Chunks are usually small enough that the default 3e5 would not bind,
% but a large chunk on a fine mesh would silently stop short of n_iters and
% look like the solve converging. See ch3_col_budget.
p = ch3_col_budget(p, n_iters, z0);

[z, out] = ch3_col_solve(p, z0);

iters_done = done + out.output.iterations; %#ok<NASGU>
save(ckpt_file, 'z', 'p', 'iters_done');

V = ch3_col_verify(z, p, false);
fprintf('[resume] chunk done: J = %.4f, max|ceq| = %.3e, verify = %.3e (ok=%d), total iters = %d\n', ...
        out.fval, out.max_ceq, V.max_dev, V.ok, done + out.output.iterations);

end
