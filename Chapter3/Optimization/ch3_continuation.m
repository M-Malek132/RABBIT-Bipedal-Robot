function [z, hist] = ch3_continuation(p, v_targets, iters_per_step, z0, logfile)
%CH3_CONTINUATION  March the walking speed, warm-starting each solve.
%
%   [z, hist] = ch3_continuation(p, v_targets)
%   [z, hist] = ch3_continuation(p, v_targets, iters_per_step, z0, logfile)
%
% Solves the gait at each speed in v_targets in turn, seeding every solve from
% the previous converged result.
%
% WHY CONTINUATION IS NOT OPTIONAL HERE.  The seed rollout walks at about
% 0.12 m/s. Asking a cold solve for 0.35 m/s makes NEC1 (L_step/T = v_des) the
% dominant residual from iteration one, and the optimizer spends its effort
% fighting that single constraint instead of shaping a gait -- observed
% behaviour was max|ceq| oscillating 2.76 -> 0.10 -> 1.05 while the cost fell
% by half, i.e. trading feasibility for objective and never settling.
%
% Marching the speed in small steps keeps every individual solve NEARLY
% FEASIBLE at its start, which is the regime SQP is good at. Each step moves
% the target by an amount the previous gait can absorb.
%
% Start v_targets at (or just above) the seed's own measured speed. Query it
% with:
%     E = ch3_col_eval(ch3_col_seed(p), p);  v_seed = E.L_step / E.T;
%
% Inputs
%   p              : parameter struct (p.v_des is overwritten per step)
%   v_targets      : vector of speeds to solve, in order
%   iters_per_step : fmincon iterations per speed (default p.max_iter)
%   z0             : optional starting decision vector (default: fresh seed)
%   logfile        : optional path; progress is appended and flushed after
%                    every step (ch3_logln), so a long march run under -batch
%                    can be watched instead of staying silent until it exits
%
% Outputs
%   z    : decision vector at the final speed reached
%   hist : struct array, one per speed, with .v_des .z .fval .exitflag
%          .max_ceq .verify_dev .verify_ok .speed
%
% A speed whose solve fails to verify stops the march: continuing from a
% solution that is not a real trajectory only propagates the problem.
%
% See also CH3_COL_SOLVE, CH3_COL_VERIFY, CH3_COL_REMESH, CH3_COL_BUDGET,
%          CH3_LOGLN.

p = ch3_upgrade_params(p);

if nargin < 3 || isempty(iters_per_step), iters_per_step = p.max_iter; end
if nargin < 4 || isempty(z0), z0 = ch3_col_seed(p); end
if nargin < 5, logfile = ''; end

hist = struct('v_des', {}, 'z', {}, 'fval', {}, 'exitflag', {}, ...
              'max_ceq', {}, 'verify_dev', {}, 'verify_ok', {}, 'speed', {});

z = z0;

for k = 1:numel(v_targets)
    pk = p;
    pk.v_des = v_targets(k);

    % Size the evaluation cap to the iteration budget, taking the mesh from z
    % rather than pk.N_nodes: a caller-supplied z0 may be finer than the p it
    % arrived with. Without this the default 3e5 caps a 61-node step at ~170
    % iterations however large iters_per_step is -- see ch3_col_budget.
    pk = ch3_col_budget(pk, iters_per_step, z);

    ch3_logln(logfile, sprintf('===== continuation step %d/%d : v_des = %.4f m/s =====', ...
                               k, numel(v_targets), pk.v_des));

    [z_new, out] = ch3_col_solve(pk, z);

    E = ch3_col_eval(z_new, pk);
    V = ch3_col_verify(z_new, pk, false);

    hist(k) = struct('v_des', pk.v_des, 'z', z_new, 'fval', out.fval, ...
                     'exitflag', out.exitflag, 'max_ceq', out.max_ceq, ...
                     'verify_dev', V.max_dev, 'verify_ok', V.ok, ...
                     'speed', E.L_step / E.T); %#ok<AGROW>

    ch3_logln(logfile, sprintf('  -> J = %.2f, max|ceq| = %.2e, speed = %.4f, verify = %.2e (ok=%d)', ...
                               out.fval, out.max_ceq, E.L_step/E.T, V.max_dev, V.ok));

    if ~V.ok
        ch3_logln(logfile, ['  STOPPING: this step did not verify as a real trajectory.' ...
                 newline '  Continuing from it would propagate a fictional gait. Refine' ...
                 newline '  the mesh (ch3_col_remesh) or take smaller speed steps.']);
        return;
    end

    z = z_new;                 % warm start the next speed
end

end
