function [z_opt, out] = ch6_lib_solve(p, L_target, z0, iters)
%CH6_LIB_SOLVE  One library gait: the Chapter-3 solve with step length pinned.
%
%   [z_opt, out] = ch6_lib_solve(p, L_target, z0, iters)
%
% Section 6.4 needs gaits at PRESCRIBED step lengths -- the thesis uses
% {0.08, 0.24, 0.40, 0.56, 0.72} m for MARLO -- because (6.22) interpolates in
% step length and needs the grid to be the thing it interpolates over.
%
% Chapter 3's transcription has no such constraint. It has a step-length FLOOR
% (p.step_len_min, which only stops the degenerate step-in-place) and NEC1, an
% average-RATE equality L_step/T = v_des, which pins a ratio rather than a
% length. So this adds one equality row,
%
%       L_step(z) = L_target,
%
% on top of ch3_col_constraints and changes nothing else. It is deliberately a
% thin wrapper rather than a new transcription: the objective, the bounds, the
% Hermite-Simpson defects, the periodicity block, every Table 3.1 limit and the
% mesh verification are all Chapter 3's, so a library gait is a Chapter-3 gait
% and can be reported by ch3_report without translation.
%
% ------------------------------------------------- why the extra row is cheap
% ch3_col_eval caches on z with a persistent, and fmincon calls the objective
% and the constraints back to back with the same z. So the ch3_col_eval call
% below is a cache hit on the one ch3_col_constraints already did -- the row
% costs a subtraction, not a second evaluation of the mesh.
%
% ------------------------------------------------------- NEC1 is turned OFF here
% Pinning L_step AND L_step/T simultaneously pins T as well, which over-
% constrains a warm start that is only close in shape. Chapter 3 already
% documents the failure mode: with the gait far from periodic, an extra
% equality becomes the dominant residual and the solve trades feasibility for
% objective without settling. The library only needs the LENGTHS to be right;
% each gait's duration is whatever it wants to be, and the resulting speed is
% recorded rather than imposed.
%
% Inputs
%   p        : parameter struct
%   L_target : step length to solve for [m]
%   z0       : warm start (required in practice -- a cold solve at a pinned
%              length is the same over-constrained problem NEC1 was)
%   iters    : fmincon iterations (default p.lib.iters)
%
% Outputs
%   z_opt : decision vector
%   out   : struct .fval .exitflag .alpha .T .L_step .speed .max_ceq .max_c
%           .verify .wall_time
%
% See also CH6_LIB_BUILD, CH3_COL_SOLVE, CH3_COL_CONSTRAINTS, CH3_COL_EVAL.

if nargin < 4 || isempty(iters), iters = p.lib.iters; end

p = ch3_upgrade_params(p);
p.enforce_nec1 = false;

if nargin < 3 || isempty(z0)
    z0 = ch3_col_seed(p);
end

% ============================================ WHY T IS BOUNDED NEAR THE SEED'S
% Chapter 3's cost is int||u||^2 PER UNIT STEP LENGTH, so asking for a SHORTER
% step raises the objective and the optimizer looks for somewhere else to put
% the loss. The one direction that is free is the step DURATION: nothing pins
% it once NEC1 is off, and a slower step costs less torque per metre.
%
% MEASURED, marching down from the 0.353 m seed with T free in Chapter 3's
% default [0.20, 1.50] s: the L* = 0.34 m solve walked T from 0.301 s to 1.431 s
% -- the upper bound -- ending at 0.239 m/s with the mesh check off by 19.2.
% It satisfied the step-length equality to 2e-03 and was not a trajectory at
% all. A library built from that is a library of fictions.
%
% Bounding T to a band around the seed's removes the escape route. It is a
% restriction on the library, not a fix to the transcription, and it is the
% right one here: (6.22) interpolates gaits that are supposed to be neighbours,
% and two gaits whose durations differ by 5x are not neighbours in any sense
% the interpolation can exploit.
if isfield(p, 'lib') && isfield(p.lib, 'T_band') && ~isempty(p.lib.T_band)
    [~, T0] = ch3_col_unpack(z0, p);
    p.T_min = max(p.T_min, T0 * (1 - p.lib.T_band));
    p.T_max = min(p.T_max, T0 * (1 + p.lib.T_band));
end

N = (numel(z0) - 1 - p.ny*p.n_ctrl) / p.nx;
[lb, ub] = ch3_col_bounds(p, N);

options = optimoptions('fmincon', ...
    'Algorithm',              'sqp', ...
    'Display',                'iter', ...
    'MaxIterations',          iters, ...
    'MaxFunctionEvaluations', p.max_fun_evals, ...
    'OptimalityTolerance',    1e-6, ...
    'ConstraintTolerance',    1e-6, ...
    'StepTolerance',          1e-10, ...
    'FiniteDifferenceType',   'central', ...
    'ScaleProblem',           false);

t0 = tic;
[z_opt, fval, exitflag] = fmincon( ...
    @(z) ch3_col_cost(z, p), z0, [], [], [], [], lb, ub, ...
    @(z) nonlcon(z, p, L_target), options);
wall = toc(t0);

[c, ceq] = nonlcon(z_opt, p, L_target);
[~, T, alpha] = ch3_col_unpack(z_opt, p);
E = ch3_col_eval(z_opt, p);

out = struct('fval', fval, 'exitflag', exitflag, 'alpha', alpha, ...
             'T', T, 'L_step', E.L_step, 'speed', E.L_step / T, ...
             'L_target', L_target, ...
             'max_ceq', max(abs(ceq)), 'max_c', max(c), ...
             'wall_time', wall, ...
             'verify', ch3_col_verify(z_opt, p, false));

fprintf(['[ch6_lib_solve] L* = %.4f -> L = %.4f (err %.2e), T = %.3f s, ' ...
         'v = %.3f m/s, J = %.4f, exit %d, %.0f s\n'], ...
        L_target, out.L_step, abs(out.L_step - L_target), T, out.speed, ...
        fval, exitflag, wall);

if ~out.verify.ok
    fprintf(['[ch6_lib_solve] *** mesh check FAILED (dev %.2e, tol %.1e). ' ...
             'This gait is not a real trajectory; do not put it in the ' ...
             'library.\n'], out.verify.max_dev, p.verify_tol);
end

end

% ---------------------------------------------------------------------------
function [c, ceq] = nonlcon(z, p, L_target)
[c, ceq] = ch3_col_constraints(z, p);
E = ch3_col_eval(z, p);            % cache hit on the call inside the line above
ceq = [ceq; E.L_step - L_target];
end
