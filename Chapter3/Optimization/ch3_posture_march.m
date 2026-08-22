function [z, hist] = ch3_posture_march(p, stages, iters_per_stage, z0, logfile)
%CH3_POSTURE_MARCH  March the torso-pitch box and hip-height band.
%
%   [z, hist] = ch3_posture_march(p, stages)
%   [z, hist] = ch3_posture_march(p, stages, iters_per_stage, z0, logfile)
%
% Solves the gait at each posture target in turn, seeding every solve from the
% previous converged result.  ch3_continuation does this for walking SPEED;
% this does it for the two POSTURE requirements, which need it for the same
% reason and one extra one.
%
% WHY MARCHING IS NOT OPTIONAL.  Both requirements are box/band constraints on
% coordinates the cost does not otherwise price, and the optimizer sits hard
% against them (see the qt_range and hip_h notes in ch3_params). Moving a bound
% the gait is resting on moves EVERY node at once: the pose at every node, the
% Bezier coefficients that generate it, and the step timing all have to shift
% together. From a converged start a small move is a nearly-feasible problem,
% which is the regime SQP handles well; a large one is a new problem with a bad
% initial guess.
%
% THE EXTRA REASON, specific to torso pitch: qt enters the problem only through
% theta = qt + q1 + q2/2, so any change in qt has to be paid for by an equal
% and opposite change in q1 + q2/2 to keep the phase clock reading the same at
% the start and end of the step. Swinging qt across zero in one move asks the
% leg angles to absorb the whole 0.4 rad at once, and the periodicity equality
% is what breaks first.
%
% Inputs
%   p               : parameter struct (its posture fields are overwritten
%                     per stage)
%   stages          : struct array; each element may set any of
%                       .qt_range   [lo hi]  torso pitch box        [rad]
%                       .hip_h      scalar   centre of height band  [m]
%                       .hip_h_tol  scalar   half-width of the band [m]
%                       .height     logical  enable/disable the band
%                     Fields a stage omits keep the value from the stage
%                     before it, so only what changes needs writing.
%   iters_per_stage : fmincon iterations per stage (default p.max_iter)
%   z0              : starting decision vector (default: fresh seed)
%   logfile         : optional path; progress is appended and flushed after
%                     every stage (ch3_logln), so a long march run under
%                     -batch can be watched instead of staying silent until it
%                     exits. ch3_lean_tall_march passes its own log down here,
%                     so a nested march interleaves into one file.
%
% Outputs
%   z    : decision vector at the final stage reached
%   hist : struct array, one per stage, with .qt_range .hip_h .hip_h_tol
%          .height .z .fval .exitflag .max_ceq .max_c .qt_lo .qt_hi
%          .hip_lo .hip_hi .verify_dev .verify_ok
%
% A stage whose solve fails to verify stops the march, exactly as in
% ch3_continuation: continuing from a solution that is not a real trajectory
% only propagates the problem.
%
% See also CH3_CONTINUATION, CH3_COL_SOLVE, CH3_COL_VERIFY, CH3_COL_BUDGET,
%          CH3_LOGLN, CH3_PARAMS.

p = ch3_upgrade_params(p);

if nargin < 3 || isempty(iters_per_stage), iters_per_stage = p.max_iter; end
if nargin < 4 || isempty(z0), z0 = ch3_col_seed(p); end
if nargin < 5, logfile = ''; end

hist = struct('qt_range', {}, 'hip_h', {}, 'hip_h_tol', {}, 'height', {}, ...
              'z', {}, 'fval', {}, 'exitflag', {}, 'max_ceq', {}, 'max_c', {}, ...
              'qt_lo', {}, 'qt_hi', {}, 'hip_lo', {}, 'hip_hi', {}, ...
              'verify_dev', {}, 'verify_ok', {});

z  = z0;
pk = p;

for k = 1:numel(stages)
    s = stages(k);

    % Only what the stage names changes; everything else carries forward.
    if isfield(s, 'qt_range')  && ~isempty(s.qt_range),  pk.qt_range = s.qt_range;  end
    if isfield(s, 'hip_h')     && ~isempty(s.hip_h),     pk.limits.hip_h = s.hip_h; end
    if isfield(s, 'hip_h_tol') && ~isempty(s.hip_h_tol), pk.limits.hip_h_tol = s.hip_h_tol; end
    if isfield(s, 'height')    && ~isempty(s.height),    pk.limits.enable.height = logical(s.height); end

    % Size the evaluation cap to the iteration budget, taking the mesh from z
    % rather than pk.N_nodes: a caller-supplied z0 may be finer than the p it
    % arrived with. Without this the default 3e5 caps a 61-node stage at ~170
    % iterations however large iters_per_stage is -- see ch3_col_budget.
    pk = ch3_col_budget(pk, iters_per_stage, z);

    ch3_logln(logfile, sprintf('===== posture stage %d/%d =====', k, numel(stages)));

    % Move the gait INTO the new pitch box before solving, rather than letting
    % fmincon clamp it there. Clamping shifts theta at every node and detonates
    % the defects (1.25e+03 measured, against 1.97e-01 for the same move done
    % as a repose); ch3_repose rotates the torso with the legs held fixed, so
    % the phase clock, both feet and the zero dynamics surface all survive
    % exactly. Set stage.repose = false to opt out.
    qt_now = ch3_col_eval(z, pk).X(3,:);
    want_repose = ~isfield(s,'repose') || isempty(s.repose) || s.repose;
    if want_repose && (min(qt_now) < pk.qt_range(1) || max(qt_now) > pk.qt_range(2))
        % Aim the LOW end of the pitch just inside the box: the optimizer
        % rides the lower bound, so that is where the gait will settle anyway.
        delta = (pk.qt_range(1) + 0.01) - min(qt_now);
        [z, ri] = ch3_repose(z, pk, delta);
        ch3_logln(logfile, sprintf('  repose  %+.4f rad : qt [%+.3f %+.3f] -> [%+.3f %+.3f]', ...
                delta, ri.qt_before(1), ri.qt_before(2), ...
                ri.qt_after(1), ri.qt_after(2)));
        [~, ceq_r] = ch3_col_constraints(z, pk);
        ch3_logln(logfile, sprintf('          warm-start residual after repose: %.3e', max(abs(ceq_r))));
    end

    ch3_logln(logfile, sprintf('  qt box   [%+.3f %+.3f] rad  (%+.1f .. %+.1f deg)', ...
            pk.qt_range(1), pk.qt_range(2), ...
            rad2deg(pk.qt_range(1)), rad2deg(pk.qt_range(2))));
    if pk.limits.enable.height
        ch3_logln(logfile, sprintf('  hip band %.3f +- %.3f m  [ENABLED]', ...
                pk.limits.hip_h, pk.limits.hip_h_tol));
    else
        ch3_logln(logfile, sprintf('  hip band %.3f +- %.3f m  [off]', ...
                pk.limits.hip_h, pk.limits.hip_h_tol));
    end

    [z_new, out] = ch3_col_solve(pk, z);

    E = ch3_col_eval(z_new, pk);
    V = ch3_col_verify(z_new, pk, false);

    qt  = E.X(3,:);
    hip = -E.X(2,:);        % pz is DOWN-positive; world height = -pz

    hist(k) = struct('qt_range', pk.qt_range, ...
                     'hip_h', pk.limits.hip_h, ...
                     'hip_h_tol', pk.limits.hip_h_tol, ...
                     'height', pk.limits.enable.height, ...
                     'z', z_new, 'fval', out.fval, 'exitflag', out.exitflag, ...
                     'max_ceq', out.max_ceq, 'max_c', out.max_c, ...
                     'qt_lo', min(qt), 'qt_hi', max(qt), ...
                     'hip_lo', min(hip), 'hip_hi', max(hip), ...
                     'verify_dev', V.max_dev, 'verify_ok', V.ok);

    ch3_logln(logfile, sprintf(['  -> J = %.2f, max|ceq| = %.2e, max c = %.2e\n' ...
             '     qt  %+.4f .. %+.4f rad (%+.1f .. %+.1f deg)\n' ...
             '     hip %.4f .. %.4f m (bob %.4f)\n' ...
             '     verify %.2e (ok=%d)'], ...
            out.fval, out.max_ceq, out.max_c, ...
            min(qt), max(qt), rad2deg(min(qt)), rad2deg(max(qt)), ...
            min(hip), max(hip), max(hip)-min(hip), V.max_dev, V.ok));

    if ~V.ok
        ch3_logln(logfile, ['  STOPPING: this stage did not verify as a real trajectory.' ...
                 newline '  Continuing from it would propagate a fictional gait. Refine' ...
                 newline '  the mesh (ch3_col_remesh) or take smaller posture steps.']);
        return;
    end

    z = z_new;                 % warm start the next stage
end

end
