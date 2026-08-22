function [z, hist, p_out] = ch3_impact_march(p, ratios, iters_per_stage, z0, logfile)
%CH3_IMPACT_MARCH  March the NEC3 impulse friction cone down to mu_s.
%
%   [z, hist]        = ch3_impact_march(p, ratios)
%   [z, hist, p_out] = ch3_impact_march(p, ratios, iters_per_stage, z0, logfile)
%
% RETURN p_out WITH z, AND USE IT.  This function enables the impact gate on its
% OWN copy of p, so the caller's p still has enable.impact = false when the
% march returns. Feeding the returned z back into ch3_col_solve with the
% caller's p therefore re-solves with the constraint RELEASED, and the optimizer
% quietly relaxes the ratio back toward the unconstrained optimum -- which looks
% exactly like the constraint failing to hold rather than never having been
% applied. Observed: a mid-march remesh-and-retry using the caller's p sprang
% the ratio from 0.761 back to 1.134 in one solve. p_out carries enable.impact
% and the mu_s_impact actually reached, so warm starts stay constrained.
%
% NEC3 asks that the impact impulse lie inside the friction cone, |Ix| <= mu_s
% Iz.  The reference gait sits at |Ix|/Iz = 2.44 against mu_s = 0.4 -- outside
% by a factor of six -- so switching the gate on in one move hands fmincon a
% start point that is not merely infeasible but infeasible in the one direction
% the gait was never shaped for.  This marches p.limits.mu_s_impact down a
% ladder instead, warm-starting each stage from the last, exactly as
% ch3_continuation does for speed and ch3_posture_march for posture.
%
% WHY THIS CONSTRAINT IS HARDER THAN THE TABLE 3.1 ONES.  Torque, GRF and
% continuous friction are all evaluated ALONG the step, so the optimizer has
% every node to spread the correction over.  The impulse is a property of ONE
% state, x_N, filtered through the impact map -- and x_N is not free: the
% periodicity equality ties it to node 1, and the guard equality pins the swing
% foot to the ground there.  So the only way to change the impulse is to change
% the whole orbit, which is precisely why it needs the smallest steps of any
% limit in the pipeline.
%
% WHAT ACTUALLY HAS TO MOVE.  The impact impulse is I = Lambda(q) v_foot, with
% Lambda = (J_sw M^-1 J_sw')^-1 the operational-space inertia at the foot.  So a
% large |Ix|/Iz is not simply "the foot is moving forward fast" -- it is the
% product of the foot's pre-impact velocity and an anisotropic inertia that a
% near-straight leg makes very stiff along its own axis.  The gait can fix the
% ratio by striking with less horizontal foot speed, by striking with the leg in
% a posture whose Lambda is less anisotropic, or by some mix; the optimizer is
% free to choose, which is why no attempt is made here to steer it.
%
% Inputs
%   p               : parameter struct. p.limits.enable.impact is forced true;
%                     everything else is left exactly as given.
%   ratios          : vector of mu_s_impact targets, descending. The first
%                     should be at or just below what the gait already
%                     measures -- read it off ch3_report.
%   iters_per_stage : fmincon iterations per stage (default p.max_iter)
%   z0              : starting decision vector (default: fresh seed)
%   logfile         : optional path; progress is appended and flushed after
%                     every stage (ch3_logln), so a long march can be watched
%                     while it runs instead of only after it finishes.
%
% Outputs
%   z     : decision vector at the last stage that converged AND verified
%   hist  : struct array, one per stage, with .mu_target .mu_measured .z .fval
%           .exitflag .max_ceq .max_c .impulse .speed .L_step .T
%           .verify_dev .verify_ok
%   p_out : the parameter struct that produced z -- impact ENABLED, with
%           mu_s_impact at the last rung actually reached (not the last rung
%           requested, which may differ if the march stopped early). Pass this,
%           not the input p, to anything that re-solves from z.
%
% A stage that fails to verify stops the march, for the reason ch3_continuation
% gives: continuing from a solution that is not a real trajectory propagates a
% fictional gait rather than fixing it.
%
% See also CH3_CONTINUATION, CH3_POSTURE_MARCH, CH3_COL_SOLVE, CH3_COL_BUDGET,
%          CH3_LOGLN, CH3_REPORT.

p = ch3_upgrade_params(p);

if nargin < 3 || isempty(iters_per_stage), iters_per_stage = p.max_iter; end
if nargin < 4 || isempty(z0), z0 = ch3_col_seed(p); end
if nargin < 5, logfile = ''; end

p.limits.enable.impact = true;

hist = struct('mu_target', {}, 'mu_measured', {}, 'z', {}, 'fval', {}, ...
              'exitflag', {}, 'max_ceq', {}, 'max_c', {}, 'impulse', {}, ...
              'speed', {}, 'L_step', {}, 'T', {}, ...
              'verify_dev', {}, 'verify_ok', {});

z  = z0;
pk = p;

% Reflects the last rung actually REACHED, so an early stop still returns a p
% that matches the z beside it.
p_out = p;
p_out.limits.enable.impact = true;

E0 = ch3_col_eval(z, pk);
ch3_logln(logfile, sprintf(['===== ch3_impact_march: %d stages =====\n' ...
                            'start: |Ix|/Iz = %.4f  (I = [%.3f %.3f] Ns)\n' ...
                            'target ladder: %s'], ...
     numel(ratios), ratio_of(E0.impulse), E0.impulse(1), E0.impulse(2), ...
     num2str(ratios(:).', '%.3f ')));

for k = 1:numel(ratios)
    pk.limits.mu_s_impact = ratios(k);

    % Size the evaluation cap to the iteration budget, taking the mesh from z
    % rather than pk.N_nodes: a caller-supplied z0 may be finer than the p it
    % arrived with. Without this the default 3e5 caps a 61-node stage at ~170
    % iterations however large iters_per_stage is -- see ch3_col_budget.
    pk = ch3_col_budget(pk, iters_per_stage, z);

    ch3_logln(logfile, sprintf('===== stage %d/%d : mu_s_impact = %.4f =====', ...
                               k, numel(ratios), ratios(k)));

    [z_new, out] = ch3_col_solve(pk, z);

    E = ch3_col_eval(z_new, pk);
    V = ch3_col_verify(z_new, pk, false);
    r = ratio_of(E.impulse);

    hist(k) = struct('mu_target', ratios(k), 'mu_measured', r, ...
                     'z', z_new, 'fval', out.fval, 'exitflag', out.exitflag, ...
                     'max_ceq', out.max_ceq, 'max_c', out.max_c, ...
                     'impulse', E.impulse, ...
                     'speed', E.L_step/E.T, 'L_step', E.L_step, 'T', E.T, ...
                     'verify_dev', V.max_dev, 'verify_ok', V.ok);

    ch3_logln(logfile, sprintf([ ...
        '  -> J = %.2f, exitflag %d, max|ceq| = %.2e, max c = %.2e\n' ...
        '     |Ix|/Iz = %.4f  (target %.4f),  I = [%.3f %.3f] Ns\n' ...
        '     speed %.4f m/s, L_step %.4f m, T %.4f s\n' ...
        '     verify %.2e (ok=%d)'], ...
        out.fval, out.exitflag, out.max_ceq, out.max_c, ...
        r, ratios(k), E.impulse(1), E.impulse(2), ...
        E.L_step/E.T, E.L_step, E.T, V.max_dev, V.ok));

    if ~V.ok
        ch3_logln(logfile, sprintf([ ...
            '  STOPPING: stage %d did not verify as a real trajectory.\n' ...
            '  Refine the mesh (ch3_col_remesh) or take smaller ratio steps.'], k));
        return;
    end

    z = z_new;                 % warm start the next stage
    p_out.limits.mu_s_impact = ratios(k);
end

end

% ---------------------------------------------------------------------------
function r = ratio_of(I)
if I(2) > 1e-9, r = abs(I(1))/I(2); else, r = Inf; end
end
