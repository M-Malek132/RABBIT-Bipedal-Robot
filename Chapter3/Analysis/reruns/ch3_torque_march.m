function ch3_torque_march(seed_file, pin_speed)
%CH3_TORQUE_MARCH  Bring a verified gait down to the declared 120 Nm box.
%
%   ch3_torque_march()                        from the 1.2 m/s gait, speed free
%   ch3_torque_march('Results/ch3_gait_posture_195.mat')
%   ch3_torque_march([], true)                keep the speed pinned while cutting
%
% WHY THIS IS NOT ch3_realizability_march. That march does exactly this job --
% 191.4 -> 120 Nm in rungs 170, 150, 142, 135, 120 -- but it ran from
% ch3_gait_fix on the 30 kg model, and both its seed and its output
% (ch3_gait_full_constrained) are in the Chapter 3 report's stale table: under
% today's 74 kg dynamics the output reads peak |u| 287 against a 120 box, mu
% 6.94, min Fz -94 N and a verify deviation of 7.9. The ladder is sound; the
% gaits it was built on are gone. This redoes it from a gait that verifies now.
%
% WHAT IS BEING ASKED. Every gait that verifies on today's model rides exactly
% 195 Nm -- posture_195 and every rung of the speed ladder down to 1.15 m/s --
% because 195 is the box they were solved with and the cost spends whatever it
% is given. The declared limit in ch3_params is 120 Nm, so this is a 38% cut,
% and sec:box120 of the Chapter 4 report shows what it costs to skip it: under
% a real 120 Nm box every law that bounds the APPLIED torque dies within three
% steps, because the reference gait's own feedforward is 1.63x the limit before
% any controller is involved.
%
% SPEED IS LET GO BY DEFAULT. Peak torque falls with walking speed, so a 38%
% cut has to be paid for somewhere: a slower gait, a shorter step, or both.
% Pinning NEC1 while cutting the box asks for both at once, which is the kind
% of jump that stalls here. The default therefore switches NEC1 OFF and lets
% the speed land where it will -- limits first, speed last, as the rest of this
% pipeline does -- and the achieved speed is reported at every rung so the
% trade is visible. pin_speed = true keeps it fixed for comparison.
%
% Rungs: 195 -> 180 -> 170 -> 160 -> 150 -> 142 -> 135 -> 128 -> 120, with the
% 142 rung kept from ch3_realizability_march, where the 150 -> 135 step stalled
% on StepTolerance and an intermediate rung fixed it.
%
% Each rung is saved, audited and verified; a rung that fails to verify stops
% the march, because a gait that is not a trajectory is not a result.
%
% Output: Results/reruns/torque_march/ (one .mat per rung, torque_march.log).
%
% See also CH3_REALIZABILITY_MARCH, CH3_SPEED_LADDER, CH3_COL_CHECK_LIMITS.

if nargin < 1 || isempty(seed_file)
    seed_file = fullfile('Results', 'reruns', 'speed_ladder', 'gait_v1200.mat');
end
if nargin < 2 || isempty(pin_speed), pin_speed = false; end

ROOT = fileparts(fileparts(fileparts(fileparts(mfilename('fullpath')))));
SPD  = fullfile(ROOT, 'Results', 'reruns', 'torque_march');
if ~exist(SPD, 'dir'), mkdir(SPD); end
logf = fullfile(SPD, 'torque_march.log');

S = load(fullfile(ROOT, seed_file));
z = S.z;
p = ch3_upgrade_params(S.p);
p.scale_problem = true;
p.limits.enable.torque = true;

E0 = ch3_col_eval(z, p);
v0 = E0.L_step / E0.T;
if pin_speed
    p.enforce_nec1 = true;  p.v_des = v0;
else
    p.enforce_nec1 = false;
end

logln(logf, '=== torque march | seed %s | speed %s | %s', seed_file, ...
      ternary(pin_speed, sprintf('pinned at %.4f', v0), 'free'), datestr(now));
report(logf, p, z, 'seed', NaN);

% The ladder is a starting list, not a fixed plan: a rung that misses is
% retried halfway between it and the last one that landed, down to MIN_STEP.
% Measured on the first run, 150 missed only on max|ceq| = 5.8e-6 against a
% 1e-6 bar while verifying at 6.5e-5 -- a near miss that a smaller step is far
% more likely to clear than abandoning the march.
MIN_STEP  = 2;
u_target  = 120;
u_done    = max([abs(E0.u(:)); abs(E0.um(:))]);   % the seed's own peak
rungs     = [180 170 160 150 142 135 128 120];
for u = rungs
    f = fullfile(SPD, sprintf('gait_u%03.0f.mat', u));
    if exist(f, 'file')
        L = load(f);
        % ONLY A LANDED RUNG IS A RUNG. Rungs are saved whether or not they
        % land, so that a long solve is never thrown away -- which means the
        % resume path has to check. The first run left a 160 Nm file with a
        % 2.9e-03 torque row and a 150 Nm file 5.8e-06 off on the equalities;
        % reloading either as a starting point would march on from a gait that
        % was never feasible.
        if isfield(L, 'landed') && L.landed
            z = L.z_try; p = L.p; u_done = u;
            report(logf, p, z, sprintf('u=%d (stored)', u), NaN);
            continue;
        end
        logln(logf, '  stored rung %d Nm did not land; re-solving it', u);
    end
    p.limits.u_max = u;
    t0 = tic;
    [z_try, out] = ch3_col_solve(p, z, struct('MaxFunctionEvaluations', 6e5, ...
                                              'MaxIterations', 600));
    secs = toc(t0);
    logln(logf, 'u_max %3d: exitflag %d | max|c| %.2e | max|ceq| %.2e | %.0f s', ...
          u, out.exitflag, out.max_c, out.max_ceq, secs);
    V = report(logf, p, z_try, sprintf('u=%d', u), secs);

    % POLISH ONLY A NEAR MISS. A second warm solve clears the residual that the
    % evaluation limit leaves behind -- the 180 Nm rung went 5.9e-06 -> 2.8e-06
    % that way. But polishing from a point that is properly infeasible is
    % dangerous: at 160 Nm, with max|c| already 2.9e-03, the polish ran away to
    % a SPURIOUS solution (peak 178.5 against its own 160 box, mu 1.327, nodes
    % 0.4 off a true rollout) and cost 272 s. The guard below keeps the
    % pre-polish gait either way, but there is no reason to start that fire.
    if out.max_c > 1e-6 && out.max_c < 1e-4 && V.ok
        logln(logf, '  polishing rung %d (max|c| %.2e)', u, out.max_c);
        t1 = tic;
        [z_p, out_p] = ch3_col_solve(p, z_try, struct('MaxFunctionEvaluations', 6e5, ...
                                                      'MaxIterations', 600));
        V_p = report(logf, p, z_p, sprintf('u=%d polished', u), toc(t1));
        if V_p.ok && out_p.max_c < out.max_c
            z_try = z_p; out = out_p; V = V_p;
        end
    end

    % SAVE EVERY RUNG, LANDED OR NOT. An earlier version saved only after the
    % verdict and lost 75 minutes of solve when the verdict went against it.
    landed = V.ok && out.max_ceq <= 1e-6 && out.max_c <= 1e-4;
    save(f, 'z_try', 'p', 'out', 'V', 'landed');

    if ~landed
        % BISECT RATHER THAN STOP, and rather than carry on regardless: an
        % earlier version computed this flag and then only broke on verify, so
        % the 160 Nm rung was carried forward with a 2.9e-03 torque row and the
        % next rung was solved from it.
        logln(logf, '  rung %d Nm missed (max|c| %.2e, max|ceq| %.2e, verify %d)', ...
              u, out.max_c, out.max_ceq, V.ok);
        [z, u_done, ok] = bisect(z, p, u_done, u, MIN_STEP, logf, SPD);
        if ~ok
            logln(logf, ['  march stops at %.0f Nm: the step is below %d Nm and the ' ...
                         'box still will not come down'], u_done, MIN_STEP);
            break;
        end
        if u_done <= u, continue; end     % bisection got past this rung
        break;
    end
    z = z_try;
    u_done = u;
end
logln(logf, '=== TORQUE_MARCH_DONE %s', datestr(now));
fprintf('TORQUE_MARCH_DONE\n');
end

% ---------------------------------------------------------------------------
function [z, u_done, ok] = bisect(z, p, u_done, u_missed, MIN_STEP, logf, SPD)
%BISECT  Halve the box step until a rung lands or the step gets too small.
step = (u_done - u_missed) / 2;
ok   = false;
while step >= MIN_STEP
    u_try = u_done - step;
    logln(logf, '  bisecting: trying %.0f Nm (step %.0f)', u_try, step);
    p.limits.u_max = u_try;
    t0 = tic;
    [z_try, out] = ch3_col_solve(p, z, struct('MaxFunctionEvaluations', 6e5, ...
                                              'MaxIterations', 600));
    V = report(logf, p, z_try, sprintf('u=%.0f', u_try), toc(t0));
    landed = V.ok && out.max_ceq <= 1e-6 && out.max_c <= 1e-4;
    save(fullfile(SPD, sprintf('gait_u%03.0f.mat', u_try)), 'z_try', 'p', 'out', 'V', 'landed');
    if landed
        z = z_try; u_done = u_try; ok = true;
        step = min(step * 1.5, u_done - u_missed);
        if u_done <= u_missed, return; end
    else
        step = step / 2;
    end
end
end

% ---------------------------------------------------------------------------
function V = report(logf, p, z, tag, secs)
E   = ch3_col_eval(z, p);
c   = ch3_col_constraints(z, p);
chk = ch3_col_check_limits(z, p);
V   = ch3_col_verify(z, p, false);
lam = [E.lam, E.lamm];
pk  = max([abs(E.u(:)); abs(E.um(:))]);
logln(logf, ['%-16s v %.4f | T %.4f | L %.4f | peak|u| %6.1f of %3.0f | ' ...
             'impulse %5.2f Ns | mu %.3f | minFz %6.1f | limits ok %d (max|c| %.2e) | ' ...
             'verify %d (%.1e) | %.0f s'], ...
      tag, E.L_step/E.T, E.T, E.L_step, pk, p.limits.u_max, norm(E.impulse), ...
      max(abs(lam(1,:))./lam(2,:)), min(lam(2,:)), chk.ok, chk.max_c, ...
      V.ok, V.max_dev, secs);
if c(4) > 1e-6
    logln(logf, '      torque row still violated: %+.3e', c(4));
end
end

function s = ternary(b, x, y)
if b, s = x; else, s = y; end
end

function logln(f, fmt, varargin)
s = sprintf(fmt, varargin{:});
fprintf('%s\n', s);
fid = fopen(f, 'a'); fprintf(fid, '%s\n', s); fclose(fid);
end
