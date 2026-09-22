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

rungs = [180 170 160 150 142 135 128 120];
for u = rungs
    f = fullfile(SPD, sprintf('gait_u%03.0f.mat', u));
    if exist(f, 'file')
        L = load(f); z = L.z; p = L.p;
        report(logf, p, z, sprintf('u=%d (stored)', u), NaN);
        continue;
    end
    p.limits.u_max = u;
    t0 = tic;
    [z_try, out] = ch3_col_solve(p, z, struct('MaxFunctionEvaluations', 6e5, ...
                                              'MaxIterations', 600));
    secs = toc(t0);
    logln(logf, 'u_max %3d: exitflag %d | max|c| %.2e | max|ceq| %.2e | %.0f s', ...
          u, out.exitflag, out.max_c, out.max_ceq, secs);
    V = report(logf, p, z_try, sprintf('u=%d', u), secs);
    if out.max_c > 1e-6 || out.max_ceq > 1e-6 || ~V.ok
        logln(logf, ['  rung %d Nm did not land (max|c| %.2e, max|ceq| %.2e, verify %d). ' ...
                     'The march stops here; the last good gait is the rung above.'], ...
              u, out.max_c, out.max_ceq, V.ok);
        break;
    end
    z = z_try;
    save(f, 'z', 'p', 'out');
end
logln(logf, '=== TORQUE_MARCH_DONE %s', datestr(now));
fprintf('TORQUE_MARCH_DONE\n');
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
