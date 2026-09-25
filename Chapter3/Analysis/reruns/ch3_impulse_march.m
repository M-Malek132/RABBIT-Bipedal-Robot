function ch3_impulse_march(seed_file, tag)
%CH3_IMPULSE_MARCH  Bring the torque-march gait under the declared 15 Ns impulse.
%
%   ch3_impulse_march()                 from the 162.5 Nm gait (torque_march/gait_u162)
%   ch3_impulse_march(seed_file, tag)   from any stored gait; tag names its files
%
% WHY. ch3_torque_march brought the box from 195 to a verified 162.5 Nm, but
% with the impulse gate (row 7 of ch3_col_constraints, ||I|| <= impulse_max)
% OFF throughout: every rung lands with 18.8-19.1 Ns against Table 3.1's 15.
% Those gaits are therefore only half realizable -- torque within its box,
% impact outside its limit -- and the Chapter 3 report says so. This is the
% second march the handoff asks for: the same gait, the torque box held where
% the first march left it, and the impulse ceiling walked down to 15.
%
% WHAT IS HELD AND WHAT IS FREE. The box stays at the seed's own u_max, so a
% landed 15 Ns rung is a gait that meets BOTH limits at 162.5 Nm, which is the
% claim worth having. The speed stays free, as in the torque march: it ran
% backwards there (1.2815 -> 1.2921 m/s as the box tightened) and pinning it
% collapsed the solve, so no second demand is put on it here.
%
% WHY RUNGS AND NOT ONE SOLVE. The impulse is a property of ONE state, x_N,
% filtered through the impact map, and x_N is tied to node 1 by periodicity --
% so the only way to change it is to change the whole orbit (ch3_impact_march
% makes the same argument for the NEC3 cone). 1 Ns rungs are ~5% cuts, the
% size of the torque march's steps that landed; a rung that misses is bisected
% down to MIN_STEP before the march gives up.
%
% KEEP THE GATE ON THE p THAT IS SAVED. Each rung file stores the p it was
% solved with, impulse gate enabled and the ceiling it was solved at.
% Re-solving a stored z with a p whose gate is off lets the impulse spring back
% toward the unconstrained optimum -- which looks exactly like the constraint
% failing (observed on the NEC3 march: 0.761 -> 1.134 in one solve).
%
% Output: Results/reruns/impulse_march/ -- one .mat per rung (<tag>_I<100*cap>),
% landed or not, with the torque march's fields (z_try, p, out, V, landed), and
% impulse_march.log. A rerun resumes from the landed rungs.
%
% See also CH3_TORQUE_MARCH, CH3_IMPACT_MARCH, CH3_COL_CONSTRAINTS.

if nargin < 1 || isempty(seed_file)
    seed_file = fullfile('Results', 'reruns', 'torque_march', 'gait_u162.mat');
end
if nargin < 2 || isempty(tag), tag = 'u162'; end

ROOT = fileparts(fileparts(fileparts(fileparts(mfilename('fullpath')))));
SPD  = fullfile(ROOT, 'Results', 'reruns', 'impulse_march');
if ~exist(SPD, 'dir'), mkdir(SPD); end
logf = fullfile(SPD, 'impulse_march.log');

src = fullfile(ROOT, seed_file);
if ~exist(src, 'file')
    logln(logf, ['=== impulse march: seed %s not found. The Chapter 3 campaign ' ...
                 'results (Results/reruns/torque_march/, speed_ladder/) are ' ...
                 'git-ignored; run this on the machine that holds them, or copy ' ...
                 'that folder here first.'], seed_file);
    error('ch3_impulse_march:noSeed', 'Seed gait not found: %s', src);
end

S = load(src);
if isfield(S, 'z_try'), z = S.z_try; else, z = S.z; end
p = ch3_upgrade_params(S.p);
p.scale_problem = true;
p.enforce_nec1  = false;                 % speed free, as in the torque march
p.limits.enable.torque  = true;          % the box the seed was solved in stays
p.limits.enable.impulse = true;          % the row this march is about

target   = 15;                           % Table 3.1 (ch3_params)
MIN_STEP = 0.25;                         % [Ns]

E0 = ch3_col_eval(z, p);
I0 = norm(E0.impulse);
p.limits.impulse_max = I0;               % the seed, reported against its own value
logln(logf, '=== impulse march | seed %s | box %.1f Nm | seed impulse %.2f Ns | %s', ...
      seed_file, p.limits.u_max, I0, datestr(now));
report(logf, p, z, 'seed', NaN);
if I0 <= target
    logln(logf, '=== the seed is already within %.0f Ns; nothing to march', target);
    logln(logf, '=== IMPULSE_MARCH_DONE %s', datestr(now));
    fprintf('IMPULSE_MARCH_DONE\n');
    return;
end

% 19, 18, 17, 16, 15 from a 19.07 Ns seed: the first rung only switches the gate
% on nearly slack, which checks that the gate alone does not disturb the solve.
rungs  = floor(I0):-1:target;
if rungs(end) > target, rungs = [rungs, target]; end
I_done = I0;

for I = rungs
    f = rung_file(SPD, tag, I);
    if exist(f, 'file')
        L = load(f);
        if isfield(L, 'landed') && L.landed
            z = L.z_try; p = L.p; I_done = I;
            report(logf, p, z, sprintf('I<=%.2f (stored)', I), NaN);
            continue;
        end
        logln(logf, '  stored rung %.2f Ns did not land; re-solving it', I);
    end

    [z_try, p_rung, out, V, landed] = solve_rung(p, z, I, logf);
    save_rung(f, z_try, p_rung, out, V, landed);

    if ~landed
        logln(logf, '  rung %.2f Ns missed (max|c| %.2e, max|ceq| %.2e, verify %d)', ...
              I, out.max_c, out.max_ceq, V.ok);
        [z, p, I_done, ok] = bisect(z, p, I_done, I, MIN_STEP, logf, SPD, tag);
        if ~ok
            logln(logf, ['  march stops at %.2f Ns: the step is below %.2f Ns and ' ...
                         'the impulse will not come down at a %.1f Nm box'], ...
                  I_done, MIN_STEP, p.limits.u_max);
            break;
        end
        if I_done <= I, continue; end    % bisection got past this rung
        break;
    end
    z = z_try; p = p_rung; I_done = I;
end

if I_done <= target + 1e-9
    logln(logf, ['=== VERDICT: a verified gait meets both the %.1f Nm box and the ' ...
                 '%.0f Ns impulse limit: %s'], p.limits.u_max, target, ...
          rung_file(SPD, tag, I_done));
else
    logln(logf, ['=== VERDICT: at a %.1f Nm box the impulse comes down only to ' ...
                 '%.2f Ns, above the declared %.0f'], p.limits.u_max, I_done, target);
end
logln(logf, '=== IMPULSE_MARCH_DONE %s', datestr(now));
fprintf('IMPULSE_MARCH_DONE\n');
end

% ---------------------------------------------------------------------------
function [z_try, p, out, V, landed] = solve_rung(p, z, I, logf)
%SOLVE_RUNG  One ceiling, with the torque march's near-miss polish.
p.limits.impulse_max = I;
opts = struct('MaxFunctionEvaluations', 6e5, 'MaxIterations', 600);
t0 = tic;
[z_try, out] = ch3_col_solve(p, z, opts);
secs = toc(t0);
logln(logf, 'I_max %5.2f: exitflag %d | max|c| %.2e | max|ceq| %.2e | %.0f s', ...
      I, out.exitflag, out.max_c, out.max_ceq, secs);
V = report(logf, p, z_try, sprintf('I<=%.2f', I), secs);

% Polish only a near miss: from a properly infeasible point the torque march's
% polish ran away to a spurious solution (ch3_torque_march, rung 160).
if out.max_c > 1e-6 && out.max_c < 1e-4 && V.ok
    logln(logf, '  polishing rung %.2f (max|c| %.2e)', I, out.max_c);
    t1 = tic;
    [z_p, out_p] = ch3_col_solve(p, z_try, opts);
    V_p = report(logf, p, z_p, sprintf('I<=%.2f polished', I), toc(t1));
    if V_p.ok && out_p.max_c < out.max_c
        z_try = z_p; out = out_p; V = V_p;
    end
end
landed = V.ok && out.max_ceq <= 1e-6 && out.max_c <= 1e-4;
end

% ---------------------------------------------------------------------------
function [z, p, I_done, ok] = bisect(z, p, I_done, I_missed, MIN_STEP, logf, SPD, tag)
%BISECT  Halve the impulse step until a rung lands or the step gets too small.
step = (I_done - I_missed) / 2;
ok   = false;
while step >= MIN_STEP
    I_try = I_done - step;
    logln(logf, '  bisecting: trying %.2f Ns (step %.2f)', I_try, step);
    [z_try, p_try, out, V, landed] = solve_rung(p, z, I_try, logf);
    save_rung(rung_file(SPD, tag, I_try), z_try, p_try, out, V, landed);
    if landed
        z = z_try; p = p_try; I_done = I_try; ok = true;
        step = min(step * 1.5, I_done - I_missed);
        if I_done <= I_missed, return; end
    else
        step = step / 2;
    end
end
end

% ---------------------------------------------------------------------------
function save_rung(f, z_try, p, out, V, landed) %#ok<INUSD>
% SAVE EVERY RUNG, LANDED OR NOT (ch3_torque_march lost 75 minutes of solve
% when it saved only after the verdict). Same fields as the torque march.
save(f, 'z_try', 'p', 'out', 'V', 'landed');
end

function f = rung_file(SPD, tag, I)
f = fullfile(SPD, sprintf('%s_I%04.0f.mat', tag, 100 * I));
end

function V = report(logf, p, z, tag, secs)
E   = ch3_col_eval(z, p);
c   = ch3_col_constraints(z, p);
chk = ch3_col_check_limits(z, p);
V   = ch3_col_verify(z, p, false);
lam = [E.lam, E.lamm];
pk  = max([abs(E.u(:)); abs(E.um(:))]);
logln(logf, ['%-20s v %.4f | T %.4f | L %.4f | peak|u| %6.1f of %5.1f | ' ...
             'impulse %5.2f of %5.2f Ns | mu %.3f | minFz %6.1f | limits ok %d ' ...
             '(max|c| %.2e) | verify %d (%.1e) | %.0f s'], ...
      tag, E.L_step/E.T, E.T, E.L_step, pk, p.limits.u_max, norm(E.impulse), ...
      p.limits.impulse_max, max(abs(lam(1,:))./lam(2,:)), min(lam(2,:)), ...
      chk.ok, chk.max_c, V.ok, V.max_dev, secs);
if c(7) > 1e-6
    logln(logf, '      impulse row still violated: %+.3e Ns', c(7));
end
if c(4) > 1e-6
    logln(logf, '      torque row violated: %+.3e', c(4));
end
end

function logln(f, fmt, varargin)
s = sprintf(fmt, varargin{:});
fprintf('%s\n', s);
fid = fopen(f, 'a'); fprintf(fid, '%s\n', s); fclose(fid);
end
