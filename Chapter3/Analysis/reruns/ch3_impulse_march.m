function ch3_impulse_march(seed_file, tag)
%CH3_IMPULSE_MARCH  Bring the torque-march gait under the declared 15 Ns impulse.
%
%   ch3_impulse_march()                 from the 162.5 Nm gait (torque_march/gait_u162)
%   ch3_impulse_march(seed_file, tag)   from any stored gait; tag names its files, e.g.
%   ch3_impulse_march(fullfile('Results','reruns','torque_march','gait_u180.mat'), 'u180')
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
% KEEP THE BEST ITERATE THAT ALREADY LANDS. The first run (2026-09-25) asked
% for 19.00 Ns from the 19.07 seed -- a 0.4% cut -- and fmincon reported
% "converged to an infeasible point" (max c 1.05e-2, in neither the torque nor
% the impulse row) after 1430 s. But its own log shows it MEETING every limit
% on the way: feasibility 2.2e-8 at iteration 19, 2.0e-7 at 25, 1.1e-6 at 85.
% The seed sits on a corner (torque, stance friction and the Fz floor all
% active), first-order optimality ran at 1e8-1e9, and SQP kept chasing the
% torque cost off the feasible set it had found. ch3_col_solve returns only the
% last iterate, so the march threw those points away and called the rung a
% miss. An OutputFcn now records the lowest-cost iterate that passes the
% landing test IN TRUE UNITS (max c <= 1e-4, max|ceq| <= 1e-6; the Feasibility
% column is in scaled units under ScaleProblem), stops the solve once PATIENCE
% iterations bring no better one, and the rung uses that iterate if the final
% one misses. It is still mesh-verified before it counts. A point that meets
% every limit and verifies is a real gait within the limits; it is not the
% torque-optimal one, and the question here is existence, not optimality.
%
% RUNGS. 1 Ns cuts below the seed (18, 17, 16, 15 from 19.07), each ~5%, the
% size of the torque march's steps that landed; a rung that misses is bisected
% down to MIN_STEP = 0.1 Ns. The first run's 0.25 Ns minimum made the 0.07 Ns
% first rung unbisectable, so one failed solve ended the march.
%
% KEEP THE GATE ON THE p THAT IS SAVED. Each rung file stores the p it was
% solved with, impulse gate enabled and the ceiling it was solved at.
% Re-solving a stored z with a p whose gate is off lets the impulse spring back
% toward the unconstrained optimum -- which looks exactly like the constraint
% failing (observed on the NEC3 march: 0.761 -> 1.134 in one solve).
%
% A STALL IS NOT A FLOOR. If the march stops above 15, that is a statement
% about continuation from this seed with this solver, as the stride march's
% negative result was; the verdict line says so.
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
MIN_STEP = 0.1;                          % [Ns]

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

rungs = floor(I0 - 0.5):-1:target;       % 18 17 16 15 from a 19.07 Ns seed
if isempty(rungs) || rungs(end) > target, rungs = [rungs, target]; end
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
            logln(logf, ['  march stops at %.2f Ns: the next cut did not land even ' ...
                         'at a %.2f Ns step'], I_done, MIN_STEP);
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
    logln(logf, ['=== VERDICT: at a %.1f Nm box the march stalls at %.2f Ns, above the ' ...
                 'declared %.0f. That is a continuation result from this seed with ' ...
                 'this solver, not a proof that no gait below it exists.'], ...
          p.limits.u_max, I_done, target);
end
logln(logf, '=== IMPULSE_MARCH_DONE %s', datestr(now));
fprintf('IMPULSE_MARCH_DONE\n');
end

% ---------------------------------------------------------------------------
function [z_try, p, out, V, landed] = solve_rung(p, z, I, logf)
%SOLVE_RUNG  One ceiling. Uses the best iterate that already lands if the
% final one does not (see the header).
p.limits.impulse_max = I;
tracker = containers.Map();
tracker('best')  = struct('z', [], 'J', inf, 'iter', NaN, 'max_c', NaN, 'max_ceq', NaN);
tracker('since') = 0;
opts = struct('MaxFunctionEvaluations', 6e5, 'MaxIterations', 600, ...
              'OutputFcn', @(zz, ov, state) track_landed(zz, ov, state, p, tracker));
t0 = tic;
[z_try, out] = ch3_col_solve(p, z, opts);
secs = toc(t0);
logln(logf, 'I_max %5.2f: exitflag %d | max|c| %.2e | max|ceq| %.2e | %.0f s', ...
      I, out.exitflag, out.max_c, out.max_ceq, secs);

best = tracker('best');
out.harvested = NaN;
if ~(out.max_ceq <= 1e-6 && out.max_c <= 1e-4) && ~isempty(best.z)
    logln(logf, ['  the final iterate misses; iterate %d already met every limit ' ...
                 '(max|c| %.2e, max|ceq| %.2e, J %.2f) and is used instead'], ...
          best.iter, best.max_c, best.max_ceq, best.J);
    z_try = best.z;
    out.max_c = best.max_c;  out.max_ceq = best.max_ceq;  out.harvested = best.iter;
end
V = report(logf, p, z_try, sprintf('I<=%.2f', I), secs);

% Polish only a near miss that was NOT harvested: from a properly infeasible
% point the torque march's polish ran away to a spurious solution (rung 160),
% and a harvested iterate already passes the landing test.
if isnan(out.harvested) && out.max_c > 1e-6 && out.max_c < 1e-4 && V.ok
    logln(logf, '  polishing rung %.2f (max|c| %.2e)', I, out.max_c);
    t1 = tic;
    [z_p, out_p] = ch3_col_solve(p, z_try, struct('MaxFunctionEvaluations', 6e5, ...
                                                  'MaxIterations', 600));
    V_p = report(logf, p, z_p, sprintf('I<=%.2f polished', I), toc(t1));
    if V_p.ok && out_p.max_c < out.max_c
        out_p.harvested = NaN;
        z_try = z_p; out = out_p; V = V_p;
    end
end
landed = V.ok && out.max_ceq <= 1e-6 && out.max_c <= 1e-4;
end

% ---------------------------------------------------------------------------
function stop = track_landed(z, ov, state, p, tracker)
%TRACK_LANDED  fmincon OutputFcn: remember the lowest-cost iterate that meets
% the landing test in true units; stop once PATIENCE iterations bring no better.
PATIENCE = 30;
stop = false;
if ~strcmp(state, 'iter'), return; end
best = tracker('best');
if ov.fval < best.J
    [c, ceq] = ch3_col_constraints(z, p);   % one evaluation; an iteration costs ~1760
    if max(c) <= 1e-4 && max(abs(ceq)) <= 1e-6
        tracker('best')  = struct('z', z, 'J', ov.fval, 'iter', ov.iteration, ...
                                  'max_c', max(c), 'max_ceq', max(abs(ceq)));
        tracker('since') = ov.iteration;
        return;
    end
end
if ~isempty(best.z) && ov.iteration - tracker('since') >= PATIENCE
    stop = true;
end
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
[cw, iw] = max(c);
logln(logf, ['%-20s v %.4f | T %.4f | L %.4f | peak|u| %6.1f of %5.1f | ' ...
             'impulse %5.2f of %5.2f Ns | mu %.4f | minFz %7.2f | limits ok %d ' ...
             '(max|c| %.2e, worst row %d) | verify %d (%.1e) | %.0f s'], ...
      tag, E.L_step/E.T, E.T, E.L_step, pk, p.limits.u_max, norm(E.impulse), ...
      p.limits.impulse_max, max(abs(lam(1,:))./lam(2,:)), min(lam(2,:)), ...
      chk.ok, cw, iw, V.ok, V.max_dev, secs);
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
