function ch3_speed_ladder(v_list, mu)
%CH3_SPEED_LADDER  Solve for a gait that walks at the speed ch3_params asks for.
%
%   ch3_speed_ladder()                 the ladder 1.45 1.35 1.25 1.20 m/s, mu 0.4
%   ch3_speed_ladder([1.45 1.35 1.25 1.20 1.10 1.00])   continue the descent
%   ch3_speed_ladder([1.45 1.35 1.25 1.20], 0.35)      and a tighter cone
%
% WHY. ch3_params sets p.v_des = 1.2 m/s with p.enforce_nec1 = true, but every
% stored gait carries enforce_nec1 = 0, so the speed equality was never imposed
% and posture_195 walks at 1.563 m/s -- 30% over the target. At that speed it
% has no margin anywhere: the peak torque sits exactly on the 195 Nm box
% (row 4 = +7.0e-6) and the impact impulse sits exactly on the mu = 0.4 friction
% cone (row 13 = +5.6e-9). Both are plausibly artefacts of solving with the
% speed left free, which lets the optimizer trade speed for cost without ever
% being told what speed the robot is meant to walk at.
%
% THE LADDER STARTS NEXT TO THE SEED, WHICH MEANS IT DESCENDS. A slow gait is
% the easy end of this problem -- less torque, a smaller impact impulse, more
% friction margin -- but the only seed available walks at 1.563 m/s, and the
% speed equality is a hard constraint on a quantity the seed misses badly.
% MEASURED: asking for 0.6 m/s directly stalled in 62 s at exitflag -2 with the
% NEC1 equality violated by 0.96, the gait never leaving 1.563. So the rungs
% step DOWN from the seed in increments of about 0.1 m/s, each warm-starting
% the next, and 1.2 -- the speed ch3_params actually asks for -- is reached
% from above. Lower speeds are then reachable by continuing the same descent.
%
% Each rung reports the margins that matter: how much of the torque box it
% uses, how much of the friction cone the landing uses, and whether the nodes
% still lie on a true trajectory (ch3_col_verify). A rung with slack in both is
% the answer to "is posture_195 simply too fast?".
%
% Output: Results/reruns/speed_ladder/ (one .mat per rung, speed_ladder.log).
% Rungs already saved are reported and skipped, so this resumes.
%
% See also CH3_MU_SENS, CH3_COL_SOLVE, CH3_COL_CHECK_LIMITS.

if nargin < 1 || isempty(v_list), v_list = [1.45 1.35 1.25 1.20]; end
if nargin < 2 || isempty(mu),     mu     = 0.40; end

ROOT = fileparts(fileparts(fileparts(fileparts(mfilename('fullpath')))));
SPD  = fullfile(ROOT, 'Results', 'reruns', 'speed_ladder');
if ~exist(SPD, 'dir'), mkdir(SPD); end
logf = fullfile(SPD, 'speed_ladder.log');

S  = load(fullfile(ROOT, 'Results', 'ch3_gait_posture_195.mat'));
z0 = S.z;
p0 = ch3_upgrade_params(S.p);
p0.limits.mu_s    = mu;
p0.enforce_nec1   = true;      % THE POINT: impose the speed equality
p0.scale_problem  = true;      % off, SQP runs away from a converged gait

% A MARCHING RUNG NEEDS FEASIBILITY, NOT OPTIMALITY. Measured on the first
% rung with a 1.5e6 budget: feasibility reached 2e-7 by iteration 121 (~215k
% evaluations) and then the solve crawled -- 136 more iterations moved the cost
% by 0.006 while first-order optimality sat at 5e7, i.e. it was never going to
% converge and would simply grind to the limit. Eight rungs of that is a day of
% compute for nothing. Each rung is only a warm start for the next, so the
% budget is set to buy a feasible point and stop.
opts = struct('MaxFunctionEvaluations', 3.0e5, 'MaxIterations', 400);

logln(logf, '=== speed ladder %s | mu %.2f | %s', mat2str(v_list), mu, datestr(now));
report(logf, p0, z0, 'baseline', NaN);

% THE FIRST RUNG IS THE SPEED THE SEED ALREADY WALKS. ch3_col_constraints says
% it outright: switch NEC1 on AT THE ACHIEVED SPEED, then march v_des. Skipping
% that step is why the first attempt at this failed -- with v_des set below the
% seed's own speed, NEC1 is violated from iteration one, SQP's first subproblem
% is infeasible, and every rung returned exitflag -2 in under 30 s without the
% gait moving at all (max|ceq| 0.11 / 0.21 / 0.31 / 0.36 at the four targets).
E0 = ch3_col_eval(z0, p0);
v_seed = E0.L_step / E0.T;
v_list = [v_seed, v_list(:).'];
logln(logf, 'first rung is the seed''s own speed, %.4f m/s (NEC1 on at the achieved speed)', v_seed);

v_target = min(v_list);        % the lowest speed asked for; the march aims here
v_done   = v_seed;             % lowest speed reached so far that actually took

z_prev = z0;
for v = v_list
    f = fullfile(SPD, sprintf('gait_v%04.0f.mat', 1000*v));
    if exist(f, 'file')
        L = load(f);
        % A SAVED RUNG IS ONLY A RUNG IF IT ACTUALLY WALKS AT ITS OWN SPEED.
        % An earlier version saved before checking, leaving a "1.10" file that
        % still walked at 1.15; on resume it was replayed as a result and the
        % march bisected from the wrong place.
        if abs(ch3_speed_of(L.z, L.p) - v) > 1e-3
            logln(logf, '  stored rung v %.3f walks at %.4f -- discarding it', ...
                  v, ch3_speed_of(L.z, L.p));
            delete(f);
        else
            z_prev = L.z; v_done = v;
            report(logf, L.p, L.z, sprintf('v=%.3f', v), L.out.wall_time);
            continue;
        end
    end
    p = p0; p.v_des = v;
    t0 = tic;
    [z, out] = ch3_col_solve(p, z_prev, opts);
    out.wall_time = toc(t0);
    save(f, 'z', 'p', 'out');
    logln(logf, 'v %.3f: exitflag %d, fval %.4f, max|c| %.2e, max|ceq| %.2e, %.0f s', ...
          v, out.exitflag, out.fval, out.max_c, out.max_ceq, out.wall_time);
    report(logf, p, z, sprintf('v=%.3f', v), out.wall_time);
    if out.max_ceq > 1e-3
        logln(logf, '  rung v %.3f did not take: max|ceq| %.2e, speed still %.4f', ...
              v, out.max_ceq, ch3_speed_of(z, p));
        delete(f);                  % a rung that never moved is not a result
        v_fail = v;
        break;
    end
    z_prev = z;
    v_done = v;
end

% ------------------------------------------------------------------ bisecting
% A FIXED STEP DIES AT A CORNER. Measured: the march ran 1.5628 -> 1.15 in 0.05
% steps and then 1.10 returned in 15 s with the speed equality untouched
% (max|ceq| = 0.05, exactly the step). At 1.15 three rows are active at once --
% torque, stance friction and the Fz floor -- so SQP's first subproblem has
% nowhere to go and it gives up rather than stepping. Halving the step walks
% around the corner; a step that has to go below MIN_STEP means the march has
% genuinely run out of room, not patience.
MIN_STEP = 0.01;
if exist('v_fail', 'var') && v_fail > v_target + MIN_STEP
    % Bisect from the speed z_prev ACTUALLY walks at, not from a bookkeeping
    % variable: on a resume every rung comes off disk, and a v_done that only
    % updates on freshly solved rungs left the march stepping from 1.5628 while
    % holding the 1.15 gait -- it then "stepped" upward to 1.281 and refused.
    v    = ch3_speed_of(z_prev, p0);
    step = max((v - v_fail) / 2, MIN_STEP);
    while v > v_target + 1e-9
        v_try = max(v_target, v - step);
        f = fullfile(SPD, sprintf('gait_v%04.0f.mat', 1000*v_try));
        p = p0; p.v_des = v_try;
        t0 = tic;
        [z, out] = ch3_col_solve(p, z_prev, opts);
        out.wall_time = toc(t0);
        took = out.max_ceq <= 1e-3;
        logln(logf, 'v %.3f (step %.3f): exitflag %d, max|c| %.2e, max|ceq| %.2e, %s, %.0f s', ...
              v_try, step, out.exitflag, out.max_c, out.max_ceq, ...
              ternary(took, 'TOOK', 'refused'), out.wall_time);
        if took
            save(f, 'z', 'p', 'out');
            report(logf, p, z, sprintf('v=%.3f', v_try), out.wall_time);
            z_prev = z; v = v_try;
            step = min(step * 1.5, 0.05);        % grow again once it is moving
        else
            step = step / 2;
            if step < MIN_STEP
                logln(logf, '  march stops at %.3f m/s: step below %.3f and still refused', ...
                      v, MIN_STEP);
                break;
            end
        end
    end
end
logln(logf, '=== SPEED_LADDER_DONE %s', datestr(now));
fprintf('SPEED_LADDER_DONE\n');
end

% ---------------------------------------------------------------------------
function report(logf, p, z, tag, secs)
E   = ch3_col_eval(z, p);
c   = ch3_col_constraints(z, p);
chk = ch3_col_check_limits(z, p);
V   = ch3_col_verify(z, p, false);
pk  = max([abs(E.u(:)); abs(E.um(:))]);
lam = [E.lam, E.lamm];
mu_step = max(abs(lam(1,:)) ./ lam(2,:));
mu_land = abs(E.impulse(1)) / E.impulse(2);
logln(logf, ['%-9s v %.4f m/s | T %.4f | L %.4f | peak|u| %7.2f Nm (%.1f%% of box) | ' ...
             'mu step %.3f land %.4f (of %.2f) | minFz %.1f N | torque row %+.2e | ' ...
             'cone row %+.2e | limits ok %d (max|c| %.2e) | verify %d (%.1e) | %.0f s'], ...
      tag, E.L_step / E.T, E.T, E.L_step, pk, 100 * pk / p.limits.u_max, ...
      mu_step, mu_land, p.limits.mu_s, min(lam(2,:)), c(4), c(13), ...
      chk.ok, chk.max_c, V.ok, V.max_dev, secs);
phys_rows(logf, c, p);
end

% ---------------------------------------------------------------------------
function phys_rows(logf, c, p)
%PHYS_ROWS  Every physical-realizability row of ch3_col_constraints, by name.
%
% c <= 0 is satisfied. A row whose gate is off is held at -1 by the constraint
% builder and is marked OFF here, because a disabled row is trivially satisfied
% and says nothing about the gait.
rows = { 4, 'peak torque |u| <= u_max      ', 'torque'
         5, 'friction cone |Fx| <= mu Fz   ', 'friction'
         6, 'min normal force Fz >= Fz_min ', 'grf'
         7, 'impact impulse ||I|| <= max   ', 'impulse'
        12, 'impulse compressive Iz >= 0   ', 'impact'
        13, 'impulse inside friction cone  ', 'impact'
         2, 'no ground penetration         ', '(always)'
         9, 'swing foot above ground       ', 'swing_clear'
        11, 'post-impact lift-off          ', 'liftoff'
        14, 'HZD fixed point exists        ', 'hzd'
        15, 'HZD fixed point stable        ', 'hzd' };
for i = 1:size(rows, 1)
    k = rows{i,1}; gate = rows{i,3};
    if strcmp(gate, '(always)'), on = true; else, on = p.limits.enable.(gate); end
    if ~on
        state = 'OFF ';
    elseif c(k) > 1e-6
        state = 'VIOL';
    elseif c(k) > -1e-4
        state = 'ACTIVE';
    else
        state = 'ok  ';
    end
    logln(logf, '    row %2d  %s %-6s %+12.4e', k, rows{i,2}, state, c(k));
end
end

function s = ternary(b, x, y)
if b, s = x; else, s = y; end
end

function v = ch3_speed_of(z, p)
E = ch3_col_eval(z, p);
v = E.L_step / E.T;
end

function logln(f, fmt, varargin)
s = sprintf(fmt, varargin{:});
fprintf('%s\n', s);
fid = fopen(f, 'a'); fprintf(fid, '%s\n', s); fclose(fid);
end
