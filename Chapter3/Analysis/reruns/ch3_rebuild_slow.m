function ch3_rebuild_slow(stages)
%CH3_REBUILD_SLOW  Build a low-speed gait from scratch at N = 81, in stages.
%
%   ch3_rebuild_slow()            every stage, in order, resuming
%   ch3_rebuild_slow('bare')      one stage
%   ch3_rebuild_slow({'grf','friction'})
%
% WHY N = 81 AND NEC1 OFF. Everything cheaper has been tried and recorded:
%
%   * Marching DOWN from posture_195 reaches 1.20 and 1.15 m/s and then walls.
%     At 1.15 the torque box, the stance friction cone and the Fz floor are all
%     active at once; no single relaxation releases it, and the pair that comes
%     closest (mu_s 0.60 with the Fz floor at 10 N) only reaches 1.10 -- i.e.
%     slower walking is bought by weakening the contact model, which is the
%     opposite of the point.
%   * ch3_stage3_from_scratch, the staged cold recipe, fails at its FIRST
%     stage under today's dynamics: exitflag -2 at ceq 3.8e-2 with no physical
%     limit enabled at all, at the default N = 41.
%   * There is no low-speed seed to borrow. ch3_gait_full_constrained, which
%     ch3_speed_march is built on, is stale under the 74 kg model (peak |u| 287
%     against a 120 box, mu 6.94, min Fz -94 N, verify deviation 7.9).
%
% So the mesh is the thing left to change. The project's own history says the
% recipe that VERIFIES from a cold seed is NEC1 off with N = 81: pinning the
% speed while also finding a periodic orbit is what makes the cold solve stall,
% and N = 41 is too coarse for these dynamics to be a real trajectory even when
% the defects look small. Speed is imposed LAST, once there is a gait to impose
% it on.
%
% STAGES, each warm-starting the last, each saved and each verified:
%
%   bare      no physical limits, NEC1 off      -- just find a periodic orbit
%   grf       + minimum normal force  (NIC1)
%   friction  + friction cone         (NIC2)    -- after grf, never before:
%                                                  |Fx|/Fz while Fz crosses zero
%                                                  is a division by ~0
%   torque    + peak torque box       (Tbl 3.1) -- as a ladder, not a step
%   impulse   + impact impulse cap    (Tbl 3.1)
%   speed     NEC1 on at the achieved speed, then march toward 0.5 m/s
%
% A stage that does not verify stops the run: ch3_col_verify failing means the
% nodes are not on a trajectory, and every number after that is fiction.
%
% Output: Results/reruns/rebuild81/ (one .mat per stage, rebuild.log).
%
% See also CH3_STAGE3_FROM_SCRATCH, CH3_SPEED_LADDER, CH3_COL_VERIFY.

all_stages = {'bare','grf','friction','torque','impulse','speed'};
if nargin < 1 || isempty(stages), stages = all_stages; end
if ischar(stages), stages = {stages}; end

ROOT = fileparts(fileparts(fileparts(fileparts(mfilename('fullpath')))));
SPD  = fullfile(ROOT, 'Results', 'reruns', 'rebuild81');
if ~exist(SPD, 'dir'), mkdir(SPD); end
logf = fullfile(SPD, 'rebuild.log');

p = ch3_params();
p.N_nodes      = 81;
p.enforce_nec1 = false;            % speed comes last
for f = {'grf','friction','torque','impulse'}
    p.limits.enable.(f{1}) = false;
end
U_TARGET = p.limits.u_max;         % captured before anything mutates p
I_TARGET = p.limits.impulse_max;

logln(logf, '=== rebuild at N = %d, NEC1 off | %s', p.N_nodes, datestr(now));

z = [];
for i = 1:numel(stages)
    st = stages{i};
    f  = fullfile(SPD, sprintf('stage_%s.mat', st));
    if exist(f, 'file')
        L = load(f); z = L.z; p = L.p;
        report(logf, p, z, [st ' (stored)'], NaN);
        continue;
    end
    if isempty(z)
        [z, p] = load_prev(SPD, st, p, logf);      % resume mid-list
    end

    t0 = tic;
    switch st
        case 'bare'
            z_seed = ch3_col_seed(p);
            E = ch3_col_eval(z_seed, p);
            % ANCHOR THE SPEED EVEN HERE. With NEC1 off AND every physical
            % limit off the problem is under-constrained, and the cost
            % (torque^2 per distance) pays for distance: measured at N = 81,
            % the bare solve ran the speed from 0.1232 to 1.1381 m/s with peak
            % torque 2335 Nm and min Fz -4949 N, ending at ceq 0.31 and a
            % verify deviation of 9.8. Anchoring NEC1 at the seed's own speed
            % costs nothing -- it is satisfied at the start -- and keeps the
            % solve near the slow gait we are actually trying to build.
            p.enforce_nec1 = true;
            p.v_des = E.L_step / E.T;
            logln(logf, 'cold seed walks at %.4f m/s, N = %d; NEC1 anchored there', ...
                  p.v_des, p.N_nodes);
            [z, out] = ch3_col_solve(p, z_seed, big_budget());
        case 'grf'
            p.limits.enable.grf = true;
            [z, out] = ch3_col_solve(p, z, big_budget());
        case 'friction'
            p.limits.enable.friction = true;
            [z, out] = ch3_col_solve(p, z, big_budget());
        case 'torque'
            [z, p, out] = torque_ladder(z, p, U_TARGET, logf);
        case 'impulse'
            p.limits.enable.impulse = true;
            p.limits.impulse_max = I_TARGET;
            [z, out] = ch3_col_solve(p, z, big_budget());
        case 'speed'
            [z, p, out] = speed_stage(z, p, logf);
        otherwise
            error('ch3_rebuild_slow:stage', 'Unknown stage "%s".', st);
    end
    secs = toc(t0);

    save(f, 'z', 'p', 'out');
    logln(logf, 'stage %-9s exitflag %d | max|c| %.2e | max|ceq| %.2e | %.0f s', ...
          st, out.exitflag, out.max_c, out.max_ceq, secs);
    V = report(logf, p, z, st, secs);
    if ~V.ok
        logln(logf, '  STOP: stage %s does not verify (dev %.2e) -- the nodes are not a trajectory', ...
              st, V.max_dev);
        break;
    end
end
logln(logf, '=== REBUILD_DONE %s', datestr(now));
fprintf('REBUILD_DONE\n');
end

% ---------------------------------------------------------------------------
function o = big_budget()
% N = 81 doubles the decision vector against the default mesh, so the cold
% stages get room; a stage that needs more than this is not short of patience.
o = struct('MaxFunctionEvaluations', 2.0e6, 'MaxIterations', 1500);
end

function [z, p, out] = torque_ladder(z, p, U_TARGET, logf)
% Walk the box down to its target in geometric rungs rather than stepping to it:
% the peak torque of an unconstrained gait is far above the box, and clamping it
% in one move is the jump that stalls.
E = ch3_col_eval(z, p);
peak = max([abs(E.u(:)); abs(E.um(:))]);
logln(logf, '  torque ladder: natural peak %.1f Nm -> target %.1f Nm', peak, U_TARGET);
p.limits.enable.torque = true;
if peak <= U_TARGET
    p.limits.u_max = U_TARGET;
    [z, out] = ch3_col_solve(p, z, big_budget());
    return;
end
n = max(3, ceil(log(peak*0.98/U_TARGET)/log(1/0.8)) + 1);
rungs = peak*0.98 * (U_TARGET/(peak*0.98)).^((0:n-1)/(n-1));
for u = rungs
    p.limits.u_max = u;
    [z, out] = ch3_col_solve(p, z, big_budget());
    logln(logf, '    rung u_max %7.1f Nm: exitflag %d, max|c| %.2e, max|ceq| %.2e', ...
          u, out.exitflag, out.max_c, out.max_ceq);
end
end

function [z, p, out] = speed_stage(z, p, logf)
% NEC1 on AT THE ACHIEVED SPEED, then march toward 0.5 -- in whichever
% direction the gait happens to sit relative to it.
E = ch3_col_eval(z, p);
v0 = E.L_step / E.T;
p.enforce_nec1 = true;
p.scale_problem = true;
logln(logf, '  speed stage: gait walks at %.4f m/s; anchoring NEC1 there', v0);
p.v_des = v0;
[z, out] = ch3_col_solve(p, z, big_budget());
step = 0.05 * sign(0.5 - v0);
v = v0;
while abs(v - 0.5) > 1e-9 && abs(step) >= 0.01
    v_try = v + step;
    if (step > 0 && v_try > 0.5) || (step < 0 && v_try < 0.5), v_try = 0.5; end
    p.v_des = v_try;
    [z_try, out_try] = ch3_col_solve(p, z, big_budget());
    if out_try.max_ceq <= 1e-3
        z = z_try; out = out_try; v = v_try;
        logln(logf, '    speed rung %.3f m/s: took (max|c| %.2e)', v, out.max_c);
    else
        step = step / 2;
        logln(logf, '    speed rung %.3f m/s: refused, halving to %.3f', v_try, abs(step));
    end
end
end

function [z, p] = load_prev(SPD, st, p, logf)
order = {'bare','grf','friction','torque','impulse','speed'};
k = find(strcmp(order, st));
for j = k-1:-1:1
    f = fullfile(SPD, sprintf('stage_%s.mat', order{j}));
    if exist(f, 'file')
        L = load(f); z = L.z; p = L.p;
        logln(logf, '  resuming %s from stored stage %s', st, order{j});
        return;
    end
end
z = [];
end

% ---------------------------------------------------------------------------
function V = report(logf, p, z, tag, secs)
E   = ch3_col_eval(z, p);
c   = ch3_col_constraints(z, p);
chk = ch3_col_check_limits(z, p);
V   = ch3_col_verify(z, p, false);
lam = [E.lam, E.lamm];
logln(logf, ['%-16s v %.4f | T %.4f | L %.4f | peak|u| %7.1f (box %.0f) | mu %.3f | ' ...
             'minFz %7.1f | limits ok %d (max|c| %.2e) | verify %d (%.1e) | %.0f s'], ...
      tag, E.L_step/E.T, E.T, E.L_step, max([abs(E.u(:));abs(E.um(:))]), ...
      p.limits.u_max, max(abs(lam(1,:))./lam(2,:)), min(lam(2,:)), ...
      chk.ok, chk.max_c, V.ok, V.max_dev, secs);
end

function logln(f, fmt, varargin)
s = sprintf(fmt, varargin{:});
fprintf('%s\n', s);
fid = fopen(f, 'a'); fprintf(fid, '%s\n', s); fclose(fid);
end
