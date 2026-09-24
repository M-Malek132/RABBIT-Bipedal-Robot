function ch3_stride_march(phase)
%CH3_STRIDE_MARCH  Shorten the stride, then try to walk slower on it.
%
%   ch3_stride_march()          both phases
%   ch3_stride_march('stride')  phase 1 only
%   ch3_stride_march('speed')   phase 2 only (needs phase 1 saved)
%
% WHY STRIDE AND NOT SPEED. Marching the speed down from posture_195 walls at
% 1.15 m/s with the torque box, the stance friction cone and the Fz floor all
% active at once. Torque is NOT what blocks it: asked for 1.0 m/s with the box
% opened to 250 Nm, the solve drew only 192.8 Nm -- 57 Nm of unused headroom --
% and still could not get below 1.1853 m/s, with the whole equality residual
% being the speed gap. What runs out is the contact: a slow gait on a long
% stride has little momentum to carry it through, so Fz sags at the end of
% stance and the friction demand climbs.
%
% The hypothesis this tests: a SHORTER stride keeps the contact healthy at low
% speed -- less leg extension at the end of stance, so Fz stays higher, and
% less horizontal redirection, so mu stays lower. It also attacks the other
% unmet limit, since a shorter step lands more slowly (the impulse on these
% gaits is 18.8-19.1 Ns against a declared 15).
%
% HOW THE STRIDE IS CONTROLLED. Row 19 of ch3_col_constraints,
% L_step <= p.limits.step_len_max, gated by p.limits.enable.step_len_max and
% off by default. It was added for this: the objective is torque-squared per
% unit distance, so a longer stride is always cheaper and nothing in the
% problem ever asks for a short one.
%
% THE INDIRECT ROUTE DOES NOT WORK, which is why the row exists. With NEC1
% pinning the speed, L = v*T, so capping T caps the stride -- but T is a BOUND
% on the decision vector, and fmincon clamps a warm start into a violated bound
% before the first iteration. Measured twice: a cap 0.065 s below the seed's T
% detonated the defects to max|c| 12.2 with the nodes 23 off a true rollout,
% and even a 0.01 s cap left them 0.19 off. A constraint row is warm-start-safe
% where a bound is not, because an infeasible start is what SQP is built for.
%
% Phase 1 walks the ceiling down with the speed FREE, so the stride shortens
% without also forcing a speed the gait cannot hold -- one demand at a time.
% Phase 2 then pins the speed on the shortest stride that landed and marches it
% down, which is the question that started this.
%
% Output: Results/reruns/stride_march/ (one .mat per rung, stride_march.log).

if nargin < 1 || isempty(phase), phase = {'stride', 'speed'}; end
if ischar(phase), phase = {phase}; end

ROOT = fileparts(fileparts(fileparts(fileparts(mfilename('fullpath')))));
SPD  = fullfile(ROOT, 'Results', 'reruns', 'stride_march');
if ~exist(SPD, 'dir'), mkdir(SPD); end
logf = fullfile(SPD, 'stride_march.log');

opts = struct('MaxFunctionEvaluations', 6e5, 'MaxIterations', 600);

%% ------------------------------------------------------- phase 1: the stride
if any(strcmp(phase, 'stride'))
    S = load(fullfile(ROOT, 'Results', 'reruns', 'speed_ladder', 'gait_v1200.mat'));
    z = S.z;
    p = ch3_upgrade_params(S.p);
    p.scale_problem = true;
    E = ch3_col_eval(z, p);
    p.enforce_nec1 = false;                      % one demand at a time
    p.limits.enable.step_len_max = true;
    % The ceiling alone is not enough: with the objective still dividing by
    % step length, SQP stalls at a step of 1e-7 with the violation frozen. The
    % stride is boxed by rows 3 and 19 here, so the normalization has nothing
    % left to protect.
    p.cost_normalize = false;
    logln(logf, '=== stride march | speed free | seed L = %.4f m | %s', E.L_step, datestr(now));
    report(logf, p, z, 'seed', NaN);

    for L_cap = round((floor(E.L_step*50)/50 - 0.02):-0.02:0.24, 4)
        f = fullfile(SPD, sprintf('stride_L%03.0f.mat', 1000*L_cap));
        if exist(f, 'file')
            L = load(f);
            if L.landed, z = L.z_try; p = L.p; report(logf, p, z, sprintf('L<=%.2f (stored)', L_cap), NaN); continue; end
        end
        p.limits.step_len_max = L_cap;
        t0 = tic;
        [z_try, out] = ch3_col_solve(p, z, opts);
        V = report(logf, p, z_try, sprintf('L<=%.2f', L_cap), toc(t0));
        landed = V.ok && out.max_ceq <= 1e-6 && out.max_c <= 1e-4;
        save(f, 'z_try', 'p', 'out', 'V', 'landed');
        if ~landed
            logln(logf, ['  stride rung L<=%.2f missed (max|c| %.2e, max|ceq| %.2e, verify %d). ' ...
                         'Phase 1 stops; the shortest landed stride is the rung above.'], ...
                  L_cap, out.max_c, out.max_ceq, V.ok);
            break;
        end
        z = z_try;
    end
end

%% -------------------------------------------------------- phase 2: the speed
if any(strcmp(phase, 'speed'))
    [z, p, tag] = shortest_landed(SPD, ROOT, logf);
    if isempty(z), logln(logf, 'phase 2 skipped: no landed stride rung'); return; end
    E = ch3_col_eval(z, p);
    logln(logf, '=== speed march on the %s stride (L = %.4f m) | %s', tag, E.L_step, datestr(now));
    % The ceiling stays on, so the stride cannot creep back out while the speed
    % comes down -- which is exactly what it would do, the objective paying for
    % distance.
    p.enforce_nec1 = true;
    p.v_des = E.L_step / E.T;
    for v = round((p.v_des - 0.05):-0.05:0.50, 4)
        f = fullfile(SPD, sprintf('short_v%04.0f.mat', 1000*v));
        if exist(f, 'file')
            L = load(f);
            if L.landed, z = L.z_try; p = L.p; report(logf, p, z, sprintf('v=%.2f (stored)', v), NaN); continue; end
        end
        p.v_des = v;
        t0 = tic;
        [z_try, out] = ch3_col_solve(p, z, opts);
        V = report(logf, p, z_try, sprintf('v=%.2f', v), toc(t0));
        landed = V.ok && out.max_ceq <= 1e-6 && out.max_c <= 1e-4;
        save(f, 'z_try', 'p', 'out', 'V', 'landed');
        if ~landed
            logln(logf, '  speed rung %.2f missed on the short stride; phase 2 stops', v);
            break;
        end
        z = z_try;
    end
end
logln(logf, '=== STRIDE_MARCH_DONE %s', datestr(now));
fprintf('STRIDE_MARCH_DONE\n');
end

% ---------------------------------------------------------------------------
function [z, p, tag] = shortest_landed(SPD, ROOT, logf) %#ok<INUSD>
z = []; p = []; tag = '';
F = dir(fullfile(SPD, 'stride_T*.mat'));
best = inf;
for k = 1:numel(F)
    L = load(fullfile(SPD, F(k).name));
    if ~L.landed, continue; end
    E = ch3_col_eval(L.z_try, L.p);
    if E.L_step < best
        best = E.L_step; z = L.z_try; p = L.p; tag = F(k).name;
    end
end
end

function V = report(logf, p, z, tag, secs)
E   = ch3_col_eval(z, p);
chk = ch3_col_check_limits(z, p);
V   = ch3_col_verify(z, p, false);
lam = [E.lam, E.lamm];
logln(logf, ['%-18s v %.4f | T %.4f | L %.4f | peak|u| %6.1f | mu %.3f | minFz %6.1f | ' ...
             'impulse %5.2f Ns | limits ok %d (max|c| %.2e) | verify %d (%.1e) | %.0f s'], ...
      tag, E.L_step/E.T, E.T, E.L_step, max([abs(E.u(:)); abs(E.um(:))]), ...
      max(abs(lam(1,:))./lam(2,:)), min(lam(2,:)), norm(E.impulse), ...
      chk.ok, chk.max_c, V.ok, V.max_dev, secs);
end

function logln(f, fmt, varargin)
s = sprintf(fmt, varargin{:});
fprintf('%s\n', s);
fid = fopen(f, 'a'); fprintf(fid, '%s\n', s); fclose(fid);
end
