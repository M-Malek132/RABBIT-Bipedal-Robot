function ch3_stride_free_theta(phase)
%CH3_STRIDE_FREE_THETA  The stride march again, with the phase endpoints free.
%
%   ch3_stride_free_theta()          both phases
%   ch3_stride_free_theta('stride')  phase 1 only
%   ch3_stride_free_theta('speed')   phase 2 only (needs phase 1 saved)
%
% ch3_stride_march could not shorten the stride by a millimetre. With
% theta_minus and theta_plus fixed, the stride is not free to begin with: both
% feet are on the ground at the strike, theta is the hip-to-foot direction, so
%
%       L_step = h_imp (tan theta_plus - tan theta_minus)
%
% (every stored gait satisfies it to ~1e-7), and the only lever left is the
% hip height at impact -- floored at 0.880 m by the hip band, i.e. a stride of
% at least 0.405 m, so the 0.40 m rung of that march was infeasible outright.
% ch3_stride_licq measures what blocks smaller cuts at the stuck iterate.
%
% This driver repeats phase 1 of that march with p.free_theta = true, so the
% sweep of the stance leg is a decision variable (inside p.theta_bounds), and
% then phase 2, the speed march on the shortest stride that lands -- the
% question the stride march was built to answer: can this robot walk slower
% on a shorter step?
%
% Same seed (the verified 1.2 m/s gait), same options as ch3_stride_march
% except the free phase: ScaleProblem on, NEC1 off in phase 1 and the torque
% cost per step (cost_normalize false) with the ceiling boxing the stride.
% Rungs of 1 cm. Each landed rung is saved with the p the solve returns
% (out.p), which carries its solved theta_minus / theta_plus.
%
% Output: Results/reruns/stride_free_theta/ (one .mat per rung, march.log).
% Runtime: minutes to tens of minutes per rung at N = 61; a few hours if it
% walks all the way down.
%
% See also CH3_STRIDE_MARCH, CH3_STRIDE_LICQ, CH3_COL_PACK.

if nargin < 1 || isempty(phase), phase = {'stride', 'speed'}; end
if ischar(phase), phase = {phase}; end

ROOT = fileparts(fileparts(fileparts(fileparts(mfilename('fullpath')))));
SPD  = fullfile(ROOT, 'Results', 'reruns', 'stride_free_theta');
if ~exist(SPD, 'dir'), mkdir(SPD); end
logf = fullfile(SPD, 'march.log');

opts = struct('MaxFunctionEvaluations', 6e5, 'MaxIterations', 600);

%% ------------------------------------------------------- phase 1: the stride
if any(strcmp(phase, 'stride'))
    S = load(fullfile(ROOT, 'Results', 'reruns', 'speed_ladder', 'gait_v1200.mat'));
    p = ch3_upgrade_params(S.p);
    p.free_theta    = true;
    p.scale_problem = true;
    p.enforce_nec1  = false;                       % one demand at a time
    p.limits.enable.step_len_max = true;
    p.cost_normalize = false;                      % the ceiling boxes the stride
    z = ch3_col_theta_augment(S.z, p);
    E = ch3_col_eval(z, p);
    logln(logf, '=== free-theta stride march | speed free | seed L = %.4f m | %s', ...
          E.L_step, datestr(now)); %#ok<TNOW1,DATST>
    report(logf, p, z, 'seed', NaN);

    for L_cap = round((floor(E.L_step*100)/100 - 0.01):-0.01:0.30, 4)
        f = fullfile(SPD, sprintf('stride_L%03.0f.mat', 1000*L_cap));
        if exist(f, 'file')
            L = load(f);
            if L.landed
                z = L.z_try; p = L.p;
                report(logf, p, z, sprintf('L<=%.2f (stored)', L_cap), NaN);
                continue;
            end
        end
        p.limits.step_len_max = L_cap;
        t0 = tic;
        [z_try, out] = ch3_col_solve(p, z, opts);
        p_try = out.p;                              % the solved theta pair
        V = report(logf, p_try, z_try, sprintf('L<=%.2f', L_cap), toc(t0));
        landed = V.ok && out.max_ceq <= 1e-6 && out.max_c <= 1e-4;
        p = p_try; %#ok<NASGU> saved below as p
        save(f, 'z_try', 'p', 'out', 'V', 'landed');
        if ~landed
            logln(logf, ['  rung L<=%.2f missed (max|c| %.2e, max|ceq| %.2e, verify %d); ' ...
                         'phase 1 stops, the shortest landed stride is the rung above.'], ...
                  L_cap, out.max_c, out.max_ceq, V.ok);
            break;
        end
        z = z_try;
    end
end

%% -------------------------------------------------------- phase 2: the speed
if any(strcmp(phase, 'speed'))
    [z, p, tag] = shortest_landed(SPD);
    if isempty(z), logln(logf, 'phase 2 skipped: no landed stride rung'); return; end
    E = ch3_col_eval(z, p);
    logln(logf, '=== speed march on the %s stride (L = %.4f m) | %s', tag, E.L_step, ...
          datestr(now)); %#ok<TNOW1,DATST>
    p.enforce_nec1 = true;
    p.v_des = E.L_step / E.T;                        % anchor at the achieved speed
    for v = round((p.v_des - 0.05):-0.05:0.50, 4)
        f = fullfile(SPD, sprintf('short_v%04.0f.mat', 1000*v));
        if exist(f, 'file')
            L = load(f);
            if L.landed
                z = L.z_try; p = L.p;
                report(logf, p, z, sprintf('v=%.2f (stored)', v), NaN);
                continue;
            end
        end
        p.v_des = v;
        t0 = tic;
        [z_try, out] = ch3_col_solve(p, z, opts);
        p_try = out.p;
        V = report(logf, p_try, z_try, sprintf('v=%.2f', v), toc(t0));
        landed = V.ok && out.max_ceq <= 1e-6 && out.max_c <= 1e-4;
        p = p_try; %#ok<NASGU>
        save(f, 'z_try', 'p', 'out', 'V', 'landed');
        if ~landed
            logln(logf, '  speed rung %.2f missed on the short stride; phase 2 stops', v);
            break;
        end
        z = z_try;
    end
end
logln(logf, '=== STRIDE_FREE_THETA_DONE %s', datestr(now)); %#ok<TNOW1,DATST>
fprintf('STRIDE_FREE_THETA_DONE\n');
end

% ---------------------------------------------------------------------------
function [z, p, tag] = shortest_landed(SPD)
z = []; p = []; tag = '';
F = dir(fullfile(SPD, 'stride_L*.mat'));
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
logln(logf, ['%-18s v %.4f | T %.4f | L %.4f | theta %+.4f..%+.4f | h_imp %.4f | ' ...
             'peak|u| %6.1f | mu %.3f | minFz %6.1f | impulse %5.2f Ns | limits ok %d ' ...
             '(max|c| %.2e) | verify %d (%.1e) | %.0f s'], ...
      tag, E.L_step/E.T, E.T, E.L_step, E.p.theta_minus, E.p.theta_plus, ...
      -E.X(2, end), max([abs(E.u(:)); abs(E.um(:))]), ...
      max(abs(lam(1,:))./lam(2,:)), min(lam(2,:)), norm(E.impulse), ...
      chk.ok, chk.max_c, V.ok, V.max_dev, secs);
end

function logln(f, fmt, varargin)
s = sprintf(fmt, varargin{:});
fprintf('%s\n', s);
fid = fopen(f, 'a'); fprintf(fid, '%s\n', s); fclose(fid);
end
