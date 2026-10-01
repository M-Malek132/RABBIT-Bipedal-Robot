function ch3_ub_ladder(seed_file, v_anchor, L_anchor, n_rungs, dv)
%CH3_UB_LADDER  March a crutch-supported gait down in speed on a shorter stride.
%
%   ch3_ub_ladder()                   w5000 seed, anchors (0.80, 0.40) (0.60, 0.36)
%                                     (0.50, 0.33) (0.40, 0.30) (0.30, 0.27) [m/s, m]
%   ch3_ub_ladder([], [], [], 1)      solve one rung and stop: the probe
%   ch3_ub_ladder(seed_file, v_anchor, L_anchor, n_rungs, dv)
%
% THE STRIDE NEEDS FREE THETA. With theta_minus, theta_plus fixed, both feet
% on the ground at the strike give
%
%       L_step = h_imp (tan theta_plus - tan theta_minus) = 0.4605 h_imp
%
% and the hip band (row 8) floors h_imp at 0.880 m, so no stride below 0.405 m
% exists: a 0.40 m ceiling is infeasible by construction and every lower one
% more so. That is the stride march's negative result (docs/HANDOFF.md,
% "Stride: not a design variable"). p.free_theta makes the leg sweep a decision
% variable, and ch3_stride_licq found a first-order direction that shortens the
% stride once it is -- so the first rung here is also the first test of that.
%
% RUNGS. The anchors are joined piecewise-linearly and walked in speed steps of
% at most dv (0.05 m/s, 1-1.5 cm of stride): the speed ladder only ever moved
% 0.05 at a time, and a 0.6 m/s target asked for directly stalled at exitflag
% -2 without the gait moving. The objective stays per-distance, so it pushes
% the stride out against the ceiling and a rung lands on its cap.
%
% A RUNG TAKES ONLY IF IT IS A GAIT: the equalities, the inequalities,
% ch3_col_verify, and the speed and stride measured on the result. Not the
% exitflag -- 0 means out of budget, feasible or not. The first miss stops the
% march and is saved as miss_v*.mat for inspection.
%
% BUDGET. 150 iterations a rung: a rung is a warm start, not a result, and the
% w50/w500 solves spent 600 iterations (5.3 h each) without converging. About
% 32 s an iteration at N = 81, so at most ~80 min a rung.
%
% Output: Results/reruns/ub_ladder_<tag>/ (git-ignored), rung_v0950.mat ...
% and ladder.log, <tag> from the seed name (ch3_ub_v100_w5000 -> w5000). A
% landed ANCHOR rung is also saved as Results/ch3_ub_ladder_<tag>_v080.mat etc.,
% with z and o like the seed. A rerun resumes from the lowest landed rung.
%
% See also CH3_STRIDE_FREE_THETA, CH3_COL_THETA_AUGMENT, CH3_UB_DEFAULTS.

if nargin < 1 || isempty(seed_file), seed_file = fullfile('Results', 'ch3_ub_v100_w5000.mat'); end
if nargin < 2 || isempty(v_anchor),  v_anchor  = [0.80 0.60 0.50 0.40 0.30]; end
if nargin < 3 || isempty(L_anchor),  L_anchor  = [0.40 0.36 0.33 0.30 0.27]; end
if nargin < 4 || isempty(n_rungs),   n_rungs   = Inf; end
if nargin < 5 || isempty(dv),        dv        = 0.05; end
if numel(v_anchor) ~= numel(L_anchor)
    error('ch3_ub_ladder:anchors', 'v_anchor and L_anchor need one entry each per anchor.');
end

ROOT = fileparts(fileparts(fileparts(mfilename('fullpath'))));
if ~isfile(seed_file), seed_file = fullfile(ROOT, seed_file); end
[~, seed_name] = fileparts(seed_file);
tag  = regexprep(seed_name, '^ch3_ub_(v\d+_)?', '');
OUT  = fullfile(ROOT, 'Results', 'reruns', ['ub_ladder_' tag]);
if ~exist(OUT, 'dir'), mkdir(OUT); end
logf = fullfile(OUT, 'ladder.log');

S = load(seed_file);
if isfield(S, 'o'), p = S.o.p; else, p = S.p; end
p = ch3_upgrade_params(p);
p.free_theta      = true;      % the stride is locked without it (header)
p.scale_problem   = true;      % unscaled SQP runs away from a converged gait
p.enforce_nec1    = true;      % without it v_des is ignored and nothing moves
p.limits.enable.step_len_max = true;
p.max_iter        = 150;
p.max_fun_evals   = 4e5;       % above 2*numel(z)*max_iter, so max_iter binds
p.checkpoint_file = '';
z = ch3_col_theta_augment(S.z, p);

% The first rung starts from the speed and stride the seed ACTUALLY walks.
E  = ch3_col_eval(z, p);
vA = [E.L_step / E.T, v_anchor(:).'];
LA = [E.L_step,       L_anchor(:).'];
rungs = zeros(0, 2);
for k = 2:numel(vA)
    n = max(1, ceil(abs(vA(k-1) - vA(k)) / dv - 1e-6));
    s = (1:n).' / n;
    rungs = [rungs; vA(k-1) + s*(vA(k) - vA(k-1)), LA(k-1) + s*(LA(k) - LA(k-1))]; %#ok<AGROW>
end
rungs = round(rungs, 4);

ch3_logln(logf, sprintf('=== ub ladder | seed %s | %d rungs to v %.2f | %s', ...
          seed_name, size(rungs, 1), rungs(end, 1), datestr(now))); %#ok<TNOW1,DATST>
report(logf, 'seed', z, p);

n_solved = 0;
for i = 1:size(rungs, 1)
    v = rungs(i, 1);  L_cap = rungs(i, 2);
    f = fullfile(OUT, sprintf('rung_v%04.0f.mat', 1000*v));
    if isfile(f)                                   % only landed rungs get this name
        R = load(f);  z = R.z;  p = R.o.p;
        report(logf, sprintf('v %.2f L<=%.3f (stored)', v, L_cap), z, p);
        continue;
    end
    if n_solved >= n_rungs, break; end

    p.v_des = v;
    p.limits.step_len_max = L_cap;
    [z_try, o] = ch3_col_solve(p, z);
    n_solved = n_solved + 1;
    E = report(logf, sprintf('v %.2f L<=%.3f', v, L_cap), z_try, o.p, o);

    landed = o.max_ceq <= 1e-6 && o.max_c <= 1e-4 && o.verify.ok ...
             && abs(E.L_step / E.T - v) <= 1e-3 && E.L_step <= L_cap + 1e-4;
    if ~landed
        save(fullfile(OUT, sprintf('miss_v%04.0f.mat', 1000*v)), 'z_try', 'o');
        ch3_logln(logf, sprintf(['  rung v %.2f missed: the march stops, the last ' ...
                  'landed gait is the rung above'], v));
        break;
    end
    z = z_try;  p = o.p;                           % o.p carries the solved theta and beta
    save(f, 'z', 'o');
    if any(abs(v - v_anchor) < 1e-6)
        save(fullfile(ROOT, 'Results', sprintf('ch3_ub_ladder_%s_v%03.0f.mat', tag, 100*v)), 'z', 'o');
    end
end
ch3_logln(logf, sprintf('=== UB_LADDER_DONE %s', datestr(now))); %#ok<TNOW1,DATST>
end

% ---------------------------------------------------------------------------
function E = report(logf, label, z, p, o)
E   = ch3_col_eval(z, p);
lam = [E.lam, E.lamm];
W   = cell2mat(arrayfun(@(s) ch3_bezier(E.p.ub.beta, s), linspace(0, 1, 41), ...
                        'UniformOutput', false));
msg = sprintf(['%-22s v %.4f | T %.4f | L %.4f | theta %+.3f..%+.3f | h_imp %.3f | ' ...
               'peak|u| %5.1f | mu %.3f | minFz %5.1f | crutch Fz mean %5.1f max %5.1f'], ...
              label, E.L_step / E.T, E.T, E.L_step, E.p.theta_minus, E.p.theta_plus, ...
              -E.X(2, end), max(abs([E.u(:); E.um(:)])), max(abs(lam(1,:)) ./ lam(2,:)), ...
              min(lam(2,:)), mean(W(2,:)), max(W(2,:)));
if nargin > 4
    msg = [msg, sprintf(' | exitflag %d, J %.1f, max|ceq| %.1e, max c %.1e, verify %d (%.1e), %.0f s', ...
                        o.exitflag, o.fval, o.max_ceq, o.max_c, o.verify.ok, ...
                        o.verify.max_dev, o.wall_time)];
end
ch3_logln(logf, msg);
end
