function ch3_dt_sweep(eps_list, dt_list)
%CH3_DT_SWEEP  The Chapter-3 controller table at several control periods.
%
%   ch3_dt_sweep()                    eps 0.5 and 0.2; dt 1, 0.5, 0.2, 0.1 ms and 0
%   ch3_dt_sweep(eps_list, dt_list)   dt = 0 means continuous control
%
% WHY. Chapter 3 explains "PD walks, the min-norm CLF-QP falls" by the CLF
% RATE: c3/eps = 0.732 at eps 0.5 cuts V by only 18% per step, and the impact
% then kicks eta back up. Chapter 4 measured something that competes with
% that explanation: holding the input for one sample, with a PERFECT model,
% leaves ||eta+|| = 0.242 after the first impact at 1 ms, 0.123 at 0.5 ms,
% 0.049 at 0.2 ms and 1.7e-4 in continuous time. So part of what Chapter 3
% attributes to the rate may be the 1 kHz hold. This separates the two: the
% same table (tab:ctrl -- posture_195, the gait's own CARE CLF, 4 steps, fall
% = stride under p.step_len_min, tracking error over 1e3, or a non-finite
% state; only clfqp_con told about the box) at each (eps, dt).
%
% It also answers whether the sampling effect is MULTIPLICATIVE (the impact
% amplifying error the step already carried) or ADDITIVE (the hold creating
% error at every impact): every step logs ||eta-|| just before the impact,
% ||eta+|| just after, and their ratio. A ratio that stays put while both
% shrink with dt is the impact amplifying an O(dt) hold error.
%
% Continuous control (dt = 0) is run for the smooth laws only (iolin_pd,
% clfqp); the constrained QP is only piecewise smooth and stalls an adaptive
% solver (the p.control_dt note in ch3_params), so its smallest dt stands in.
% Validity there is scored from the contact force recomputed at every output
% point under the law's own u(x).
%
% Output: Results/reruns/ch3_dt_sweep/ (sweep.log, one .mat per run).
% Runtime: roughly 30-60 minutes for the default grid; the 0.1 ms runs are
% the slow ones.
%
% See also CH3_TABLE_RERUN, CH3_VALIDITY, CH3_STEP.

if nargin < 1 || isempty(eps_list), eps_list = [0.5 0.2]; end
if nargin < 2 || isempty(dt_list),  dt_list  = [1e-3 5e-4 2e-4 1e-4 0]; end

ROOT = fileparts(fileparts(fileparts(fileparts(mfilename('fullpath')))));
SPD  = fullfile(ROOT, 'Results', 'reruns', 'ch3_dt_sweep');
if ~exist(SPD, 'dir'), mkdir(SPD); end
logf = fullfile(SPD, 'sweep.log');

S = load(fullfile(ROOT, 'Results', 'ch3_gait_posture_195.mat'));
z = S.z; p0 = ch3_upgrade_params(S.p);
E = ch3_col_eval(z, p0);
p0 = E.p;
[X, ~, alpha] = ch3_col_unpack(z, p0);
peak_ff = max([abs(E.u(:)); abs(E.um(:))]);
p0.T_max = 0.5;
n_steps = 4;

configs = {'iolin_pd', 0.6; 'clfqp', 0.6; 'clfqp_con', 0.6; 'clfqp_con', 1.0};
logln(logf, '=== ch3_dt_sweep | %s | peak ff %.1f Nm | mu_s %.2f', ...
      datestr(now), peak_ff, p0.limits.mu_s); %#ok<TNOW1,DATST>

for e = eps_list
    for dt = dt_list
        for c = 1:size(configs, 1)
            [name, frac] = configs{c, :};
            if dt == 0 && strcmp(name, 'clfqp_con'), continue; end
            f = fullfile(SPD, sprintf('%s_box%03.0f_eps%03.0f_dt%04.0fus.mat', ...
                                      name, 100*frac, 100*e, 1e6*dt));
            if exist(f, 'file')
                L = load(f); r = L.r;
            else
                pc = p0;
                pc.eps = e; pc.control_dt = dt;
                pc.controller = name;
                pc.limits.u_max = frac * peak_ff;
                pc.limits.enable.torque = strcmp(name, 'clfqp_con');
                t0 = tic;
                r = run_steps(X(:, 1), alpha, pc, n_steps);
                r.seconds = toc(t0);
                r.name = name; r.box = pc.limits.u_max; r.eps = e; r.dt = dt;
                save(f, 'r');
            end
            logln(logf, ['eps %.2f dt %6.1f us %-9s box %5.1f | steps %d valid %d | ' ...
                         'max|eta| %s | |eta-| %s | |eta+| %s | +/- %s | peak %.1f | %.0f s'], ...
                  r.eps, 1e6*r.dt, r.name, r.box, r.steps, r.valid, vec(r.eta_max), ...
                  vec(r.eta_minus), vec(r.eta_plus), vec(r.eta_plus ./ r.eta_minus), ...
                  r.peak, r.seconds);
        end
    end
end
logln(logf, '=== DT_SWEEP_DONE %s', datestr(now)); %#ok<TNOW1,DATST>
fprintf('DT_SWEEP_DONE\n');
end

% ---------------------------------------------------------------------------
function r = run_steps(x, alpha, p, n_steps)
r = struct('steps', 0, 'valid', 0, 'eta_max', [], 'eta_minus', [], ...
           'eta_plus', [], 'peak', 0, 'fell', false, 'why', '');
still_valid = true;
for k = 1:n_steps
    if ~all(isfinite(x)), r.fell = true; r.why = 'non-finite state'; break; end
    s = ch3_step(x, alpha, p);
    if p.control_dt == 0
        s = post_lambda(s, alpha, p);          % contact force under u(x)
    end
    eta = zeros(1, numel(s.t));
    for j = 1:numel(s.t)
        [y, yd] = ch3_outputs(s.x(:, j), alpha, p);
        eta(j) = norm([y; yd]);
    end
    [ym, ydm] = ch3_outputs(s.x_end, alpha, p);
    [yp, ydp] = ch3_outputs(s.x_next, alpha, p);
    r.eta_max(end+1)   = max(eta);
    r.eta_minus(end+1) = norm([ym; ydm]);
    r.eta_plus(end+1)  = norm([yp; ydp]);
    if ~isempty(s.u)
        r.peak = max(r.peak, max(abs(s.u(:))));
    else
        for j = 1:numel(s.t)
            r.peak = max(r.peak, max(abs(ch3_control(s.x(:, j), alpha, p))));
        end
    end
    fell = ~s.ok || s.L_step < p.step_len_min || ~all(isfinite(s.x_next)) || ~(max(eta) < 1e3);
    V = ch3_validity(struct('steps', s), p);
    if still_valid && ~fell && V.available && V.valid_steps == 1
        r.valid = r.valid + 1;
    else
        still_valid = false;
    end
    if fell
        r.fell = true;
        r.why = sprintf('step %d: ok %d, L %.3f, max|eta| %.2e', k, s.ok, s.L_step, max(eta));
        break;
    end
    r.steps = k;
    x = s.x_next;
end
end

% ---------------------------------------------------------------------------
function s = post_lambda(s, alpha, p)
%POST_LAMBDA  Contact force at every output point of a continuous-control step.
lam = zeros(2, numel(s.t));
for j = 1:numel(s.t)
    u = ch3_control(s.x(:, j), alpha, p);
    [~, ~, aux] = ch3_control_affine(s.x(:, j), p);
    lam(:, j) = aux.lam_drift + aux.lam_in * u;
end
s.lambda = lam;
s.u = zeros(p.nu, 0);          % no held torque in continuous time
end

function s = vec(v)
s = ['[' sprintf('%.3g ', v) ']'];
end

function logln(f, fmt, varargin)
s = sprintf(fmt, varargin{:});
fprintf('%s\n', s);
fid = fopen(f, 'a'); fprintf(fid, '%s\n', s); fclose(fid);
end
