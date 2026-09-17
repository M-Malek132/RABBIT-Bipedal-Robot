function ch4_eps020_diagnostics(stage)
%CH4_EPS020_DIAGNOSTICS  The three Chapter 4 diagnostics the report still
% quotes at eps = 0.35, rerun at the chapter's eps = 0.20 under the actuator
% rating, scored for physical validity.
%
%   ch4_eps020_diagnostics()             every stage, in order
%   ch4_eps020_diagnostics('kappa')      one stage
%
% Stages
%   'kappa'      tab:kappa. rclfqp_con, 3 steps, boundary layer
%                kappa in {0 .33 .66 .99 1.32 3.30} x Cases I-III (scale 1,
%                1.5, 0.7): median per-sample change of the held torque
%                (inf-norm), max ||eta|| over the recorded points, valid steps.
%   'growth'     "growth and plateau". rclfqp_con at scale 1, 25 steps: the
%                default law, then with the boundary layer off, the contact rows
%                off, the box off, D2 = 0 and D1 = 0. Per-step peak V, whether
%                it grows, whether the run falls, valid steps. (At eps 0.35 the
%                growth felled Case I in step 21, survived every ablation but
%                D1 = 0.)
%   'predictor'  the thesis predictor's theta_hat overshoot. l1 with
%                p.l1.predictor = 'thesis', no normalization, 25 steps, scales
%                0.7 and 1.5 (and 1 as reference): max ||theta_hat|| / max
%                ||theta_true||, and the median of the pointwise ratio. The
%                plant predictor (the default) alongside for contrast.
%   'summary'    Results/reruns/ch4_eps020/summary.log from the saved runs.
%
% Every run is saved on its own and skipped when present, so a crashed session
% resumes. ONE MATLAB SESSION AT A TIME (see ch4_rerun_all).
%
% See also CH4_RERUN_ALL, CH4_LONG_RERUNS, CH4_RUN_ENTRY.

if nargin < 1 || isempty(stage)
    for s = {'kappa', 'growth', 'predictor', 'summary'}
        ch4_eps020_diagnostics(s{1});
    end
    return;
end

ROOT = [fileparts(fileparts(fileparts(fileparts(mfilename('fullpath'))))) filesep];
SPD  = [fullfile(ROOT, 'Results', 'reruns', 'ch4_eps020') filesep];
if ~exist(SPD, 'dir'), mkdir(SPD); end
logf = [SPD 'progress.log'];

[x0, alpha, p] = ch4_load_gait();
if abs(p.eps - 0.20) > 1e-12
    error('ch4_eps020_diagnostics:eps', 'ch4_params has eps %.2f, expected 0.20.', p.eps);
end
logline(logf, 'stage %s start %s | eps %.2f | rating %.0f', stage, datestr(now), ...
        p.eps, p.box.rating);

switch stage
    case 'kappa',     stage_kappa(SPD, logf, x0, alpha, p);
    case 'growth',    stage_growth(SPD, logf, x0, alpha, p);
    case 'predictor', stage_predictor(SPD, logf, x0, alpha, p);
    case 'summary',   stage_summary(SPD);
    otherwise
        error('ch4_eps020_diagnostics:stage', ...
              'Unknown stage "%s" (kappa|growth|predictor|summary).', stage);
end
logline(logf, 'STAGE %s DONE', stage);
fprintf('STAGE_DONE %s\n', stage);
end

% ---------------------------------------------------------------------------
function stage_kappa(SPD, logf, x0, alpha, p)
for kap = [0 0.33 0.66 0.99 1.32 3.30]
    for s = [1 1.5 0.7]
        f = sprintf('%skappa_k%03.0f_s%03.0f.mat', SPD, 100*kap, 100*s);
        if exist(f, 'file'), continue; end
        pc = ch4_run_params(p, 'rclfqp_con', s, p.box.rating);
        pc.rclf.boundary_layer = kap;
        R = run_one(x0, alpha, pc, 3, false); %#ok<NASGU>
        info = struct('stage', 'kappa', 'kappa', kap, 'scale', s, 'N', 3, ...
                      'box', p.box.rating); %#ok<NASGU>
        save(f, 'info', 'R');
        logline(logf, 'kappa %.2f scale %.2f: %d/3 valid %d | med|du| %.1f Nm | max|eta| %.3f | %s', ...
                kap, s, R.n, R.valid, R.du_med, R.max_eta, R.reason);
    end
end
end

% ---------------------------------------------------------------------------
function stage_growth(SPD, logf, x0, alpha, p)
N = 25;
V = {
  'default',  @(q) q
  'nolayer',  @(q) setf(q, 'rclf.boundary_layer', 0)
  'norows',   @(q) setf(setf(q, 'limits.enable.friction', false), 'limits.enable.grf', false)
  'nobox',    @(q) setf(q, 'limits.enable.torque', false)
  'd2zero',   @(q) setf(q, 'rclf.delta2_max', 0)
  'd1zero',   @(q) setf(q, 'rclf.delta1_max', 0)
};
for iv = 1:size(V, 1)
    [tag, fset] = V{iv, :};
    f = sprintf('%sgrowth_%s.mat', SPD, tag);
    if exist(f, 'file'), continue; end
    pc = fset(ch4_run_params(p, 'rclfqp_con', 1, p.box.rating));
    R = run_one(x0, alpha, pc, N, false); %#ok<NASGU>
    info = struct('stage', 'growth', 'variant', tag, 'scale', 1, 'N', N, ...
                  'box', p.box.rating); %#ok<NASGU>
    save(f, 'info', 'R');
    logline(logf, 'growth %-8s: %d/%d valid %d | Vpk step1 %.3g last %.3g | max|eta| %.3f | %s', ...
            tag, R.n, N, R.valid, first(R.vpk), last(R.vpk), R.max_eta, R.reason);
end
end

% ---------------------------------------------------------------------------
function stage_predictor(SPD, logf, x0, alpha, p)
N = 25;
for pred = {'thesis', 'plant'}
    for s = [0.7 1.5 1]
        f = sprintf('%spredictor_%s_s%03.0f.mat', SPD, pred{1}, 100*s);
        if exist(f, 'file'), continue; end
        pc = ch4_run_params(p, 'l1', s, p.box.rating);
        pc.l1.predictor = pred{1};
        pc.l1.normalized_rate = 0;
        R = run_one(x0, alpha, pc, N, true); %#ok<NASGU>
        info = struct('stage', 'predictor', 'predictor', pred{1}, 'scale', s, ...
                      'N', N, 'box', p.box.rating); %#ok<NASGU>
        save(f, 'info', 'R');
        logline(logf, ['predictor %-6s scale %.2f: %d/%d valid %d | max|theta_hat| %.0f ' ...
                       'max|theta| %.0f ratio %.2f median ratio %.2f | %s'], ...
                pred{1}, s, R.n, N, R.valid, R.th_hat_max, R.th_true_max, ...
                R.th_ratio, R.th_ratio_med, R.reason);
    end
end
end

% ---------------------------------------------------------------------------
function R = run_one(x0, alpha, pc, N, is_l1)
t0 = tic;
sim = ch4_simulate(x0, alpha, pc, N);
Vd  = ch4_validity(sim, pc);
R = struct('n', sim.n_ok, 'reason', sim.reason, 'valid', Vd.valid_steps, ...
           'validity', Vd, 'du_med', NaN, 'du_med_step', nan(1, sim.n_ok), ...
           'peak_u', NaN, 'max_eta', NaN, 'vpk', nan(1, sim.n_ok), ...
           'epk', nan(1, sim.n_ok), 'th_hat_max', NaN, 'th_true_max', NaN, ...
           'th_ratio', NaN, 'th_ratio_med', NaN, 'Fz_min', NaN, 'secs', NaN);
if sim.n_ok == 0, R.secs = toc(t0); return; end

% held torque: change from one control period to the next, within each step
du = [];
for k = 1:sim.n_ok
    U = sim.steps(k).u;
    d = max(abs(diff(U, 1, 2)), [], 1);
    if ~isempty(d), R.du_med_step(k) = median(d); end
    du = [du, d]; %#ok<AGROW>
end
R.du_med = median(du);
R.peak_u = max(arrayfun(@(s) max(abs(s.u(:))), sim.steps(1:sim.n_ok)));

% V, eta and (L1) theta per step, from ch4_forces at the run-entry density
th_hat = []; th_true = [];
t_off = 0;
for k = 1:sim.n_ok
    st = sim.steps(k);
    F = ch4_forces(t_off + st.t, st.x, alpha, pc, st.xi, t_off + st.t_xi, 700);
    t_off = t_off + st.T;
    en = vecnorm([F.y; F.ydot], 2, 1);
    R.epk(k) = max(en);
    R.vpk(k) = max(F.V);
    R.Fz_min = min([R.Fz_min, F.Fz_min], [], 'omitnan');
    if is_l1
        th_hat  = [th_hat,  vecnorm(F.theta_hat, 2, 1)]; %#ok<AGROW>
        th_true = [th_true, vecnorm(F.theta_true, 2, 1)]; %#ok<AGROW>
    end
end
R.max_eta = max(R.epk);
if is_l1 && ~isempty(th_hat)
    R.th_hat_max  = max(th_hat);
    R.th_true_max = max(th_true);
    % With a perfect model (mass scale 1) the true theta is zero to rounding,
    % and any ratio to it is noise over noise -- it came out as 1e13. The
    % overshoot this stage measures only means something under a perturbation,
    % so leave the ratios NaN when there is no model error to estimate.
    if R.th_true_max > 1
        R.th_ratio = R.th_hat_max / R.th_true_max;
        big = th_true > 0.1 * R.th_true_max;      % skip near-zero denominators
        if any(big), R.th_ratio_med = median(th_hat(big) ./ th_true(big)); end
    end
end
R.secs = toc(t0);
end

% ---------------------------------------------------------------------------
function stage_summary(SPD)
fid = fopen([SPD 'summary.log'], 'w');
F = dir([SPD '*.mat']);
for f = F(:).'
    S = load([SPD f.name]); I = S.info; R = S.R;
    switch I.stage
        case 'kappa'
            fprintf(fid, 'KAPPA %.2f scale %.2f | %d/3 valid %d | med|du| %7.1f Nm | max|eta| %.3f | peak|u| %.0f | minFz %.0f | %s\n', ...
                    I.kappa, I.scale, R.n, R.valid, R.du_med, R.max_eta, R.peak_u, R.Fz_min, R.reason);
        case 'growth'
            fprintf(fid, 'GROWTH %-8s | %2d/%d valid %2d | Vpk per step %s| max|eta| %.3f | minFz %.0f | %s\n', ...
                    I.variant, R.n, I.N, R.valid, sprintf('%.3g ', R.vpk), R.max_eta, R.Fz_min, R.reason);
        case 'predictor'
            fprintf(fid, 'PRED %-6s scale %.2f | %2d/%d valid %2d | max|th_hat| %.0f max|th| %.0f ratio %.2f med %.2f | max|eta| %.3f | %s\n', ...
                    I.predictor, I.scale, R.n, I.N, R.valid, R.th_hat_max, R.th_true_max, ...
                    R.th_ratio, R.th_ratio_med, R.max_eta, R.reason);
    end
end
fprintf(fid, 'DONE\n');
fclose(fid);
end

% ---------------------------------------------------------------------------
function q = setf(q, name, val)
parts = strsplit(name, '.');
q = setfield(q, parts{:}, val); %#ok<SFLD>
end

function v = first(x)
if isempty(x), v = NaN; else, v = x(1); end
end

function v = last(x)
if isempty(x), v = NaN; else, v = x(end); end
end

function logline(logf, fmt, varargin)
fid = fopen(logf, 'a');
fprintf(fid, [fmt '\n'], varargin{:});
fclose(fid);
end
