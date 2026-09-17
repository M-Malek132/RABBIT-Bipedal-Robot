function ch4_long_reruns(stage)
%CH4_LONG_RERUNS  Rerun the Chapter-4 rate sweep and long-horizon tables at the
% single actuator rating (p.box.rule = 'rating'), scoring every run for
% physical validity (ch4_validity). Every run is saved as soon as it completes,
% so a crashed MATLAB session resumes where it stopped.
%
%   ch4_long_reruns(stage)
%
% Output: Results/reruns/ch4/ (one .mat per run, progress.log, summary.log).
% The 'robust' stage needs the rating-rule robust and Case IV results that
% ch4_main writes (ch4_rerun_all runs those first). The run lists reproduce the
% report's tables; the 'nine' runs and variants are the ones in
% CH4_UNCERTAINTY.md 5a.
%
%   stage: 'robust'  rclfqp_con past 25 steps (60; 120 at scale 3)  -> tab:plateau
%          'sweep'   predictor rate at x1.5, with and without normalization -> tab:rate
%          'nine'    the nine long runs x the variants still in the code -> tab:mitig, tab:window
%          'oos'     six fresh random-load sequences (seeds 4-9)          -> tab:oos
%          'summary' write summary.log from the saved runs

root = fileparts(fileparts(fileparts(fileparts(mfilename('fullpath')))));
RES  = [fullfile(root, 'Results') filesep];
SPD  = [fullfile(RES, 'reruns', 'ch4') filesep];
if ~exist(SPD, 'dir'), mkdir(SPD); end
logf = [SPD 'progress.log'];

switch stage
    case 'robust',  stage_robust(SPD, RES, logf);
    case 'sweep',   stage_sweep(SPD, logf);
    case 'nine',    stage_long(SPD, logf, nine_runs(), nine_variants());
    case 'oos',     stage_long(SPD, logf, oos_runs(), oos_variants());
    case 'summary', stage_summary(SPD);
    otherwise, error('ch4_long_reruns:stage', 'Unknown stage "%s" (robust|sweep|nine|oos|summary).', stage);
end
logline(logf, 'STAGE %s DONE', stage);
fprintf('STAGE_DONE %s\n', stage);
end

% ---------------------------------------------------------------------------
function stage_robust(SPD, RES, logf)
Rr = newest_result(RES, 'robust');
R4 = newest_result(RES, 'case4');
x0 = Rr.x0; alpha = Rr.alpha;
runs = {'rob100', Rr.p, 1, 60; 'rob070', Rr.p, 0.7, 60; ...
        'rob150', Rr.p, 1.5, 60; 'rob300', R4.p_case4, 3, 120};
for r = 1:size(runs, 1)
    [tag, p, s, N] = runs{r, :};
    outf = [SPD tag '.mat'];
    if exist(outf, 'file'), continue; end
    if s == 3, box = p.box.rating_case4; else, box = p.box.rating; end
    pc = ch4_run_params(p, 'rclfqp_con', s, box);
    t0 = tic;
    sim = ch4_simulate(x0, alpha, pc, N);
    M = metrics(sim, alpha, pc, false); %#ok<NASGU>
    info = struct('tag', tag, 'ctrl', 'rclfqp_con', 'scale', s, 'box', box, 'N', N, ...
                  'd1', p.rclf.delta1_max, 'd2', p.rclf.delta2_max, ...
                  'source', Rr.file); %#ok<NASGU>
    save(outf, 'info', 'M');
    logline(logf, '%s: %d/%d (%s) box %.0f [%.0fs]', tag, sim.n_ok, N, sim.reason, box, toc(t0));
end
end

% ---------------------------------------------------------------------------
function stage_sweep(SPD, logf)
[x0, alpha, p] = ch4_load_gait();
for nrm = [0.75 0]
    for a = [500 632 700 800 900]
        f = sprintf('%ssweep_n%03.0f_a%04.0f.mat', SPD, nrm*100, a);
        if exist(f, 'file'), continue; end
        pa = p; pa.l1.predictor_rate = a; pa.l1.normalized_rate = nrm;
        C = ch4_compare_controllers(x0, alpha, pa, 'l1', struct('scales', 1.5, ...
                'controllers', {{'l1', 'l1_con'}}, 'n_steps', 25, ...
                'store_traj', false, 'verbose', false));
        rows = arrayfun(@(e) struct('name', e.name, 'steps', e.steps_completed, ...
                                    'max_eta', e.max_eta, 'Fz_min', e.Fz_min, ...
                                    'box', e.u_box, 'reason', e.reason, ...
                                    'valid', e.valid_steps, 'mu_max', e.mu_max), C); %#ok<NASGU>
        save(f, 'a', 'nrm', 'rows');
        for e = C(:).'
            logline(logf, 'sweep nrm %.2f a %d %-6s %d steps max|eta| %.2f box %.0f %s', ...
                    nrm, a, e.name, e.steps_completed, e.max_eta, e.u_box, e.reason);
        end
    end
end
end

% ---------------------------------------------------------------------------
function stage_long(SPD, logf, base, variants)
[x0, alpha, p] = ch4_load_gait();
box = p.box.rating;
for iv = 1:size(variants, 1)
    [vtag, fset, only] = variants{iv, :};
    for c = 1:size(base, 1)
        [tag, ctrl, s, mL, rngL, sd, N] = base{c, :};
        if ~isempty(only) && ~any(strcmp(tag, only)), continue; end
        name = sprintf('%s_%s', vtag, tag);
        outf = [SPD name '.mat'];
        if exist(outf, 'file'), continue; end
        pc = ch4_run_params(p, ctrl, s, box);
        pc.uncertainty.load_mass = mL;
        pc.load_random_range     = rngL;
        pc.load_seed             = sd;
        pc = fset(pc);
        t0 = tic;
        sim = ch4_simulate(x0, alpha, pc, N);
        M = metrics(sim, alpha, pc, true); %#ok<NASGU>
        info = struct('name', name, 'variant', vtag, 'tag', tag, 'ctrl', ctrl, ...
                      'scale', s, 'load', mL, 'seed', sd, 'box', box, 'N', N); %#ok<NASGU>
        save(outf, 'info', 'M');
        logline(logf, '%s: %d/%d (%s) [%.0fs]', name, sim.n_ok, N, sim.reason, toc(t0));
    end
end
end

function base = nine_runs()
% tag, controller, scale, load, range, seed, steps
base = {
  's100',   'l1_con', 1,   0,  [],     11, 60
  's070',   'l1_con', 0.7, 0,  [],     11, 60
  's150',   'l1_con', 1.5, 0,  [],     11, 60
  's150l1', 'l1',     1.5, 0,  [],     11, 120
  'rnd11',  'l1_con', 1,   0,  [0 70], 11, 60
  'rnd1',   'l1_con', 1,   0,  [0 70], 1,  60
  'rnd2',   'l1_con', 1,   0,  [0 70], 2,  60
  'rnd3',   'l1_con', 1,   0,  [0 70], 3,  60
  'kg46',   'l1_con', 1,   46, [],     11, 60
};
end

function base = oos_runs()
base = cell(6, 7);
for k = 1:6
    base(k, :) = {sprintf('rnd%d', k+3), 'l1_con', 1, 0, [0 70], k+3, 60};
end
end

function V = nine_variants()
% tag, parameter change, runs it applies to ({} = all)
dt6 = {'s150', 's150l1', 'rnd11', 'rnd1', 'rnd2', 'rnd3'};
V = {
  'none',     @(q) setf(q, 'l1.normalized_rate', 0),                               {}
  'leak002',  @(q) setf(setf(q, 'l1.normalized_rate', 0), 'l1.alpha_leak', 2),         {}
  'leak010',  @(q) setf(setf(q, 'l1.normalized_rate', 0), 'l1.alpha_leak', 10),        {}
  'leak050',  @(q) setf(setf(q, 'l1.normalized_rate', 0), 'l1.alpha_leak', 50),        {}
  'continuous', @(q) setf(setf(q, 'l1.normalized_rate', 0), 'l1.impact_estimate', 'continuous'), {}
  'fold',     @(q) setf(setf(q, 'l1.normalized_rate', 0), 'l1.impact_estimate', 'fold'), {}
  'cap150',   @(q) setf(setf(q, 'l1.normalized_rate', 0), 'l1.alpha_regressor_rate', 1.5), {}
  'cap100',   @(q) setf(setf(q, 'l1.normalized_rate', 0), 'l1.alpha_regressor_rate', 1),   {}
  'cap050',   @(q) setf(setf(q, 'l1.normalized_rate', 0), 'l1.alpha_regressor_rate', 0.5), {}
  'dt0500us', @(q) setf(setf(q, 'l1.normalized_rate', 0), 'control_dt', 5e-4),           dt6
  'nrm010',   @(q) setf(q, 'l1.normalized_rate', 0.1),                             {}
  'nrm050',   @(q) setf(q, 'l1.normalized_rate', 0.5),                             {}
  'nrm075',   @(q) setf(q, 'l1.normalized_rate', 0.75),                            {}
  'nrm100',   @(q) setf(q, 'l1.normalized_rate', 1),                               {}
  'nrm150',   @(q) setf(q, 'l1.normalized_rate', 1.5),                             {}
  'nrm200',   @(q) setf(q, 'l1.normalized_rate', 2),                               {}
};
end

function V = oos_variants()
V = {
  'none',     @(q) setf(q, 'l1.normalized_rate', 0),                               {}
  'cap050',   @(q) setf(setf(q, 'l1.normalized_rate', 0), 'l1.alpha_regressor_rate', 0.5), {}
  'nrm075',   @(q) setf(q, 'l1.normalized_rate', 0.75),                            {}
  'nrm100',   @(q) setf(q, 'l1.normalized_rate', 1),                               {}
  'dt0500us', @(q) setf(setf(q, 'l1.normalized_rate', 0), 'control_dt', 5e-4),           {}
};
end

function q = setf(q, name, val)
parts = strsplit(name, '.');
q = setfield(q, parts{:}, val); %#ok<SFLD>
end

% ---------------------------------------------------------------------------
function M = metrics(sim, alpha, pc, is_l1)
% Per-step peaks of V and ||eta||, step duration/length/speed, loads, and for
% L1 the estimate sizes and the median alpha/beta cosine.
clf = ch3_res_clf(pc);
n = sim.n_ok;
M = struct('n', n, 'reason', sim.reason, 'T', [], 'L', [], 'loads', [], ...
           'vpk', nan(1, n), 'epk', nan(1, n), 'apk', nan(1, n), ...
           'bpk', nan(1, n), 'cosab', nan(1, n));
% physical validity at full solver resolution (ch4_step records the true force)
M.validity = ch4_validity(sim, pc);
if n > 0
    M.T = [sim.steps.T]; M.L = [sim.steps.L_step];
end
if isfield(sim, 'loads'), M.loads = sim.loads; end
for k = 1:n
    st = sim.steps(k);
    idx = 1:3:numel(st.t);
    Vs = zeros(1, numel(idx)); es = Vs;
    for j = 1:numel(idx)
        [~, ~, o] = ch3_outputs(st.x(:, idx(j)), alpha, pc);
        Vs(j) = ch3_clf_eval(o.eta, clf, pc.eps); es(j) = norm(o.eta);
    end
    M.vpk(k) = max(Vs); M.epk(k) = max(es);
    if is_l1 && isfield(st, 'xi') && ~isempty(st.xi)
        a = 0; b = 0; cs = [];
        for j = 1:size(st.xi, 2)
            ss = ch4_l1_state('unpack', pc, st.xi(:, j));
            na = norm(ss.alpha_hat); nb = norm(ss.beta_hat);
            a = max(a, na); b = max(b, nb);
            if na > 1 && nb > 1
                cs(end+1) = (ss.alpha_hat.' * ss.beta_hat) / (na * nb); %#ok<AGROW>
            end
        end
        M.apk(k) = a; M.bpk(k) = b;
        if ~isempty(cs), M.cosab(k) = median(cs); end
    end
end
end

% ---------------------------------------------------------------------------
function stage_summary(SPD)
fid = fopen([SPD 'summary.log'], 'w');
F = dir([SPD '*.mat']);
for f = F(:).'
    S = load([SPD f.name]);
    if isfield(S, 'rows')
        for r = S.rows(:).'
            fprintf(fid, 'SWEEP nrm %.2f a %4d %-6s %2d/25 valid %2d max|eta| %6.2f minFz %6.0f maxmu %.2f box %4.0f %s\n', ...
                    S.nrm, S.a, r.name, r.steps, r.valid, r.max_eta, r.Fz_min, r.mu_max, r.box, r.reason);
        end
        continue;
    end
    I = S.info; M = S.M; n = M.n;
    if n == 0
        fprintf(fid, '%-22s 0/%d %s\n', f.name, I.N, M.reason); continue;
    end
    blk = arrayfun(@(b0) max(M.epk(b0:min(b0+9, n))), 1:10:n);
    vbl = arrayfun(@(b0) mean(M.vpk(b0:min(b0+9, n))), 1:10:n);
    spd = arrayfun(@(b0) mean(M.L(b0:min(b0+9, n)) ./ M.T(b0:min(b0+9, n))), 1:10:n);
    Vd = M.validity;
    fprintf(fid, '%-22s VALID %3d/%3d first %s in step %g | bad %.3f%% | minFz %.0f maxmu %.2f\n', ...
            f.name, Vd.valid_steps, n, Vd.first_kind, Vd.first_step, 100 * Vd.frac_invalid, ...
            Vd.Fz_min, Vd.mu_max);
    fprintf(fid, ['%-22s %3d/%3d | median step max|eta| %.2f | max %.2f | |eta| per 10 %s | ' ...
                  'Vpk mean per 10 %s | speed per 10 %s | max|a| %.0f max|b| %.0f cos med %.2f | %s\n'], ...
            f.name, n, I.N, median(M.epk), max(M.epk), sprintf('%.2f ', blk), ...
            sprintf('%.3g ', vbl), sprintf('%.3f ', spd), max(M.apk), max(M.bpk), ...
            median(M.cosab, 'omitnan'), M.reason);
end
fprintf(fid, 'DONE\n');
fclose(fid);
end

% ---------------------------------------------------------------------------
function R = newest_result(RES, field)
% Newest ch4_result .mat with a non-empty FIELD, run under the rating rule.
F = dir([RES 'ch4_result_*.mat']);
[~, order] = sort({F.name});
for f = fliplr(F(order).')
    w = whos('-file', [RES f.name]);
    if ~any(strcmp({w.name}, field)), continue; end
    R = load([RES f.name]);
    if isempty(R.(field)), continue; end
    pp = R.p;
    if ~strcmp(ch4_box_rule(pp), 'rating'), continue; end
    R.file = f.name;
    return;
end
error('no rating-rule result with a %s sweep in %s', field, RES);
end

function logline(logf, fmt, varargin)
fid = fopen(logf, 'a');
fprintf(fid, [fmt '\n'], varargin{:});
fclose(fid);
end
