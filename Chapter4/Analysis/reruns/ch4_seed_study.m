function ch4_seed_study(stage, seeds)
%CH4_SEED_STUDY  Fall rates under random load with enough seeds to compare them.
%
%   ch4_seed_study()                   'runs' then 'summary', seeds 101-150
%   ch4_seed_study('summary')
%   ch4_seed_study('runs', seeds)
%
% WHY. Chapter 4 chose the normalization rate rho = 0.75 on nine runs and
% checked it on six fresh random-load sequences: 0 of 6 falls against 4 of 6
% without normalization. The comparison was right to make and too small to
% settle anything -- the two-sided Fisher exact p is 0.061, and the 95% Wilson
% intervals, [0, 0.39] and [0.30, 0.90], overlap. Fifty sequences per variant
% make the same question answerable.
%
% THE RUNS. l1_con, 1 kHz unless the variant says otherwise, a hip load redrawn
% every step from 0-70 kg (p.load_seed = the seed), 60 steps, the 556 Nm
% rating, contact validity scored (ch4_validity). Variants:
%   none      the first three fixes only (no normalization)
%   nrm075    normalized adaptation, rho = 0.75 (the chapter's choice)
%   dt0500us  no normalization, 0.5 ms control period
%   pwc       the piecewise-constant adaptation law (p.l1.adaptation)
% Seeds 101-150 are disjoint from every seed the chapter has used (1-11), so
% nothing here was selected on.
%
% THE SUMMARY. Per variant: falls (a run that ends before 60 steps), runs
% contact-valid for all 60 steps, and the median valid step count, with Wilson
% 95% intervals on both rates; then the two-sided Fisher exact p of every pair
% of variants, for falls and for fully valid runs.
%
% Output: Results/reruns/ch4_seeds/ (one .mat per run, progress.log,
% summary.log). Runtime: about 50 x (20 + 20 + 36 + 20) s, 80 minutes; runs
% resume where they stopped.
%
% See also CH4_LONG_RERUNS, CH4_WILSON, CH4_FISHER_EXACT.

if nargin < 1 || isempty(stage), stage = {'runs', 'summary'}; end
if ischar(stage), stage = {stage}; end
if nargin < 2 || isempty(seeds), seeds = 101:150; end

ROOT = fileparts(fileparts(fileparts(fileparts(mfilename('fullpath')))));
SPD  = fullfile(ROOT, 'Results', 'reruns', 'ch4_seeds');
if ~exist(SPD, 'dir'), mkdir(SPD); end
logf = fullfile(SPD, 'progress.log');
N_STEPS = 60;

variants = {
  'none',     @(q) setf(q, 'l1.normalized_rate', 0)
  'nrm075',   @(q) setf(q, 'l1.normalized_rate', 0.75)
  'dt0500us', @(q) setf(setf(q, 'l1.normalized_rate', 0), 'control_dt', 5e-4)
  'pwc',      @(q) setf(setf(q, 'l1.adaptation', 'pwc'), 'l1.pwc_rate', 0)
};

if any(strcmp(stage, 'runs'))
    [x0, alpha, p] = ch4_load_gait();
    for sd = seeds
        for iv = 1:size(variants, 1)
            [vtag, fset] = variants{iv, :};
            f = fullfile(SPD, sprintf('%s_s%03d.mat', vtag, sd));
            if exist(f, 'file'), continue; end
            pc = ch4_run_params(p, 'l1_con', 1, p.box.rating);
            pc.load_random_range = [0 70];
            pc.load_seed         = sd;
            pc = fset(pc);
            t0 = tic;
            try
                sim = ch4_simulate(x0, alpha, pc, N_STEPS);
                V   = ch4_validity(sim, pc);
                r   = struct('n', sim.n_ok, 'valid', V.valid_steps, ...
                             'reason', sim.reason, 'first_kind', V.first_kind);
            catch err
                r = struct('n', 0, 'valid', 0, 'reason', ['error: ' err.message], ...
                           'first_kind', '');
            end
            r.variant = vtag; r.seed = sd; r.N = N_STEPS; r.seconds = toc(t0);
            save(f, 'r');
            logline(logf, '%-9s seed %3d: %2d/%d valid %2d | %s [%.0f s]', vtag, sd, ...
                    r.n, N_STEPS, r.valid, r.reason, r.seconds);
        end
    end
end

if any(strcmp(stage, 'summary'))
    fid = fopen(fullfile(SPD, 'summary.log'), 'w');
    names = variants(:, 1).';
    fall = zeros(1, numel(names)); full = fall; nrun = fall; medv = nan(1, numel(names));
    for iv = 1:numel(names)
        F = dir(fullfile(SPD, sprintf('%s_s*.mat', names{iv})));
        n = []; v = []; NN = [];
        for k = 1:numel(F)
            L = load(fullfile(SPD, F(k).name));
            n(end+1) = L.r.n; v(end+1) = L.r.valid; NN(end+1) = L.r.N; %#ok<AGROW>
        end
        nrun(iv) = numel(n);
        fall(iv) = sum(n < NN);
        full(iv) = sum(v == NN);
        if ~isempty(v), medv(iv) = median(v); end
        [fl, fh] = ch4_wilson(fall(iv), nrun(iv));
        [vl, vh] = ch4_wilson(full(iv), nrun(iv));
        fprintf(fid, ['%-9s runs %2d | falls %2d (%.2f, 95%% [%.2f, %.2f]) | fully valid %2d ' ...
                      '(%.2f, [%.2f, %.2f]) | median valid steps %.1f\n'], names{iv}, nrun(iv), ...
                fall(iv), fall(iv) / max(nrun(iv), 1), fl, fh, full(iv), ...
                full(iv) / max(nrun(iv), 1), vl, vh, medv(iv));
    end
    fprintf(fid, 'pairwise two-sided Fisher exact p (falls | fully valid):\n');
    for i = 1:numel(names)
        for j = i+1:numel(names)
            if nrun(i) == 0 || nrun(j) == 0, continue; end
            fprintf(fid, '  %-9s vs %-9s  %.4f | %.4f\n', names{i}, names{j}, ...
                    ch4_fisher_exact(fall(i), nrun(i), fall(j), nrun(j)), ...
                    ch4_fisher_exact(full(i), nrun(i), full(j), nrun(j)));
        end
    end
    fprintf(fid, 'DONE\n');
    fclose(fid);
    type(fullfile(SPD, 'summary.log'));
end
fprintf('SEED_STUDY_DONE\n');
end

% ---------------------------------------------------------------------------
function q = setf(q, name, val)
parts = strsplit(name, '.');
q = setfield(q, parts{:}, val); %#ok<SFLD>
end

function logline(logf, fmt, varargin)
s = sprintf(fmt, varargin{:});
fprintf('%s\n', s);
fid = fopen(logf, 'a'); fprintf(fid, '%s\n', s); fclose(fid);
end
