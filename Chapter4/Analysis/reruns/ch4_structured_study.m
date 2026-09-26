function ch4_structured_study(stage, n_models)
%CH4_STRUCTURED_STUDY  The Chapter-4 controllers against a structured prior set.
%
%   ch4_structured_study()               'bounds', then 'runs', then 'summary'
%   ch4_structured_study(stage)          one of them
%   ch4_structured_study('runs', n)      n test draws (default 20)
%
% The chapter's mass-scale study tests the easiest uncertainty there is (a
% uniform scale: an input-gain error plus a proportional drift, with the zero
% dynamics unchanged), and its robust law is given bounds measured from the
% very perturbation it faces. This replaces both halves.
%
%   'bounds'   D1, D2 over the PRIOR set (ch4_uncertainty_set's defaults:
%              per-link mass 10%, COM 2 cm, inertia 20%, reflected rotor
%              inertia to 0.2 kg m^2, joint friction, torque bias), from 200
%              draws with seed 1 plus the two corners, along the gait's nodes
%              and a jittered cloud, times 1.2 -- fixed before any run.
%   'runs'     n FRESH draws (seed 2, never used for the bounds). On each, the
%              baselines (clfqp, clfqp_con), the robust law with the prior
%              bounds (scalar Delta2 as the chapter runs it, and the 'matrix'
%              form, which is the literal worst case for a non-scalar Delta2),
%              and L1 (the gradient law at nrm 0.75, and the piecewise-constant
%              law), 25 steps each, all constrained laws at the 556 Nm rating.
%              Two variants: 'model' (link and actuator errors only) and
%              'full' (plus measurement noise and a one-sample actuation
%              delay). Contact validity is scored as everywhere else.
%   'summary'  per controller and variant: runs completed, runs valid for all
%              25 steps (with Wilson 95% intervals), median valid steps and
%              median max||eta||.
%
% Output: Results/reruns/ch4_structured/ (bounds.mat, one .mat per run,
% progress.log, summary.log). Runtime: 'bounds' a few minutes; 'runs' about
% n x 12 runs x 10-20 s (40-80 minutes at n = 20).
%
% See also CH4_UNCERTAINTY_SET, CH4_DELTA_BOUNDS_PRIOR, CH4_LINK_DYNAMICS.

if nargin < 1 || isempty(stage), stage = {'bounds', 'runs', 'summary'}; end
if ischar(stage), stage = {stage}; end
if nargin < 2 || isempty(n_models), n_models = 20; end

ROOT = fileparts(fileparts(fileparts(fileparts(mfilename('fullpath')))));
SPD  = fullfile(ROOT, 'Results', 'reruns', 'ch4_structured');
if ~exist(SPD, 'dir'), mkdir(SPD); end
logf = fullfile(SPD, 'progress.log');

[x0, alpha, p, meta] = ch4_load_gait();
S  = load(meta.file);
pg = ch3_upgrade_params(S.p);
X  = ch3_col_unpack(S.z, ch3_col_effective_params(S.z, pg));
spec = ch4_uncertainty_set();

%% ---------------------------------------------------------------- bounds
fb = fullfile(SPD, 'bounds.mat');
if any(strcmp(stage, 'bounds')) && ~exist(fb, 'file')
    t0 = tic;
    Ub = ch4_uncertainty_set(200, 1, spec);
    Bo = ch4_delta_bounds_prior(X, alpha, p, Ub, struct('spec', spec));
    Bj = ch4_delta_bounds_prior(X, alpha, p, Ub(1:40), ...
                                struct('spec', spec, 'n_jitter', 3, 'jitter', 0.05));
    B = struct('orbit', Bo, 'cloud', Bj, 'spec', spec); %#ok<NASGU>
    save(fb, 'B');
    logline(logf, ['BOUNDS orbit: max|Delta1| %.1f  max|Delta2| %.4f (scalar part %.4f) ' ...
                   '-> D1 %.1f D2 %.4f | cloud sd 0.05: max|Delta1| %.1f max|Delta2| %.4f ' ...
                   '| feasible %d | %.0f s'], Bo.n1_max, Bo.n2_max, Bo.n2_scalar_max, ...
            Bo.delta1_max, Bo.delta2_max, Bj.n1_max, Bj.n2_max, Bo.feasible, toc(t0));
end

%% ------------------------------------------------------------------ runs
if any(strcmp(stage, 'runs'))
    if ~exist(fb, 'file')
        error('ch4_structured_study:bounds', 'run the bounds stage first');
    end
    L  = load(fb);
    Bo = L.B.orbit;
    Ut = ch4_uncertainty_set(n_models, 2, spec);
    box = p.box.rating;
    ctrls = { ...
      'clfqp',        'clfqp',      @(q) q
      'clfqp_con',    'clfqp_con',  @(q) q
      'rclf_scalar',  'rclfqp_con', @(q) setf(setf(setf(q, 'rclf.delta1_max', Bo.delta1_max), ...
                                             'rclf.delta2_max', Bo.delta2_max), 'rclf.delta2_model', 'scalar')
      'rclf_matrix',  'rclfqp_con', @(q) setf(setf(setf(q, 'rclf.delta1_max', Bo.delta1_max), ...
                                             'rclf.delta2_max', Bo.delta2_max), 'rclf.delta2_model', 'matrix')
      'l1_nrm075',    'l1_con',     @(q) setf(q, 'l1.normalized_rate', 0.75)
      'l1_pwc',       'l1_con',     @(q) setf(setf(q, 'l1.adaptation', 'pwc'), 'l1.pwc_rate', 0)
    };
    variants = {'model', 'full'};
    for iv = 1:numel(variants)
        for m = 1:n_models
            u = Ut(m);
            if strcmp(variants{iv}, 'model')
                u.noise = struct('q', 0, 'dq', 0, 'seed', u.noise.seed);
                u.delay = 0;
            end
            for c = 1:size(ctrls, 1)
                [tag, ctrl, fset] = ctrls{c, :};
                f = fullfile(SPD, sprintf('%s_%s_m%02d.mat', variants{iv}, tag, m));
                if exist(f, 'file'), continue; end
                pc = ch4_run_params(p, ctrl, 1, box);
                pc.uncertainty = u;
                pc = fset(pc);
                t0 = tic;
                try
                    sim = ch4_simulate(x0, alpha, pc, 25);
                    V   = ch4_validity(sim, pc);
                    r   = summarize_run(sim, V, alpha, pc);
                catch err
                    r = struct('n', 0, 'reason', ['error: ' err.message], ...
                               'valid', 0, 'eta_max', NaN, 'u_peak', NaN, ...
                               'Fz_min', NaN, 'mu_max', NaN);
                end
                r.tag = tag; r.variant = variants{iv}; r.model = m; r.seconds = toc(t0);
                save(f, 'r');
                logline(logf, '%-6s %-12s m%02d: %2d/25 valid %2d | max|eta| %7.2f | peak|u| %6.1f | %s [%.0f s]', ...
                        r.variant, r.tag, m, r.n, r.valid, r.eta_max, r.u_peak, r.reason, r.seconds);
            end
        end
    end
end

%% --------------------------------------------------------------- summary
if any(strcmp(stage, 'summary'))
    F = dir(fullfile(SPD, '*_m*.mat'));
    R = struct('tag', {}, 'variant', {}, 'n', {}, 'valid', {}, 'eta_max', {});
    for k = 1:numel(F)
        L = load(fullfile(SPD, F(k).name));
        R(end+1) = struct('tag', L.r.tag, 'variant', L.r.variant, 'n', L.r.n, ...
                          'valid', L.r.valid, 'eta_max', L.r.eta_max); %#ok<AGROW>
    end
    fid = fopen(fullfile(SPD, 'summary.log'), 'w');
    if exist(fb, 'file')
        LB = load(fb);
        fprintf(fid, 'PRIOR BOUNDS D1 %.1f D2 %.4f (orbit, x%.1f); cloud max|Delta1| %.1f max|Delta2| %.4f\n', ...
                LB.B.orbit.delta1_max, LB.B.orbit.delta2_max, LB.B.orbit.safety, ...
                LB.B.cloud.n1_max, LB.B.cloud.n2_max);
    end
    keys = unique(strcat({R.variant}, '|', {R.tag}), 'stable');
    for k = 1:numel(keys)
        sel = R(strcmp(strcat({R.variant}, '|', {R.tag}), keys{k}));
        n  = numel(sel);
        nc = sum([sel.n] == 25);
        nv = sum([sel.valid] == 25);
        [lo_c, hi_c] = ch4_wilson(nc, n);
        [lo_v, hi_v] = ch4_wilson(nv, n);
        fprintf(fid, ['%-20s draws %2d | completed %2d [%.2f, %.2f] | valid 25/25 %2d [%.2f, %.2f] ' ...
                      '| median valid %4.1f | median max|eta| %6.2f\n'], keys{k}, n, nc, lo_c, hi_c, ...
                nv, lo_v, hi_v, median([sel.valid]), median([sel.eta_max], 'omitnan'));
    end
    fprintf(fid, 'DONE\n');
    fclose(fid);
end
fprintf('STRUCTURED_STUDY_DONE\n');
end

% ---------------------------------------------------------------------------
function r = summarize_run(sim, V, alpha, pc)
r = struct('n', sim.n_ok, 'reason', sim.reason, 'valid', V.valid_steps, ...
           'eta_max', NaN, 'u_peak', NaN, 'Fz_min', V.Fz_min, 'mu_max', V.mu_max);
if sim.n_ok > 0
    e = 0;
    for j = 1:3:size(sim.x, 2)
        [y, yd] = ch3_outputs(sim.x(:, j), alpha, pc);
        e = max(e, norm([y; yd]));
    end
    r.eta_max = e;
    r.u_peak  = max(abs([sim.steps.u]), [], 'all');
end
end

function q = setf(q, name, val)
parts = strsplit(name, '.');
q = setfield(q, parts{:}, val); %#ok<SFLD>
end

function logline(logf, fmt, varargin)
s = sprintf(fmt, varargin{:});
fprintf('%s\n', s);
fid = fopen(logf, 'a'); fprintf(fid, '%s\n', s); fclose(fid);
end
