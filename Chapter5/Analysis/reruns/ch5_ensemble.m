function out = ch5_ensemble(stage)
%CH5_ENSEMBLE  Replace Chapter 5's single-realization numbers with ensembles.
%
%   ch5_ensemble('pendulum')   ensemble over solvers x initial conditions
%   ch5_ensemble('dtscale')    O(dt) study on the LINEAR plant
%   ch5_ensemble()             both
%
% WHY. Section "the trajectory is not reproducible" measures rk4 vs ode45
% diverging by 5.6e-1 on the pendulum and then correctly says no claim may
% rest on one pendulum trajectory -- but tab:pendulum still reports p99|tau|,
% max|tau|, % active, |theta1(T)-pi| and the theta2 excursion, every one of
% which is a property of exactly one realization. This measures them over an
% ensemble so they can be reported as ranges.
%
% Likewise the O(dt) table has three points on the NONLINEAR plant and the
% third changes sign, so it does not show O(dt) at all. The scaling is a
% statement about sampling, not about the pendulum, so 'dtscale' runs it on
% the spring-mass, where rk4 and ode45 agree to 3.2e-14 and a trajectory IS
% reproducible.
%
% Outputs go to Results/reruns/ch5/ and each run is skipped if already saved.

if nargin < 1 || isempty(stage), stage = {'pendulum','dtscale'}; end
if ischar(stage), stage = {stage}; end

root = fileparts(fileparts(fileparts(fileparts(mfilename('fullpath')))));
odir = fullfile(root, 'Results', 'reruns', 'ch5');
if ~exist(odir, 'dir'), mkdir(odir); end
out = struct();

for i = 1:numel(stage)
    switch stage{i}
        case 'pendulum', out.pendulum = pendulum_ensemble(odir);
        case 'dtscale',  out.dtscale  = dt_scaling(odir);
        otherwise, error('ch5_ensemble:stage', 'Unknown stage "%s".', stage{i});
    end
    fprintf('CH5_ENSEMBLE_DONE %s\n', stage{i});
end
end

% ===========================================================================
function R = pendulum_ensemble(odir)
%PENDULUM_ENSEMBLE  Both solvers x perturbed x0, for the three controllers.
%
% The perturbation is applied to theta0 only, so every member still starts at
% rest with unstretched springs; it moves the arm off the exact unstable
% equilibrium by up to 2 degrees, which is the scale of a real initial-pose
% error and is far larger than the 5.6e-1 rad solver divergence is small.

f    = fullfile(odir, 'pendulum_ensemble.mat');
ckpt = fullfile(odir, 'pendulum_ckpt.mat');
if exist(f, 'file'), R = load_field(f, 'R'); fprintf(' already done: %s\n', f); return; end

rows_done = struct([]);
if exist(ckpt, 'file')
    rows_done = load_field(ckpt, 'rows');
    fprintf(' resuming from checkpoint: %d runs already done\n', numel(rows_done));
end

solvers = {'rk4', 'ode45'};
% 1 nominal + 8 perturbed poses, deterministic.
rng(5);
dth = [zeros(2,1), (2*pi/180) * (2*rand(2,8) - 1)];
runs = struct('controller', {'clfqp','ecbfclfqp','ecbfclfqp'}, 'level', {[], -1.0, -0.5});

rows = rows_done;
for r = 1:numel(runs)
    for s = 1:numel(solvers)
        for m = 1:size(dth,2)
            % Everything must go through ch5_params: it calls ch5_system with
            % the list of names the caller set, and a SECOND bare ch5_system(p)
            % would silently restore the pendulum's own defaults for
            % integrator and control_dt (ch5_system lines 184, 200) -- which
            % would run every "ode45" member as rk4.
            args = {'system', 'pendulum', 'controller', runs(r).controller, ...
                    'integrator', solvers{s}};
            if ~isempty(runs(r).level)
                args = [args, {'constraint', runs(r).level}]; %#ok<AGROW>
            end
            p = ch5_params(args{:});
            assert(strcmp(p.integrator, solvers{s}), 'integrator was overridden');

            % WALL-CLOCK BUDGET. A perturbed x0 whose trajectory diverges makes
            % ode45 shrink its step without bound: the run never returns, while
            % rk4 on the same member fails fast. A run that cannot finish in
            % BUDGET seconds is recorded as skipped rather than stalling the
            % ensemble. Completed runs here take 20-40 s, so the budget is
            % loose enough not to truncate a healthy run.
            budget = 180;
            deadline = tic;
            p.ode_opts.OutputFcn = @(t, y, flag) deadline_stop(deadline, budget);
            x0 = ch5_x0(p);
            x0(1:2) = x0(1:2) + dth(:,m);
            if ~isempty(rows) && any(arrayfun(@(q) strcmp(q.controller, runs(r).controller) ...
                    && isequal(q.level, runs(r).level) && strcmp(q.solver, solvers{s}) ...
                    && q.member == m, rows))
                continue;                % already in the checkpoint
            end
            fprintf('  %-9s %-5s level %s  member %d ... ', runs(r).controller, ...
                    solvers{s}, mat2str(runs(r).level), m);

            rows(end+1).controller = runs(r).controller; %#ok<AGROW>
            rows(end).level  = runs(r).level;
            rows(end).solver = solvers{s};
            rows(end).member = m;

            % Corollary 5.2 at THIS x0, before simulating. The poles were
            % chosen at the nominal pose; whether they remain admissible a
            % couple of degrees away is the chapter's own stated limitation,
            % so it is measured per member rather than assumed.
            rows(end).adm_ok = NaN; rows(end).adm_slack = NaN;
            if ~strcmpi(runs(r).controller, 'clfqp')
                adm = ch5_ecbf_admissible(x0, p);
                rows(end).adm_ok    = adm.ok;
                rows(end).adm_slack = min(adm.slack);
            end

            % An inadmissible x0 can diverge and hand the QP a non-finite
            % value. That is a result about the construction, not a crash to
            % hide, so it is recorded and the ensemble continues.
            t0 = tic;
            try
                sim = ch5_simulate(p, 'x0', x0);
                % A run that ABORTS early still returns an h_min, over a
                % truncated window. Reporting that next to a completed run's
                % h_min would compare different experiments, so completion is
                % recorded and the summary splits on it.
                rows(end).ran     = true;
                rows(end).ok      = sim.ok;
                rows(end).reason  = sim.reason;
                % sim.t is the FULL preallocated grid even when the loop
                % breaks early, so sim.t(end) is always T and says nothing
                % about how far the run got. sim.n counts the samples that
                % were actually logged; use that.
                rows(end).n_samp  = sim.n;
                rows(end).t_end   = sim.t(min(max(sim.n,1), numel(sim.t)));
                rows(end).t_frac  = sim.n / numel(sim.t);
                rows(end).h_min   = sim.h_min;
                rows(end).tau_p99 = prctile(max(abs(sim.u), [], 1), 99);
                rows(end).tau_max = max(abs(sim.u(:)));
                rows(end).pct_act = 100 * mean(sim.cbf_active);
                rows(end).th1_err = abs(sim.x(1,end) - pi);
                rows(end).th2_min = min(sim.x(2,:));
                rows(end).feas    = all(sim.feasible);
                rows(end).err     = '';
                fprintf('h_min %+.4f  th2_min %+.3f  adm %d  ok %d  t %.2f/%.0f  %.0fs\n', ...
                        sim.h_min, rows(end).th2_min, rows(end).adm_ok, ...
                        sim.ok, sim.t(end), p.T, toc(t0));
            catch ME
                if strcmp(ME.identifier, 'ch5_ensemble:budget') || ...
                        contains(ME.message, 'ch5_ensemble:budget')
                    rows(end).err = 'skipped:walltime';
                else
                    rows(end).err = ME.identifier;
                end
                [rows(end).h_min, rows(end).tau_p99, rows(end).tau_max, ...
                 rows(end).pct_act, rows(end).th1_err, rows(end).th2_min] = deal(NaN);
                rows(end).ran    = false;
                rows(end).ok     = false;
                rows(end).reason = rows(end).err;
                rows(end).n_samp = NaN;
                rows(end).t_end  = NaN;
                rows(end).t_frac = NaN;
                rows(end).feas   = false;
                fprintf('FAILED (adm %d, slack %+.3g): %s\n', ...
                        rows(end).adm_ok, rows(end).adm_slack, ME.message(1:min(60,end)));
            end
            rows(end).cpu = toc(t0);
            save(ckpt, 'rows');          % checkpoint: a kill costs one run
        end
    end
end

R = summarize(rows);
save(f, 'R', 'rows');
print_summary(R);
end

% ---------------------------------------------------------------------------
function R = summarize(rows)
%SUMMARIZE  Range (min..max) of each metric over the ensemble, per controller.
key = arrayfun(@(r) sprintf('%s|%s', r.controller, mat2str(r.level)), rows, ...
               'UniformOutput', false);
[u, ~, g] = unique(key, 'stable');
flds = {'h_min','tau_p99','tau_max','pct_act','th1_err','th2_min'};
R = struct([]);
for k = 1:numel(u)
    sel = rows(g == k);
    R(k).name = u{k};
    R(k).n    = numel(sel);
    R(k).n_ran        = sum([sel.ran]);
    R(k).n_complete   = sum([sel.ok] == 1);
    R(k).t_frac_lo    = min([sel.t_frac]);
    R(k).n_admissible = sum([sel.adm_ok] == 1);
    R(k).adm_slack_lo = min([sel.adm_slack]);
    % Metrics are summarized over runs that COMPLETED; an aborted run's
    % h_min describes a shorter experiment and is counted, not averaged in.
    ran = sel([sel.ok] == 1);
    R(k).all_feasible = ~isempty(ran) && all([ran.feas]);
    R(k).all_safe     = ~isempty(ran) && all([ran.h_min] >= 0);
    for j = 1:numel(flds)
        v = [ran.(flds{j})];
        if isempty(v), v = NaN; end
        R(k).([flds{j} '_lo']) = min(v);
        R(k).([flds{j} '_hi']) = max(v);
    end
end
end

function print_summary(R)
fprintf('\n--- pendulum ensemble (ranges over solvers x initial poses) ---\n');
for k = 1:numel(R)
    fprintf(['%-20s n=%2d ran=%2d complete=%2d adm=%2d (min slack %+.3g)\n' ...
             '    h_min [%+.4f, %+.4f]  p99|tau| [%.1f, %.1f]  max|tau| [%.1f, %.1f]\n' ...
             '    %%act [%.1f, %.1f]  |th1(T)-pi| [%.4f, %.4f]  min th2 [%+.3f, %+.3f]' ...
             '  safe=%d feas=%d\n'], ...
        R(k).name, R(k).n, R(k).n_ran, R(k).n_complete, R(k).n_admissible, R(k).adm_slack_lo, ...
        R(k).h_min_lo, R(k).h_min_hi, ...
        R(k).tau_p99_lo, R(k).tau_p99_hi, R(k).tau_max_lo, R(k).tau_max_hi, ...
        R(k).pct_act_lo, R(k).pct_act_hi, R(k).th1_err_lo, R(k).th1_err_hi, ...
        R(k).th2_min_lo, R(k).th2_min_hi, R(k).all_safe, R(k).all_feasible);
end
end

% ===========================================================================
function D = dt_scaling(odir)
%DT_SCALING  |min(h,0)| vs control period on the LINEAR spring-mass plant.
%
% The spring-mass is the right plant for this: rk4 and ode45 agree to 3.2e-14
% there, so a change in the excursion is the sampling and not a different
% trajectory. Six periods over a decade and a half give a log-log slope.

f = fullfile(odir, 'dt_scaling.mat');
if exist(f, 'file'), D = load_field(f, 'D'); fprintf(' already done: %s\n', f); return; end

dts = [4e-3 2e-3 1e-3 5e-4 2.5e-4 1.25e-4];
D = struct('dt', num2cell(dts), 'h_min', [], 'exc', [], 'hdot_max', [], 'bound', []);
for i = 1:numel(dts)
    p = ch5_params('system', 'springmass', 'controller', 'ecbfclfqp', ...
                   'integrator', 'rk4', 'control_dt', dts(i), 'constraint', 3.00);
    assert(p.control_dt == dts(i), 'control_dt was overridden');
    fprintf('  dt = %.3e ... ', dts(i));
    sim = ch5_simulate(p);
    hh = sim.h(~isnan(sim.h)); tt = sim.t(~isnan(sim.h));
    hdm = max(abs(diff(hh) ./ diff(tt)));
    D(i).h_min    = sim.h_min;
    D(i).exc      = abs(min(sim.h_min, 0));
    D(i).hdot_max = hdm;
    D(i).bound    = 10 * dts(i) * hdm;
    fprintf('h_min %+.4e  excursion %.4e  max|hdot| %.3f\n', ...
            sim.h_min, D(i).exc, hdm);
end

e = [D.exc]; d = [dts];
k = e > 0;
if nnz(k) >= 2
    slope = polyfit(log(d(k)), log(e(k)), 1);
    fprintf('\n--- dt scaling (spring-mass, linear) ---\n');
    fprintf('log-log slope of excursion vs dt over %d positive points: %.3f\n', ...
            nnz(k), slope(1));
else
    fprintf('\n--- dt scaling: %d of %d periods are strictly safe ---\n', ...
            nnz(~k), numel(d));
end
save(f, 'D');
end

% ---------------------------------------------------------------------------
function v = load_field(f, name)
S = load(f, name); v = S.(name);
end

% ---------------------------------------------------------------------------
function status = deadline_stop(tstart, budget)
%DEADLINE_STOP  ode45 OutputFcn that aborts a run past its wall-clock budget.
if toc(tstart) > budget
    error('ch5_ensemble:budget', ...
          'run exceeded its %g s wall-clock budget', budget);
end
status = 0;
end
