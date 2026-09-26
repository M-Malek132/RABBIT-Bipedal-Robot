function ch5_ensemble_probe(stage)
%CH5_ENSEMBLE_PROBE  Why the pendulum ECBF runs diverge, and how far the
% completed ones leave the set, measured rather than argued.
%
%   ch5_ensemble_probe()            'runs' then 'summary'
%   ch5_ensemble_probe('summary')   re-read the saved runs
%
% WHAT IS BEING TESTED. ch5_ensemble ran 36 ECBF members on the pendulum (two
% levels x two solvers x nine initial poses) and only 10 completed; 11 went
% non-finite, seven of them between 0.226 and 0.256 s. The baseline CLF-QP's
% one-sample torque spike -- the min-norm law dividing by ||LgV|| as it passes
% near zero -- sits at 0.224 s. Two explanations are open, and they predict
% different things:
%
%   (a) the CLF row's singularity (LgV -> 0) meets an ACTIVE barrier row. Then
%       ||LgV|| is at a minimum, ||mu|| or the slack delta at a maximum, and
%       the barrier row active, in the samples just before each failure --
%       whatever the integrator.
%   (b) the variable-step ode45 chasing the discontinuity of the held control.
%       Then the failures are ode45's: its internal step count per control
%       period explodes before them (13 of the 26 non-completions were runs
%       that exceeded the 180 s wall-clock budget, which is that signature),
%       and a FIXED-step RK4 at the control period does not fail there.
%
% So every ECBF member is rerun three ways -- the ensemble's rk4 (10 substeps)
% and ode45, plus rk4 with ONE step per control period -- with ch5_simulate's
% 'log' on: ||LgV||, psi, the QP exit flag, ||Lb||, ||mu||, delta, y_rb, the
% barrier row's activity and, for ode45, the steps per period. The baseline
% is rerun once per pose to place its own LgV minimum in time.
%
% THE DIMENSIONAL BOUND (the report's |min h| <= c dt max|hdot|, c = 10,
% applied so far only to the spring-mass): for every COMPLETED ECBF member,
% the excursion |min(h, 0)|, max|hdot| over the run (hdot = eta_b(2)), and the
% effective c = |min h| / (dt max|hdot|). The report's -0.0017 is then either
% inside the sampling envelope or not -- a number, not "of that order".
%
% The poses are ch5_ensemble's exactly: rng(5), then 2 x 8 perturbations of up
% to 2 degrees on theta(0), member 1 unperturbed.
%
% Output: Results/reruns/ch5_probe/ (one .mat per run, keeping the last 400
% logged samples of a failed run; probe.log; summary.log).
% Runtime: 54 ECBF runs + 9 baseline runs; completed runs take 20-60 s, a
% failing ode45 run up to the same 180 s budget: about an hour.
%
% See also CH5_ENSEMBLE, CH5_SIMULATE, CH5_CTRL_ECBF_CLF_QP.

if nargin < 1 || isempty(stage), stage = {'runs', 'summary'}; end
if ischar(stage), stage = {stage}; end

ROOT = fileparts(fileparts(fileparts(fileparts(mfilename('fullpath')))));
OUT  = fullfile(ROOT, 'Results', 'reruns', 'ch5_probe');
if ~exist(OUT, 'dir'), mkdir(OUT); end
logf = fullfile(OUT, 'probe.log');

% ch5_ensemble's poses, bit for bit
rng(5);
dth = [zeros(2, 1), (2*pi/180) * (2*rand(2, 8) - 1)];
BUDGET = 180;
KEEP   = 400;                                  % samples kept before a failure

solvers = struct('name', {'rk4', 'ode45', 'rk4_1'}, ...
                 'integrator', {'rk4', 'ode45', 'rk4'}, 'substeps', {10, 10, 1});
runs = struct('controller', {'clfqp', 'ecbfclfqp', 'ecbfclfqp'}, 'level', {[], -1.0, -0.5});

if any(strcmp(stage, 'runs'))
    for r = 1:numel(runs)
        for s = 1:numel(solvers)
            % the baseline only to locate its LgV minimum: one solver is enough
            if isempty(runs(r).level) && ~strcmp(solvers(s).name, 'rk4'), continue; end
            for m = 1:size(dth, 2)
                tag = sprintf('%s_%s_%s_m%d', runs(r).controller, lvl(runs(r).level), ...
                              solvers(s).name, m);
                f = fullfile(OUT, [tag '.mat']);
                if exist(f, 'file'), continue; end

                args = {'system', 'pendulum', 'controller', runs(r).controller, ...
                        'integrator', solvers(s).integrator, 'n_substeps', solvers(s).substeps};
                if ~isempty(runs(r).level)
                    args = [args, {'constraint', runs(r).level}]; %#ok<AGROW>
                end
                p = ch5_params(args{:});
                assert(strcmp(p.integrator, solvers(s).integrator) && ...
                       p.n_substeps == solvers(s).substeps, 'solver settings were overridden');
                deadline = tic;
                p.ode_opts.OutputFcn = @(t, y, flag) deadline_stop(deadline, BUDGET);
                x0 = ch5_x0(p);
                x0(1:2) = x0(1:2) + dth(:, m);

                R = struct('tag', tag, 'controller', runs(r).controller, ...
                           'level', runs(r).level, 'solver', solvers(s).name, ...
                           'member', m, 'ok', false, 'reason', '', 't_end', NaN, ...
                           'h_min', NaN, 'hdot_max', NaN, 'dt', p.control_dt, ...
                           'c_eff', NaN, 'within_c10', NaN, 'win', [], ...
                           't_LgV_min', NaN, 'LgV_min', NaN, 'cpu', NaN);
                t0 = tic;
                try
                    sim = ch5_simulate(p, 'x0', x0, 'log', true);
                    R = fill_from_sim(R, sim, KEEP);
                catch err
                    R.reason = err.message;
                    if contains(err.message, 'budget'), R.reason = 'skipped: wall-clock budget'; end
                end
                R.cpu = toc(t0);
                save(f, 'R');
                logline(logf, ['%-34s ok %d | t_end %6.3f | h_min %+.4f | c_eff %7.2f | ' ...
                               'LgV min %.2e at %.3f s | %s [%.0f s]'], tag, R.ok, R.t_end, ...
                        R.h_min, R.c_eff, R.LgV_min, R.t_LgV_min, R.reason, R.cpu);
            end
        end
    end
end

if any(strcmp(stage, 'summary'))
    F = dir(fullfile(OUT, '*_m*.mat'));
    A = [];
    for k = 1:numel(F)
        L = load(fullfile(OUT, F(k).name));
        A = [A, L.R]; %#ok<AGROW>
    end
    fid = fopen(fullfile(OUT, 'summary.log'), 'w');
    base = A(strcmp({A.controller}, 'clfqp'));
    if ~isempty(base)
        fprintf(fid, 'BASELINE LgV minimum over [0, 0.4] s, per pose: %s s\n', ...
                sprintf('%.3f ', [base.t_LgV_min]));
    end
    for lv = [-1.0 -0.5]
        for s = {'rk4', 'ode45', 'rk4_1'}
            sel = A(strcmp({A.solver}, s{1}) & arrayfun(@(a) isequal(a.level, lv), A));
            if isempty(sel), continue; end
            done = sel([sel.ok]);
            fail = sel(~[sel.ok]);
            fprintf(fid, ['ECBF %+.1f %-6s: %d of %d complete | failure times %s | ' ...
                          'worst completed h_min %+.4f, c_eff max %.2f, within c = 10: %d of %d\n'], ...
                    lv, s{1}, numel(done), numel(sel), sprintf('%.3f ', [fail.t_end]), ...
                    min([done.h_min, Inf]), max([done.c_eff, -Inf]), ...
                    sum([done.within_c10] == 1), numel(done));
            for q = fail(:).'
                w = q.win;
                if isempty(w)
                    fprintf(fid, '    m%d: %s\n', q.member, q.reason);
                    continue;
                end
                tl = max(1, numel(w.t) - 19):numel(w.t);
                fprintf(fid, ['    m%d t_end %.3f: last %d samples -- min ||LgV|| %.2e, ' ...
                              'max ||mu|| %.2e, max delta %.2e, barrier active %d, ' ...
                              'QP exit flags %s, max integrator steps/period %s | %s\n'], ...
                        q.member, q.t_end, numel(tl), min(w.LgV(tl)), max(w.mu(tl)), ...
                        max(w.delta(tl)), sum(w.act(tl)), ...
                        mat2str(unique(w.exitflag(tl))), num2str(max(w.ode_steps)), ...
                        q.reason);
            end
        end
    end
    fprintf(fid, 'DONE\n');
    fclose(fid);
    type(fullfile(OUT, 'summary.log'));
end
fprintf('CH5_PROBE_DONE\n');
end

% ---------------------------------------------------------------------------
function R = fill_from_sim(R, sim, KEEP)
n = sim.n;
R.ok     = sim.ok;
R.reason = sim.reason;
R.t_end  = sim.t(max(n, 1));
R.h_min  = sim.h_min;
if size(sim.eta_b, 1) >= 2
    hd = sim.eta_b(2, 1:n);
    R.hdot_max = max(abs(hd(isfinite(hd))));
end
if sim.ok && isfinite(R.hdot_max) && R.hdot_max > 0
    exc = abs(min(sim.h_min, 0));
    R.c_eff      = exc / (R.dt * R.hdot_max);
    R.within_c10 = exc <= 10 * R.dt * R.hdot_max;
end
lg = sim.log;
t  = sim.t(1:n);
early = t <= 0.4;
if any(early) && any(isfinite(lg.LgV(early)))
    [R.LgV_min, i] = min(lg.LgV(early));
    R.t_LgV_min = t(i);
end
if ~sim.ok
    k = max(1, n - KEEP + 1):n;
    R.win = struct('t', t(k), 'h', sim.h(k), 'LgV', lg.LgV(k), 'psi', lg.psi(k), ...
                   'exitflag', lg.exitflag(k), 'Lb', lg.Lb(k), 'mu', lg.mu(k), ...
                   'delta', sim.delta(k), 'y_rb', sim.y_rb(k), ...
                   'act', sim.cbf_active(k), 'ode_steps', lg.ode_steps(k));
end
end

function s = lvl(v)
if isempty(v), s = 'base'; else, s = sprintf('p%03.0f', -100 * v); end
end

function status = deadline_stop(tstart, budget)
if toc(tstart) > budget
    error('ch5_probe:budget', 'run exceeded its %g s wall-clock budget', budget);
end
status = 0;
end

function logline(logf, fmt, varargin)
s = sprintf(fmt, varargin{:});
fprintf('%s\n', s);
fid = fopen(logf, 'a'); fprintf(fid, '%s\n', s); fclose(fid);
end
