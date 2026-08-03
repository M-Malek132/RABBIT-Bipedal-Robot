function out = ch5_main(varargin)
%CH5_MAIN  Run the Chapter-5 study end to end.
%
%   out = ch5_main()
%   out = ch5_main('studies', {'pendulum'}, 'animate', false)
%
% Chapter 5 has no optimization and no robot. It validates one construction --
% the Exponential Control Barrier Function -- on two plants chosen so that
% between them they cover what the construction claims:
%
%   (1) SERIAL SPRING-MASS      linear, 1 input, RELATIVE DEGREE 6   Fig 5.3
%   (2) 2-LINK PENDULUM         nonlinear, 2 inputs, REL. DEGREE 4   Fig 5.4
%       with elastic actuators
%
% Each study runs the same three controllers over the same trajectory:
%
%   a. CLF-QP                stability only -- violates the constraint
%   b. ECBF-CLF-QP           at the chapter's first constraint level
%   c. ECBF-CLF-QP           at a tighter level, poles UNCHANGED
%
% and (c) is the one that carries the argument. The poles are held fixed
% between (b) and (c) precisely so that any difference in the response is
% attributable to the constraint rather than to retuning -- which is what lets
% Fig. 5.3's caption claim that peak forces and response speed are unchanged.
%
% ---------------------------------------- the third study, which is not a figure
% 'relative_degree' runs the Section 5.1 reciprocal CBF on both a
% relative-degree-1 constraint and a relative-degree-6 one, on the SAME plant.
% It works on the first and cannot even be written down on the second. That
% contrast is the chapter's entire motivation and it is the one claim the two
% figures above cannot make, since both of their constraints are high relative
% degree and Section 5.1 simply has nothing to say about either.
%
% Options (name/value)
%   'studies'  cell of {'springmass','pendulum','relative_degree'} (default all)
%   'plot'     draw and save figures (default true)
%   'save'     write a .mat of everything (default true)
%   'animate'  write the pendulum comparison GIF (default true)
%   'fast'     rk4 instead of ode45, for a quick pass (default false)
%   'resume'   path to the result .mat: names the output file, and loads it if
%              it already exists, so a run can be split across invocations
%   'only'     indices of the runs to attempt within each study (default all).
%              Runs already present in the resumed file are never redone.
%   any ch5_params field, applied to EVERY study, including dotted names
%
% Output
%   out : struct .springmass .pendulum .relative_degree .figs .gif .file
%
% See also CH5_PARAMS, CH5_SIMULATE, CH5_REPORT, CH5_PLOT_SPRINGMASS.

own = {'studies','plot','save','animate','fast','resume','only'};
o   = struct('studies', {{'springmass','pendulum','relative_degree'}}, ...
             'plot', true, 'save', true, 'animate', true, 'fast', false, ...
             'resume', '', 'only', []);

pv = {};
for k = 1:2:numel(varargin)
    if any(strcmpi(varargin{k}, own))
        o.(lower(varargin{k})) = varargin{k+1};
    else
        pv(end+1:end+2) = varargin(k:k+1); %#ok<AGROW>
    end
end
if o.fast
    pv(end+1:end+2) = {'integrator', 'rk4'};
end

results_dir = fullfile(fileparts(mfilename('fullpath')), '..', 'Results');
if ~exist(results_dir, 'dir'), mkdir(results_dir); end
stamp   = datestr(now, 'yyyy-mm-dd_HH-MM-SS'); %#ok<TNOW1,DATST>
out_dir = fullfile(results_dir, sprintf('ch5_%s', stamp));

out = struct('springmass', [], 'pendulum', [], 'relative_degree', [], ...
             'figs', gobjects(0), 'gif', '', 'file', '', 'stamp', stamp);

if o.save
    out.file = fullfile(results_dir, sprintf('ch5_result_%s.mat', stamp));
end

%% ---------------------------------------------------------------- resume
% Load a previous run and extend it, so the studies can be done in separate
% invocations. Combined with the per-study checkpointing below, this makes the
% whole chapter restartable at study granularity.
%
% That is not a nicety on this machine. MATLAB here is killed at random by a
% heap race in its bundled libcurl (the add-on registry, on a background
% thread), and the rate got bad enough that six consecutive attempts at the
% full run died partway through the pendulum study. Splitting a 6-minute run
% into three 1-minute ones is the difference between finishing and not.
%
%   ch5_main('studies', {'springmass'})
%   ch5_main('studies', {'pendulum'}, 'resume', <that file>, 'plot', false)
%   ch5_main('studies', {}, 'resume', <that file>)          figures only
% 'resume' names the output file AND loads it if it is already there, so the
% same argument works for the first invocation of a split run and for every
% one after it.
if ~isempty(o.resume)
    out.file = o.resume;

    % The stamp keys the figure directory, so take it from the filename to keep
    % a split run's outputs together.
    tok = regexp(o.resume, 'ch5_result_(.+)\.mat$', 'tokens', 'once');
    if ~isempty(tok), stamp = tok{1}; end

    if exist(o.resume, 'file')
        R = load(o.resume);
        for f = {'springmass','pendulum','relative_degree'}
            if isfield(R, f{1}) && ~isempty(R.(f{1}))
                out.(f{1}) = R.(f{1});
            end
        end
        if isfield(R, 'stamp'), stamp = R.stamp; end
        fprintf(' resuming %s (have: %s)\n', o.resume, strjoin(have(out), ', '));
    else
        fprintf(' starting %s\n', o.resume);
    end

    out.stamp = stamp;
    out_dir   = fullfile(results_dir, sprintf('ch5_%s', stamp));
end

fprintf('\n================ CHAPTER 5 ================\n');
fprintf(' Exponential Control Barrier Functions\n');

%% =================================================== (1) serial spring mass
if any(strcmpi('springmass', o.studies))
    fprintf('\n--- Fig 5.3: serial spring-mass, relative degree 6 ---\n');

    specs = { ...
        'clfqp',     3.15, 'CLF-QP (no barrier)', 'x3max'; ...
        'ecbfclfqp', 3.15, 'ECBF-CLF-QP, x_3^{max}-x_{3d} = 15 cm', 'x3max'; ...
        'ecbfclfqp', 3.00, 'ECBF-CLF-QP, x_3^{max}-x_{3d} = 0 cm',  'x3max'};

    out.springmass = run_specs('springmass', specs, out.springmass, o, pv);
    if all_done(out.springmass, specs)
        ch5_report(out.springmass, 'title', ...
                   'Serial spring-mass (relative degree 6)');
    end
    checkpoint(out, o);
end

%% ============================================ (2) pendulum, elastic actuators
if any(strcmpi('pendulum', o.studies))
    fprintf('\n--- Fig 5.4: 2-link pendulum with elastic actuators, rel. degree 4 ---\n');

    specs = { ...
        'clfqp',     -1.0, 'CLF-QP (no barrier)',          'p2min'; ...
        'ecbfclfqp', -1.0, 'ECBF-CLF-QP, p_{2min} = -1.0 m', 'p2min'; ...
        'ecbfclfqp', -0.5, 'ECBF-CLF-QP, p_{2min} = -0.5 m', 'p2min'};

    out.pendulum = run_specs('pendulum', specs, out.pendulum, o, pv);
    if all_done(out.pendulum, specs)
        ch5_report(out.pendulum, 'title', ...
                   '2-link pendulum with elastic actuators (relative degree 4)');
    end
    checkpoint(out, o);
end

%% ======================== (3) why Section 5.1 is not enough, on one plant
if any(strcmpi('relative_degree', o.studies))
    fprintf('\n--- Section 5.1 vs 5.2: the same plant, two relative degrees ---\n');
    out.relative_degree = ch5_relative_degree_study(pv);
    checkpoint(out, o);
end

%% -------------------------------------------------------------------- plots
if o.plot
    figs = gobjects(0);
    if ~isempty(out.springmass)
        figs = [figs, ch5_plot_springmass(out.springmass, out_dir)];
    end
    if ~isempty(out.pendulum)
        figs = [figs, ch5_plot_pendulum(out.pendulum, out_dir)];
    end
    out.figs = figs;
    if ~isempty(figs), fprintf('\n Figures: %s\n', out_dir); end
end

%% ------------------------------------------- save BEFORE the animation
% Same ordering, and the same reason, as ch4_main: the sweeps above cost
% minutes and the GIF below renders hundreds of frames through getframe, which
% can take the whole MATLAB process down rather than raising something
% catchable. Make the numbers durable first.
%
% checkpoint() has already written after each study, for a reason specific to
% this machine: MATLAB here aborts at random inside a background thread of its
% bundled libcurl (the add-on registry), several minutes into a run and with no
% warning. A single save at the end lost a completed spring-mass study to a
% crash during the pendulum one. Saving per study makes each result durable as
% soon as it exists.
if o.save
    checkpoint(out, o);
    fprintf(' Saved: %s\n', out.file);
end

%% ---------------------------------------------------------------- animation
if o.animate && ~isempty(out.pendulum)
    if ~exist(out_dir, 'dir'), mkdir(out_dir); end
    gif = fullfile(out_dir, 'ch5_pendulum_compare.gif');
    try
        out.gif = ch5_animate_pendulum(out.pendulum, gif);
        fprintf(' GIF: %s\n', out.gif);
        if o.save
            gifpath = out.gif; %#ok<NASGU>
            save(out.file, 'gifpath', '-append');
        end
    catch err
        fprintf(' ch5_animate_pendulum skipped: %s\n', err.message);
    end
end

fprintf('==========================================\n');
fprintf('CH5_MAIN_DONE\n\n');

end

% ---------------------------------------------------------------------------
function runs = run_specs(system, specs, runs, o, pv)
%RUN_SPECS  Simulate the specs not already present, keeping the table order.
%
% Each row is {controller, constraint level, label, level-name}. A run already
% in `runs` (matched by name) is kept as-is, so an invocation can add one row
% at a time and the table fills in across invocations.
%
% o.only restricts which rows are attempted -- ch5_main('only', 3) runs the
% third alone. That granularity exists for the same reason the resume option
% does: on this machine MATLAB is killed at random within a couple of minutes,
% and a study of three 30-second runs is far more likely to finish as three
% invocations than as one.

if isempty(runs)
    runs = struct('name', {}, 'label', {}, 'sim', {});
end

want = o.only;
if isempty(want), want = 1:size(specs,1); end

for i = 1:size(specs,1)
    nm = spec_name(specs, i);

    if any(strcmp({runs.name}, nm))
        continue;                       % already have it
    end
    if ~ismember(i, want)
        continue;                       % not requested this time
    end

    p = ch5_params('system', system, 'controller', specs{i,1}, ...
                   'constraint', specs{i,2}, pv{:});
    fprintf('  %-11s %s = %+.2f ... ', specs{i,1}, specs{i,4}, specs{i,2});
    s = ch5_simulate(p);
    fprintf('h_min = %+.3e (%.0f s)\n', s.h_min, s.cpu);

    runs(end+1) = struct('name', nm, 'label', specs{i,3}, 'sim', s); %#ok<AGROW>
end

% Restore the table's order regardless of the order they were produced in --
% the figures put one column per spec and the columns must not shuffle.
order = zeros(1, numel(runs));
for k = 1:numel(runs)
    for i = 1:size(specs,1)
        if strcmp(runs(k).name, spec_name(specs, i)), order(k) = i; end
    end
end
[~, ix] = sort(order);
runs = runs(ix);

end

function nm = spec_name(specs, i)
nm = sprintf('%s %s=%+.2f', specs{i,1}, specs{i,4}, specs{i,2});
end

function tf = all_done(runs, specs)
tf = ~isempty(runs) && numel(runs) == size(specs,1);
end

% ---------------------------------------------------------------------------
function names = have(out)
%HAVE  Which studies the loaded result already contains.
names = {};
for f = {'springmass','pendulum','relative_degree'}
    if ~isempty(out.(f{1})), names{end+1} = f{1}; end %#ok<AGROW>
end
if isempty(names), names = {'nothing'}; end
end

% ---------------------------------------------------------------------------
function checkpoint(out, o)
%CHECKPOINT  Write everything finished so far to the run's .mat.
%
% Cheap relative to the sweeps, and it means an interrupted run still leaves
% behind the studies that did complete.
if ~o.save || isempty(out.file), return; end
R = rmfield(out, 'figs'); %#ok<NASGU>
save(out.file, '-struct', 'R');
end

% ---------------------------------------------------------------------------
function R = ch5_relative_degree_study(pv)
%CH5_RELATIVE_DEGREE_STUDY  The motivating comparison, on the springmass.
%
% Two constraints on the SAME plant, from the same initial condition:
%
%   xdot1 <= 0.8      relative degree 1  -- Section 5.1 applies
%   x3    <= 3.00     relative degree 6  -- Section 5.1 does not
%
% Holding the plant fixed and varying only the relative degree is what makes
% this a controlled comparison. If the reciprocal CBF failed on a different
% system one could blame the system.

R = struct('rd1', [], 'rd6', [], 'rd1_base', []);

% --- relative degree 1: the reciprocal CBF of Section 5.1 works
c1 = struct('type', 'v1_max', 'value', 0.8);

p = ch5_params('system','springmass','controller','cbfclfqp', ...
               'constraint', c1, 'ecbf.poles', 2.0, pv{:});
R.rd1 = ch5_simulate(p);

p.controller = 'clfqp';
R.rd1_base = ch5_simulate(p);

fprintf('  rel. degree 1  (xdot1 <= 0.8)\n');
fprintf('    CLF-QP        max xdot1 = %.4f   h_min = %+.4f\n', ...
        max(R.rd1_base.x(4,:)), R.rd1_base.h_min);
fprintf('    CBF-CLF-QP    max xdot1 = %.4f   h_min = %+.4f   <- Section 5.1 works\n', ...
        max(R.rd1.x(4,:)), R.rd1.h_min);

% --- relative degree 6: the same controller has no barrier row to write
q = ch5_params('system','springmass','controller','cbfclfqp', ...
               'constraint', 3.00, pv{:});
b = ch5_barrier(ch5_x0(q), q);
R.rd6 = ch5_simulate(q);

fprintf('  rel. degree 6  (x3 <= 3.00)\n');
fprintf('    L_g h = %s  ->  the barrier row contains no control at all\n', ...
        mat2str(b.Lgh));
fprintf('    CBF-CLF-QP    max x3    = %.4f   h_min = %+.4f   <- violated\n', ...
        max(R.rd6.x(3,:)), R.rd6.h_min);

% BOTH failure modes are possible and which one occurs is a property of the
% run, not something to assert in advance. With L_g h = 0 the row reads
% 0 >= -gamma h^3 - L_f h: vacuously true while the right-hand side is
% negative, and unsatisfiable the moment it turns positive -- which happens as
% h -> 0 with the cart still moving, since then gamma h^3 no longer covers
% L_f h. So the controller is inert for most of the run and then infeasible
% near the boundary. Report the measured split rather than describing one.
n_inf = sum(~R.rd6.feasible);
fprintf('    QP infeasible on %d/%d samples (%.1f%%)\n', ...
        n_inf, R.rd6.n, 100*n_inf/R.rd6.n);
if n_inf == 0
    fprintf(['    The row was VACUOUS throughout: the QP always solved and the\n' ...
             '    controller silently reduced to the CLF-QP. The quiet failure.\n']);
else
    fprintf(['    Both failure modes appeared: the row is vacuous while\n' ...
             '    -gamma h^3 - L_f h <= 0, then UNSATISFIABLE once the cart\n' ...
             '    approaches the limit with speed. The chapter''s "high control\n' ...
             '    inputs or become infeasible", in that order.\n']);
end
fprintf('    Either way the constraint is not enforced. This is what\n');
fprintf('    Section 5.2 exists to fix.\n');

end
