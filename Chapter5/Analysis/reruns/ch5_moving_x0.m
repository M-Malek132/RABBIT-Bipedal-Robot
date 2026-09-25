function out = ch5_moving_x0(stage)
%CH5_MOVING_X0  Exercise Corollary 5.2 from a MOVING initial state.
%
%   ch5_moving_x0()              both plants
%   ch5_moving_x0('springmass')  the linear plant only
%   ch5_moving_x0('pendulum')    the nonlinear plant only
%
% WHY. Every Chapter 5 run starts at rest, so hdot(0) = 0, every lower bound of
% Corollary 5.2 is zero and each pole's margin is the pole itself (tab:adm).
% The corollary's dependence on x0 -- the limitation Section 5.2.3 closes with
% -- is therefore never exercised. This gives the same four configurations an
% initial velocity and measures three things:
%
%   1. WHERE THE CHAIN BREAKS. s_adm, the speed at which x0 leaves the chain
%      C_1..C_(rb-1) -- the part ch5_ecbf_admissible can check exactly, and the
%      part Theorem 5.1 needs, since y_rb >= 0 IS the barrier row the QP
%      enforces from the first sample. s_top is where the drift-evaluated y_rb
%      goes negative too (the QP is then active from t = 0). The reports' rule
%      of thumb, -hdot(0) <= p_1 h(0), is only the FIRST rung; it is logged as
%      s_1 so the two can be compared.
%   2. WHETHER LEAVING IT MATTERS. Closed-loop runs at multiples of s_adm with
%      the gate off, recording h_min. Corollary 5.2 is sufficient, not
%      necessary, so an inadmissible x0 need not violate; how far past s_adm
%      the barrier actually holds is the corollary's conservatism, measured.
%   3. THE REMEDY. At the inadmissible speeds, the smallest uniform scale c on
%      the poles that restores the chain, and a run with poles 1.05*c*p_b.
%
% THE DIRECTIONS. spring-mass: a rigid translation toward the boundary,
% xdot1 = xdot2 = xdot3 = s. The springs stay unstretched, so every derivative
% of h above the first is zero at x0 and the chain has a closed form,
%
%       y_i(x0) = e_i h0 - e_(i-1) s  >= 0   <=>   s <= h0 / sum_(j<=i) 1/p_j
%
% (e_i the elementary symmetric polynomials of p_1..p_i). The script checks the
% bisection against it, which tests ch5_ecbf_admissible as a side effect. With
% p_b = 0.12*(10:15) the first rung allows p_1 h0 = 3.60 m/s at h0 = 3 but the
% fifth only 0.85 -- the rule of thumb overstates the tolerance about 4x.
%
% pendulum: link 1 and its motor rotating together, thetadot1 = thetadotm1 = s,
% from the upright start. There hdot(0) = 0 for every s -- the end effector is
% at the top of its circle -- so the first rung never binds and the chain
% breaks, if it does, through the centripetal terms in the higher rungs.
%
% Output: Results/reruns/ch5_moving_x0/ -- one .mat per run, <tag>_summary.mat
% per configuration, and moving_x0.log (the file to read; stdout is buffered
% under -batch). Every run is skipped if its file exists, so a killed session
% resumes where it stopped. Prints CH5_MOVING_X0_DONE at the very end.

if nargin < 1 || isempty(stage), stage = {'springmass', 'pendulum'}; end
if ischar(stage), stage = {stage}; end

root = fileparts(fileparts(fileparts(fileparts(mfilename('fullpath')))));
odir = fullfile(root, 'Results', 'reruns', 'ch5_moving_x0');
if ~exist(odir, 'dir'), mkdir(odir); end
logf = fullfile(odir, 'moving_x0.log');

% The four configurations of tab:adm, and the scan range for each plant.
cfg = struct( ...
    'system', {'springmass', 'springmass', 'pendulum', 'pendulum'}, ...
    'level',  {3.15, 3.00, -1.0, -0.5}, ...
    'tag',    {'sm315', 'sm300', 'pd100', 'pd050'}, ...
    's_max',  {5, 5, 20, 20}, ...
    'unit',   {'m/s', 'm/s', 'rad/s', 'rad/s'});

out = struct('tag', {}, 'B', {}, 'runs', {});
for c = 1:numel(cfg)
    if ~any(strcmpi(stage, cfg(c).system)), continue; end
    out(end+1) = run_config(cfg(c), odir, logf); %#ok<AGROW>
end

logln(logf, 'CH5_MOVING_X0_DONE %s | %s', strjoin(stage, ','), datestr(now));
end

% ===========================================================================
function R = run_config(cf, odir, logf)
p = ch5_params('system', cf.system, 'controller', 'ecbfclfqp', ...
               'constraint', cf.level);
% The gate is what is under test, so it is evaluated here, per run, and the
% simulator is told not to act on it.
p.ecbf.admissibility = 'off';

x_rest = ch5_x0(p);
d      = direction(p);
b0     = ch5_barrier(x_rest, p);
poles0 = p.ecbf.poles(:).';

logln(logf, '=== %s | %s, level %.2f, h0 %.3f, rb %d, poles %s | %s', cf.tag, ...
      cf.system, cf.level, b0.h, b0.rb, mat2str(poles0, 4), datestr(now));

%% ------------------------------------------------ 1. where the chain breaks
B = boundary(p, x_rest, d, cf.s_max);
logln(logf, ['  boundary: first rung s_1 %s | exact chain s_adm %s | with the ' ...
             'drift rung s_top %s  (%s)'], fmt_s(B.s_1), fmt_s(B.s_adm), ...
      fmt_s(B.s_top), cf.unit);
logln(logf, '  just past s_adm: %s', B.msg_past);

if strcmpi(cf.system, 'springmass')
    s_cf = b0.h ./ cumsum(1 ./ poles0);       % rung i allows s <= s_cf(i)
    B.s_closed = s_cf;
    dev = max(abs([B.s_1 - s_cf(1), B.s_adm - s_cf(b0.rb-1), ...
                   B.s_top - s_cf(b0.rb)]));
    logln(logf, '  closed form h0/sum(1/p): %s  -> bisection agrees to %.1e', ...
          mat2str(s_cf, 4), dev);
end

if ~isfinite(B.s_adm)
    logln(logf, '  the chain holds over the whole scan (s <= %g); nothing to run', ...
          cf.s_max);
    R = struct('tag', cf.tag, 'B', B, 'runs', struct([]));
    save(fullfile(odir, [cf.tag '_summary.mat']), 'R');
    return;
end

%% --------------------------------------- 2. closed loop at multiples of s_adm
mult   = [0 0.5 0.9 1.1 1.5 2 3];
remedy = [1.1 1.5 3];
runs   = struct([]);
for k = 1:numel(mult)
    s  = mult(k) * B.s_adm;
    x0 = x_rest + s * d;
    f  = fullfile(odir, sprintf('%s_x%03.0f.mat', cf.tag, 100*mult(k)));
    r  = one_run(p, x0, f);
    r.mult = mult(k); r.s = s; r.c = 1;
    report(logf, r, cf.unit);
    runs = [runs, r]; %#ok<AGROW>

    %% ------------------------------------------------------- 3. the remedy
    if any(abs(mult(k) - remedy) < 1e-12) && ~r.adm_ok
        c = pole_scale(p, x0, poles0);
        if ~isfinite(c)
            logln(logf, '    remedy: no uniform pole scale up to 1e4 restores the chain');
            continue;
        end
        q  = p;
        q.ecbf.poles = 1.05 * c * poles0;
        f  = fullfile(odir, sprintf('%s_x%03.0f_fast.mat', cf.tag, 100*mult(k)));
        r  = one_run(q, x0, f);
        r.mult = mult(k); r.s = s; r.c = 1.05 * c;
        logln(logf, '    remedy: chain restored at c = %.3f; run with 1.05c = %.3f', ...
              c, 1.05 * c);
        report(logf, r, cf.unit);
        runs = [runs, r]; %#ok<AGROW>
    end
end

R = struct('tag', cf.tag, 'B', B, 'runs', runs);
save(fullfile(odir, [cf.tag '_summary.mat']), 'R');
end

% ---------------------------------------------------------------------------
function d = direction(p)
%DIRECTION  Unit initial velocity, toward the boundary.
switch lower(p.system)
    case 'springmass'
        d = [0 0 0 1 1 1].';          % rigid translation: springs unstretched
    case 'pendulum'
        d = [0 0 0 0 1 0 1 0].';      % thetadot1 = thetadotm1: spring 1 unstretched
    otherwise
        error('ch5_moving_x0:system', 'Unknown system "%s".', p.system);
end
end

% ---------------------------------------------------------------------------
function B = boundary(p, x_rest, d, s_max)
%BOUNDARY  First speed at which each condition fails: grid scan, then bisection.
b = ch5_barrier(x_rest, p);
e = ch5_ecbf_gain(p, b.rb);
adm_at = @(s) ch5_ecbf_admissible(x_rest + s * d, p, e);

tests = struct( ...
    'name', {'s_1', 's_adm', 's_top'}, ...
    'fn',   {@(a) a.y(2) >= 0, @(a) a.ok, @(a) a.ok && a.y(end) >= 0});

grid = linspace(0, s_max, 401);
A = arrayfun(adm_at, grid, 'UniformOutput', false);
A = [A{:}];

B = struct('s_max', s_max);
for t = 1:numel(tests)
    pass = arrayfun(tests(t).fn, A);
    k = find(~pass, 1);
    if isempty(k)
        B.(tests(t).name) = inf;
        continue;
    end
    if k == 1
        B.(tests(t).name) = 0;       % fails at rest: the configuration itself is bad
        continue;
    end
    lo = grid(k-1); hi = grid(k);
    for it = 1:60
        mid = (lo + hi) / 2;
        if tests(t).fn(adm_at(mid)), lo = mid; else, hi = mid; end
    end
    B.(tests(t).name) = lo;
end

B.msg_past = '(chain holds over the scan)';
if isfinite(B.s_adm)
    a = adm_at(B.s_adm * (1 + 1e-6) + 1e-9);
    B.msg_past = a.msg;
end
end

% ---------------------------------------------------------------------------
function c = pole_scale(p, x0, poles0)
%POLE_SCALE  Smallest uniform scale c >= 1 on the poles admitting x0.
if chain_ok(p, x0, poles0), c = 1; return; end
lo = 1; hi = 2;
while ~chain_ok(p, x0, hi * poles0)
    lo = hi; hi = 2 * hi;
    if hi > 1e4, c = NaN; return; end
end
for it = 1:60
    mid = (lo + hi) / 2;
    if chain_ok(p, x0, mid * poles0), hi = mid; else, lo = mid; end
end
c = hi;
end

function ok = chain_ok(p, x0, poles)
p.ecbf.poles = poles;
b  = ch5_barrier(x0, p);
a  = ch5_ecbf_admissible(x0, p, ch5_ecbf_gain(p, b.rb));
ok = a.ok;
end

% ---------------------------------------------------------------------------
function r = one_run(p, x0, f)
%ONE_RUN  One closed-loop run from x0, checkpointed to f.
if exist(f, 'file')
    S = load(f, 'r');
    r = S.r;
    r.stored = true;
    return;
end

b   = ch5_barrier(x0, p);
adm = ch5_ecbf_admissible(x0, p, ch5_ecbf_gain(p, b.rb));

r = struct('adm_ok', adm.ok, 'adm_slack', min(adm.slack(1:end-1)), ...
           'adm_msg', adm.msg, 'y_top', adm.y(end), 'poles', p.ecbf.poles(:).', ...
           'ran', false, 'ok', false, 'reason', '', 'h_min', NaN, ...
           'h_min_t', NaN, 'violated', NaN, 'u_max', NaN, 'pct_act', NaN, ...
           'feas', false, 't_frac', NaN, 'cpu', NaN, 'stored', false);

% A diverging run makes ode45 shrink its step without bound; the budget turns
% that into a recorded result instead of a stalled session (as ch5_ensemble).
budget = 300;
t0 = tic;
if strcmpi(p.integrator, 'ode45')
    p.ode_opts.OutputFcn = @(t, y, flag) deadline_stop(t0, budget);
end
try
    sim = ch5_simulate(p, 'x0', x0);
    n = max(sim.n, 1);
    r.ran      = true;
    r.ok       = sim.ok;
    r.reason   = sim.reason;
    r.h_min    = sim.h_min;
    r.h_min_t  = sim.h_min_t;
    r.violated = sim.h_min < 0;
    r.u_max    = max(abs(sim.u(:, 1:n)), [], 'all');
    r.pct_act  = 100 * mean(sim.cbf_active(1:n));
    r.feas     = all(sim.feasible(1:n));
    r.t_frac   = sim.n / numel(sim.t);
catch ME
    r.reason = ME.message;
end
r.cpu = toc(t0);
save(f, 'r');
end

% ---------------------------------------------------------------------------
function report(logf, r, unit)
stored = '';
if r.stored, stored = ' (stored)'; end
logln(logf, ['  x%-4.2f s %7.4f %-5s c %5.3f | adm %d slack %+9.3g y_rb %+9.3g | ' ...
             'h_min %+.4e at t %6.2f | violated %d | max|u| %8.2f | act %5.1f%% | ' ...
             'feas %d | done %.2f | %s | %.0f s%s'], ...
      r.mult, r.s, unit, r.c, r.adm_ok, r.adm_slack, r.y_top, r.h_min, ...
      r.h_min_t, r.violated, r.u_max, r.pct_act, r.feas, r.t_frac, ...
      r.reason, r.cpu, stored);
end

function s = fmt_s(v)
if isinf(v), s = 'none'; else, s = sprintf('%.4f', v); end
end

function status = deadline_stop(tstart, budget)
if toc(tstart) > budget
    error('ch5_moving_x0:budget', 'run exceeded its %g s wall-clock budget', budget);
end
status = 0;
end

function logln(f, fmt, varargin)
s = sprintf(fmt, varargin{:});
fprintf('%s\n', s);
fid = fopen(f, 'a'); fprintf(fid, '%s\n', s); fclose(fid);
end
