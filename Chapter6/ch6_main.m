function out = ch6_main(varargin)
%CH6_MAIN  Run the Chapter-6 studies end to end.
%
%   out = ch6_main()
%   out = ch6_main('studies', {'stones'}, 'plot', false)
%
% Chapter 6 is Chapter 3's controller with barrier rows added to its QP, so
% unlike Chapter 5 it has no plants of its own -- it runs on RABBIT and on the
% Chapter-3 reference gait. What it adds is a set of position constraints and
% the machinery to enforce them.
%
% The studies, in the chapter's order:
%
%   'obstacle'  Section 6.1.2, Fig. 6.4 -- avoid a low ceiling (6.2), then an
%               overhead obstacle at a specific location (6.5)
%   'stones'    Sections 6.1.3-6.1.4, Figs. 6.5-6.7 -- walk a set of randomly
%               generated discrete footholds
%   'moving'    Section 6.2, Fig. 6.10 -- stepping stones that move with time
%   'width'     Section 6.3, Figs. 6.13-6.15 -- step-length and step-width
%               constraints together, on the swing-foot surrogate (DURUS is not
%               modelled here; see ch6_foot3d)
%   'range'     Section 6.1.4's headline claim -- the range of step length one
%               nominal gait plus a CBF can deliver, measured rather than quoted
%   'table'     Section 6.4, Table 6.1 -- the three controllers of (6.27),
%               reduced grid (see ch6_table61)
%
% Options (name/value)
%   'studies'  cell of the names above          (default all but 'table')
%   'plot'     draw and save figures            (default true)
%   'save'     write a .mat of everything       (default true)
%   'params'   a ch6_params struct to start from
%
% 'table' is out of the default set because it is the expensive one: even the
% reduced grid is a few hundred multi-step walking simulations. Ask for it
% explicitly, and read ch6_table61's header on what a reduced grid does and
% does not say.
%
% Output
%   out : struct with one field per study, plus .p and .dir
%
% See also CH6_PARAMS, CH6_SIMULATE, CH6_REPORT, CH6_TABLE61, CH6_TEST_ALL.

opt = struct('studies', {{'obstacle','stones','moving','width','range'}}, ...
             'plot', true, 'save', true, 'params', []);
for k = 1:2:numel(varargin), opt.(varargin{k}) = varargin{k+1}; end

p = opt.params;
if isempty(p), p = ch6_params(); end

root = fileparts(fileparts(mfilename('fullpath')));
stamp = datestr(now, 'yyyy-mm-dd_HH-MM-SS'); %#ok<TNOW1,DATST>
dir_out = fullfile(root, 'Results', ['ch6_' stamp]);
if opt.plot || opt.save
    if ~exist(dir_out, 'dir'), mkdir(dir_out); end
end

[x0, alpha] = load_reference(root);

fprintf('\n=============== CHAPTER 6 ===============\n');
fprintf(' gait: L_step %.4f m, controller %s, barrier %s\n', ...
        0.3533, p.controller, p.cbf.form);
fprintf(' poles (%.4g, %.4g), torque box %s at %.0f Nm\n', ...
        p.cbf.gamma_b, p.cbf.gamma, onoff(p.limits.enable.torque), ...
        p.limits.u_max);
fprintf(' output: %s\n', dir_out);

out = struct('p', p, 'dir', dir_out);
want = @(s) any(strcmpi(s, opt.studies));

%% ================================================== 6.1.2 overhead obstacles
if want('obstacle')
    fprintf('\n---------- 6.1.2  avoiding overhead obstacles ----------\n');
    ob = struct();
    for typ = {'ceiling', 'circle'}
        q = p;
        q.cbf.problem   = 'obstacle';
        q.obstacle.type = typ{1};
        q.n_steps       = 3;

        s = ch6_simulate(x0, alpha, q);
        R = ch6_report(s, q);
        ob.(typ{1}) = struct('sim', s, 'report', R);

        if opt.plot && s.n_ok > 0
            ch6_plot_barriers(s, q, fullfile(dir_out, ...
                sprintf('ch6_obstacle_%s.png', typ{1})));
            close;
        end
    end
    out.obstacle = ob;
end

%% ================================================ 6.1.3-6.1.4 discrete footholds
if want('stones')
    fprintf('\n---------- 6.1.3  walking over discrete footholds ----------\n');
    q = p;
    q.cbf.problem = 'stones';
    terr = ch6_terrain(p.n_steps, p.demo_band, p.mc.stone_sz, p.mc.seed);

    s = ch6_simulate(x0, alpha, q, terr);
    R = ch6_report(s, q);
    out.stones = struct('sim', s, 'report', R, 'terrain', terr);

    if opt.plot && s.n_ok > 0
        ch6_plot_stones(s,   fullfile(dir_out, 'ch6_stones_placement.png')); close;
        ch6_plot_barriers(s, q, fullfile(dir_out, 'ch6_stones_barriers.png')); close;
        ch6_plot_forces(s,   q, fullfile(dir_out, 'ch6_stones_forces.png')); close;
        ch6_animate(s,       q, fullfile(dir_out, 'ch6_stones_stick.png')); close;
        ch6_phase_plot(s,    q, alpha, fullfile(dir_out, 'ch6_stones_phase.png')); close;
    end
end

%% ==================================================== 6.2 moving stepping stones
if want('moving')
    fprintf('\n---------- 6.2  time-varying stepping stones ----------\n');
    mv = struct();
    for mot = {'linear', 'sinusoidal'}
        q = p;
        q.cbf.problem    = 'stones';
        q.stones.motion  = mot{1};
        q.stones.v_stone = 0.10;
        q.stones.amp     = 0.04;
        q.stones.freq    = 1.0;
        q.n_steps        = 6;

        terr = ch6_terrain(q.n_steps, p.demo_band, p.mc.stone_sz, p.mc.seed+7);
        s = ch6_simulate(x0, alpha, q, terr);
        R = ch6_report(s, q);
        mv.(mot{1}) = struct('sim', s, 'report', R);

        if opt.plot && s.n_ok > 0
            ch6_plot_stones(s, fullfile(dir_out, ...
                sprintf('ch6_moving_%s.png', mot{1})), ...
                sprintf('Time-varying stepping stones (%s)', mot{1}));
            close;
        end
    end
    out.moving = mv;
end

%% =============================================== 6.3 step length and step width
if want('width')
    fprintf('\n---------- 6.3  3D stepping stones (surrogate) ----------\n');
    spec = struct('T', 0.45, 'l_nom', 0.366, 'w_nom', 0.233, 'w0', 0.233, ...
                  'h_apex', 0.10, 'a_max', 60, 'stone_sz', 0.05);
    q = p;  q.cbf.gamma_b = 30;  q.cbf.gamma = 30;

    cases = {{'length', 0.28, 0.233}, ...
             {'width',  0.366, 0.28}, ...
             {'both',   0.30,  0.27}};
    wd = struct();
    for i = 1:3
        s = spec;
        s.cases = cases{i}{1};  s.l_d = cases{i}{2};  s.w_d = cases{i}{3};
        r = ch6_foot3d(q, s);
        wd.(sprintf('case%d', i)) = r;
        fprintf('  Case %d (%-6s): l_s = %.4f%s   w_s = %.4f%s   %d/%d feasible\n', ...
                i, s.cases, r.l_s, mark(r.in_l, strcmp(s.cases,'width')), ...
                r.w_s, mark(r.in_w, strcmp(s.cases,'length')), ...
                sum(r.log.feasible), numel(r.log.feasible));
    end
    out.width = wd;
end

%% ============================================ 6.1.4 the achievable step range
if want('range')
    fprintf('\n---------- 6.1.4  achievable step-length range ----------\n');
    out.range = ch6_step_range(x0, alpha, p);
end

%% =========================================================== 6.4 Table 6.1
if want('table')
    fprintf('\n---------- 6.4  Table 6.1 ----------\n');
    lib = [];
    try
        lib = ch6_lib_load(p);
        fprintf('  gait library: %d gaits, L in [%.3f, %.3f] m\n', ...
                numel(lib.L), lib.L(1), lib.L(end));
    catch ME
        fprintf('  no gait library (%s)\n  -> controllers I and III are skipped\n', ...
                ME.message);
    end
    out.table = ch6_table61(p, alpha, lib);
end

%% ------------------------------------------------------------------- save
if opt.save
    f = fullfile(dir_out, 'ch6_results.mat');
    save(f, '-struct', 'out');
    fprintf('\n Saved: %s\n', f);
end

fprintf('=========================================\n');

end

% ---------------------------------------------------------------------------
function [x0, alpha] = load_reference(root)
%LOAD_REFERENCE  The Chapter-3 reference gait: periodic, mesh-verified, stable.
f = fullfile(root, 'Results', 'ch3_reference_gait.mat');
if ~exist(f, 'file')
    error('ch6_main:noGait', ...
          ['Chapter 6 runs on the Chapter-3 reference gait and it is not at\n' ...
           '  %s\n' ...
           'Solve one with ch3_main, or point ch6_main at another result.'], f);
end
S = load(f);
[X, ~, alpha] = ch3_col_unpack(S.z_opt, S.p);
x0 = X(:,1);
end

function s = onoff(b)
if b, s = 'on'; else, s = 'off'; end
end

function s = mark(in_window, not_constrained)
if not_constrained, s = ' (free)';
elseif in_window,   s = ' [on stone]';
else,               s = ' *** MISSED ***';
end
end
