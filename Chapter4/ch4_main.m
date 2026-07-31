function out = ch4_main(varargin)
%CH4_MAIN  Run the Chapter-4 study end to end.
%
%   out = ch4_main()
%   out = ch4_main('n_steps', 5, 'presets', {'l1'}, 'plot', false)
%
% Chapter 4 does not solve an optimization. It takes the Chapter-3 gait, breaks
% the controller's model, and measures what each control scheme does about it.
% The pipeline is therefore short and reads in the order the chapter argues:
%
%   (0) load the nominal gait               ch4_load_gait
%   (1) MEASURE the uncertainty             ch4_delta_bounds        (4.4), (4.10)
%   (2) robust CLF-QP sweep                 ch4_compare_controllers Section 4.1.4
%   (3) L1 adaptive sweep                   ch4_compare_controllers Section 4.2.4
%   (4) figures                             ch4_plot_uncertainty
%
% STEP 1 IS NOT OPTIONAL AND IT COMES FIRST FOR A REASON. The robust
% controller's guarantee is conditional on Delta1max, Delta2max actually
% bounding the uncertainty; a robust controller run outside its own bound is
% not a robust controller, it is an aggressive one. So the bounds are measured
% against the gait before anything is run with them, and the measurement is
% reported next to the numbers in p.rclf so a mismatch is visible rather than
% buried.
%
% Options (name/value)
%   'gait'      path to a Chapter-3 result .mat (default the upright gait)
%   'presets'   cell of {'robust','l1'} (default both)
%   'n_steps'   steps per run (default 3, matching the chapter's figures)
%   'plot'      draw and save figures (default true)
%   'save'      write a .mat of everything (default true)
%   any ch4_params field, including dotted nested names, e.g.
%   ch4_main('l1.omega_c', 100, 'rclf.delta2_model', 'matrix')
%
% Output
%   out : struct .p .x0 .alpha .meta .bounds .robust .l1 .figs .file
%
% See also CH4_PARAMS, CH4_LOAD_GAIT, CH4_COMPARE_CONTROLLERS, CH4_REPORT.

%% --- split our own options from ch4_params overrides --------------------
own = {'gait','presets','n_steps','plot','save'};
o   = struct('gait', '', 'presets', {{'robust','l1'}}, 'n_steps', 3, ...
             'plot', true, 'save', true);

pv = {};
for k = 1:2:numel(varargin)
    if any(strcmpi(varargin{k}, own))
        o.(lower(varargin{k})) = varargin{k+1};
    else
        pv(end+1:end+2) = varargin(k:k+1); %#ok<AGROW>
    end
end

[x0, alpha, p, meta] = ch4_load_gait(o.gait, pv{:});

results_dir = fullfile(fileparts(mfilename('fullpath')), '..', 'Results');
if ~exist(results_dir, 'dir'), mkdir(results_dir); end
stamp = datestr(now, 'yyyy-mm-dd_HH-MM-SS'); %#ok<TNOW1,DATST>

fprintf('\n================ CHAPTER 4 ================\n');
fprintf(' gait   %s\n', meta.file);
fprintf(' T = %.4f s, L = %.4f m, v = %.4f m/s\n', meta.T, meta.L_step, meta.v_avg);
fprintf(' control at %.0f Hz, eps %.2f, CLF via %s\n', ...
        1/p.control_dt, p.eps, p.clf_construction);

%% --- (1) measure the uncertainty before designing against it ------------
% Sample along an actual rollout of the baseline rather than the collocation
% nodes: 15 nodes is too coarse to bound anything, and the controller runs on
% the trajectory, not on the mesh.
pb = p; pb.controller = 'clfqp';
pb.uncertainty = struct('mass_scale', 1, 'load_mass', 0);
sim_b = ch4_simulate(x0, alpha, pb, 2);

if sim_b.n_ok == 0
    error('ch4_main:baseline', ...
          ['The nominal baseline did not complete a step (%s), so there is ' ...
           'no trajectory to measure uncertainty along.'], sim_b.reason);
end

B = ch4_delta_bounds(sim_b.x(:, 1:4:end), alpha, p, [1 0.7 1.5 3], ...
                     struct('n_jitter', 0, 'safety', 1.2));

fprintf(' p.rclf in use: delta1_max %.1f, delta2_max %.3f (%s)\n', ...
        p.rclf.delta1_max, p.rclf.delta2_max, p.rclf.delta2_model);
covered_1 = p.rclf.delta1_max >= max([B.per_scale(2:3).n1]);
covered_2 = p.rclf.delta2_max >= max([B.per_scale(2:3).n2]);
fprintf(' covers Cases I-III on the orbit: Delta1 %s, Delta2 %s\n', ...
        yn(covered_1), yn(covered_2));
if ~(covered_1 && covered_2)
    fprintf([' NOTE: a bound below the measured value means the robust\n' ...
             ' guarantee does not cover these cases. ch4_report flags each\n' ...
             ' run that leaves the set.\n']);
end

out = struct('p', p, 'x0', x0, 'alpha', alpha, 'meta', meta, 'bounds', B, ...
             'robust', [], 'l1', [], 'figs', [], 'file', '');

%% --- (2)(3) the sweeps ---------------------------------------------------
copts = struct('n_steps', o.n_steps, 'store_traj', true, 'verbose', true);

if any(strcmpi('robust', o.presets))
    out.robust = ch4_compare_controllers(x0, alpha, p, 'robust', copts);
end
if any(strcmpi('l1', o.presets))
    out.l1 = ch4_compare_controllers(x0, alpha, p, 'l1', copts);
end

%% --- (4) figures ---------------------------------------------------------
if o.plot
    figs = gobjects(0);
    if ~isempty(out.robust)
        d = fullfile(results_dir, sprintf('ch4_robust_%s', stamp));
        figs = [figs, ch4_plot_uncertainty(out.robust, p, d)];
    end
    if ~isempty(out.l1)
        d = fullfile(results_dir, sprintf('ch4_l1_%s', stamp));
        figs = [figs, ch4_plot_uncertainty(out.l1, p, d)];
    end
    out.figs = figs;
end

%% --- save ----------------------------------------------------------------
if o.save
    out.file = fullfile(results_dir, sprintf('ch4_result_%s.mat', stamp));
    R = rmfield(out, 'figs'); %#ok<NASGU>
    save(out.file, '-struct', 'R');
    fprintf(' Saved: %s\n', out.file);
end

fprintf('==========================================\n\n');

end

% ---------------------------------------------------------------------------
function s = yn(b)
if b, s = 'yes'; else, s = 'NO'; end
end
