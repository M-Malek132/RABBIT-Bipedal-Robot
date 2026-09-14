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
%   (2) robust CLF-QP, Cases I-III          ch4_compare_controllers Section 4.1.4
%  (2b) robust CLF-QP, Case IV (scale 3)    ch4_compare_controllers Section 4.1.4
%   (3) L1 adaptive sweep                   ch4_compare_controllers Section 4.2.4
%   (4) figures                             ch4_plot_uncertainty
%   (5) animation                           ch4_animate
%
% STEP 0 REFUSES A GAIT THAT IS NOT AN ORBIT OF THE CURRENT DYNAMICS. Every
% number below measures a controller AGAINST the gait, so a gait that is no
% longer periodic on the model on the path turns every table and figure into a
% measurement of the reference instead. ch4_load_gait re-evaluates it
% (meta.orbit); this stops rather than producing a full, plausible-looking
% result set on it -- which is what happened for ten days after the
% 2026-09-02 dynamics regeneration.
%
% STEP 1 IS NOT OPTIONAL AND IT COMES FIRST FOR A REASON. The robust
% controller's guarantee is conditional on Delta1max, Delta2max actually
% bounding the uncertainty; a robust controller run outside its own bound is
% not a robust controller, it is an aggressive one. So the bounds are measured
% against the gait before anything is run with them, and -- unless the caller
% set rclf.delta1_max / rclf.delta2_max explicitly -- the measured bounds over
% the sweep's Cases I-III ARE the bounds used, and Case IV runs on its own,
% measured the same way at scale 3. Fixed defaults are fitted to
% one particular gait and need not cover another; a sweep that only printed
% "covers: NO" and ran anyway would certify nothing.
%
% THE SWEEPS RUN 25 STEPS, NOT THE THREE THE CHAPTER'S FIGURES SHOW. Three
% steps is long enough to separate a controller that falls at once from one
% that does not, and too short to separate one that converges from one that
% drifts. At the old eps = 0.35 the robust law's three-step rows looked
% converged in every case, while over 25 steps it fell in step 21 of Case I --
% with a perfect model. 25 is also the horizon of the chapter's own Fig. 4.6.
% Pass 'n_steps', 3 for the thesis-style three-step figures.
%
% Options (name/value)
%   'gait'      path to a Chapter-3 result .mat (default ch4_load_gait's)
%   'presets'   cell of {'robust','case4','l1'} (default all three)
%   'n_steps'   steps per run (default 25; see above)
%   'plot'      draw and save figures (default true)
%   'save'      write a .mat of everything (default true)
%   'animate'   write comparison GIFs (default true)
%   'anim_scales'      mass scales to animate (default [1.5 0.7])
%   'anim_controllers' who to race (default clfqp / rclfqp_con / l1)
%   'anim_steps'       steps per GIF (default 4). Separate from n_steps
%                      because the GIF is a fixed 150 frames: 26 steps in
%                      it would be six frames a step.
%   any ch4_params field, including dotted nested names, e.g.
%   ch4_main('l1.omega_c', 100, 'rclf.delta2_model', 'matrix')
%
% Output
%   out : struct .p .x0 .alpha .meta .bounds .robust .case4 .p_case4 .l1
%                .figs .gifs .file   (p_case4: the parameters Case IV ran with)
%
% See also CH4_PARAMS, CH4_LOAD_GAIT, CH4_COMPARE_CONTROLLERS, CH4_REPORT.

%% --- split our own options from ch4_params overrides --------------------
own = {'gait','presets','n_steps','plot','save','animate', ...
       'anim_scales','anim_controllers','anim_steps'};
o   = struct('gait', '', 'presets', {{'robust','case4','l1'}}, 'n_steps', 25, ...
             'plot', true, 'save', true, 'animate', true, ...
             'anim_scales', [1.5 0.7], ...
             'anim_controllers', {{'clfqp','rclfqp_con','l1'}}, ...
             'anim_steps', 4);

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
fprintf(' T = %.4f s, L = %.4f m, v = %.4f m/s, peak |u| %.1f Nm\n', ...
        meta.T, meta.L_step, meta.v_avg, p.gait_u_peak);
fprintf(' orbit on current dynamics: defect %.1e, periodicity %.1e, ||eta+|| %.1e (tol %.0e)\n', ...
        meta.orbit.defect, meta.orbit.periodicity, meta.orbit.eta_post, meta.orbit.tol);
fprintf(' control at %.0f Hz, eps %.2f, CLF via %s\n', ...
        1/p.control_dt, p.eps, p.clf_construction);

%% --- (0) the gait must be an orbit of the robot on the path --------------
if ~meta.orbit.ok
    error('ch4_main:notAnOrbit', ...
          ['"%s" is not a periodic orbit of the current dynamics (see the ' ...
           'residuals above), so every Chapter-4 number would measure the ' ...
           'reference rather than the controllers. Re-solve the gait in ' ...
           'Chapter 3, or pass ''gait'' pointing at one that verifies.'], ...
          meta.file);
end

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

SAFETY = 1.2;
B = ch4_delta_bounds(sim_b.x(:, 1:4:end), alpha, p, [1 0.7 1.5 3], ...
                     struct('n_jitter', 0, 'safety', SAFETY));

% Adopt what was measured over Cases I-III for their sweep. Case IV (scale 3)
% is not folded in: that would size the bounds -- and so the aggressiveness of
% the robust law -- for a perturbation Cases I-III never face. It gets its own
% bounds below instead. A bound the caller set explicitly is kept as given.
in_sweep = ismember([B.per_scale.mass_scale], [1 0.7 1.5]);
need_1   = max([B.per_scale(in_sweep).n1]);
need_2   = max([B.per_scale(in_sweep).n2]);

set_by_caller = lower(pv(1:2:end));
src = {'caller', 'caller'};
if ~any(ismember(set_by_caller, {'rclf', 'rclf.delta1_max'}))
    p.rclf.delta1_max = SAFETY * need_1;
    src{1} = sprintf('measured x %.1f', SAFETY);
end
if ~any(ismember(set_by_caller, {'rclf', 'rclf.delta2_max'}))
    p.rclf.delta2_max = SAFETY * need_2;
    src{2} = sprintf('measured x %.1f', SAFETY);
end

fprintf(' p.rclf in use: delta1_max %.1f (%s), delta2_max %.3f (%s), %s model\n', ...
        p.rclf.delta1_max, src{1}, p.rclf.delta2_max, src{2}, p.rclf.delta2_model);
covered_1 = p.rclf.delta1_max >= need_1;
covered_2 = p.rclf.delta2_max >= need_2;
fprintf(' covers Cases I-III on the orbit: Delta1 %s, Delta2 %s\n', ...
        yn(covered_1), yn(covered_2));
if ~(covered_1 && covered_2)
    fprintf([' NOTE: a bound below the measured value means the robust\n' ...
             ' guarantee does not cover these cases. ch4_report flags each\n' ...
             ' run that leaves the set.\n']);
end
if p.rclf.delta2_max >= 1
    fprintf([' NOTE: delta2_max >= 1 -- the robust CLF-QP is not pointwise\n' ...
             ' feasible and will report every sample infeasible.\n']);
end

% Case IV's own bounds: the same rule, at scale 3 alone. A pure mass scale fixes
% ||Delta2|| = |1/s - 1| exactly -- 2/3 at s = 3, against 3/7 at worst in Cases
% I-III -- so the Cases I-III bounds would run Case IV outside the robust law's
% own hypothesis, and 2/3 x SAFETY = 0.8 is still below the feasibility limit
% of 1.
p_case4 = p;
if any(strcmpi('case4', o.presets))
    k4 = find([B.per_scale.mass_scale] == 3, 1);
    if ~any(ismember(set_by_caller, {'rclf', 'rclf.delta1_max'}))
        p_case4.rclf.delta1_max = SAFETY * B.per_scale(k4).n1;
    end
    if ~any(ismember(set_by_caller, {'rclf', 'rclf.delta2_max'}))
        p_case4.rclf.delta2_max = SAFETY * B.per_scale(k4).n2;
    end
    fprintf(' Case IV bounds: delta1_max %.1f, delta2_max %.3f (covers scale 3: Delta1 %s, Delta2 %s)\n', ...
            p_case4.rclf.delta1_max, p_case4.rclf.delta2_max, ...
            yn(p_case4.rclf.delta1_max >= B.per_scale(k4).n1), ...
            yn(p_case4.rclf.delta2_max >= B.per_scale(k4).n2));
end

out = struct('p', p, 'x0', x0, 'alpha', alpha, 'meta', meta, 'bounds', B, ...
             'robust', [], 'case4', [], 'p_case4', [], 'l1', [], ...
             'figs', [], 'gifs', {{}}, 'file', '');

%% --- (2)(2b)(3) the sweeps -----------------------------------------------
copts = struct('n_steps', o.n_steps, 'store_traj', true, 'verbose', true);

if any(strcmpi('robust', o.presets))
    out.robust = ch4_compare_controllers(x0, alpha, p, 'robust', copts);
end
if any(strcmpi('case4', o.presets))
    out.case4   = ch4_compare_controllers(x0, alpha, p_case4, 'case4', copts);
    out.p_case4 = p_case4;
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
    if ~isempty(out.case4)
        d = fullfile(results_dir, sprintf('ch4_case4_%s', stamp));
        figs = [figs, ch4_plot_uncertainty(out.case4, out.p_case4, d, {'IV'})];
    end
    if ~isempty(out.l1)
        d = fullfile(results_dir, sprintf('ch4_l1_%s', stamp));
        figs = [figs, ch4_plot_uncertainty(out.l1, p, d)];
    end
    out.figs = figs;
end

%% --- save, BEFORE the animation -----------------------------------------
% Ordering is deliberate. The sweeps above cost minutes; the animation below
% renders hundreds of frames through getframe, and a renderer that dies takes
% the whole MATLAB process with it rather than raising something the try/catch
% could hold. Saving first makes the expensive numerical results durable no
% matter what the graphics stack does -- measured the hard way, on a run that
% completed every sweep and every figure and then died mid-GIF with nothing
% written.
if o.save
    out.file = fullfile(results_dir, sprintf('ch4_result_%s.mat', stamp));
    R = rmfield(out, 'figs'); %#ok<NASGU>
    save(out.file, '-struct', 'R');
    fprintf(' Saved: %s\n', out.file);
end

%% --- (5) animation -------------------------------------------------------
% One GIF per perturbation, each showing every controller on the SAME robot at
% the same instant. The tables say the baseline "completed 1 step"; this is
% where you see what that means.
if o.animate
    out.gifs = {};
    for s = o.anim_scales
        pa = p;
        pa.uncertainty.mass_scale = s;

        % The boxes the sweeps used at this scale, so each panel is the
        % controller its table row scored rather than one at a default box.
        k = [];
        if ~isempty(out.robust), k = find([out.robust.mass_scale] == s, 1); end
        if ~isempty(k)
            pa.limits.u_max = out.robust(k).u_box;
        elseif any(ismember(o.anim_controllers, {'clfqp_con', 'rclfqp_con'}))
            fprintf([' ch4_animate: no robust sweep at scale %.2f, so the ' ...
                     'constrained robust laws animate at limits.u_max %.0f Nm\n'], ...
                    s, pa.limits.u_max);
        end
        if ~isempty(out.l1)
            k = find([out.l1.mass_scale] == s, 1);
            if ~isempty(k), pa.l1.u_max = out.l1(k).u_box; end
        end

        gif = fullfile(results_dir, ...
                       sprintf('ch4_walk_%s_scale%03.0f.gif', stamp, s*100));
        try
            ch4_animate(x0, alpha, pa, o.anim_controllers, o.anim_steps, gif);
            out.gifs{end+1} = gif;
        catch err
            fprintf(' ch4_animate skipped for scale %.2f: %s\n', s, err.message);
        end
    end

    % Record which GIFs actually got written. Cheap append rather than a second
    % full write of the trajectory data.
    if o.save && ~isempty(out.gifs)
        gifs = out.gifs; %#ok<NASGU>
        save(out.file, 'gifs', '-append');
    end
end

fprintf('==========================================\n\n');

end

% ---------------------------------------------------------------------------
function s = yn(b)
if b, s = 'yes'; else, s = 'NO'; end
end
