function C = ch4_compare_controllers(x0, alpha, p, preset, opts)
%CH4_COMPARE_CONTROLLERS  The Chapter-4 experiment: same gait, wrong model.
%
%   C = ch4_compare_controllers(x0, alpha, p)
%   C = ch4_compare_controllers(x0, alpha, p, preset)
%   C = ch4_compare_controllers(x0, alpha, p, preset, opts)
%
% Runs one optimized gait under several controllers, across several model
% perturbations that none of them are told about. This is the numerical
% validation of both halves of the chapter -- Section 4.1.4 for the robust
% controller and Section 4.2.4 for the adaptive one -- and the table it prints
% is what Figures 4.2 and 4.8 plot.
%
% ------------------------------------------------------------------ presets
% 'robust'  Section 4.1.4.  Controllers A/B/C = min-norm CLF-QP, CLF-QP with
%           torque saturation, robust CLF-QP with torque saturation, over
%           Cases I-III (mass scale 1, 1.5, 0.7).
%
%           THE TORQUE BOX IS PER-CASE, AND THAT IS THE POINT. Section 4.1.4
%           sets it "slightly below the maximum torque that controller A uses"
%           for each perturbation -- 60, 80, 150 Nm for scales 0.7, 1, 1.5.
%           A fixed box across cases would be an unfair comparison: a 1.5x
%           robot genuinely needs more torque, and the interesting question is
%           whether the robust controller keeps its convergence rate when it is
%           held to what the baseline actually drew, not whether it survives an
%           arbitrary limit.
%
%           SO THE BOX IS MEASURED, NOT COPIED. Those three numbers belong to
%           the thesis's 32 kg robot. On the current 74 kg model controller A
%           drew 196 / 779 / 507 Nm on posture_195 (eps 0.35, three steps --
%           the measurements in this note all date from that setting), so
%           the copied boxes sat
%           at 12-41% of the baseline's peak, starved both constrained laws off
%           their feet in every case -- the perfect-model case included --
%           and the sweep measured nothing. The thesis's RULE is applied
%           instead: controller A runs first at each scale, and
%
%               box = max( opts.box_frac * peak|u| of A ,  scale * p.gait_u_peak )
%
%           with box_frac = 0.8. The floor is the true robot's own feedforward
%           peak -- s times the mass needs s times the torque to follow the
%           same motion -- and it binds in exactly one place, Case I. With a
%           perfect model A draws barely more than the feedforward itself
%           (196.4 vs 195.0 Nm), so any box "slightly below" it is below what
%           the gait needs: at 0.8x (157 Nm) clfqp_con fell in step 3 and
%           rclfqp_con in step 1, while at the 195 Nm floor both walk all
%           three. Pass opts.u_box to fix the boxes by hand instead.
%
% 'case4'   Section 4.1.4's Case IV: the same three controllers at mass scale 3,
%           boxes measured by the same rule. It is a separate preset because it
%           needs separate BOUNDS. For a uniform scale s, ||Delta2|| is exactly
%           |1/s - 1| -- 2/3 at s = 3, against 3/7 for the worst of Cases I-III
%           -- so bounds sized for Cases I-III leave Case IV outside the robust
%           law's own hypothesis, while sizing every case for Case IV would make
%           Cases I-III measure a law tuned for a perturbation they never face.
%           ch4_main runs it with bounds measured at scale 3.
%
% READ Vend/Vmx LOOSELY FOR 'rclfqp_con' WITH p.rclf.boundary_layer = 0. The
% exact law's worst-case term chatters at the sample rate (see
% ch4_ctrl_rclf_qp), so V at the final sample lands anywhere in an
% order-of-magnitude band: the same Case I run gave 0.021 at a 195.0 Nm box
% and 0.19 at 195 Nm plus a floating-point hair. Steps, max||eta|| and min Fz
% are the columns that survive that. The default boundary layer removes the
% chatter, and with it this caveat.
%
% 'l1'      Section 4.2.4.  Controllers A/B/C = CLF-QP, L1 + CLF-QP, L1 +
%           CLF-QP with torque saturation, over mass scales 1, 0.7, 1.5. The
%           thesis's 65 Nm is replaced by p.l1.u_max, which ch4_load_gait
%           sizes to the gait (see the note there). The box is the same at
%           every scale in the thesis form, and grows with a heavier robot
%           under p.l1.constrain_applied (see below).
%
% ------------------------------------------------------- what to read off it
% The headline column is max|eta| -- the tracking error the chapter plots. The
% claim being tested is not "the robust/adaptive controller is better on
% average"; it is the sharper one that its convergence behaviour is UNCHANGED
% across the perturbations, while the baseline's degrades. So compare each
% controller's row DOWN the scales, not across controllers within a scale.
%
% Vend/Vmx makes that explicit: the CLF left at the end of the run as a
% fraction of the largest excursion it reached. Every impact throws the outputs
% off and the controller has one step to recover; this is how much of the worst
% throw was still outstanding when the run ended. A controller holding its rate
% keeps this small in EVERY case; one losing the fight has it climb toward 1 as
% the perturbation grows.
%
% Every run is configured by ch4_run_params -- a box only for the laws whose
% formulation has one, contact rows as p has them -- so ch4_animate races
% exactly the controllers this table scores.
%
% "min Fz" is the TRUE model's normal force at the stance foot. The stance
% contact is integrated as a pin, so a controller can demand a NEGATIVE normal
% force -- the ground pulling the foot down -- and keep walking in simulation.
% A negative entry means that run's steps are not physically realizable,
% whatever the tracking columns say.
%
% Inputs
%   x0, alpha : an optimized gait (ch3_col_unpack of a solved z)
%   p         : parameter struct
%   preset    : 'robust' (default), 'case4' or 'l1'
%   opts      : struct overriding .controllers .scales .u_box .box_frac
%               .n_steps .verbose .store_traj. u_box = [] (the 'robust' and
%               'case4' default) measures the boxes as above; a scalar or one
%               value per scale fixes them.
%
% Output
%   C : struct array, one entry per (controller, scale), with
%         .name .mass_scale .u_box .steps_completed .peak_torque
%         .box_violation .max_eta .final_eta .V0 .Vend .V_ratio
%         .delta_max .qp_infeasible .int_u2 .Fz_min .grf_pred_error .reason
%         .step_T .step_L   duration [s] and length [m] of each completed step,
%                           so whether a run slows down can be read directly
%         .traj (t, y, u, V, theta_hat) when opts.store_traj
%
% See also CH4_SIMULATE, CH4_FORCES, CH4_PLOT_UNCERTAINTY, CH3_COMPARE_CONTROLLERS.

if nargin < 4 || isempty(preset), preset = 'robust'; end
if nargin < 5, opts = struct(); end

switch lower(preset)
    case 'robust'
        % u_box empty = measured per case from controller A, see the header.
        def = struct('controllers', {{'clfqp', 'clfqp_con', 'rclfqp_con'}}, ...
                     'scales',      [1 1.5 0.7], ...
                     'u_box',       []);
    case 'case4'
        % Case IV alone. p must carry bounds that cover it -- see the header.
        def = struct('controllers', {{'clfqp', 'clfqp_con', 'rclfqp_con'}}, ...
                     'scales',      3, ...
                     'u_box',       []);
    case 'l1'
        % The box comes from p, not from a constant here: ch4_load_gait has
        % already raised it to the gait's own peak torque if the Section
        % 4.2.4 figure of 65 Nm could not execute this gait. Hard-coding 65
        % would undo that on every run and starve the inner QP -- see the
        % note in ch4_load_gait for what that does to the adaptation.
        % u_box = [] here means "from p", resolved below once the scales
        % are known -- not "measured", which is the robust preset's meaning.
        def = struct('controllers', {{'clfqp', 'l1', 'l1_con'}}, ...
                     'scales',      [1 0.7 1.5], ...
                     'u_box',       []);
    otherwise
        error('ch4_compare_controllers:preset', ...
              'Unknown preset "%s" (expected robust|case4|l1).', preset);
end

def.box_frac   = 0.8;
def.n_steps    = 25;                   % see ch4_main: three steps hide drift
def.verbose    = true;
def.store_traj = true;
opts = fill_defaults(opts, def);

names  = opts.controllers;
scales = opts.scales;
boxes  = opts.u_box;

% THE L1 BOX. Under the thesis form it bounds mu1 alone, so one box from p
% serves every scale. Under p.l1.constrain_applied it bounds the TOTAL torque,
% and a box that cannot carry the heaviest robot makes the sweep measure the
% box: s times the mass needs s times the torque for the same motion, and at
% 1.5x a 244 Nm box is below the 293 Nm the true robot's feedforward alone
% needs. Measured: 'l1_con' fell in step 2-3 at 1.5x in every such variant,
% and walked all 25 steps once the box was scaled. So the box grows with the
% robot, and never shrinks below the nominal gait's box, since the inner QP
% delivers the nominal feedforward before the adaptation has learned anything.
if strcmpi(preset, 'l1') && isempty(boxes)
    l1o = ch4_l1_opts(p);
    if l1o.constrain_applied
        boxes = p.l1.u_max * max(1, scales);
    else
        boxes = repmat(p.l1.u_max, size(scales));
    end
end

measure_box = isempty(boxes);
if ~measure_box
    if numel(boxes) == 1, boxes = repmat(boxes, size(scales)); end
    if numel(boxes) ~= numel(scales)
        error('ch4_compare_controllers:uBox', ...
              'opts.u_box must be scalar or match opts.scales (%d entries).', ...
              numel(scales));
    end
end

u_ff_peak = p.gait_u_peak;
if isempty(u_ff_peak), u_ff_peak = 0; end

C = struct('name', {}, 'mass_scale', {}, 'u_box', {}, 'steps_completed', {}, ...
           'peak_torque', {}, 'box_violation', {}, 'max_eta', {}, ...
           'final_eta', {}, 'V0', {}, 'Vend', {}, 'V_ratio', {}, ...
           'delta_max', {}, 'qp_infeasible', {}, 'int_u2', {}, ...
           'Fz_min', {}, 'grf_pred_error', {}, 'reason', {}, ...
           'step_T', {}, 'step_L', {}, 'traj', {});

if opts.verbose
    fprintf('\n%s\n CHAPTER 4 -- preset "%s", %d steps per run\n', ...
            repmat('=',1,105), preset, opts.n_steps);
    fprintf(' robust bounds: delta1_max %.1f, delta2_max %.3f (%s model)\n', ...
            p.rclf.delta1_max, p.rclf.delta2_max, p.rclf.delta2_model);
    % A pure mass scale fixes ||Delta2|| = |1/s - 1| exactly, so a robust law
    % about to run outside its own bound can be named before it runs.
    if any(strncmpi(names, 'rclfqp', 6)) && p.uncertainty.load_mass == 0
        [d2_need, iw] = max(abs(1 ./ scales - 1));
        if p.rclf.delta2_max < d2_need
            fprintf([' NOTE: delta2_max %.3f is below ||Delta2|| = %.3f at ' ...
                     'scale %.2g: the robust law runs there outside its own ' ...
                     'hypothesis.\n'], p.rclf.delta2_max, d2_need, scales(iw));
        end
    end
    fprintf(' L1: Gamma %.0e, filter %.0f rad/s, sample %.0f Hz\n', ...
            p.l1.Gamma, p.l1.omega_c, 1/p.control_dt);
    if measure_box
        fprintf(' box per case: max(%.2f x clfqp peak, scale x %.1f Nm feedforward peak)\n', ...
                opts.box_frac, u_ff_peak);
    end
    fprintf('%s\n', repmat('-',1,105));
    fprintf(' %-11s %6s %7s %6s %9s %9s %10s %10s %9s %7s %8s\n', ...
            'controller', 'scale', 'box', 'steps', 'peak|u|', 'over box', ...
            'max|eta|', 'fin|eta|', 'Vend/Vmx', 'QPfail', 'min Fz');
    fprintf('%s\n', repmat('-',1,105));
end

for is = 1:numel(scales)

    runs = cell(1, numel(names));

    if measure_box
        % Controller A first: what it draws at THIS perturbation sets the box.
        % max() skips a NaN peak, so a baseline that never completed a step
        % leaves the box on the feedforward floor rather than on NaN.
        eA = run_one(x0, alpha, p, 'clfqp', scales(is), [], opts);
        iA = find(strcmpi(names, 'clfqp'), 1);
        if ~isempty(iA), runs{iA} = eA; end

        box = max([opts.box_frac * eA.peak_torque, scales(is) * u_ff_peak]);
        if ~(isfinite(box) && box > 0)
            error('ch4_compare_controllers:noBox', ...
                  ['Cannot size the box at mass scale %.2f: controller A ' ...
                   'completed no step and p.gait_u_peak is empty. Pass ' ...
                   'opts.u_box.'], scales(is));
        end
    else
        box = boxes(is);
    end

    for ic = 1:numel(names)
        if isempty(runs{ic})
            runs{ic} = run_one(x0, alpha, p, names{ic}, scales(is), box, opts);
        end
        entry = runs{ic};

        % Recorded for every row, the baselines included: "over box" for a law
        % that was not told the limit is the point of the column.
        entry.u_box = box;
        if ~isnan(entry.peak_torque)
            entry.box_violation = max(entry.peak_torque - box, 0);
        end

        C(end+1) = entry; %#ok<AGROW>

        if opts.verbose
            fprintf(' %-11s %6.2f %7.0f %6d %9.1f %9.1f %10.3e %10.3e %9.3f %7d %8.0f\n', ...
                    entry.name, entry.mass_scale, entry.u_box, ...
                    entry.steps_completed, entry.peak_torque, ...
                    entry.box_violation, entry.max_eta, entry.final_eta, ...
                    entry.V_ratio, entry.qp_infeasible, entry.Fz_min);
        end
    end
    if opts.verbose && is < numel(scales)
        fprintf('%s\n', repmat('-',1,105));
    end
end

if opts.verbose
    fprintf('%s\n', repmat('=',1,105));
    fprintf([' Read DOWN each controller across scales: the claim is that the\n' ...
             ' robust / adaptive laws hold their convergence behaviour while the\n' ...
             ' baseline degrades. "over box" is how far the COMMANDED torque left\n' ...
             ' the limit -- nonzero only for laws not told the limit exists.\n' ...
             ' "min Fz" is the TRUE normal force [N]; below zero the ground had\n' ...
             ' to pull the foot down, so those steps are not realizable.\n']);
    fprintf('%s\n\n', repmat('=',1,105));
end

end

% ---------------------------------------------------------------------------
function entry = run_one(x0, alpha, p, name, scale, box, opts)
%RUN_ONE  Simulate one (controller, scale) pair and summarize it.
% box = [] keeps p's boxes; only the laws that carry a box read it.
pc  = ch4_run_params(p, name, scale, box);
sim = ch4_simulate(x0, alpha, pc, opts.n_steps);

entry = struct('name', name, 'mass_scale', scale, ...
               'u_box', NaN, 'steps_completed', sim.n_ok, ...
               'peak_torque', NaN, 'box_violation', NaN, ...
               'max_eta', NaN, 'final_eta', NaN, ...
               'V0', NaN, 'Vend', NaN, 'V_ratio', NaN, ...
               'delta_max', NaN, 'qp_infeasible', NaN, ...
               'int_u2', NaN, 'Fz_min', NaN, 'grf_pred_error', NaN, ...
               'reason', sim.reason, ...
               'step_T', [sim.steps.T], 'step_L', [sim.steps.L_step], ...
               'traj', []);

if sim.n_ok == 0, return; end

% ch4_forces decimates to a fixed number of points, 2000 by default -- about
% 670 a step at the three-step horizon it was set for. Keep that density per
% step rather than per run, or a 25-step run is analysed on a grid eight times
% coarser and its min Fz can step over a contact-force dip the short run would
% have caught.
F = ch4_forces(sim.t, sim.x, alpha, pc, sim.xi, sim.t_xi, ...
               ceil(2000 * max(opts.n_steps, 3) / 3));

eta_n = vecnorm([F.y; F.ydot], 2, 1);

entry.peak_torque    = F.torque_max;
entry.max_eta        = max(eta_n);
entry.final_eta      = eta_n(end);
entry.delta_max      = F.delta_max;
entry.qp_infeasible  = F.qp_infeasible;
entry.int_u2         = F.int_u2;
entry.Fz_min         = F.Fz_min;
entry.grf_pred_error = F.grf_pred_error;

% V is NaN for the pure-feedforward and PD laws; fall back to
% ||eta||^2 so the column still means something for them.
Vs = F.V;
if all(isnan(Vs)), Vs = eta_n.^2; end
entry.V0   = Vs(1);
entry.Vend = Vs(end);

% RECOVERED FRACTION, not Vend/V0.
%
% The run starts ON the periodic orbit, so V0 is essentially zero
% and Vend/V0 is a ratio of a real number to numerical noise -- it
% came out as 1e12 and said nothing. What the chapter's claim is
% actually about is whether the controller RECOVERS from the
% excursion each impact creates, so measure exactly that: the
% residual left at the end as a fraction of the largest excursion
% reached. Bounded in [0,1]; near 0 means converged, near 1 means
% the controller never got the error back.
entry.V_ratio = Vs(end) / max(max(Vs), realmin);

if opts.store_traj
    % F.t / F.x, not sim.t / sim.x: ch4_forces decimates, and every
    % other field here lives on ITS grid.
    entry.traj = struct('t', F.t, 'y', F.y, 'u', F.u, ...
                        'V', F.V, 'eta_n', eta_n, ...
                        'theta_hat', F.theta_hat, ...
                        'theta_true', F.theta_true, ...
                        'x', F.x, 'lambda', F.lambda);
end
end

% ---------------------------------------------------------------------------
function s = fill_defaults(s, d)
f = fieldnames(d);
for i = 1:numel(f)
    if ~isfield(s, f{i}) || isempty(s.(f{i}))
        s.(f{i}) = d.(f{i});
    end
end
end
