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
% 'l1'      Section 4.2.4.  Controllers A/B/C = CLF-QP, L1 + CLF-QP, L1 +
%           CLF-QP with torque saturation at a FIXED 65 Nm, over mass scales
%           1, 0.7, 1.5.
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
% Inputs
%   x0, alpha : an optimized gait (ch3_col_unpack of a solved z)
%   p         : parameter struct
%   preset    : 'robust' (default) or 'l1'
%   opts      : struct overriding .controllers .scales .u_box .n_steps
%               .verbose .store_traj
%
% Output
%   C : struct array, one entry per (controller, scale), with
%         .name .mass_scale .u_box .steps_completed .peak_torque
%         .box_violation .max_eta .final_eta .V0 .Vend .V_ratio
%         .delta_max .qp_infeasible .int_u2 .Fz_min .grf_pred_error .reason
%         .traj (t, y, u, V, theta_hat) when opts.store_traj
%
% See also CH4_SIMULATE, CH4_FORCES, CH4_PLOT_UNCERTAINTY, CH3_COMPARE_CONTROLLERS.

if nargin < 4 || isempty(preset), preset = 'robust'; end
if nargin < 5, opts = struct(); end

switch lower(preset)
    case 'robust'
        def = struct('controllers', {{'clfqp', 'clfqp_con', 'rclfqp_con'}}, ...
                     'scales',      [1 1.5 0.7], ...
                     'u_box',       [80 150 60]);
    case 'l1'
        def = struct('controllers', {{'clfqp', 'l1', 'l1_con'}}, ...
                     'scales',      [1 0.7 1.5], ...
                     'u_box',       [65 65 65]);
    otherwise
        error('ch4_compare_controllers:preset', ...
              'Unknown preset "%s" (expected robust|l1).', preset);
end

def.n_steps    = 3;
def.verbose    = true;
def.store_traj = true;
opts = fill_defaults(opts, def);

names  = opts.controllers;
scales = opts.scales;
boxes  = opts.u_box;

if numel(boxes) == 1, boxes = repmat(boxes, size(scales)); end
if numel(boxes) ~= numel(scales)
    error('ch4_compare_controllers:uBox', ...
          'opts.u_box must be scalar or match opts.scales (%d entries).', ...
          numel(scales));
end

C = struct('name', {}, 'mass_scale', {}, 'u_box', {}, 'steps_completed', {}, ...
           'peak_torque', {}, 'box_violation', {}, 'max_eta', {}, ...
           'final_eta', {}, 'V0', {}, 'Vend', {}, 'V_ratio', {}, ...
           'delta_max', {}, 'qp_infeasible', {}, 'int_u2', {}, ...
           'Fz_min', {}, 'grf_pred_error', {}, 'reason', {}, 'traj', {});

if opts.verbose
    fprintf('\n%s\n CHAPTER 4 -- preset "%s", %d steps per run\n', ...
            repmat('=',1,96), preset, opts.n_steps);
    fprintf(' robust bounds: delta1_max %.1f, delta2_max %.3f (%s model)\n', ...
            p.rclf.delta1_max, p.rclf.delta2_max, p.rclf.delta2_model);
    fprintf(' L1: Gamma %.0e, filter %.0f rad/s, sample %.0f Hz\n', ...
            p.l1.Gamma, p.l1.omega_c, 1/p.control_dt);
    fprintf('%s\n', repmat('-',1,96));
    fprintf(' %-11s %6s %7s %6s %9s %9s %10s %10s %9s %7s\n', ...
            'controller', 'scale', 'box', 'steps', 'peak|u|', 'over box', ...
            'max|eta|', 'fin|eta|', 'Vend/Vmx', 'QPfail');
    fprintf('%s\n', repmat('-',1,96));
end

for is = 1:numel(scales)
    for ic = 1:numel(names)

        pc = p;
        pc.controller             = names{ic};
        pc.uncertainty.mass_scale = scales(is);
        pc.limits.u_max           = boxes(is);
        pc.l1.u_max               = boxes(is);

        % Only the laws that carry a box in their formulation are TOLD about
        % it. The baselines are not -- that asymmetry is the experiment.
        pc.limits.enable.torque = any(strcmpi(names{ic}, ...
                                       {'clfqp_con', 'rclfqp_con'}));

        sim = ch4_simulate(x0, alpha, pc, opts.n_steps);

        entry = struct('name', names{ic}, 'mass_scale', scales(is), ...
                       'u_box', boxes(is), 'steps_completed', sim.n_ok, ...
                       'peak_torque', NaN, 'box_violation', NaN, ...
                       'max_eta', NaN, 'final_eta', NaN, ...
                       'V0', NaN, 'Vend', NaN, 'V_ratio', NaN, ...
                       'delta_max', NaN, 'qp_infeasible', NaN, ...
                       'int_u2', NaN, 'Fz_min', NaN, 'grf_pred_error', NaN, ...
                       'reason', sim.reason, 'traj', []);

        if sim.n_ok > 0
            F = ch4_forces(sim.t, sim.x, alpha, pc, sim.xi, sim.t_xi);

            eta_n = vecnorm([F.y; F.ydot], 2, 1);

            entry.peak_torque    = F.torque_max;
            entry.box_violation  = max(F.torque_max - boxes(is), 0);
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

        C(end+1) = entry; %#ok<AGROW>

        if opts.verbose
            fprintf(' %-11s %6.2f %7.0f %6d %9.1f %9.1f %10.3e %10.3e %9.3f %7d\n', ...
                    entry.name, entry.mass_scale, entry.u_box, ...
                    entry.steps_completed, entry.peak_torque, ...
                    entry.box_violation, entry.max_eta, entry.final_eta, ...
                    entry.V_ratio, entry.qp_infeasible);
        end
    end
    if opts.verbose && is < numel(scales)
        fprintf('%s\n', repmat('-',1,96));
    end
end

if opts.verbose
    fprintf('%s\n', repmat('=',1,96));
    fprintf([' Read DOWN each controller across scales: the claim is that the\n' ...
             ' robust / adaptive laws hold their convergence behaviour while the\n' ...
             ' baseline degrades. "over box" is how far the COMMANDED torque left\n' ...
             ' the limit -- nonzero only for laws not told the limit exists.\n']);
    fprintf('%s\n\n', repmat('=',1,96));
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
