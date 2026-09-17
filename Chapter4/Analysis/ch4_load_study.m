function S = ch4_load_study(x0, alpha, p, opts)
%CH4_LOAD_STUDY  Section 4.2.4's unknown-load experiment (Fig. 4.11).
%
%   S = ch4_load_study(x0, alpha, p)
%   S = ch4_load_study(x0, alpha, p, opts)
%
% The chapter closes its L1 section on a harder perturbation than a mass scale:
% the robot carries a mass on its torso that the controller is never told about,
% redrawn at random on every step (Fig. 4.11a) or held fixed (Fig. 4.11b, the
% norm of the torques). The claim is that L1 carries up to 94% of the robot's
% own weight.
%
% ------------------------------------------------------------------ the loads
% SCALED TO THE ROBOT, NOT COPIED. The thesis's 0-30 kg random range and its
% 10 / 15 / 20 kg fixed loads belong to a 32 kg robot. The same shares of this
% 74 kg model are 0-70 kg and 23 / 35 / 46 kg (x 74/32), which are the defaults;
% the thesis's kilograms would test a 41% load and call it 94%.
%
% The load is p.uncertainty.load_mass: a point mass at the torso base, i.e. at
% the hip. It adds to M and G but brings no rotational inertia and does not
% move the torso's centre of mass (ch4_control_affine).
%
% ------------------------------------------------ what a load does, measured
% A load is NOT a small mass scale, and the difference drives this study.
% Along the nominal orbit of posture_195:
%
%   load [kg]                  23      35      46      70
%   ||Delta2||               1.99    2.30    2.47    2.72
%   its isotropic part       0.19    0.23    0.26    0.30
%   min eig(I + Delta2)      0.49    0.41    0.37    0.30
%
% A mass scale has Delta2 = (1/s - 1) I: isotropic, and ||Delta2|| < 1 for any
% s above one half. A hip load changes the decoupling matrix far more in some
% output directions than in others. Every eigenvalue of I + Delta2 stays real
% and positive, so no direction's input gain reverses, but the matrix is
% strongly non-normal (singular values 0.11-2.9 at 70 kg). ||Delta2|| > 1 is
% also outside anything the robust CLF-QP's bound can cover (ch4_ctrl_rclf_qp),
% which is why this study, like the thesis's, is an L1 one. The verbose table
% prints these measurements for the loads it runs.
%
% ------------------------------------------------------------------- the box
% UNDER p.box.rule = 'rating' (the default, ch4_params) every case runs at the
% one actuator rating p.box.rating, which was sized once from the heaviest load
% in the design envelope (70 kg) and not from the case at hand. The controllers
% are never told the load, and now neither is their box. The feedforward peak
% below is still measured and printed, as a diagnostic of how close each load
% runs to the rating.
%
% UNDER 'thesis' (the rule before 2026-09-18) each case gets its own box.
% 'l1_con' bounds the torque it applies (p.l1.constrain_applied), so its box
% must be able to carry the loaded robot, and a hip load does not raise the
% torque demand in proportion to mass. The box is therefore the gait box's
% headroom over its feedforward, p.l1.u_max / p.gait_u_peak (1.25), times the
% LOADED robot's own feedforward peak along the nominal orbit, at the heaviest
% load the case applies:
%
%   load [kg]                   0      23      35      46      70
%   true feedforward peak     195     288     330     366     445   Nm
%   box                       244     360     412     458     556   Nm
%
% For a uniform scale s the same rule gives s times the gait box, which is the
% 'l1' preset's rule in ch4_compare_controllers. Sized by mass ratio instead
% (320 / 359 / 395 / 474 Nm), 'l1_con' fell in step 2 at 70 kg and in step 3
% at 46 kg; at these boxes it walked all 25 steps of both.
%
% ------------------------------------------------------- the fair baseline
% The thesis compares nothing here. This study runs the unconstrained CLF-QP of
% Section 4.2.4 AND 'clfqp_con' at the same box and contact rows as 'l1_con'.
% The second is not decoration: under load the unconstrained baseline's peak
% torque is 2-3.5 times 'l1_con's box (729-1797 Nm), so a comparison against it
% alone cannot say whether the adaptation or the torque is doing the work.
% Measured over 25 steps: 'clfqp_con' fell in every case, within 3-10 steps,
% while 'l1_con' at the same box walked every step of every case.
%
% Inputs
%   x0, alpha : an optimized gait
%   p         : parameter struct (p.gait_u_peak filled by ch4_load_gait)
%   opts      : struct overriding
%                 .controllers  default {'clfqp','clfqp_con','l1','l1_con'}
%                 .range        random load range [kg], default [0 70]; [] skips
%                 .loads        fixed loads [kg], default [23 35 46]
%                 .X_orbit      nominal-orbit states to size the boxes along
%                               (default: every other state of a 2-step
%                               nominal rollout)
%                 .n_steps .verbose .store_traj   as ch4_compare_controllers
%
% Output
%   S : struct array of ch4_run_entry rows, one per (case, controller), with
%       .u_box and .box_violation set and .case_label naming the load. Random
%       cases come first and carry load_mass = NaN.
%
% See also CH4_RUN_ENTRY, CH4_PLOT_LOAD, CH4_SIMULATE, CH4_CONTROL_AFFINE.

if nargin < 4, opts = struct(); end
opts = fill_defaults(opts, struct('controllers', {{'clfqp', 'clfqp_con', 'l1', 'l1_con'}}, ...
                                  'range', [0 70], 'loads', [23 35 46], ...
                                  'X_orbit', [], 'n_steps', 25, ...
                                  'verbose', true, 'store_traj', true));

if isempty(p.gait_u_peak)
    error('ch4_load_study:noGait', ...
          ['p.gait_u_peak is empty. Build p with ch4_load_gait, which sizes ' ...
           'the L1 box to the gait this study scales it from.']);
end
headroom = p.l1.u_max / p.gait_u_peak;
rule     = ch4_box_rule(p);

X = opts.X_orbit;
if isempty(X)
    pb = p;
    pb.controller        = 'clfqp';
    pb.uncertainty       = struct('mass_scale', 1, 'load_mass', 0);
    pb.load_random_range = [];
    sim_b = ch4_simulate(x0, alpha, pb, 2);
    X = sim_b.x(:, 1:2:end);
end

%% --- the cases, each with its box ----------------------------------------
cases = struct('label', {}, 'load_mass', {}, 'range', {}, 'box', {}, ...
               'ff_peak', {}, 'n2', {}, 'n2_iso', {}, 'eig_min', {});
if ~isempty(opts.range)
    m = measure(X, alpha, p, max(opts.range));
    cases(end+1) = struct('label', sprintf('random %g-%g kg', opts.range), ...
                          'load_mass', 0, 'range', opts.range, ...
                          'box', case_box(rule, p, headroom, m), 'ff_peak', m.ff_peak, ...
                          'n2', m.n2, 'n2_iso', m.n2_iso, 'eig_min', m.eig_min);
end
for mL = opts.loads(:).'
    m = measure(X, alpha, p, mL);
    cases(end+1) = struct('label', sprintf('%g kg', mL), ...
                          'load_mass', mL, 'range', [], ...
                          'box', case_box(rule, p, headroom, m), 'ff_peak', m.ff_peak, ...
                          'n2', m.n2, 'n2_iso', m.n2_iso, 'eig_min', m.eig_min); %#ok<AGROW>
end

if opts.verbose
    l1o = ch4_l1_opts(p);
    fprintf('\n%s\n CHAPTER 4 -- unknown load (Fig. 4.11), %d steps per run\n', ...
            repmat('=',1,105), opts.n_steps);
    if isfinite(l1o.loop_gain_max)
        nrm = sprintf('adaptation normalized above %.2f rad/sample', ...
                      sqrt(l1o.loop_gain_max) * p.control_dt);
    else
        nrm = 'adaptation not normalized';
    end
    fprintf(' L1: Gamma %.0e, filter %.0f rad/s, predictor %s, %s\n', ...
            p.l1.Gamma, p.l1.omega_c, l1o.predictor, nrm);
    fprintf(' box rule: %s\n', rule);
    fprintf(' %-16s %14s %10s %10s %14s %9s\n', 'case', 'feedfwd peak', ...
            '||Delta2||', 'isotropic', 'min eig(I+D2)', 'box');
    for c = cases
        fprintf(' %-16s %11.1f Nm %10.3f %10.3f %14.3f %6.0f Nm\n', c.label, ...
                c.ff_peak, c.n2, c.n2_iso, c.eig_min, c.box);
    end
    fprintf('%s\n', repmat('-',1,105));
    fprintf(' %-9s %-16s %6s %6s %9s %10s %8s %7s %7s  %s\n', 'controller', ...
            'load', 'box', 'steps', 'peak|u|', 'max|eta|', 'min Fz', 'Fz<0 %', ...
            'QPfail', 'reason');
    fprintf('%s\n', repmat('-',1,105));
end

%% --- the runs -------------------------------------------------------------
S = [];
for c = cases
    for k = 1:numel(opts.controllers)
        pc = ch4_run_params(p, opts.controllers{k}, 1, c.box);
        pc.uncertainty.load_mass = c.load_mass;
        pc.load_random_range     = c.range;

        e = ch4_run_entry(x0, alpha, pc, opts);
        e.u_box = c.box;
        if ~isnan(e.peak_torque)
            e.box_violation = max(e.peak_torque - c.box, 0);
        end
        e.case_label = c.label;
        S = [S, e]; %#ok<AGROW>

        if opts.verbose
            fz_neg = NaN;
            if ~isempty(e.traj), fz_neg = 100 * mean(e.traj.lambda(2, :) < 0); end
            fprintf(' %-9s %-16s %6.0f %6d %9.1f %10.3f %8.0f %7.1f %7d  %s\n', ...
                    e.name, c.label, e.u_box, e.steps_completed, e.peak_torque, ...
                    e.max_eta, e.Fz_min, fz_neg, e.qp_infeasible, e.reason);
        end
    end
    if opts.verbose, fprintf('%s\n', repmat('-',1,105)); end
end

if opts.verbose
    fprintf([' "box" binds on clfqp_con and l1_con only. "min Fz" and "Fz<0" are\n' ...
             ' the TRUE normal force under the load each step actually carried.\n']);
    fprintf('%s\n\n', repmat('=',1,105));
end

end

% ---------------------------------------------------------------------------
function b = case_box(rule, p, headroom, m)
%CASE_BOX  The box one load case runs at, under the rule in force.
if strcmp(rule, 'rating')
    b = p.box.rating;
else
    b = headroom * m.ff_peak;
end
end

% ---------------------------------------------------------------------------
function m = measure(X, alpha, p, mL)
%MEASURE  What carrying mL does along the orbit: feedforward peak and Delta2.
unc = struct('mass_scale', 1, 'load_mass', mL);
m = struct('ff_peak', 0, 'n2', 0, 'n2_iso', 0, 'eig_min', inf);
for k = 1:size(X, 2)
    D = ch4_uncertainty(X(:, k), alpha, p, unc);
    m.ff_peak = max(m.ff_peak, max(abs(D.LgLfy_true \ D.Lf2y_true)));
    m.n2      = max(m.n2, D.n2);
    m.n2_iso  = max(m.n2_iso, D.n2_scalar);
    m.eig_min = min(m.eig_min, min(real(eig(eye(p.ny) + D.Delta2))));
end
end

function s = fill_defaults(s, d)
% Fills MISSING fields only -- unlike ch4_compare_controllers', which also
% replaces empty ones -- so an explicit opts.range = [] skips the random case.
f = fieldnames(d);
for i = 1:numel(f)
    if ~isfield(s, f{i})
        s.(f{i}) = d.(f{i});
    end
end
end
