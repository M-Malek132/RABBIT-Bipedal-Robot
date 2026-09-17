function entry = ch4_run_entry(x0, alpha, pc, opts)
%CH4_RUN_ENTRY  Simulate one configured Chapter-4 run and score it as a table row.
%
%   entry = ch4_run_entry(x0, alpha, pc, opts)
%
% The one place a comparison row is computed, so the mass-scale sweeps
% (ch4_compare_controllers) and the load study (ch4_load_study) score their runs
% the same way. pc arrives fully configured -- law, perturbation, boxes, as
% ch4_run_params and the caller set them -- and nothing here changes it.
%
% ------------------------------------------- forces under a load that changes
% ch4_forces evaluates the TRUE contact force and the true uncertainty with ONE
% model, pc.uncertainty, for the whole trajectory it is handed. When
% pc.load_random_range redraws the carried mass at every step, that model is
% the right one for no step in particular, and min Fz -- the column that says
% whether a run is walking at all -- would be computed for a robot that was
% never simulated. So a randomized run is analysed step by step, each step
% under the load it was simulated with (sim.loads), and the pieces are joined.
%
% Inputs
%   x0, alpha : an optimized gait
%   pc        : parameter struct for this run
%   opts      : struct with .n_steps .store_traj
%
% Output
%   entry : .name .mass_scale .u_box .steps_completed .peak_torque
%           .box_violation .max_eta .final_eta .V0 .Vend .V_ratio .delta_max
%           .qp_infeasible .int_u2 .Fz_min .grf_pred_error .reason
%           .step_T .step_L .load_mass .loads .traj
%           .valid_steps .invalid_step .invalid_kind .frac_invalid .mu_max
%           the physical-validity score (ch4_validity): steps completed
%           before the first sample at which the true contact force lifts
%           off or exceeds the friction coefficient, where that sample is,
%           and how much of the run violated the contact at all.
%           u_box and box_violation stay NaN for the caller, which knows the
%           box. load_mass is NaN when the load was redrawn every step, and
%           loads lists the carried mass on every attempted step.
%
% See also CH4_COMPARE_CONTROLLERS, CH4_LOAD_STUDY, CH4_FORCES, CH4_SIMULATE.

sim = ch4_simulate(x0, alpha, pc, opts.n_steps);

randomized = ~isempty(pc.load_random_range);
if randomized
    load_mass = NaN;
else
    load_mass = pc.uncertainty.load_mass;
end

entry = struct('name', pc.controller, 'mass_scale', pc.uncertainty.mass_scale, ...
               'u_box', NaN, 'steps_completed', sim.n_ok, ...
               'peak_torque', NaN, 'box_violation', NaN, ...
               'max_eta', NaN, 'final_eta', NaN, ...
               'V0', NaN, 'Vend', NaN, 'V_ratio', NaN, ...
               'delta_max', NaN, 'qp_infeasible', NaN, ...
               'int_u2', NaN, 'Fz_min', NaN, 'grf_pred_error', NaN, ...
               'reason', sim.reason, ...
               'step_T', [sim.steps.T], 'step_L', [sim.steps.L_step], ...
               'load_mass', load_mass, 'loads', sim.loads, ...
               'traj', [], ...
               'valid_steps', NaN, 'invalid_step', NaN, 'invalid_kind', '', ...
               'frac_invalid', NaN, 'mu_max', NaN);

% At full solver resolution and under each step's own load: ch4_step records
% the true contact force as it integrates, so this needs no re-analysis.
Vd = ch4_validity(sim, pc);
entry.valid_steps  = Vd.valid_steps;
entry.invalid_step = Vd.first_step;
entry.invalid_kind = Vd.first_kind;
entry.frac_invalid = Vd.frac_invalid;
entry.mu_max       = Vd.mu_max;

if sim.n_ok == 0, return; end

% ch4_forces decimates to a fixed number of points, 2000 by default -- about
% 670 a step at the three-step horizon it was set for. Keep that density per
% step rather than per run, or a 25-step run is analysed on a grid eight times
% coarser and its min Fz can step over a contact-force dip the short run would
% have caught.
if randomized
    F = forces_by_step(sim, alpha, pc, ceil(2000 / 3));
else
    F = ch4_forces(sim.t, sim.x, alpha, pc, sim.xi, sim.t_xi, ...
                   ceil(2000 * max(opts.n_steps, 3) / 3));
end

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
function F = forces_by_step(sim, alpha, pc, per_step)
%FORCES_BY_STEP  ch4_forces over each completed step under its own load, joined.
parts = cell(1, sim.n_ok);
t_off = 0;
for k = 1:sim.n_ok
    st = sim.steps(k);
    pk = pc;
    pk.uncertainty.load_mass = sim.loads(k);
    pk.load_random_range     = [];
    parts{k} = ch4_forces(t_off + st.t, st.x, alpha, pk, ...
                          st.xi, t_off + st.t_xi, per_step);
    t_off = t_off + st.T;
end
P = [parts{:}];
F = struct('t', [P.t], 'x', [P.x], 'u', [P.u], 'y', [P.y], 'ydot', [P.ydot], ...
           'V', [P.V], 'lambda', [P.lambda], ...
           'theta_hat', [P.theta_hat], 'theta_true', [P.theta_true], ...
           'torque_max', max([P.torque_max]), ...
           'delta_max', max([P.delta_max]), ...
           'qp_infeasible', sum([P.qp_infeasible]), ...
           'int_u2', sum([P.int_u2]), ...
           'Fz_min', min([P.Fz_min]), ...
           'grf_pred_error', max([P.grf_pred_error]));
end
