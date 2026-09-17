function p = ch4_params(varargin)
%CH4_PARAMS  Single source of truth for the Chapter-4 pipeline.
%
%   p = ch4_params()                 Chapter-3 defaults plus Chapter-4 fields
%   p = ch4_params('name',value,...) with overrides (any ch3 or ch4 field)
%
% Chapter 4 changes ONE thing about Chapter 3 and then builds two controllers
% around that change: the controller no longer knows the plant.
%
%   f, g          TRUE model   -- what the robot actually is
%   ftil, gtil    NOMINAL model -- what the controller was designed against
%
% Everything in Chapter 3 assumed these coincide. Here they do not, and
% p.uncertainty is what separates them: the simulation integrates the TRUE
% model while every controller is built from the NOMINAL one. Set the
% uncertainty to zero and Chapter 4 collapses back onto Chapter 3 exactly --
% that identity is asserted in ch4_test_model rather than merely claimed.
%
% THE TWO ANSWERS.
%
%   'rclfqp'   Section 4.1.  Assume the uncertainty is bounded, and satisfy the
%              RES-CLF condition for the WORST case inside that bound. Buys a
%              guarantee; pays for it by being aggressive even when the model
%              happens to be perfect.
%   'l1'       Section 4.2.  Estimate the uncertainty online and cancel it.
%              Pays nothing when the model is perfect (it estimates zero), but
%              offers a bounded-error guarantee rather than a worst-case one.
%
% See also CH3_PARAMS, CH4_CTRL_RCLF_QP, CH4_CTRL_L1, CH4_MAIN.

p = ch3_params();

%% ---------------------------------------------------- CLF convergence rate
% Chapter 3 ships eps = 0.5, chosen there so the PD baseline stays inside
% RABBIT's torque envelope. Chapter 4 tightens it twice, for reasons specific
% to what this chapter measures.
%
% eps sets the required convergence rate (c3/eps): it has to be fast enough
% that the outputs re-converge between impacts, since each footstrike expands
% eta and the controller gets exactly one step to beat that expansion. At
% eps = 0.5 the Chapter-3 min-norm CLF-QP does NOT manage it even with a
% PERFECT model. Measured on posture_195 at 1 kHz over four steps, the
% post-impact error ||eta+|| reads 0.34, 0.32, 0.19, 0.39 at eps = 0.5 --
% never settling -- against 0.24, 0.09, 0.19, 0.01 at eps = 0.35. (Some of
% that per-impact kick is the 1 kHz sample-and-hold itself, not the model:
% ||eta+|| after the first step is 0.24 / 0.12 / 0.05 at a 1 / 0.5 / 0.2 ms
% period, with a perfect model and continuous control giving 1.7e-4.)
%
% That matters here more than it did in Chapter 3, because the L1 controller's
% REFERENCE MODEL IS THAT CONTROLLER (Section 4.2.2). L1 promises to make the
% perturbed system behave like the reference model; if the reference model is
% itself impact-marginal, L1 faithfully reproduces marginal behaviour and the
% Section 4.2.4 comparison measures the reference model rather than the
% adaptation.
%
% 0.35 WAS STILL TOO SLOW FOR THE ROBUST LAW, AND ONLY A LONG RUN SHOWS IT.
% Over three steps 'rclfqp_con' looked converged. Over 25 it drifted in every
% case, and in Case I -- a perfect model -- it fell in step 21 (max||eta|| 17 /
% 8.3 / 11.4 in Cases I-III). The drift is the Delta1 term against the slow
% rate: it survived the boundary layer, the contact rows, the torque box and
% D2 = 0, and vanished with D1 = 0. At eps = 0.20 (ch4_main, 25 steps, kappa =
% 1) the same law walks all 25 steps in every case, max||eta|| 1.13 / 4.11 /
% 4.38, and the nominal CLF-QP stays on the orbit (max||eta|| 0.25). The drift
% is slower there, not gone: the per-step CLF peak still grows ~70-500x over
% the run. The baselines change character too: under perturbation they no
% longer fall within 25 steps but track at max||eta|| 14-17 with erratic step
% times -- degraded walking instead of a fall, which is still the contrast
% Remark 4.6 draws.
%
% Raise it back to 0.5 to reproduce the Chapter-3 defaults exactly; expect the
% baselines to fail earlier if you do.
p.eps = 0.20;

%% ------------------------------------------------------------- controller
% Extends the Chapter-3 list (ff | iolin_pd | clfqp | clfqp_con) with:
%
%   'rclfqp'      Section 4.1.2, robust CLF-QP, eq (4.12)
%   'rclfqp_con'  Section 4.1.3, robust CLF-QP + constraints, eq (4.13)
%   'l1'          Section 4.2.2, L1 adaptive with a CLF-QP reference model
%   'l1_con'      Section 4.2.3, the same with torque saturation, eq (4.36-37)
%
% The Chapter-3 names still work and still mean the same thing; under nonzero
% uncertainty they are the BASELINES the chapter measures against (controller
% A in both Section 4.1.4 and Section 4.2.4).
p.controller = 'clfqp';

%% ------------------------------------------------- the true-vs-nominal gap
% How the TRUE plant differs from the NOMINAL model the controller holds.
% The nominal model is untouched either way -- see ch4_control_affine.
%
%   mass_scale  Multiply every link mass and inertia by this factor. Sections
%               4.1.4 and 4.2.4 both sweep it: 1 (no uncertainty), 1.5, 0.7,
%               and 3 for the extreme Case IV. M, the Coriolis vector and
%               gravity are all LINEAR in the mass parameters, so this scales
%               all three by the same factor and nothing else. Each scale is
%               nonetheless re-derived from scratch rather than applied as sM
%               (ch4_case_dynamics; registered: 0.5, 0.7, 1.5, 3), and
%               ch4_test_model's KKT-split check is what confirms the two agree.
%
%   load_mass   A point mass rigidly attached at the torso base [kg], for the
%               "carrying unknown mass on the torso" study in Section 4.2.4
%               (load randomised in 0-30 kg per step, up to 94% of body mass).
%               Unlike mass_scale this is NOT a uniform scaling: it changes the
%               mass matrix non-uniformly, so it also changes the impact map.
p.uncertainty = struct('mass_scale', 1.0, ...
                       'load_mass',  0.0);

% Section 4.2.4's load experiment redraws the carried mass at every step from
% this range. Empty disables the redraw and load_mass stays fixed.
p.load_random_range = [];        % e.g. [0 30] for the Fig. 4.11a experiment

% Seed of that redraw. Every controller in a comparison sees the same sequence;
% changing the seed is how a result is checked against more than one sequence.
p.load_seed = 11;

%% --------------------------------------------- robust CLF-QP (Section 4.1)
% Bounds (4.10) on the induced uncertainty (4.4).  These are bounds on the
% I/O-LINEARIZED uncertainty, not on the mass error: Delta1 has units of
% ydd (rad/s^2) and Delta2 is dimensionless.
%
% HOW TO CHOOSE THEM.  Do not guess. ch4_uncertainty measures the actual
% Delta1, Delta2 along a gait for a given mass_scale, and ch4_delta_bounds
% turns that measurement into these two numbers. The defaults below are that
% measurement for the default gait (posture_195), along a nominal CLF-QP
% rollout, over Cases I-III (mass scale 1, 1.5, 0.7), with a 1.2 safety
% factor: max ||Delta1|| 232.6 and ||Delta2|| 0.4286, both from Case III.
% They matter only when p is used standalone -- ch4_main re-measures on
% whatever gait it loads and adopts the result unless these are set
% explicitly, since a bound fitted to one gait need not fit another (the
% previous defaults, 250 and 0.45, were fitted to the old upright gait and
% leave this one 7% and 5% of margin instead of the intended 20%).
%
% WHAT THE MEASUREMENT REVEALS, and it is worth knowing before tuning these.
% For a uniform mass/inertia scale s, the constrained dynamics split exactly
% (derived and asserted in ch4_test_model):
%
%       Delta2 = (1/s - 1) I,        Delta1 = -(1/s - 1) Lftil^2 y.
%
% Three consequences:
%   * delta2_max is not an estimate at all. For Case III (s = 0.7) it is
%     |1/0.7 - 1| = 0.4286 exactly, and for Case IV (s = 3) it is 0.6667.
%   * ||Delta2|| < 1 -- the feasibility condition of the robust QP -- holds iff
%     s > 0.5. Worst-case robustness cannot reach below half the nominal mass,
%     however the bounds are tuned.
%   * Delta1 is proportional to the OUTPUT DRIFT Lftil^2 y, which grows away
%     from the orbit. A single constant delta1_max valid over a neighborhood is
%     therefore several times larger than the value valid on the orbit (measured
%     on posture_195: 232.6 along the rollout, 1140.6 with the sampled states
%     jittered by sd 0.05 -- 4.9x). Since the
%     commanded ||mu|| scales as delta1_max/(1 - delta2_max), that difference is
%     the whole of the "unnecessarily aggressive" limitation Section 4.1.4
%     closes on -- in a form you can put a number to.
p.rclf = struct();
p.rclf.delta1_max = 279.1;       % ||Delta1|| bound          [rad/s^2]
p.rclf.delta2_max = 0.5143;      % ||Delta2|| bound          [-]

% HOW THE MAX IN (4.11) IS TAKEN OVER THE Delta2 BALL.
%
%   'scalar'  Delta2 = d2 * I with |d2| <= delta2_max.  The worst case of the
%             term LgV*Delta2*mu is then delta2_max*|LgV*mu|, and a bound on an
%             absolute value is TWO LINEAR INEQUALITIES -- one per sign. This
%             is the reduction Chapter 4 relies on (Remark 4.5) to keep the
%             min-max problem an honest quadratic program.
%
%   'matrix'  The unstructured ball ||Delta2||_2 <= delta2_max, which is the
%             literal reading of (4.10).  The worst case is then
%             delta2_max*||LgV||*||mu||, a SECOND-ORDER CONE constraint, not a
%             linear one. ch4_ctrl_rclf_qp solves the unconstrained case of
%             this EXACTLY in closed form; the constrained case falls back on
%             ||mu||_2 <= ||mu||_1, which is linear-representable and
%             conservative by at most sqrt(ny) = 2.
%
% 'scalar' is the default because it is what makes (4.12) a QP. Use 'matrix'
% to check how much of the guarantee rests on that structural assumption.
p.rclf.delta2_model = 'scalar';

% BOUNDARY LAYER on the robust term (kappa in ch4_ctrl_rclf_qp).
%
% Near the orbit the exact robust law applies a correction of fixed magnitude
% D1/(1 - D2) along -LgV'/||LgV||, a unit vector that flips sign whenever LgV
% crosses zero: sliding-mode control, which a 1 kHz sample-and-hold turns into
% chatter at the sample rate. The layer saturates that unit vector over a band
% of ||LgV||, measured in units of the thickness at which one held sample of the
% correction, on the worst-case input gain 1 + D2, lands exactly on LgV = 0.
% One sample then multiplies LgV by at worst 1 - 1/kappa, so
%
%   kappa = 0        the exact law of (4.12)/(4.13), chatter included
%   kappa > 1/2      the sampled loop along LgV converges for every d2 in bound
%   kappa >= 1       ... and without overshoot, i.e. without chatter
%
% The price is the guarantee: (4.11) still holds exactly outside the layer, and
% inside it the worst case may exceed it by at most D1*phi/4, which bounds the
% tracking error instead of driving it to zero. control_dt = 0 disables the
% layer, since continuous control has no sample to overshoot.
p.rclf.boundary_layer = 1;

% CONTACT ROWS: p.limits.enable.friction / .grf stay ON, as ch3_params sets
% them. In the constrained CLF-QPs they add the friction cone and the
% normal-force floor of (4.13), written on the NOMINAL model (Remark 4.4). They
% are not optional decoration: without the floor, the exact robust law's
% chattering worst-case term pulled the true stance foot into the ground for
% 17-42% of the samples in each case, and even with the boundary layer above,
% 25 steps at mass scale 0.7 still went negative on 1.8% -- see ch4_run_params.

%% ------------------------------------------ L1 adaptive control (Sect. 4.2)
p.l1 = struct();

% STATE PREDICTOR.  Section 4.2.2's predictor (4.19) is driven by its own copy
% of the reference model, eta_hat_dot = F eta_hat + G mu1_hat + G(mu2 +
% theta_hat), and the error dynamics (4.28) it relies on do not follow from it:
% the min-norm CLF-QP is nonlinear, so mu1(eta_hat) - mu1(eta) reaches the
% prediction error and the adaptation reads it as model error (ch4_l1_deriv
% has the derivation).
%
%   'thesis'  (4.19) as written, with its sampled-data advance as written.
%             Reproduces the Section 4.2 results in Results/ stamped before
%             this option existed, with Gamma = 1e4, Gamma_alpha = [] and
%             constrain_applied = false as well.
%   'plant'   the standard L1 state predictor, driven by what the robot
%             actually received:
%                 eta_hat_dot = F eta + G(mu + theta_hat) - a(eta_hat - eta)
%             so eta_tilde_dot = -a eta_tilde + G theta_tilde, exactly. Under
%             sampled control it holds mu and reads eta at both ends of each
%             period, as the plant does (ch4_l1_advance).
p.l1.predictor = 'plant';

% PREDICTOR RATE a [rad/s], 'plant' only.  With the plant predictor the beta
% channel of the estimator closes the loop s^2 + a s + Gamma, so a sets its
% damping: a/(2 sqrt(Gamma)); 2 sqrt(1e5) = 632 is critical.
%
% 800 (damping 1.26), because critical damping sits on an edge at mass scale
% 1.5. Measured over 25 steps (eps 0.20): 'l1' fell in step 16 at a = 500 and
% at a = 632, and walked all 25 steps at a = 700 / 800 / 900 with max||eta||
% 2.2 / 2.2 / 2.3. 632.46 instead of 632 was already enough to flip it. At
% a = 800, scales 1 and 0.7 read 0.14 and 4.8, as good as at 632.
% 'l1_con' at 1.5x stays fragile at every rate: it walked at 500, 632, 800 and
% 900 but fell in step 12 at 700 -- see CH4_UNCERTAINTY.md §5a.
p.l1.predictor_rate = 800;

% ALPHA'S OWN ADAPTATION GAIN; [] uses Gamma.  alpha_hat's regressor is
% ||eta||, so its loop gain is Gamma_alpha*||eta||^2, which grows exactly when
% the tracking degrades. 0 freezes alpha_hat at zero and leaves beta_hat to
% carry theta, which (4.17) says it can. Both estimates stay on by default:
% under the plant predictor, beta alone tracked worse in every case measured
% (25 steps, eps 0.20) -- at Gamma = 1e4 'l1' fell at 1.5x in step 10 with
% beta alone and walked all 25 steps at max||eta|| 2.5 with both.
p.l1.Gamma_alpha = [];

% LEAKAGE ON alpha_hat [1/s], a sigma-modification: alpha_hat_dot gains
% -alpha_leak * alpha_hat. 0 is the adaptation law as (4.26) writes it. See
% ch4_l1_deriv for why alpha_hat, and only alpha_hat, would need it.
%
% OFF BY DEFAULT: over the long runs (CH4_UNCERTAINTY.md §5a) a leak of 2, 10
% or 50 /s did not prevent the falls, and at 1.5x l1_con fell sooner at every
% rate -- alpha_hat does real work within a step there. Kept so that
% measurement can be repeated.
p.l1.alpha_leak = 0;

% WHAT THE ESTIMATES DO AT A FOOTSTRIKE: 'carry' | 'continuous' | 'fold'
% (ch4_l1_state). ||eta|| jumps at impact, and carrying alpha_hat and beta_hat
% unchanged makes theta_hat jump by alpha_hat times that jump. 'carry' is the
% law as Section 4.2 writes it and the default: keeping theta_hat continuous
% made the long runs fall sooner, not later. Kept for the same reason as the
% leak.
p.l1.impact_estimate = 'carry';

% CAP ON alpha's REGRESSOR, as the fastest the estimator loop may run [rad per
% control sample]; 0 or [] for no cap, the law as written. See ch4_l1_deriv.
%
% OFF BY DEFAULT BECAUSE IT TRADES, measured over 60-120 steps at 1 kHz (plain
% l1 at 1.5x; l1_con at 1, 0.7, 1.5, 46 kg and four random 0-70 kg sequences):
% runs falling, of nine -- 5 uncapped, 3 at 1.5, 2 at 1, 1 at 0.5 -- while at
% 0.5 l1_con's error at 1.5x rises from 2.0-4.5 to 6.0-7.4 per 10-step block
% and the random-load runs that survive pass through excursions to
% max||eta|| 18-48. CH4_UNCERTAINTY.md §5a has the rest.
p.l1.alpha_regressor_rate = 0;

% NORMALIZED ADAPTATION, the other limit on the estimator loop's speed [rad per
% control sample]; 0 or [] for none, the law as written. Both adaptation laws
% are divided by m^2 = max(1, (Gamma + Gamma_alpha ||eta||^2) / (rate/dt)^2),
% so the loop never runs faster than the rate while theta_hat keeps
% alpha_hat*||eta|| in full. At or below sqrt(Gamma)*dt (0.32 here) it is the
% textbook form m^2 = 1 + (Gamma_alpha/Gamma) ||eta||^2. See ch4_l1_deriv.
%
% 0.75, THE BEST 1 kHz REMEDY MEASURED, over the cap's nine long runs plus six
% more random 0-70 kg sequences (seeds 4-9): runs falling, of 15 -- 9 as
% written, 4 with the cap at 0.5, 1 at 0.75, 2 at 1. It is a window: at 0.5,
% 1.5 and 2 two of the first nine fall again, and the textbook form loses
% seven of them. At 0.75 'l1_con' at 1.5x keeps its typical step (median
% per-step max||eta|| 2.15, against 1.85 as written and 4.44 capped), but
% random-load runs can pass through excursions to max||eta|| ~50. A 0.5 ms
% control period without it fell in 1 of 12 and tracked flatter.
% CH4_UNCERTAINTY.md §5a has the tables.
p.l1.normalized_rate = 0.75;

% WHAT 'l1_con' CONSTRAINS.  false: Section 4.2.3 as written -- the torque box
% on mu1 alone, no contact rows, so mu2 is applied outside every constraint.
% true: mu2 enters the QP as a known offset, and the box and the friction and
% normal-force rows (as p.limits.enable has them) bound the TOTAL torque the
% robot receives. See ch4_ctrl_l1.
p.l1.constrain_applied = true;

% ADAPTATION GAIN Gamma in (4.26).  The bound (4.34)-(4.35) on the estimation
% error shrinks like 1/||Gamma||, so "sufficiently large" is the whole design
% rule -- the estimate is allowed to be as fast as the integrator can follow.
% Decoupling estimation speed from control smoothness is exactly what the
% low-pass filter below is for, so this can be large without putting
% high-frequency content into the torque.
%
% 1e5, not the 1e4 the thesis form ran at. Under the plant predictor the beta
% channel's bandwidth is sqrt(Gamma): 100 rad/s at 1e4, below the 150 rad/s
% filter it feeds, 316 rad/s at 1e5. Measured over 25 steps at eps 0.20, with
% the box scaled to the robot and the predictor critically damped, 'l1_con'
% at scales 0.7 / 1.5 went from max||eta|| 6.2 / 35.5 at 1e4 (the 1.5x
% estimate running away) to 5.1 / 7.6 at 1e5, true normal force positive
% throughout. Set it back to 1e4 to reproduce the thesis form.
p.l1.Gamma = 1e5;

% LOW-PASS FILTER C(s) in (4.23), first order with unit DC gain:
%   mu2 = -C(s) theta_hat,   C(s) = omega_c / (s + omega_c).
% Section 4.2.4 uses 150 rad/s (~23 Hz).
%
% This filter is not a detail. theta_hat is deliberately fast and therefore
% ragged; feeding it straight into the torque would put high-frequency content
% on the ground reaction force and break the unilateral contact constraint --
% the robot would chatter its foot off the ground. The filter is what lets the
% estimator be fast AND the control be smooth.
p.l1.omega_c = 150;

% PROJECTION BOUNDS for (4.26).  The projection operator keeps the estimates
% inside these balls, which is what makes alpha_tilde, beta_tilde bounded in
% (4.32) and hence the whole error bound (4.35) finite. Same units as Delta1.
p.l1.alpha_max = 200;            % ||alpha_hat||_2 <= alpha_max
p.l1.beta_max  = 400;            % ||beta_hat||_2  <= beta_max
p.l1.proj_eps  = 0.1;            % smoothing band of the projection, in (0,1]

% TORQUE SATURATION for 'l1_con', eq (4.36)-(4.37). Section 4.2.4 uses 65 Nm
% on all four joints.
%
% Note the scope, stated in Section 4.2.3: the saturation is imposed on the
% CLF-QP component mu1 only, NOT on the adaptive component mu2. The realized
% torque can therefore leave the box by whatever mu2 contributes; ch4_forces
% reports that overshoot rather than hiding it. That is the thesis form, with
% p.l1.constrain_applied = false; the default bounds the applied torque.
%
% 65 Nm is below the peak torque of every gait in Results/ (195-465 Nm), so
% ch4_load_gait raises it to 1.25x the gait's own peak (243.8 Nm on
% posture_195); see the note there for what a starved box does to the
% adaptation, and why the headroom. This value stands for a gait that fits
% inside it.
p.l1.u_max = 65;

% Peak joint torque of the loaded gait, filled in by ch4_load_gait. Empty
% when p is built standalone, since a bare parameter struct has no gait.
p.gait_u_peak = [];

% TORQUE BOX OF THE CONSTRAINED LAWS IN THE SWEEPS.  How ch4_compare_controllers
% and ch4_load_study size the box for clfqp_con, rclfqp_con and l1_con.
%
%   'rating'  ONE ACTUATOR RATING for every controller, mass scale and load, as
%             a hardware torque limit would be. It is sized once, from the
%             envelope the robot is designed for, and never from the case being
%             run: 1.25x (the L1 headroom, see ch4_load_gait) the heaviest
%             feedforward peak in that envelope on posture_195 -- mass scales
%             up to 1.5 (293 Nm) and hip loads up to 70 kg (445 Nm) -- so
%             1.25 x 445 = 556 Nm. Case IV (scale 3, feedforward 585 Nm) lies
%             outside that envelope and runs as its own experiment at
%             1.25 x 585 = 731 Nm. Nothing about a case, including what
%             another controller drew under it, reaches its box.
%   'thesis'  the per-case rules this chapter used before 2026-09-18, which
%             size each box FROM the perturbation: robust / case4
%             max(0.8 x clfqp's peak under that case, s x feedforward peak);
%             l1 244 x max(1, s); load 1.25 x the loaded robot's feedforward.
%             Kept to reproduce the earlier tables.
%
% The two ratings are for posture_195. Another gait needs its own envelope.
p.box.rule         = 'rating';
p.box.rating       = 556;
p.box.rating_case4 = 731;

% PREDICTOR RESET AT IMPACT.  eta jumps discontinuously at every footstrike,
% so eta_hat must be told about it or eta_tilde = eta_hat - eta would register
% the impact as a huge phantom uncertainty and the adaptation would chase it.
%   true   re-seed eta_hat = eta+ at the start of each step (eta_tilde = 0).
%   false  carry eta_hat through unchanged.
% The estimates alpha_hat, beta_hat are ALWAYS carried across the impact: they
% describe a property of the robot, which the footstrike does not change.
p.l1.reset_predictor = true;

%% ------------------------------------------------------------ integration
% Chapter 4 defaults to sampled-data control, unlike Chapter 3.
%
% Two independent reasons, either one sufficient. First, the same argument as
% p.control_dt in ch3_params: the constrained QPs are only piecewise smooth in
% x and stall an adaptive solver. Second, and specific to this chapter, the L1
% controller HAS INTERNAL STATE -- predictor, estimates, filter -- integrated
% alongside the plant, and with Gamma ~ 1e4 that state is far stiffer than the
% robot. Holding the control over a fixed period keeps the plant integration
% from being dragged down to the adaptation timescale.
%
% 1 kHz matches the chapter's framing of a QP solved "well above 1 kHz".
p.control_dt = 1e-3;

%% -------------------------------------------------------------- overrides
% Applied here rather than deferred to ch3_params so that ch4-only fields can
% be overridden too. Nested fields are addressable with dots, e.g.
%   ch4_params('uncertainty.mass_scale', 1.5, 'l1.omega_c', 100)
for k = 1:2:numel(varargin)
    p = set_field(p, varargin{k}, varargin{k+1});
end

p.n_ctrl = p.bez_deg + 1;
p.ny     = numel(p.iact);

end

% ---------------------------------------------------------------------------
function s = set_field(s, name, value)
parts = strsplit(name, '.');
if ~isfield(s, parts{1})
    error('ch4_params:unknownField', 'Unknown parameter "%s".', name);
end
if numel(parts) == 1
    s.(parts{1}) = value;
else
    s.(parts{1}) = set_field(s.(parts{1}), strjoin(parts(2:end), '.'), value);
end
end
