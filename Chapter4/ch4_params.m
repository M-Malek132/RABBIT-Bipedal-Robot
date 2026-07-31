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
% RABBIT's torque envelope. Chapter 4 tightens it, and the reason is specific
% to what this chapter measures.
%
% eps sets the required convergence rate (c3/eps): it has to be fast enough
% that the outputs re-converge between impacts, since each footstrike expands
% eta and the controller gets exactly one step to beat that expansion. At
% eps = 0.5 the Chapter-3 min-norm CLF-QP does NOT manage it on this gait even
% with a PERFECT model -- measured step times 0.370, 0.383, 0.445 s against a
% nominal 0.368, i.e. drifting away rather than settling.
%
% That matters here more than it did in Chapter 3, because the L1 controller's
% REFERENCE MODEL IS THAT CONTROLLER (Section 4.2.2). L1 promises to make the
% perturbed system behave like the reference model; if the reference model is
% itself impact-marginal, L1 faithfully reproduces marginal behaviour and the
% Section 4.2.4 comparison measures the reference model rather than the
% adaptation. At eps = 0.35 the same controller holds 0.370, 0.380, 0.418 and
% the robust controller holds the nominal 0.368 across every perturbation --
% which is the regime the chapter's figures depict.
%
% Raise it back to 0.5 to reproduce the Chapter-3 defaults exactly; expect the
% baselines to fail earlier if you do.
p.eps = 0.35;

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
% Both perturbations are exactly representable, so no dynamics are regenerated
% and the nominal model is untouched -- see ch4_control_affine.
%
%   mass_scale  Multiply every link mass and inertia by this factor. Sections
%               4.1.4 and 4.2.4 both sweep it: 1 (no uncertainty), 1.5, 0.7,
%               and 3 for the extreme Case IV. Because M, the Coriolis vector
%               and gravity are all LINEAR in the mass parameters, this scales
%               all three by the same factor and nothing else -- which is why
%               the perturbation can be applied without re-deriving anything.
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

%% --------------------------------------------- robust CLF-QP (Section 4.1)
% Bounds (4.10) on the induced uncertainty (4.4).  These are bounds on the
% I/O-LINEARIZED uncertainty, not on the mass error: Delta1 has units of
% ydd (rad/s^2) and Delta2 is dimensionless.
%
% HOW TO CHOOSE THEM.  Do not guess. ch4_uncertainty measures the actual
% Delta1, Delta2 along a gait for a given mass_scale, and ch4_delta_bounds
% turns that measurement into these two numbers. The defaults below are the
% measured values along the reference gait for Cases I-III (mass scale 1, 1.5,
% 0.7) with a 1.2 safety factor; run ch4_delta_bounds for any other operating
% point.
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
%     here: 226 on the gait, ~914 with the state jittered by 0.05). Since the
%     commanded ||mu|| scales as delta1_max/(1 - delta2_max), that difference is
%     the whole of the "unnecessarily aggressive" limitation Section 4.1.4
%     closes on -- in a form you can put a number to.
p.rclf = struct();
p.rclf.delta1_max = 250;         % ||Delta1|| bound          [rad/s^2]
p.rclf.delta2_max = 0.45;        % ||Delta2|| bound          [-]

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

%% ------------------------------------------ L1 adaptive control (Sect. 4.2)
p.l1 = struct();

% ADAPTATION GAIN Gamma in (4.26).  The bound (4.34)-(4.35) on the estimation
% error shrinks like 1/||Gamma||, so "sufficiently large" is the whole design
% rule -- the estimate is allowed to be as fast as the integrator can follow.
% Decoupling estimation speed from control smoothness is exactly what the
% low-pass filter below is for, so this can be large without putting
% high-frequency content into the torque.
p.l1.Gamma = 1e4;

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
% reports that overshoot rather than hiding it.
p.l1.u_max = 65;

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
