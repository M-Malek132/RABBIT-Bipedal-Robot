function p = ch3_params(varargin)
%CH3_PARAMS  Single source of truth for the Chapter-3 pipeline.
%
%   p = ch3_params()                 default parameter struct
%   p = ch3_params('name',value,...) defaults with overrides
%
% The Chapter-3 pipeline is:
%
%   (1) hybrid model            xdot = f(x)+g(x)u,  x+ = Delta(x-)
%   (2) virtual constraints     y = y0(q) - yd(s(q),alpha)
%   (3) offline optimization    min (1/L_step) int ||u||^2   ->  alpha
%   (4) I/O linearization       ydd = Lf^2 y + LgLf y u
%   (5) PD baseline             mu = -Kp/eps^2 y - Kd/eps ydot
%   (6) RES-CLF                 Vdot_eps + (c3/eps) V_eps <= 0
%   (7) CLF-QP                  min mu'mu  s.t. CLF inequality
%   (8) constrained CLF-QP      + torque box / friction / GRF, slack delta
%
% Every stage reads its knobs from here so that a single edit changes the
% whole pipeline consistently.
%
% NOTE ON UNITS AND LIMITS.  The thesis Table 3.1 numbers are ATRIAS numbers
% (63 kg, 50:1 harmonic drives, so its "|u| <= 5 Nm" is MOTOR torque = 250 Nm
% at the joint).  RABBIT is ~30 kg and DIRECT DRIVE -- u here IS joint torque.
% The limits below are therefore RABBIT-scaled placeholders; the intended
% workflow is to solve once with p.limits.enable all false, read the measured
% ranges off ch3_report, and only then tighten.  See ch3_limits_from_report.

p = struct();

%% ------------------------------------------------------------------ model
% Generalized coordinates q = [px y qt q1 q2 q3 q4]'
%   px,y   torso base position   (y is DOWN-positive; world height = -y)
%   qt     torso pitch
%   q1,q2  stance hip, knee      (actuated)
%   q3,q4  swing  hip, knee      (actuated)
p.nq   = 7;                     % generalized coordinates
p.nu   = 4;                     % actuators
p.nx   = 2*p.nq;                % full state dimension (14)
p.iact = (4:7)';                % actuated coordinate indices
p.ny   = numel(p.iact);         % number of outputs (4)
p.g0   = 9.8062;                % gravity used in the symbolic derivation
p.mass = 30;                    % total mass [kg] (10 torso + 4 x 5 leg links).
                                % Matches Mass_Properties in the symbolic
                                % derivation, and M(1,1) -- which for a
                                % translational floating-base coordinate IS the
                                % total mass. Nothing reads this field; it is
                                % here so the number is stated somewhere.

% y0(q) = H*q.  For RABBIT the controlled variables are simply the four
% actuated joint angles (hip and knee of each leg), so H selects them.
p.H = zeros(p.ny, p.nq);
p.H(:, p.iact) = eye(p.ny);

%% --------------------------------------------------- phase variable theta
% theta = qt + q1 + q2/2 = c*q, the ABSOLUTE STANCE LEG ANGLE.
%
% Linear in q on purpose.  The geometric hip->foot angle atan2(dx,dz) is
% EXACTLY this expression for every physical pose (both links have length
% 1/2), but atan2 wraps at +-pi and a wrap mid-solve makes the optimizer
% chase a phantom discontinuity.  The linear form is the same function with
% the branch cut removed, and its gradient is the CONSTANT row c -- which is
% what makes ds/dq constant and the Lie derivatives in ch3_io_lin cheap.
% This is Westervelt's hypothesis HH6, on which the HZD decoupling rests.
p.c_theta = [0 0 1 1 0.5 0 0];

% theta sweeps monotonically from theta_minus (start of step) to theta_plus
% (swing-foot strike).  s = (theta - theta_minus)/(theta_plus - theta_minus).
p.theta_minus = -0.15;
p.theta_plus  =  0.30;

%% ------------------------------------------ virtual constraint parametrization
% yd(s,alpha) is a vector of 4 curves in the phase s.
%   'bezier'  Chapter-3 form: degree p.bez_deg Bernstein polynomials,
%             analytic first and second derivatives (ch3_bezier).
%   'bspline' clamped B-spline via the repo's existing BSpline.m /
%             BSpline_derivative.m, kept as a cross-check.
%
% These two agree EXACTLY when n_ctrl = bez_deg+1 and the B-spline degree is
% bez_deg: a clamped B-spline with degree+1 control points IS the Bezier
% curve of that degree.  ch3_test_basis exploits that as a unit test.
%
% THE DEFAULT IS NOT THAT CASE, AND IT COSTS ORDER -- MEASURED. bsp_deg = 3
% with 6 control points is a piecewise cubic, C^2 at its interior knots and no
% smoother, and ch3_yd central-differences its d2yd/ds2. Hermite-Simpson needs
% more smoothness than that, so the seed-rollout defect converges at
%
%       bsp_deg = 3   order 1.61,  max|defect| = 1.675e-03 at N = 49
%       bsp_deg = 5   order 4.51,  max|defect| = 3.874e-07 at N = 49
%       'bezier'      order 4.51,  max|defect| = 3.874e-07 at N = 49
%
% about 4300x in transcription accuracy. ch3_test_collocation's 4th-order
% assertion is therefore left FAILING ON PURPOSE -- an honest report of the
% tradeoff, not a defect to tune away. Use 'bezier' or bsp_deg = 5 to recover.
%
% Anything deriving coefficients from endpoint identities must ask the basis
% rather than assume Bezier; see ch3_seed. The reference gait in Results/ was
% solved under 'bezier', which is why two such bugs survived there unseen.
p.basis   = 'bspline';
p.bez_deg = 5;                  % Bezier degree M
p.n_ctrl  = p.bez_deg + 1;      % alpha columns per output (M+1)
p.bsp_deg = 3;                  % degree used only when basis = 'bspline'

%% ------------------------------------------------------------- controller
% Selects the swing-phase feedback law used to RUN a gait -- forward
% simulation, animation, and the force/torque recovery that follows one --
% because those all route through ch3_control.
%
% IT DOES NOT AFFECT THE STAGE-3 SOLVE. Nothing in Optimization/ calls
% ch3_control at all: ch3_col_dynamics goes straight to ch3_io_lin and runs on
% the FEEDFORWARD alone. The gait is designed ON the zero dynamics surface --
% node 1 is pinned to y = ydot = 0 and u_ff renders that surface invariant, so
% eta stays zero across the whole step and any PD or CLF term would be
% multiplying zero. Changing this field runs the SAME alpha under a different
% law; it does not produce a different gait. See ch3_col_dynamics for the
% derivation and ch3_main for where the two stages divide.
%   'iolin_pd'  stage 5: I/O linearization + PD on the linearized system
%   'clfqp'     stage 7: unconstrained CLF-QP
%   'clfqp_con' stage 8: CLF-QP + torque box (+ friction / GRF if enabled)
%   'ff'        pure feedforward u_ff, no output feedback (diagnostic only)
p.controller = 'clfqp';

% Stage 5 gains.  mu = -(1/eps^2) Kp y - (1/eps) Kd ydot.
%
% Kp = 100, Kd = 20 places the nominal transverse poles at a critically damped
% 10 rad/s; eps then scales that rate by 1/eps, so the effective bandwidth is
% 10/eps rad/s and the effective proportional gain is 100/eps^2.
%
% eps must be SMALL enough that the outputs re-converge between impacts -- each
% impact expands eta and the controller has exactly one step to beat that
% expansion -- and LARGE enough that the torque stays inside the actuators.
% That tension is the entire motivation for stages 6-8: at eps = 0.1 the
% effective gain is 1e4 and a 0.3 rad output error alone commands ~6.7 kNm,
% which RABBIT (direct drive, tens of Nm) cannot produce, and PD has no way to
% back off. The constrained CLF-QP does. Default here is a value RABBIT can
% actually track; lower it to see stage 8 start earning its keep.
p.eps = 0.50;
p.Kp  = eye(p.ny) * 100;
p.Kd  = eye(p.ny) *  20;

%% ---------------------------------------------------------------- RES-CLF
% V_eps(eta) = eta' P_eps eta with eta = [y; ydot] (8x1) and
% P_eps = I_eps P I_eps,  I_eps = blkdiag((1/eps) I, I).
%
%   'care'  F'P + PF - P G G' P + Q = 0   (control ARE)
%           The standard RES-CLF construction.  Guarantees that a feasible mu
%           exists pointwise, which is what makes the stage-7 QP always
%           solvable in the unconstrained case.
%   'lyap'  A'P + PA = -Q with A = [0 I; -Kp -Kd]
%           The literal form written in Chapter 3.  Also a valid CLF (the
%           stage-5 mu always satisfies the decrease condition) but tied to
%           the PD gains.
p.clf_construction = 'lyap';
p.Q_clf = eye(2*p.ny);          % Q in the (C)ARE / Lyapunov equation
p.clf_slack_penalty = 1e6;      % "p" multiplying delta^2 in stage 8

%% -------------------------------------------------- gait / speed targets
% The seed rollout walks at ~0.12 m/s with a 0.084 m step; the known reference
% gait for this robot is ~0.355 m/s. Asking for much more than that in a cold
% solve makes NEC1 the dominant residual and the solve fights it from the
% first iteration. Use ch3_continuation to march v_des beyond this.
p.v_des        = 1.2;          % NEC1 average walking rate [m/s]
p.step_len_min = 0.15;          % floor on step length, kills "step in place"

% TORSO PITCH BOX [rad].  Unlike the rails in ch3_col_bounds this one is a
% design requirement, not a guard against bad arithmetic.
%
% qt is UNACTUATED (p.iact = 4:7), so pitching the torso costs the objective
% nothing directly, and leaning shifts gravity load off the actuators. Nothing
% else in the problem mentions qt either: it enters only through
% theta = qt + q1 + q2/2, which pins one LINEAR COMBINATION of the pose, not
% the pose -- a large negative qt is freely offset by a large positive q2/2.
% Left at +-1 rad the optimizer duly spends that freedom: the first reference
% gait solved here walks with the torso pitched 46 deg BACKWARD for the whole
% step (qt in [-0.870, -0.806], a 3.6 deg variation), knees at 73-108 deg of
% flexion, at 1.17 m/s. It is periodic, mesh-verified and stable -- and it
% looks nothing like walking.
%
% Tighten this to force an upright torso, and tighten it INCREMENTALLY with
% warm starts, exactly like the Table 3.1 limits: clamping straight to upright
% from a leaning solution moves every node at once. Expect J to rise; an
% upright gait genuinely costs more torque, which is the honest answer rather
% than a regression.
%
% SIGN: qt < 0 leans the torso BACKWARD, qt > 0 leans it FORWARD.  The box is
% what picks the lean, not a preference term in the cost -- with the load
% argument above pushing qt monotonically negative, the optimizer parks on
% whichever bound is lower and stays there. So the LOWER bound is the design
% knob: a solve given [-0.30 0.45] settles at -0.300, dead on the bound.
% Asking for a forward lean therefore means a box that is entirely positive,
% e.g. [0.08 0.25] for the ~5 deg lean of human walking; a box that merely
% straddles zero will still come back leaning backward.
%
% The default stays wide because a COLD solve needs it wide -- same reasoning
% as enforce_nec1 below. Reach a leaning target with ch3_posture_march.
p.qt_range     = [-1.0 1.0];    % [qt_min qt_max]  (see SIGN above)

% Gate on the NEC1 speed equality. ON by default here, but turn it OFF for a
% COLD solve: getting a periodic gait at all is the hard part, and pinning the
% speed simultaneously is what makes a cold solve stall (see
% ch3_col_constraints). Then turn it back on at the achieved speed and march
% with ch3_continuation.
p.enforce_nec1 = true;

%% ---------------------------------- Table 3.1 physical realizability limits
% Each limit is individually gated, and EVERY GATE IS ON below -- full fidelity
% rather than staging up to it. That is right for re-solving a gait that already
% exists and wrong for finding a new one.
%
% A disabled limit is still evaluated and reported, just not enforced, so the
% staged workflow stays available and is the fallback when a cold solve stalls:
% gates off, measure off ch3_report, then re-enable GRF first (get Fz > 0), then
% friction, then torque, then impulse. The order matters -- enabling friction
% while Fz still crosses zero makes |Fx/Fz| blow up and fmincon chases a
% divide-by-zero.
p.limits = struct();
p.limits.u_max       = 120;     % |u_i| <= u_max                    [Nm]
p.limits.impulse_max = 15;      % ||impact impulse||_2 <= impulse_max [Ns]
p.limits.mu_s        = 0.4;     % |Fx| <= mu_s * Fz     NIC2        [-]
p.limits.Fz_min      = 50;      % Fz >= Fz_min          NIC1        [N]
p.limits.clearance   = 0.05;    % swing-foot height at mid-step     [m]

%% -------------------------------- Section 6.3.4 constraint set (NIC / NEC)
% The book's numbered constraints, in its own order. ch3_col_constraints maps
% each to a row and says which hypothesis it discharges. Everything here is
% gated for the same measure-then-tighten reason the Table 3.1 limits are.
%
% NIC1 = grf (Fz_min above), NIC2 = friction (mu_s above) -- both already had
% homes, so they keep their existing knobs rather than gaining duplicates.

% NIC3 -- "swing leg end height to ensure S intersects Z (ONLY) at the end of
% the step". Two separate requirements hide in that sentence:
%
%   (a) the swing foot is STRICTLY above the ground during the step, so the
%       guard cannot fire early. A margin, not >= 0: at exactly zero the guard
%       is grazed and the step can terminate anywhere.
%   (b) the crossing at the end is TRANSVERSAL -- the foot is genuinely moving
%       down when it strikes. Without this, "height = 0" at node N is satisfied
%       by a foot that touches and rises, which is a tangency, not a strike,
%       and the Poincare map is not even locally well defined there.
%
% The margin cannot be applied at the two ENDS of the step: the foot is at
% height zero at node 1 by construction (it just impacted) and at node N by the
% guard equality, so a margin there is infeasible by definition. See
% ch3_col_constraints for exactly which nodes and midpoints it covers.
% MEASURED: the reference gait's worst interior clearance is 0.0029 m and it
% strikes at -3.88 m/s. Note how small that clearance is next to the 0.19 m it
% reaches at mid-step -- the binding points are right next to the endpoints,
% which is why the margin here is millimetres and p.limits.clearance, the
% mid-step style constraint, is centimetres. They are not the same knob.
p.limits.sw_clear_min   = 1e-3; % strict swing-foot clearance, interior  [m]
p.limits.sw_strike_rate = 0.05; % foot must be descending at strike    [m/s]

% NEC2 -- vertical component of the post-impact swing-leg velocity is positive.
% After relabeling the new swing leg is the OLD STANCE leg, so this is the
% condition that the trailing foot actually lifts off.
%
% MEASURED: the reference gait lifts off at 0.0495 m/s -- positive, but only
% just, and two orders below its own 3.88 m/s strike rate. The default is set
% an order of magnitude below that so it reads as what it is, a strict-
% positivity margin, rather than as a hidden style requirement that the one
% gait in the repo happens to fail by 1%. Raise it if you want a gait that
% picks its trailing foot up decisively.
p.limits.liftoff_rate   = 0.01; % d/dt(new swing-foot height) at t=0+  [m/s]

% NEC3 -- validity of the impact. The rigid plastic impact model is only
% meaningful if the impulse it predicts is one the ground can actually apply:
% compressive (the floor cannot pull the foot down) and inside the friction
% cone (no slip during the collision). Same mu_s as NIC2.
%
% MEASURED, AND IT FAILS: the reference gait's impulse is mostly HORIZONTAL
% where a walking impact should be mostly vertical, giving |Ix|/Iz = 1.35
% against mu_s = 0.4. The impact map imposes J_sw dq+ = 0, "the foot sticks";
% at that ratio it would not stick, it would skid. The continuous-phase
% friction (NIC2) is a comfortable 0.194, so this is invisible unless the
% impulse is checked separately -- exactly why the book lists NEC3 apart from
% NIC2. Enabling this gate changes the gait; it is not a cosmetic limit.
%
% READ THAT NUMBER ON A FINE MESH.  At N = 21 the same gait measures 2.44 with
% Iz = 4.86 Ns; at N = 41 it measures 1.35 with Iz = 8.34 Ns. The N = 21 gait
% PASSES ch3_col_verify (6.75e-04 against tol 1e-03), so this is not a spurious
% discrete solution -- it is a real trajectory whose IMPULSE is nonetheless 80%
% wrong. ch3_col_verify bounds max|X_node - X_true| over the step; the impulse
% is Lambda(q_N) v_foot(x_N), evaluated at ONE endpoint, and Iz is small enough
% that a 7e-04 state error swamps it. Every Table 3.1 quantity is an extremum
% over the whole step and survives a coarse mesh; NEC3 is the one constraint in
% the set that does not. Refine before believing it, and before laddering to
% it -- a ladder calibrated on the coarse number never becomes active.
p.limits.impulse_z_min  = 0;    % Iz >= this (compressive)              [Ns]

% The impulse cone gets its OWN coefficient, empty meaning "use mu_s". Not
% because the physics differ -- it is the same floor -- but because a gait that
% starts 6x outside the cone has to be MARCHED in, and marching p.limits.mu_s
% would silently drag the continuous-phase constraint (NIC2) along with it.
% ch3_impact_march steps this down from whatever the gait measures to mu_s.
p.limits.mu_s_impact    = [];   % |Ix| <= this * Iz; [] -> p.limits.mu_s   [-]

% HH6 -- theta strictly monotonic. Also what makes s a legitimate clock and
% carries the forward-progression content of HGW6.
% MEASURED: the reference gait's slowest thetadot over the step is 1.25 rad/s.
p.limits.thetadot_min   = 0.10; % thetadot >= this over the step      [rad/s]

% HH2 -- decoupling matrix invertible on Z, as a floor on its SMALLEST SINGULAR
% VALUE rather than on rcond or a determinant. det is scale-blind (it can stay
% large while one direction collapses) and rcond is a norm ratio; sigma_min is
% the actual distance to a singular matrix, which is the quantity the
% hypothesis is about.
% MEASURED: the reference gait's worst sigma_min over the step is 0.114.
p.limits.dec_min        = 1e-3; % sigma_min(LgLf y) >= this

% HIP-HEIGHT BAND [m].  Like qt_range, a design requirement rather than a
% hardware limit -- and gated for the same reason the others are.
%
% The hip is the torso base, so its world height is -pz, and NOTHING in the
% problem otherwise mentions pz: it is not in theta, not actuated, and a lower
% hip buys knee torque cheaply (a bent knee has a shorter gravity moment arm).
% Left free the optimizer duly crouches -- the upright-torso gait solved here
% walks with the hip at 0.855 m and the stance knee at 54-64 deg of flexion,
% i.e. a permanent half-squat that is periodic, stable, and does not look like
% walking.
%
% The band's CENTRE is the standing height, its WIDTH the tolerated bob.  Legs
% are two 0.5 m links, so 1.0 m is full extension and hip_h is best read as a
% fraction of that: 0.90 leaves 10% of knee bend, which is about where human
% walking sits. Do not push it past ~0.95 -- the leg straightens, the stance
% knee approaches its singularity, and the KKT solve in ch3_col_dynamics loses
% conditioning.
%
% Enable it the way the Table 3.1 limits are enabled: warm-started, one stage
% at a time, marching the centre from whatever the current gait measures.
% ch3_posture_march does this.
p.limits.hip_h       = 0.90;    % centre of the hip-height band      [m]
p.limits.hip_h_tol   = 0.02;    % half-width (so 4 cm of bob total)  [m]

p.limits.enable = struct('torque',  true, ...
                         'impulse', true, ...
                         'friction',true, ...   % NIC2
                         'grf',     true, ...   % NIC1
                         'clearance', true, ...
                         'height',  true, ...
                         'swing_clear', true, ...  % NIC3
                         'liftoff',     true, ...  % NEC2
                         'impact',      true, ...  % NEC3
                         'hzd',         true, ...  % NEC4 + NEC5
                         'phase_mono',  true, ...  % HH6
                         'decoupling',  true);     % HH2

%% -------------------------------------------- hybrid zero dynamics (NEC4/5)
% Grid for the quadratures in ch3_zero_dynamics. Two sizes on purpose:
%
%   hzd_grid        used for REPORTING. Accuracy matters, cost does not.
%   hzd_grid_solve  used INSIDE the constraints when enable.hzd is on.
%
% The quadratures converge at O(h^2): delta_zero^2 on the reference gait reads
% 0.76156 / 0.76028 / 0.75996 / 0.759875 at 41 / 81 / 161 / 321 points. So 41 is
% already three digits -- ample for a constraint whose job is to keep delta^2
% away from 1 -- while the report can afford 161.
%
% WHAT enable.hzd COSTS -- MEASURED. fmincon finite-differences the
% constraints over every decision variable, and each evaluation of
% ch3_zero_dynamics costs n_grid points of about four 7x7 solves each. On the
% reference gait (319 variables) one constraint evaluation goes from 0.003 s to
% 0.241 s when this gate is on, so one gradient goes from about 1 s to about
% 77 s -- a 80x tax on every fmincon iteration.
%
% It is ON in the shipped defaults anyway, along with every other gate: this
% configuration solves at full fidelity rather than staging. Budget for that
% 80x when timing a solve, and if a COLD solve stalls, the first thing to try
% is the staged workflow the Table 3.1 limits describe -- turn this off, read
% the measured delta_zero^2 and zeta*_2 off ch3_report, and re-enable once the
% gait exists. hzd_grid_solve is the other dial: it is what this gate actually
% evaluates, and 41 points is already three digits of delta_zero^2.
%
% AND NOTE WHAT ALREADY ENFORCES IT. This transcription imposes periodicity
% through Delta DIRECTLY, so a converged solve is already at the fixed point --
% NEC4 and NEC5 are then a CHECK on the fixed point it found, not the mechanism
% that finds one. The book needs them as constraints because its optimization
% parametrizes alpha alone and never propagates a full state; here they earn
% their keep as a stability certificate, and as constraints only when you want
% to steer a solve away from a marginally stable orbit.
p.hzd_grid       = 161;
p.hzd_grid_solve = 41;

% Margins turning the book's two STRICT inequalities into closed constraints.
% An optimizer cannot satisfy a strict inequality; without a margin it parks on
% delta^2 = 1 exactly, which is neutral stability reported as success.
p.hzd_tol = struct('delta_margin', 1e-3, ...   % delta^2 <= 1 - this
                   'delta_min',    1e-8, ...   % delta^2 >= this
                   'zeta_margin',  1e-6);      % zeta* clears V_max/delta^2 by this

%% ------------------------------------------------------- direct collocation
% Node count is the main cost/accuracy dial. Every fmincon gradient costs
% roughly (n_vars x 2 x 3(N-1)) evaluations of the closed-loop dynamics, and
% n_vars itself grows as 14N -- so the work scales like N^2. 15 nodes keeps a
% solve to minutes while leaving Hermite-Simpson (3rd order) plenty accurate
% for a step this smooth.
p.N_nodes  = 41;                % Hermite-Simpson nodes per step
p.T_min    = 0.20;              % step duration bounds [s]
p.T_max    = 1.50;
p.dq_max   = 20;                % box on joint velocities, keeps fmincon sane

p.max_iter      = 300;
p.max_fun_evals = 3e5;

% Tolerance for ch3_col_verify: how far the node states may sit from a tight
% ODE rollout before the solution is rejected as a spurious discrete solution.
% Small defects do NOT imply a real trajectory -- see ch3_col_verify.
p.verify_tol = 1e-3;

% Periodic checkpoint of the decision vector during a solve. Empty disables.
% Set it to a path and a crashed or interrupted solve can be resumed from the
% last checkpoint instead of restarted.
p.checkpoint_file  = '';
p.checkpoint_every = 10;

% Mid-step swing-knee flexion baked into the seed profile [rad]. Without it
% the seed does not walk -- the swing foot scuffs down onto the stance foot
% and the step length is ~0. See ch3_seed for the sign derivation.
% Measured sweep (see ch3_test_collocation): lift 0 -> L_step 0.002 m, phase
% only 52% complete; lift 1.2 -> L_step 0.084 m, phase 90% complete. Returns
% diminish past ~1.2.
p.seed_knee_lift = 1.2;

%% ------------------------------------------------------------ integration
p.ode_reltol = 1e-8;
p.ode_abstol = 1e-9;
p.ode_maxstep = 5e-3;
p.guard_min_time = 0.05;        % ignore guard crossings before this [s], so
                                % the step does not terminate on the swing
                                % foot still being at z=0 at t=0

% CONTROL SAMPLE PERIOD [s].  0 = evaluate the feedback continuously, i.e. let
% ode45 call the controller at every stage of every RK step.
%
% Set this NONZERO for the constrained CLF-QP.  That controller is only
% piecewise smooth in x: the QP's active set changes as torque bounds engage
% and disengage, and u(x) kinks at every switch. An adaptive explicit solver
% treats each kink as a failed error test and shrinks h to resolve it, which is
% not a slow simulation but a stalled one -- measured here at 51908 RHS
% evaluations to advance 0.0015 s of a 0.3009 s step, with h collapsed to ~3e-8.
%
% Holding u over a fixed period removes the kinks from the integrand entirely
% (within a period the RHS is the smooth f + g u with u a constant) and is also
% what the hardware actually does -- the chapter's own framing is a QP solved
% "well above 1 kHz", not a continuous-time control law. So this is the more
% faithful model as well as the tractable one.
p.control_dt = 0;

%% -------------------------------------------------------------- overrides
for k = 1:2:numel(varargin)
    field = varargin{k};
    if ~isfield(p, field)
        error('ch3_params:unknownField', 'Unknown parameter "%s".', field);
    end
    p.(field) = varargin{k+1};
end

% Keep derived quantities consistent if the caller changed bez_deg.
p.n_ctrl = p.bez_deg + 1;
p.ny     = numel(p.iact);

end
