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
% Generalized coordinates q = [px pz qt q1 q2 q3 q4]'
%   px,pz  torso base position   (pz is DOWN-positive; world height = -pz)
%   qt     torso pitch
%   q1,q2  stance hip, knee      (actuated)
%   q3,q4  swing  hip, knee      (actuated)
p.nq   = 7;                     % generalized coordinates
p.nu   = 4;                     % actuators
p.nx   = 2*p.nq;                % full state dimension (14)
p.iact = (4:7)';                % actuated coordinate indices
p.ny   = numel(p.iact);         % number of outputs (4)
p.g0   = 9.8062;                % gravity used in the symbolic derivation
p.mass = 32;                    % total mass [kg] (10 torso + 4 x 5 leg links + 2)

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
p.basis   = 'bezier';
p.bez_deg = 5;                  % Bezier degree M
p.n_ctrl  = p.bez_deg + 1;      % alpha columns per output (M+1)
p.bsp_deg = 3;                  % degree used only when basis = 'bspline'

%% ------------------------------------------------------------- controller
% Selects the swing-phase feedback law used EVERYWHERE (simulation, cost,
% constraints, force recovery) because they all route through ch3_control.
%   'iolin_pd'  stage 5: I/O linearization + PD on the linearized system
%   'clfqp'     stage 7: unconstrained CLF-QP
%   'clfqp_con' stage 8: CLF-QP + torque box (+ friction / GRF if enabled)
%   'ff'        pure feedforward u_ff, no output feedback (diagnostic only)
p.controller = 'iolin_pd';

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
p.clf_construction = 'care';
p.Q_clf = eye(2*p.ny);          % Q in the (C)ARE / Lyapunov equation
p.clf_slack_penalty = 1e6;      % "p" multiplying delta^2 in stage 8

%% -------------------------------------------------- gait / speed targets
% The seed rollout walks at ~0.12 m/s with a 0.084 m step; the known reference
% gait for this robot is ~0.355 m/s. Asking for much more than that in a cold
% solve makes NEC1 the dominant residual and the solve fights it from the
% first iteration. Use ch3_continuation to march v_des beyond this.
p.v_des        = 0.35;          % NEC1 average walking rate [m/s]
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
p.qt_range     = [-1.0 1.0];    % [qt_min qt_max]

% Gate on the NEC1 speed equality. Turn it OFF for the first solve: getting a
% periodic gait at all is the hard part, and pinning the speed simultaneously
% is what makes a cold solve stall (see ch3_col_constraints). Then turn it on
% at the achieved speed and march with ch3_continuation.
p.enforce_nec1 = true;

%% ---------------------------------- Table 3.1 physical realizability limits
% Each limit is individually gated.  A DISABLED limit is still evaluated and
% reported -- it is just not enforced -- so you can measure before tightening.
% Enable them in this order: GRF first (get Fz > 0), then friction, then
% torque, then impulse.  Enabling friction while Fz still crosses zero makes
% |Fx/Fz| blow up and fmincon chases a divide-by-zero.
p.limits = struct();
p.limits.u_max       = 120;     % |u_i| <= u_max                    [Nm]
p.limits.impulse_max = 15;      % ||impact impulse||_2 <= impulse_max [Ns]
p.limits.mu_s        = 0.4;     % |Fx| <= mu_s * Fz                 [-]
p.limits.Fz_min      = 50;      % Fz >= Fz_min (foot stays loaded)  [N]
p.limits.clearance   = 0.05;    % swing-foot height at mid-step     [m]

p.limits.enable = struct('torque',  false, ...
                         'impulse', false, ...
                         'friction',false, ...
                         'grf',     false, ...
                         'clearance', true);

%% ------------------------------------------------------- direct collocation
% Node count is the main cost/accuracy dial. Every fmincon gradient costs
% roughly (n_vars x 2 x 3(N-1)) evaluations of the closed-loop dynamics, and
% n_vars itself grows as 14N -- so the work scales like N^2. 15 nodes keeps a
% solve to minutes while leaving Hermite-Simpson (3rd order) plenty accurate
% for a step this smooth.
p.N_nodes  = 15;                % Hermite-Simpson nodes per step
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
