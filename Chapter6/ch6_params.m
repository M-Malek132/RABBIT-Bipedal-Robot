function p = ch6_params(varargin)
%CH6_PARAMS  Single source of truth for the Chapter-6 stepping-stone pipeline.
%
%   p = ch6_params()                 defaults
%   p = ch6_params('name',value,...) defaults with overrides
%
% Chapter 6 is Chapter 3's controller with extra rows in the QP. It therefore
% STARTS from ch3_params and adds only what is new, so that a change to the
% robot, the phase variable, the RES-CLF or the Table 3.1 limits propagates
% here with no second edit. Every ch3 field is still reachable by name.
%
% What Chapter 6 adds:
%
%   p.cbf        the barrier construction    (Sections 6.1.1, 6.2)
%   p.stones     the discrete-foothold geometry, R1/R2   (6.1.3, 6.9)
%   p.obstacle   the overhead-obstacle geometry          (6.1.2, 6.2/6.5)
%   p.width      the 3D step-width geometry, R3/R4       (6.3.2, 6.19/6.20)
%   p.lib        the gait library                        (6.4, 6.22-6.24)
%   p.mc         the Table 6.1 Monte-Carlo settings
%
% ------------------------------------------------------------- CONTROLLER NAMES
% p.controller extends ch3's list with one new value:
%
%   'ff' | 'iolin_pd' | 'clfqp' | 'clfqp_con'   exactly as in Chapter 3
%   'cbf_clf_qp'                                Chapter 6: (6.12) / (6.26)
%
% The Chapter-3 names are kept working because Table 6.1 needs them as the
% BASELINES it is measuring against -- controller I ("gait library", no CBF)
% is 'clfqp_con' with p.lib.enable true, and controller II ("CBF, one nominal
% gait") is 'cbf_clf_qp' with p.lib.enable false.
%
% See also CH3_PARAMS, CH6_CONTROL, CH6_BARRIER, CH6_MAIN.

p = ch3_params();

%% ------------------------------------------------------------------ controller
% Chapter 6's QP is Chapter 3's stage-8 QP with barrier rows appended, so it
% inherits the sampled-data requirement for exactly the same reason: the active
% set changes as barriers and torque bounds engage, u(x) kinks at every switch,
% and an adaptive explicit solver stalls on the kinks. See ch3_params on
% control_dt -- measured there at 51908 RHS evaluations for 0.5% of a step.
%
% Chapter 6 makes it WORSE than Chapter 3 did, because the barrier rows switch
% more often than the torque box does, so control_dt is nonzero by DEFAULT here
% rather than being something the caller has to remember. 1 kHz is the rate the
% chapter itself quotes ("solved in under 1 ms").
p.controller = 'cbf_clf_qp';
p.control_dt = 1e-3;

%% ------------------------------------------------ Section 6.1.1: the barrier
% THE ONE IDEA OF 6.1.1, AND WHAT IT REALLY IS.
%
% Everything Chapter 6 constrains is a POSITION: a foot placement, a head
% height, a step width. Position constraints g_b(q) >= 0 have relative degree
% TWO -- u does not appear in gdot -- so the Section 5.1 reciprocal CBF cannot
% be written down for them (its row is 0*u >= ..., see ch5_barrier). Section
% 6.1.1's fix is to build a relative-degree-ONE function out of g_b,
%
%       h_CBF(q,qdot) = gamma_b g_b(q) + gdot_b(q,qdot) >= 0            (6.1)
%
% and put the reciprocal barrier B = 1/h_CBF on THAT.
%
% h_CBF >= 0 keeps g_b >= 0 by the argument under (6.1): at g_b = 0 with the
% constraint about to be violated we would need gdot_b < 0, which makes
% h_CBF < 0 first. So (6.1) is a strictly stronger condition that is also
% enforceable, which is the trade.
%
% IT IS THE rb = 2 EXPONENTIAL CBF, WRITTEN OUT BY HAND. Take (6.1) and ask
% for h_CBF to decay no faster than exponentially, hdot_CBF >= -gamma2 h_CBF:
%
%       gamma_b gdot + gddot >= -gamma2 (gamma_b g + gdot)
%   <=> gddot >= -(gamma_b gamma2) g - (gamma_b + gamma2) gdot
%   <=> gddot >= -Kb eta_b,   eta_b = [g; gdot],  Kb = [g1 g2, g1+g2]
%
% which is Definition 5.1 with poles (gamma_b, gamma2) -- exactly what
% ch5_ecbf_gain builds. ch6_test_barrier asserts the two rows are identical.
% So Section 6.1 is not a different construction from Section 5.2; it is the
% rb = 2 case of it, and the ONLY substantive difference between 6.1 and 6.2
% is the class-K function on the right-hand side:
%
%   'reciprocal'   hdot_CBF >= -gamma * h_CBF^3    Section 6.1, from B = 1/h
%   'exponential'  hdot_CBF >= -gamma2 * h_CBF     Section 6.2/6.4, (6.12)
%
% The cubic is what B = 1/h gives after substituting into Bdot <= gamma/B; see
% ch5_barrier. It is MUCH weaker than the linear one when h is small (h^3 << h)
% and much stronger when h is large. Near the boundary -- which is the only
% place a barrier is doing anything -- the cubic therefore permits a faster
% approach, which is why the exponential form is the one Section 6.4 reaches
% for when it also has to respect a friction cone and a torque box.
% =========================================================== CHOOSING THE POLES
% This is the decision that governs everything in Chapter 6 and the chapter does
% not discuss it. The two requirements on gamma_b pull in opposite directions.
%
%   FROM BELOW.  At touchdown the barrier is nearly tight BY CONSTRUCTION -- the
%   foot is landing on the edge of a small stone -- while it is still moving at
%   metres per second. h_CBF = gamma_b g + gdot >= 0 therefore needs
%
%           gamma_b >= |gdot| / g
%
%   and on the reference gait (3.88 m/s strike, +-2.5 cm stones) that ratio
%   reaches about 170 rad/s for g_ST2.
%
%   FROM ABOVE.  The row demands gddot >= -Kb eta_b, whose dominant term near
%   h_CBF = 0 is -(gamma_b + gamma) gdot ~ gamma_b |gdot|. At gamma_b = 170 and
%   |gdot| ~ 2 m/s that is 340 m/s^2 of swing-foot acceleration -- which no
%   torque box on this robot can supply, so the QP goes infeasible instead.
%
% THEY ARE INCOMPATIBLE FOR THIS GAIT. What the tuning below buys is that the
% infeasible samples end up confined to the last few milliseconds of the step,
% AFTER the placement is determined: the foot arrives on the stone, and the
% certificate lapses at the terminal instant. ch6_report counts those samples
% and ch6_step_range separates "landed" from "landed with the guarantee intact"
% for exactly this reason.
%
% MEASURED, one step from the reference fixed point (which lands at 0.3537 m),
% sweeping the desired foothold over 0.15:0.025:0.55 m, torque box 300 Nm,
% friction and GRF rows on, PD reference. "Landed" = the foot came down inside
% the window:
%
%     poles    stone +-2.5cm     stone +-5cm       stone +-7.5cm
%     20       0.150 - 0.200     0.150 - 0.250     0.150 - 0.250
%     40       0.150 - 0.325     0.150 - 0.350     0.150 - 0.375   <- widest
%     80       ragged            ragged            ragged
%
% Below 20 the barrier fights the nominal gait everywhere; above 40 it does
% nothing until it demands the impossible, and the ranges break into islands.
% Bigger stones help monotonically, which is the geometry talking: g at
% touchdown grows with the distance past l_min, and the required pole is
% |gdot|/g.
p.cbf = struct();
p.cbf.problem = 'stones';       % 'stones' | 'obstacle' | 'none'
p.cbf.form    = 'exponential';  % 'reciprocal' (6.1) | 'exponential' (6.2/6.4)
p.cbf.gamma_b = 40;             % gamma_b in (6.1); the FIRST ECBF pole
p.cbf.gamma   = 40;             % gamma2 (exponential) or gamma (reciprocal)

% ------------------------------------------------- what the QP's cost pulls to
% (6.26) is written as min ||mu||^2 -- the MINIMUM-NORM control that certifies
% the CLF rate. That is 'ff', and it is not the default.
%
% MEASURED: 10 steps from the reference fixed point, foothold window so wide it
% can never bind, and the BARRIER ROWS REMOVED ENTIRELY:
%
%     PD                                10/10 steps, l_s = 0.3537 every step
%     min-norm CLF-QP (= 'ff', nb = 0)   4/10, l_s drifts 0.358 -> 0.371 ->
%                                        0.418 -> falls
%
% The Chapter-6 QP with zero barriers reproduces the second line exactly, and
% adding the barriers back changes nothing. So the instability is in the COST,
% not the barrier -- which is worth stating plainly, because every symptom of it
% (drifting step length, saturated torque, infeasible samples) looks like a
% barrier-tuning problem and no amount of pole tuning touches it. Chapter 3 had
% already measured the same thing about its stage-7 law: "minimum-norm is not
% the same as well-behaved".
%
% 'pd' makes the program a SAFETY FILTER on the Chapter-3 PD law: identical to
% PD when no barrier is active, and the smallest departure from PD that
% satisfies one when it is. Only the objective changes -- the CLF row, the
% barrier rows, the contact rows and the box are untouched, so nothing (6.26)
% guarantees is weakened.
p.cbf.reference = 'pd';         % 'pd' (safety filter) | 'ff' (literal (6.26))

% ------------------------------------------------- the CLF slack penalty, p1
% Chapter 3 defaults this to 1e6. Lowered here mostly for readability of delta:
% MEASURED over 10 steps with the PD reference, p1 in {1e0, 1e2, 1e4} makes NO
% difference at all -- 10/10 steps at l_s = 0.3537 and 109.6 Nm in every case --
% because with mu_ref = mu_PD the QP already sits where it wants to be and delta
% is whatever the PD law happens to leave. 1e4 keeps delta on a scale where it
% still reads as "how much convergence rate was given up", which is the quantity
% Remark 3.2 is about.
%
% It DOES matter with p.cbf.reference = 'ff', where the CLF row and the cost are
% the same objective and p1 sets the trade -- but that path destabilises this
% gait for reasons p1 cannot fix (see p.cbf.reference).
p.clf_slack_penalty = 1e4;

% Barrier rows are HARD -- they never get a slack variable. Only the CLF row
% does. When the two conflict, safety wins and tracking waits; reversing that
% gives a controller that is occasionally unsafe and always convergent, which
% is the wrong trade and not what (6.12)/(6.26) say. See ch5_ctrl_ecbf_clf_qp.
%
% The consequence is that the QP CAN be infeasible, and this is where it
% happens: barrier rows plus a torque box. ch6_ctrl_cbf_clf_qp reports it
% (qp.feasible) instead of silently returning something; ch6_simulate counts
% infeasible samples per step and ch6_report prints them.
p.cbf.min_LgLf = 1e-10;         % |LgLf g| below this = row has no grip on u

%% ------------------------------------------- Section 6.1.3: discrete footholds
% The two circles of Fig. 6.3, in swing-foot coordinates r = (l_f, h_f)
% measured FROM THE STANCE FOOT:
%
%   O1 at (-R1, 0),        radius R1 + l_max,   foot INSIDE     -> l_s <= l_max
%   O2 at (l_min/2, -R2),  radius sqrt(R2^2 + (l_min/2)^2),
%                          foot OUTSIDE                         -> l_s >= l_min
%
% At h_f = 0 the first gives |R1 + l_s| <= R1 + l_max, i.e. l_s <= l_max; the
% second gives (l_s - l_min/2)^2 >= (l_min/2)^2, i.e. l_s <= 0 or l_s >= l_min,
% and the foot is in front so l_s >= l_min. That is (6.8), and ch6_test_barrier
% checks it by sampling rather than by trusting the algebra.
%
% R2 IS THE SCUFFING KNOB (Remark 6.2). O2 passes through (0,0) and (l_min,0)
% and bulges UPWARD in between, so "outside O2" is simultaneously the lower
% step-length bound and a swing-foot clearance requirement. The bulge height at
% mid-step is
%
%       c = sqrt(R2^2 + (l_min/2)^2) - R2      ->      R2 = (a^2 - c^2)/(2c)
%
% with a = l_min/2. So R2 is not a free tuning constant: pick the clearance you
% want and it follows. That inversion is ch6_R2_from_clearance, and it is why
% the default below is a CLEARANCE in metres rather than a radius.
%
% NOTE the clearance cannot exceed a = l_min/2 (the formula returns R2 <= 0 at
% c = a, a semicircle). Short stones therefore buy less clearance, which is
% geometrically honest: a 10 cm step has no room for a 6 cm arch.
p.stones = struct();
p.stones.R1          = 0.50;    % [m] radius offset of circle O1
p.stones.R2_mode     = 'clearance';   % 'clearance' -> from R2_clearance
                                      % 'fixed'     -> use R2_fixed
p.stones.R2_clearance = 0.05;   % [m] desired mid-step swing-foot clearance
p.stones.R2_fixed     = 0.13;   % [m] used when R2_mode = 'fixed'

% Nominal foothold window used when no terrain is supplied.
p.stones.l_min = 0.25;          % [m]
p.stones.l_max = 0.60;          % [m]

% Section 6.2: the stones MOVE. l_min(t), l_max(t) are affine or sinusoidal in
% the time since the start of the current step, and they FREEZE the moment the
% foot lands ("these stepping stones move with time, stopping once a foot is
% placed on it"). ch6_stone_level evaluates them together with the first and
% second time derivatives the barrier needs -- the time dependence is explicit,
% so gdot picks up a dg/dt term and gddot picks up d2g/dt2 and a cross term.
% Getting those wrong shows up as a barrier that is satisfied at every sample
% and violated between them; ch6_test_barrier finite-differences them.
p.stones.motion  = 'static';    % 'static' | 'linear' | 'sinusoidal'
p.stones.v_stone = 0.0;         % [m/s] drift rate for 'linear'
p.stones.amp     = 0.0;         % [m]   amplitude for 'sinusoidal'
p.stones.freq    = 0.0;         % [Hz]  frequency for 'sinusoidal'

% THE RESOLVED STONE FOR THE CURRENT STEP. p.stones is the terrain description;
% p.stone is the single foothold the barrier is looking at right now, with R1
% and R2 already resolved to numbers. ch6_simulate rewrites it at every impact
% and ch6_resolve_stone builds it; the default here just makes a bare
% ch6_barrier call work without a simulation around it.
p.stone = ch6_resolve_stone(p.stones, p.stones.l_min, p.stones.l_max);

%% ------------------------------------------- Section 6.1.2: overhead obstacles
% Two constraints on the head (torso top), in coordinates measured from the
% stance foot:
%
%   'ceiling'   g_C = h_r - h_H >= 0                                    (6.2)
%   'circle'    g_O = |(l_H,h_H) - (l_m, h_m + R1o)| - R1o >= 0         (6.5)
%
% The second keeps the head OUTSIDE a disc of radius R1o sitting on top of the
% obstacle at (l_m, h_m) -- the green circle of Fig. 6.1 -- so unlike the
% ceiling it only costs head height where the obstacle actually is.
%
% THE NUMBERS BELOW ARE MEASURED OFF THE REFERENCE GAIT, not chosen. Running it
% with no obstacle, the torso top travels
%
%       h_H in [1.265, 1.285] m       l_H in [-0.659, -0.305] m
%
% both relative to the stance foot. Two things follow, and neither is guessable:
%
%   * A ceiling has to sit just under 1.265 m to be a constraint at all, and not
%     much under it to be a satisfiable one. At h_r = 1.05 -- a plausible-looking
%     round number -- the head starts 22 cm inside the obstacle and the QP is
%     infeasible from the first sample; the run then reports a friction
%     violation, which is true and completely misleading about the cause.
%
%   * The head is BEHIND the stance foot for the whole step. The reference gait
%     walks with the torso pitched back about 50 degrees (the repo README
%     documents this), so 0.75 m of torso puts the top 0.3-0.66 m behind the
%     hip. An obstacle placed at a positive l_m is one the head never reaches,
%     and the 'circle' constraint would sit inactive while appearing to work.
p.obstacle = struct();
p.obstacle.type = 'ceiling';    % 'ceiling' | 'circle'
p.obstacle.h_r  = 1.25;         % [m] ceiling height above the stance foot
p.obstacle.l_m  = -0.45;        % [m] obstacle position, in the head's own path
p.obstacle.h_m  = 1.24;         % [m] obstacle height above the stance foot
p.obstacle.R1o  = 0.20;         % [m] radius of the keep-out disc, R1 in (6.4)

%% -------------------------------------------- Section 6.3.2: step-width circles
% Top view, swing-foot coordinates (l_f, w_f):
%
%   O3 tangent to w = w_max FROM BELOW, foot INSIDE   -> w_f <= w_max
%   O4 tangent to w = w_min FROM ABOVE, foot INSIDE   -> w_f >= w_min
%
% both passing through the initial swing-foot position, with
%
%       R = (l3^2 + dw^2) / (2 dw),   dw = |w_max - w0| or |w0 - w_min|   (6.19)
%
% l3 is the only free parameter: the horizontal offset between the initial foot
% position and the circle centre. It is a CONSERVATISM dial. R grows like l3^2,
% so a large l3 gives a huge circle that is locally almost the half-plane
% w_f <= w_max and barely constrains l_f; a small l3 gives a tight circle that
% couples step width to step length. Section 6.3's Case 3 needs both bounds at
% once, so l3 wants to be large enough that O3 and O4 do not fight over l_f.
%
% ON THE SENSE OF O4. The thesis writes the lower bound as "O4F >= R4", i.e.
% the foot OUTSIDE the disc. That cannot bound w_f from below: outside a finite
% disc is not contained in {w_f >= w_min} in any direction, and at the tangency
% abscissa the outside condition selects w_f <= w_min, which is the wrong side.
% Being INSIDE a disc tangent to the line from above does give w_f >= w_min,
% and it is what "the same principle can be applied" actually means -- O3's
% principle is containment. It is implemented that way, and ch6_test_barrier
% asserts BOTH: that the inside form implies the bound, and that the outside
% form does not. See ch6_bar_width.
p.width = struct();
p.width.l3    = 0.30;           % [m] centre offset for both circles
p.width.w_min = 0.18;           % [m]
p.width.w_max = 0.28;           % [m]
p.width.w0    = 0.233;          % [m] previous step width (w0 in 6.19)

%% ------------------------------------------------ Section 6.4: the gait library
% A = { alpha(L_step) : L_min <= L_step <= L_max }                      (6.24)
%
% built by linearly interpolating the Bezier coefficients of a handful of gaits
% solved at fixed step lengths, (6.22)-(6.23), with linear EXTRAPOLATION past
% the ends. Remark 6.6 blames part of the residual failure rate in Table 6.1 on
% exactly that extrapolation, so ch6_lib_alpha reports when it is extrapolating
% and ch6_report prints how often it happened.
%
% The thesis library is MARLO's, at {0.08 0.24 0.40 0.56 0.72} m. RABBIT is a
% different robot with a different nominal gait (the repo's reference gait
% steps 0.353 m), so the grid below is RABBIT-scaled and the ACHIEVED grid is
% whatever ch6_lib_build converges -- it is stored in the library file and read
% back, never assumed.
% ------------------------------------------------------------ building it
% p.lib.targets IS THE MARCH. Each solve warm-starts from its neighbour, so the
% spacing has to be small enough that a solve begins nearly feasible -- and
% every gait it produces is kept, since a denser library only helps (6.22).
%
% MEASURED, and this is why the spacing is 3 cm rather than the 8 cm a five-gait
% grid would suggest: marching 0.353 -> 0.400 in one step converged to the
% target to 3e-15 and still failed the mesh check at 3.19e-02 (tol 1e-03), and
% marching 0.353 -> 0.340 in one step went infeasible with the duration at its
% upper bound. Chapter 3 documents the same thing about speed continuation.
%
% T_band and N_nodes are the two guards that make a pinned step length behave;
% ch6_lib_solve and ch6_lib_build explain each at the point of use.
p.lib = struct();
p.lib.enable  = false;          % false = one nominal gait (controller II)
p.lib.file    = '';             % '' -> Results/ch6_gait_library.mat
p.lib.targets = 0.26 : 0.03 : 0.44;   % [m] the march, and the library
p.lib.iters   = 250;            % fmincon iterations per library gait
p.lib.T_band  = 0.25;           % T within +-25% of the seed's, see ch6_lib_solve
p.lib.N_nodes = 31;             % remesh the seed before marching

%% ------------------------------------------------------- Table 6.1 Monte Carlo
% The thesis runs 100 problem sets x 10 stones for each of 7 ranges and 3
% controllers -- 21000 simulated walking steps. Each step here is an
% event-terminated ode45 rollout with a QP at 1 kHz, so the full grid is days.
% ch6_table61 therefore takes the counts from here and PRINTS what it ran, so a
% reduced table is never mistaken for the full one.
p.mc = struct();
p.mc.n_trials = 20;             % problem sets per (controller, range) cell
p.mc.n_stones = 10;             % stones per problem set
p.mc.stone_sz = 0.05;           % [m] l_max - l_min, the stone size
p.mc.ranges   = {[0.28 0.34], [0.24 0.34], [0.20 0.34], [0.16 0.34]};
p.mc.seed     = 6;              % rng seed, so the table is reproducible

% THE DEMONSTRATION TERRAIN, and why it is not centred on the nominal gait.
% ch6_step_range measures what one nominal gait plus a CBF can actually deliver
% on this robot: roughly [0.15, 0.33] m against a 0.3537 m nominal -- a wide
% shortening range and almost no lengthening (see p.cbf.reference and
% docs/CH6_STEPPING_STONES.md for why that asymmetry is physical rather than a
% tuning artifact). A demo terrain drawn from [0.28, 0.40], which looks like the
% natural choice, spends most of its stones OUTSIDE that envelope and the run
% fails on the first one -- correctly, but it demonstrates the limit rather than
% the method.
%
% So the demo band sits inside the measured envelope, and the envelope itself is
% what the 'range' study reports. The two are separate on purpose: one shows the
% method working, the other shows where it stops.
p.demo_band = [0.20 0.32];      % [m] desired step lengths for the ch6_main runs

%% -------------------------------------------- (6.25) and input saturation
% (6.26) carries three rows Chapter 3's stage-8 QP only optionally carried, and
% Chapter 6 turns them ON by default, because the whole claim of Section 6.4 is
% that footstep placement is achieved WITHOUT violating them. Measuring them
% while not enforcing them would be measuring a different controller.
%
% The numbers are RABBIT-scaled versions of the thesis's, which are MARLO's:
%
%   thesis            here            why
%   Fv_st >= 150 N    Fz >= 50 N      MARLO is ~63 kg, RABBIT ~32 kg, and the
%                                     reference gait's own minimum is 165 N
%   |Fh/Fv| <= 0.6    same            a friction coefficient is not robot-scaled
%   |u| <= 5 Nm       |u| <= 300 Nm   MARLO's 5 Nm is MOTOR torque through a
%                                     50:1 drive, i.e. 250 Nm at the joint;
%                                     RABBIT here is DIRECT DRIVE, so u IS joint
%                                     torque
%
% ON THE 300 Nm. The reference gait itself needs 110 Nm, so the box is not what
% makes the gait possible -- it is what limits how far the barrier can bend it.
% MEASURED: at 150 Nm the achievable foothold range collapses to a couple of
% centimetres either side of the nominal. So this number is a determinant of
% Section 6.1.4's headline range, and quoting a range without it would be
% quoting nothing.
p.limits.enable.torque   = true;
p.limits.enable.friction = true;
p.limits.enable.grf      = true;
p.limits.u_max  = 300;          % [Nm]  joint torque, direct drive
p.limits.mu_s   = 0.6;          % [-]   kf in (6.25)
p.limits.Fz_min = 50;           % [N]   delta_N in (6.25)

%% ---------------------------------------------------------------- integration
% A stepping-stone run is many steps, and a step that never strikes must not
% burn the whole budget before the run notices.
p.n_steps = 10;                 % steps per stepping-stone run

%% -------------------------------------------------------------- overrides
% Nested structs are addressed with dots: ch6_params('stones.l_max', 0.55).
for k = 1:2:numel(varargin)
    p = setfield_dotted(p, varargin{k}, varargin{k+1});
end

p.n_ctrl = p.bez_deg + 1;
p.ny     = numel(p.iact);

end

% ---------------------------------------------------------------------------
function s = setfield_dotted(s, name, value)
%SETFIELD_DOTTED  Assign p.a.b.c = value from the string 'a.b.c'.
%
% ch3_params rejects unknown field names, and that check is worth keeping --
% a typo'd override that silently does nothing is the worst kind of parameter
% bug. subsref/subsasgn would accept any path, so the existence of every level
% is verified explicitly here first.
parts = strsplit(name, '.');

cur = s;
for i = 1:numel(parts)
    if ~isstruct(cur) || ~isfield(cur, parts{i})
        error('ch6_params:unknownField', 'Unknown parameter "%s".', name);
    end
    cur = cur.(parts{i});
end

sub = struct('type', repmat({'.'}, 1, numel(parts)), 'subs', parts);
s   = subsasgn(s, sub, value);
end
