function [c, ceq] = ch3_col_constraints(z, p)
%CH3_COL_CONSTRAINTS  Stage 3 constraints: periodicity, dynamics, Table 3.1.
%
%   [c, ceq] = ch3_col_constraints(z, p)
%
% EQUALITIES (ceq), in order:
%
%   Node 1 -- put the step on the zero dynamics surface and fix the gauge (13)
%       y(x_1)       = 0        (4)   start on Z
%       ydot(x_1)    = 0        (4)   ... and with zero transverse velocity
%       theta(x_1)   = theta_-  (1)   the phase clock reads 0 at the start
%       P_st(q_1)    = [0;0]    (2)   stance foot at the origin (gauge choice)
%       J_st(q_1)dq_1= 0        (2)   stance foot not moving
%
%   WHY ONLY AT NODE 1.  y = 0, ydot = 0 and the contact conditions are all
%   INVARIANT under the collocation dynamics: u_ff renders ydd = 0 exactly, so
%   y evolves linearly in t, and Hermite-Simpson is exact for cubics -- it
%   reproduces a linear y with no truncation error at all. Imposing them at
%   every node instead would add 12(N-1) equations that the defects already
%   imply, over-determining the system (532 equations against 319 unknowns at
%   N = 21) and leaving fmincon to fight a rank-deficient KKT matrix.
%   ch3_report measures the residual drift across the step to confirm this
%   holds in practice rather than only in theory.
%
%   Dynamics -- Hermite-Simpson defects                     (nx*(N-1))
%
%   Node N -- land on the switching surface S                       (2)
%       theta(x_N) = theta_+    (1)
%       swing-foot height = 0   (1)   x_N is genuinely ON the guard
%
%   Periodicity -- the HZD condition                               (13)
%       Delta(x_N) - x_1 = 0, EXCLUDING px. px is the translation gauge and
%       must advance; requiring it to repeat would demand the robot end where
%       it started, i.e. not walk.
%
%   NEC1 -- average walking rate                                    (1)
%       L_step / T = v_des
%
% INEQUALITIES (c), always length 17 regardless of what is enabled:
%
%       row  constraint                                  gate           source
%      ----  ------------------------------------------  -------------  ------
%        1   mid-step swing-foot clearance               clearance      style
%        2   no ground penetration at interior nodes     (always)       —
%        3   step-length floor                           (always)       —
%        4   peak torque       |u|  <= u_max             torque         Tbl 3.1
%        5   friction cone     |Fx| <= mu_s Fz           friction       NIC2
%        6   minimum normal force  Fz >= Fz_min          grf            NIC1
%        7   impact impulse   ||I|| <= impulse_max       impulse        Tbl 3.1
%        8   hip-height band                             height         style
%        9   swing foot strictly above ground            swing_clear    NIC3
%       10   transversal strike (foot descending at N)   swing_clear    NIC3
%       11   post-impact swing-leg lift-off              liftoff        NEC2
%       12   impact impulse compressive   Iz >= 0        impact         NEC3
%       13   impact impulse inside the friction cone     impact         NEC3
%       14   existence of the fixed point                hzd            NEC4
%       15   stability of the fixed point                hzd            NEC5
%       16   theta strictly monotonic                    phase_mono     HH6
%       17   decoupling matrix invertible on Z           decoupling     HH2
%
%   GATED CONSTRAINTS ARE HELD AT -1, NOT REMOVED.  A disabled inequality is
%   trivially satisfied but still present, so c has a constant length and
%   fmincon's problem dimensions never change between runs. That is what makes
%   "enable them one at a time, warm-starting each phase" a safe workflow
%   instead of a new problem each time.
%
% -------------------------------------------------------------------------
% SECTION 6.3.4 AND WHAT IT ASKS FOR
%
% The book's constraint set is chosen so that a feasible point satisfies the
% stability conditions (5.79)/(5.80), the gait hypotheses HGW2, HI3, HGW6, the
% output hypotheses HH2, HH4, HH5, HH6, and a desired average walking speed.
% Those hypotheses are the PURPOSE of the constraints, not extra rows -- except
% where this transcription would otherwise leave one unchecked, which is why
% HH2 and HH6 get rows of their own (16, 17).
%
%   HGW2 (single support, with double support instantaneous)
%        rows 6 and 9: the stance foot stays loaded and the swing foot stays
%        strictly clear, so the robot is in single support throughout and
%        touches down only at the end.
%   HI3  (the impact is a valid rigid plastic impact)
%        rows 11, 12, 13: the trailing foot leaves, the impulse is compressive,
%        and it does not slip.
%   HGW6 (forward progression)
%        rows 3 and 16: a positive step length and a strictly increasing theta.
%   HH2  (decoupling matrix invertible) -> row 17.
%   HH6  (theta a valid phase variable) -> row 16.
%   HH4/HH5 (hybrid invariance, Delta(S ∩ Z) ⊂ Z)
%        NOT a row. In this transcription it is IMPLIED: node 1 is pinned to
%        y = ydot = 0 and the periodicity block equates Delta(x_N) to node 1,
%        so Delta(x_N) lands on Z whenever both are satisfied. Adding it again
%        would duplicate rows the periodicity block already spans and hand
%        fmincon a rank-deficient KKT matrix -- the same reason node 1 is the
%        only node where y = 0 is imposed. E.eta_post measures it instead, and
%        ch3_report prints it, so the implication is verified rather than
%        assumed.
%   (5.79)/(5.80) -> rows 14 and 15, via ch3_zero_dynamics.
%   desired average walking speed -> NEC1, in the EQUALITY block above.
%
% WHERE THE BOOK'S "NEC" LABEL IS LOOSE.  NEC2-NEC5 are filed under nonlinear
% EQUALITY constraints but every one of them is written as an inequality
% ("is positive", "zeta*_2 > ...", "0 < delta^2 < 1"), and they are imposed
% that way here. NEC1 is the only member of the list that is genuinely an
% equality, and it is the only one in ceq.
%
% See also CH3_COL_EVAL, CH3_ZERO_DYNAMICS, CH3_COL_COST, CH3_COL_SOLVE.

E  = ch3_col_eval(z, p);
nq = p.nq;
N  = E.N;

X     = E.X;
alpha = E.alpha;

%% ============================== EQUALITIES ==============================

% --- node 1: on Z, phase at zero, stance foot pinned at the origin --------
[y1, yd1] = ch3_outputs(X(:,1), alpha, p);
theta1    = p.c_theta(:).' * X(1:nq, 1);
foot1     = P_st(X(1:nq, 1));
footvel1  = J_st(X(1:nq, 1)) * X(nq+1:2*nq, 1);

ceq_start = [ y1; ...
              yd1; ...
              theta1 - p.theta_minus; ...
              foot1; ...
              footvel1 ];

% --- dynamics -------------------------------------------------------------
ceq_dyn = E.defect(:);

% --- node N: on the switching surface S -----------------------------------
thetaN = p.c_theta(:).' * X(1:nq, N);
ceq_end = [ thetaN - p.theta_plus; ...
            E.sw_h(N) ];

% --- periodicity through Delta, excluding px ------------------------------
ceq_per = E.x_next(2:end) - X(2:end, 1);

% --- NEC1: average walking rate -------------------------------------------
% GATED, and worth gating.  Finding a PERIODIC gait is the hard part of this
% problem; pinning the speed at the same time is what makes a cold solve
% stall. Measured: from the hand seed, periodicity alone starts at a residual
% of 0.78 while NEC1 starts at 0 (the seed trivially walks at its own speed) --
% so demanding a DIFFERENT speed immediately fights the constraint that was
% already the binding one, and feasibility oscillates instead of settling.
%
% The workflow this enables: solve with NEC1 off to get any periodic gait,
% then switch it on at the achieved speed and march v_des with
% ch3_continuation. Held at 0 when disabled so ceq keeps a constant length.
if ~isfield(p,'enforce_nec1') || p.enforce_nec1
    ceq_rate = E.L_step / E.T - p.v_des;
else
    ceq_rate = 0;
end

ceq = [ceq_start; ceq_dyn; ceq_end; ceq_per; ceq_rate];

%% ============================= INEQUALITIES =============================
c = -ones(17, 1);

% 1. swing-foot clearance at mid-step (NIC3): height >= clearance
k_mid = max(2, min(N-1, round((N+1)/2)));
c(1)  = p.limits.clearance - E.sw_h(k_mid);
if ~p.limits.enable.clearance, c(1) = -1; end

% 2. no ground penetration at interior nodes (endpoints are on the ground by
%    construction, so they are excluded)
if N > 2
    c(2) = -min(E.sw_h(2:N-1));
else
    c(2) = -1;
end

% 3. step-length floor
c(3) = p.step_len_min - E.L_step;

% 4. peak torque, over nodes AND midpoints (a midpoint can exceed both of its
%    neighbours, so checking nodes only would under-report the peak)
if p.limits.enable.torque
    c(4) = max([abs(E.u(:)); abs(E.um(:))]) - p.limits.u_max;
end

% 5. friction cone: |Fx| - mu_s Fz <= 0
if p.limits.enable.friction
    lam_all = [E.lam, E.lamm];
    c(5) = max(abs(lam_all(1,:)) - p.limits.mu_s * lam_all(2,:));
end

% 6. minimum normal force: Fz_min - Fz <= 0
if p.limits.enable.grf
    lam_all = [E.lam, E.lamm];
    c(6) = p.limits.Fz_min - min(lam_all(2,:));
end

% 7. impact impulse magnitude
if p.limits.enable.impulse
    c(7) = norm(E.impulse) - p.limits.impulse_max;
end

% 8. hip-height band: |h_hip - hip_h| <= hip_h_tol, over nodes AND midpoints.
%    pz is DOWN-positive, so world height is -pz and the torso base IS the hip
%    (ch3_body_points takes it from the same transform) -- the height reads
%    straight off row 2, no forward kinematics. Midpoints are included for the
%    reason the torque check includes them: a Hermite-Simpson midpoint is a
%    real point of the trajectory and can leave a band both its neighbours
%    stay inside.
%
%    ONE BAND DOES TWO JOBS. Its CENTRE sets how tall the robot stands, which
%    is what takes the crouch out; its WIDTH caps the vertical travel of the
%    hip, which is what stops the gait bouncing. Pinning only the mean would
%    leave the bob free; capping only the bob would let the optimizer keep the
%    crouch it already prefers, since a low hip buys knee torque cheaply.
%
%    Like qt_range this is a DESIGN REQUIREMENT, not a physical limit, and it
%    fights the objective the same way: expect J to rise. Tighten it the way
%    the Table 3.1 limits are tightened -- incrementally, warm-starting each
%    stage (ch3_posture_march), never straight to the target.
if p.limits.enable.height
    h_hip = -[X(2,:), E.xm(2,:)];
    c(8)  = max(abs(h_hip - p.limits.hip_h)) - p.limits.hip_h_tol;
end

%% ------------------------- Section 6.3.4: NIC3, NEC2-NEC5, HH2, HH6 ------

% 9. NIC3(a). The swing foot is STRICTLY above the ground during the step, so
%    that S meets Z only at the end.  Row 2 already forbids penetration; this
%    is the strict version, and the difference matters: at height exactly zero
%    the guard is grazed and the step can terminate anywhere in the interior.
%
%    WHICH POINTS IT COVERS, AND WHY NOT ALL OF THEM.  The foot is at height
%    zero at node 1 (it has just impacted) and at node N (the guard equality
%    puts it there), so demanding a positive margin at either is infeasible by
%    construction.  The two ADJACENT MIDPOINTS are the subtle case: midpoint 1
%    sits half an interval after lift-off and midpoint N-1 half an interval
%    before strike, so at a typical h = 0.02 s and a lift-off rate of 0.3 m/s
%    the foot is only ~3 mm up there. A 1 cm margin would reject a perfectly
%    good gait for being sampled too close to its own endpoints. They are held
%    to non-negativity instead, which row 2 already covers, and the margin
%    applies to interior nodes 2..N-1 and interior midpoints 2..N-2.
if p.limits.enable.swing_clear
    h_int = E.sw_h(2:N-1);
    if N > 3
        h_int = [h_int, E.sw_hm(2:N-2)];
    end
    if isempty(h_int)
        c(9) = -1;
    else
        c(9) = p.limits.sw_clear_min - min(h_int);
    end
end

% 10. NIC3(b). The crossing at node N is TRANSVERSAL -- the swing foot is
%     genuinely moving downward when it strikes.  Without this, "height = 0 at
%     node N" is equally satisfied by a foot that touches and lifts, which is a
%     tangency: the guard is not crossed, the impact never fires, and the
%     Poincare map is not locally well defined.  sw_hd is d/dt of the swing-foot
%     height, so descending means sw_hd <= -sw_strike_rate.
if p.limits.enable.swing_clear
    c(10) = E.sw_hd(N) + p.limits.sw_strike_rate;
end

% 11. NEC2. The vertical component of the post-impact swing-leg velocity is
%     positive.  Delta relabels, so the new swing leg is the leg that was
%     STANCE during the step just finished: this is the condition that the
%     trailing foot leaves the ground rather than scuffing along it.  Fail it
%     and the "instantaneous double support" of the gait hypotheses is not
%     instantaneous -- the model claims single support while two feet are on
%     the floor.
if p.limits.enable.liftoff
    c(11) = p.limits.liftoff_rate - E.sw_hd_post;
end

% 12-13. NEC3. Validity of the impact.  The rigid plastic impact model returns
%     whatever impulse solves its momentum balance; it does not check that the
%     ground could apply it.  Two things have to hold for the model to mean
%     anything:
%
%       12. COMPRESSIVE.  Iz >= 0. A negative normal impulse is the floor
%           pulling the foot down, which is not a contact force.
%       13. NO SLIP.  |Ix| <= mu_s Iz, the friction cone at impulse level. The
%           impact map imposes J_sw dq+ = 0, i.e. the foot sticks; that
%           assumption is only self-consistent if the tangential impulse
%           required to make it stick is one friction can supply.
%
%     Same mu_s as NIC2, and written as |Ix| - mu_s Iz <= 0 rather than as a
%     ratio so it stays finite as Iz -> 0.
if p.limits.enable.impact
    mu_imp = p.limits.mu_s_impact;
    if isempty(mu_imp), mu_imp = p.limits.mu_s; end
    c(12) = p.limits.impulse_z_min - E.impulse(2);
    c(13) = abs(E.impulse(1)) - mu_imp * E.impulse(2);
end

% 14-15. NEC4 and NEC5, the stability conditions (5.79) and (5.80).
%
%     These are the only rows that do not read off the trajectory: they are
%     properties of the RESTRICTED POINCARE MAP of the hybrid zero dynamics,
%     which ch3_zero_dynamics builds from alpha alone by quadrature.
%
%       14. EXISTENCE   zeta*_2 > V_zero^MAX / delta_zero^2
%           Multiply by delta^2 and the left side is the kinetic coordinate at
%           the START of the step; the condition then says the robot never
%           stalls before the swing foot strikes.
%       15. STABILITY   0 < delta_zero^2 < 1
%           delta_zero^2 is the slope of a scalar affine map, hence its
%           eigenvalue.
%
%     COST.  One evaluation is p.hzd_grid_solve points of roughly four 7x7
%     solves each, and fmincon finite-differences it once per decision
%     variable. That is why this gate defaults to false and why the grid used
%     here is the coarse one -- see the note in ch3_params.
if p.limits.enable.hzd
    Zd = ch3_zero_dynamics(alpha, p, p.hzd_grid_solve);
    c(14) = Zd.nec4;
    c(15) = Zd.nec5;
end

% 16. HH6: theta is a valid phase variable, i.e. STRICTLY MONOTONIC along the
%     step.  This is what lets s replace time. If thetadot reaches zero the
%     phase stalls, s stops advancing, and yd(s) freezes mid-step; if it
%     changes sign the gait runs backwards through its own profile. Checked at
%     midpoints too, since thetadot can dip between two healthy nodes.
%
%     It also carries the computable content of HGW6: thetadot > 0 throughout,
%     together with the step-length floor in row 3, is forward progression.
if p.limits.enable.phase_mono
    c(16) = p.limits.thetadot_min - min([E.thd, E.thdm]);
end

% 17. HH2: the decoupling matrix LgLf y is invertible on Z.  Everything from
%     stage 4 onwards divides by it -- u_ff, the PD law, both QPs -- so this is
%     the hypothesis the entire control design rests on, and a gait that drifts
%     toward a singular decoupling matrix produces enormous torques near the
%     singularity while looking perfectly feasible in position.
%
%     Bounded below by the SMALLEST SINGULAR VALUE, which is literally the
%     distance to the nearest singular matrix. A determinant floor would not do:
%     det is a product, so it stays comfortably large while one singular value
%     collapses and another blows up.
if p.limits.enable.decoupling
    c(17) = p.limits.dec_min - E.dec_min;
end

end
