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
% INEQUALITIES (c), always length 8 regardless of what is enabled:
%
%   1  swing-foot clearance at mid-step        (NIC3)
%   2  no ground penetration at interior nodes
%   3  step-length floor
%   4  peak torque            |u|  <= u_max          (gated)
%   5  friction cone          |Fx| <= mu_s Fz        (gated)
%   6  minimum normal force    Fz  >= Fz_min         (gated)
%   7  impact impulse        ||I|| <= impulse_max    (gated)
%   8  reserved / spare, held at -1
%
%   GATED CONSTRAINTS ARE HELD AT -1, NOT REMOVED.  A disabled inequality is
%   trivially satisfied but still present, so c has a constant length and
%   fmincon's problem dimensions never change between runs. That is what makes
%   "enable them one at a time, warm-starting each phase" a safe workflow
%   instead of a new problem each time.
%
% See also CH3_COL_EVAL, CH3_COL_COST, CH3_COL_SOLVE.

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
c = -ones(8, 1);

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

end
