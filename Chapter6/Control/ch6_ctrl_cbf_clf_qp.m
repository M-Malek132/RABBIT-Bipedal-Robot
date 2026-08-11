function [u, qp] = ch6_ctrl_cbf_clf_qp(Lf2y, LgLfy, u_ff, info, B, p, u_ref)
%CH6_CTRL_CBF_CLF_QP  The Chapter-6 controller, (6.12) and (6.26).
%
%   [u, qp] = ch6_ctrl_cbf_clf_qp(Lf2y, LgLfy, u_ff, info, B, p)
%   [u, qp] = ch6_ctrl_cbf_clf_qp(Lf2y, LgLfy, u_ff, info, B, p, u_ref)
%
%       min_{u,d1}  ||mu - mu_ref||^2 + p1 d1^2                       (6.26)
%       s.t.  Vdot(x,mu) + lambda V(eta) <= d1        (CLF, RELAXED)
%             Bdot_i(x,u) + gamma B_i(x) >= 0         (ECBF, HARD)  i = 1..nb
%             Fv_st(x,u) >= delta_N > 0               (Normal force)
%             |Fh_st(x,u)| <= kf Fv_st(x,u)           (Friction cone)
%             u_min <= u <= u_max                     (Input saturation)
%
% with mu = LgLf y (u - u_ff), so every row is affine in u and this is a QP --
% Remark 6.4. Chapter 3's stage-8 QP is this program with mu_ref = 0 and the
% barrier rows deleted; ch6_test_qp asserts that equivalence.
%
% ======================================================== WHAT mu_ref IS FOR
% (6.26) is written with mu_ref = 0: the MINIMUM-NORM control that certifies the
% CLF rate at each instant. That is the literal reading and it is available
% (p.cbf.reference = 'ff'). It is not the default, and the reason is measured.
%
% Chapter 3 already recorded that the min-norm CLF-QP is the WORST of its three
% controllers on this robot -- "minimum-norm is not the same as well-behaved",
% peak torque 203 Nm against PD's 110, max|eta| 2.77 against 0.018 -- because
% asking for the least mu that certifies the rate leaves no margin, and under
% sampling that error compounds.
%
% MEASURED HERE, 10 steps from the reference fixed point with NO BARRIERS AT ALL
% and a foothold window so wide it can never bind:
%
%     PD                        10/10 steps, l_s = 0.3537 every step, 110 Nm
%     min-norm CLF-QP            4/10 steps, l_s drifts 0.3577 -> 0.3705 ->
%                                0.4178 -> the robot falls, 156 Nm
%     Chapter-6 QP, nb = 0       4/10, identical to the line above
%     Chapter-6 QP + barriers    4/10, identical again
%
% The last two lines are the point. The Chapter-6 QP with zero barrier rows
% destabilises the gait exactly as the min-norm CLF-QP does, and adding the
% barriers changes nothing -- so the failure is INHERITED FROM THE COST, and
% attributing it to the barrier (or to pole tuning, or to the stone size) would
% have been wrong in a way no amount of retuning would have fixed.
%
% Setting mu_ref to the PD control makes the program a SAFETY FILTER: with no
% barrier active it returns the PD law exactly, and when one is active it makes
% the smallest change to PD that satisfies it. The CLF row, the barrier rows,
% the contact rows and the box are all unchanged -- this is the objective, not
% the constraints, so nothing that (6.26) guarantees is weakened.
%
% It is also the closer reading of what Section 6.4 says the controller does:
% "the ECBF-CLF-QP controller TRACKS THE OUTPUTS corresponding to this gait by
% solving a quadratic program in real-time to find the control input that
% follows this new gait while maintaining all above constraints".
%
% ============================================================ WHAT BENDS AND WHAT DOESN'T
% ONLY THE CLF ROW IS SLACKED. The barrier rows are hard. When they conflict --
% and on a short stone at speed they genuinely do -- safety wins and tracking
% waits. That is what (6.26) says, and reversing it gives a controller that is
% occasionally unsafe and always convergent, which is the wrong trade.
%
% The consequence is that this QP CAN be infeasible, and that is not a defect
% to be smoothed over: a hard barrier row plus a torque box is a statement about
% the robot's actuators, and when it has no solution the honest report is that
% the requested foothold is not reachable with this gait and these actuators.
% Remark 6.6 attributes part of the residual failure rate in Table 6.1 to
% exactly this. So the infeasible branch below does not quietly relax anything;
% it saturates the feedforward, sets qp.feasible = false, and lets ch6_simulate
% count it and ch6_report print it.
%
% =========================================================== WORKING IN u, NOT mu
% Chapter 5's QPs use mu as the decision variable; this one uses u. The reason
% is the three rows Chapter 5 did not have: the torque box is a plain BOUND in
% u (and a dense linear constraint in mu), and the contact wrench is
% lambda = lam_drift + lam_in u straight out of the same KKT solve that produced
% f and g -- so the friction cone and the minimum normal force are EXACT affine
% rows in u and would have to be composed through LgLfy^-1 in mu. Same feasible
% set either way; one of them is better conditioned and shorter.
%
% The cost is then ||LgLfy (u - u_ref)||^2 rather than ||mu - mu_ref||^2 written
% out -- the same number, since mu - mu_ref = LgLfy (u - u_ref).
%
% ==================================================================== SCALING
% p.clf_slack_penalty is a RELATIVE weight. LgLfy carries the inverse inertia,
% so the mu-block of the Hessian lands around 1e-3 for RABBIT; pairing that with
% a literal 1e6 gives quadprog a condition number near 1e9 and a meaningless
% delta. The penalty is scaled by the mu-block's own magnitude, exactly as in
% ch3_ctrl_clf_qp, so the ratio is p.clf_slack_penalty independent of the
% robot's units.
%
% Barrier rows get row-scaled for the same reason and a sharper one: Lf2g
% contains rdot' Hg rdot, and at the 3.9 m/s the swing foot reaches through a
% circle of radius 0.13 m that term is O(100) while the CLF row is O(1). An
% interior-point solver reads the resulting row-norm spread as ill-conditioning
% and returns a different point.
%
% Inputs
%   Lf2y, LgLfy, u_ff, info : from ch3_io_lin
%   B     : 1 x nb struct array of barriers from ch6_barrier (may be empty)
%   p     : parameter struct
%   u_ref : the torque the cost pulls towards. Omitted or empty means u_ff,
%           i.e. mu_ref = 0, i.e. the literal (6.26). ch6_control supplies the
%           PD control here when p.cbf.reference = 'pd'.
%
% Outputs
%   u  : nu x 1 joint torque
%   qp : struct
%          .V .LfV .LgV .psi .delta .active .exitflag .feasible
%          .h        1 x nb   the barrier values g_i
%          .h_cbf    1 x nb   gamma_b g_i + gdot_i                     (6.1)
%          .margin   1 x nb   hdot_CBF + kappa(h_CBF), >= 0 is the guarantee
%          .cbf_active 1 x nb which barrier rows are tight
%          .n_dead   number of rows with no grip on u
%          .mu       ny x 1
%
% See also CH6_CBF_ROW, CH6_BARRIER, CH3_CTRL_CLF_QP, CH5_CTRL_ECBF_CLF_QP.

nu  = p.nu;
nb  = numel(B);

if nargin < 7 || isempty(u_ref), u_ref = u_ff; end

clf = ch3_res_clf(p);
eta = info.eta;

[V, LfV, LgV] = ch3_clf_eval(eta, clf, p.eps);
psi = LfV + (clf.c3 / p.eps) * V;

qp = struct('V', V, 'LfV', LfV, 'LgV', LgV, 'psi', psi, ...
            'delta', 0, 'active', psi > 0, 'exitflag', 1, 'feasible', true, ...
            'h', nan(1,nb), 'h_cbf', nan(1,nb), 'margin', nan(1,nb), ...
            'cbf_active', false(1,nb), 'n_dead', 0, 'viol', 0, ...
            'mu', zeros(p.ny,1));

A   = LgLfy;
AtA = A.' * A;

%% ------------------------------------------------------------- barrier rows
Ab = zeros(0, nu);
bb = zeros(0, 1);
rows = repmat(ch6_cbf_row(dummy_barrier(nu), p), 1, max(nb,1));
rows = rows(1:nb);

infeasible_row = false;

for i = 1:nb
    r = ch6_cbf_row(B(i), p);
    rows(i) = r;

    qp.h(i)     = B(i).g;
    qp.h_cbf(i) = r.h_cbf;

    if r.live
        [ra, rbs] = scale_row(r.A, r.b);
        Ab = [Ab; ra];      %#ok<AGROW>
        bb = [bb; rbs];     %#ok<AGROW>
    else
        qp.n_dead = qp.n_dead + 1;
        if ~r.vacuous
            % 0*u <= negative. No control can satisfy it; saying so is the
            % point of Section 5.1's failure mode being reported rather than
            % hidden. Do not hand this row to quadprog -- it would report
            % infeasibility with no indication of which row caused it.
            infeasible_row = true;
        end
    end
end

%% ----------------------------------------------------------------- the QP
% z = [u; delta].  Cost ||A(u - u_ref)||^2 + p1 delta^2, expanded; the constant
% u_ref' AtA u_ref is dropped because it does not move the argmin.
H  = 2 * blkdiag(AtA, slack_penalty(AtA, p));
H  = (H + H.')/2;
fq = [-2 * (AtA * u_ref); 0];

%% -------------------------------------------------- contact rows, (6.25)
% Kept SEPARATE from the barrier rows, because they are separate in kind. The
% friction cone and the minimum normal force are statements about the physics
% of the contact: violating them does not degrade the guarantee, it invalidates
% the model -- the stance foot slips or leaves the ground and the single-support
% dynamics being integrated stop describing the robot. The barrier rows are
% statements about a design goal. When the two cannot both hold, the physics has
% to win, and that ordering is only expressible if they are built apart.
Ap = zeros(0, nu);
bp = zeros(0, 1);

lam_d = info.aux.lam_drift;
lam_i = info.aux.lam_in;

if p.limits.enable.friction
    ms = p.limits.mu_s;
    Ap = [Ap; ...
          [ lam_i(1,:) - ms*lam_i(2,:)]; ...
          [-lam_i(1,:) - ms*lam_i(2,:)]];
    bp = [bp; ...
          -( lam_d(1) - ms*lam_d(2)); ...
          -(-lam_d(1) - ms*lam_d(2))];
end

if p.limits.enable.grf
    Ap = [Ap; -lam_i(2,:)];
    bp = [bp;  lam_d(2) - p.limits.Fz_min];
end

%% ---------------------------------------------------------------- assemble
% relaxed CLF row:  (LgV A) u - delta <= -psi + LgV A u_ff
Aineq = [LgV * A, -1];
bineq = -psi + LgV * A * u_ff;

if ~isempty(Ab)
    Aineq = [Aineq; Ab, zeros(size(Ab,1), 1)];
    bineq = [bineq; bb];
end
if ~isempty(Ap)
    Aineq = [Aineq; Ap, zeros(size(Ap,1), 1)];
    bineq = [bineq; bp];
end

if p.limits.enable.torque
    lb = [-p.limits.u_max * ones(nu,1); 0];
    ub = [ p.limits.u_max * ones(nu,1); inf];
else
    lb = [-inf(nu,1); 0];
    ub = [ inf(nu,1); inf];
end

if infeasible_row
    exitflag = -2;
    z = [];
else
    [z, ~, exitflag] = quadprog(H, fq, Aineq, bineq, [], [], lb, ub, ...
                                [u_ref; 0], qp_options());
end

%% ------------------------------------------------------------------ report
if exitflag <= 0 || isempty(z)
    % ===================================================== THE INFEASIBLE BRANCH
    % The hard barrier rows and the input box have no common solution. This is a
    % real statement about the robot -- Remark 6.6 is about exactly it -- and it
    % is REPORTED (qp.feasible = false) whatever happens next.
    %
    % What happens next still matters, because the simulation has to continue in
    % order to say what went wrong. The rows are relaxed IN ORDER OF WHAT THEY
    % MEAN, which is not the same as in order of how much they hurt:
    %
    %   input box            KEPT. An actuator limit cannot be exceeded; a
    %                        controller that "relaxes" it is reporting torques
    %                        the robot cannot produce.
    %   friction cone, GRF   KEPT. Violating these does not degrade a guarantee,
    %                        it invalidates the MODEL -- the stance foot slips or
    %                        unloads and the single-support dynamics being
    %                        integrated stop describing the robot. A run that
    %                        relaxes them reports a trajectory that does not
    %                        exist. (Measured, before this was fixed: |Fx/Fz|
    %                        reached 8.6 against a 0.6 cone during the infeasible
    %                        samples, and every number after that point was
    %                        fiction.)
    %   barrier rows         RELAXED, by one common slack s, minimised.
    %   CLF row              DROPPED. It is the row that was allowed to bend in
    %                        the first place, and keeping it lets a stability
    %                        requirement compete with a safety one at the moment
    %                        safety is already losing.
    %
    % The alternative -- hold the saturated feedforward -- discards every row
    % including the satisfiable ones, and measured, that turns a marginal step
    % into a fall: the run then reports "the robot fell" when the truth is "the
    % QP was a few Nm short for three samples".
    %
    % None of this is a relaxation of (6.26). The program above is solved first
    % and is the one that carries the guarantee; this is the answer to a program
    % that has no solution, not a softer version of one that does. qp.viol says
    % how far short it fell, so "infeasible by 1e-9" and "infeasible by
    % 400 m/s^2" are distinguishable.
    [u, qp.viol] = least_violation(Ab, bb, Ap, bp, lb(1:nu), ub(1:nu), ...
                                   u_ref, AtA);
    qp.exitflag = exitflag;
    qp.feasible = false;
else
    u = z(1:nu);
    qp.delta    = z(end);
    qp.exitflag = exitflag;
    qp.viol     = 0;
end

qp.mu = A * (u - u_ff);

% The residual the guarantee is about, recomputed from the u actually applied
% so it reports the closed loop rather than the solver's opinion of it.
for i = 1:nb
    qp.margin(i)     = rows(i).b - rows(i).A * u;
    qp.cbf_active(i) = rows(i).live && ...
                       (qp.margin(i) <= 1e-7 * max(1, abs(rows(i).b)));
end

end

% ---------------------------------------------------------------------------
function [u, viol] = least_violation(Ab, bb, Ap, bp, lb, ub, u_ref, AtA)
%LEAST_VIOLATION  The best u available when the hard rows cannot all be met.
%
%       min_{u,s}  s + w ||u - u_ref||^2_{AtA}
%       s.t.  Ab u <= bb + s     (barriers, relaxed by one common slack)
%             Ap u <= bp         (friction cone and GRF, still HARD)
%             s >= 0,  lb <= u <= ub
%
% One common slack across the barrier rows, minimised. The tiny quadratic term
% is a tie-breaker, not an objective: without it every u on the optimal face is
% equally good, the solver returns an arbitrary one, and consecutive samples
% jump between them -- which the sampled-data integration sees as a control that
% chatters.
%
% Ap stays hard. See the caller for why: those rows are the physics, and a
% trajectory that violates them is not a trajectory of the model being
% integrated.
nu = numel(u_ref);
w  = 1e-6 * max(1, mean(diag(AtA)));

u_sat = min(max(u_ref, lb), ub);

if isempty(Ab)
    % No barrier rows to relax, so the box (plus possibly the contact rows) was
    % the whole problem and there is nothing left to trade.
    u = u_sat;  viol = 0;
    return;
end

H  = 2 * blkdiag(w * AtA, 0);
H  = (H + H.')/2 + 1e-12 * eye(nu+1);      % PSD -> PD for the solver
fq = [-2 * w * (AtA * u_ref); 1];

A = [Ab, -ones(size(Ab,1), 1)];
b = bb;
if ~isempty(Ap)
    A = [A; Ap, zeros(size(Ap,1), 1)];
    b = [b; bp];
end

[z, ~, ef] = quadprog(H, fq, A, b, [], [], [lb; 0], [ub; inf], ...
                      [u_sat; 0], qp_options());

if ef <= 0 || isempty(z)
    % With s free this can still fail, because Ap and the box are BOTH hard: a
    % pose where no admissible torque keeps the foot inside the friction cone
    % has no answer at this level and the right report is that there was none.
    u    = u_sat;
    viol = inf;
else
    u    = z(1:nu);
    viol = max(0, z(end));
end
end

% ---------------------------------------------------------------------------
function s = slack_penalty(AtA, p)
mu_scale = mean(diag(AtA));
if ~isfinite(mu_scale) || mu_scale <= 0, mu_scale = 1; end
s = p.clf_slack_penalty * mu_scale;
end

% ---------------------------------------------------------------------------
function [a, b] = scale_row(a, b)
%SCALE_ROW  Normalize a row so quadprog sees comparable row norms.
%
% Dividing both sides by the same positive number does not move the halfspace,
% so this is a pure change of units and cannot change the feasible set.
s = max(1, norm([a, b], inf));
a = a / s;
b = b / s;
end

% ---------------------------------------------------------------------------
function B = dummy_barrier(nu)
%DUMMY_BARRIER  A zero barrier, used only to give `rows` its struct type.
%
% MATLAB has no way to preallocate a struct array of a type defined elsewhere
% without an instance of it, and building `rows` by growth instead would make
% the field list of ch6_cbf_row an implicit contract of this file.
B = struct('label', '', 'g', 0, 'gdot', 0, 'eta_b', [0;0], ...
           'Lf2g', 0, 'LgLfg', zeros(1,nu), ...
           'controllable', false, 'singular', false);
end

% ---------------------------------------------------------------------------
function opts = qp_options()
%QP_OPTIONS  Built once. See ch3_ctrl_clf_qp for the measurement that motivates
% the persistent: constructing the options object cost 2.80 ms against a 0.98 ms
% solve, i.e. 74% of the controller's time, and this sits inside an ODE
% right-hand side.
persistent o
if isempty(o)
    o = optimoptions('quadprog', 'Display', 'off', ...
                     'Algorithm', 'interior-point-convex');
end
opts = o;
end
