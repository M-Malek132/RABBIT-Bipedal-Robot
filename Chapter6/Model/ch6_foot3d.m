function out = ch6_foot3d(p, spec)
%CH6_FOOT3D  Section 6.3 on a 3D swing-foot surrogate, not on DURUS.
%
%   out = ch6_foot3d(p, spec)
%
% ============================================================== WHAT THIS IS FOR
% Section 6.3 validates simultaneous step-length and step-width constraints on
% DURUS, a 23-DoF 3D humanoid. THERE IS NO DURUS MODEL IN THIS REPOSITORY, and
% inventing one would make every number in Section 6.3.3 unfalsifiable -- a
% 23-DoF model nobody can check is worse than no model, because it produces
% plots that look like results.
%
% So Section 6.3 is split into the part that does not need the robot and the
% part that does.
%
%   WHAT IS CHECKED EXACTLY, with no robot at all (ch6_test_barrier):
%     * (6.19): the circles O3/O4 are tangent to w = w_max / w = w_min and pass
%       through the initial foot position.
%     * containment implies the bound, on sampled points, for both O3 and O4.
%     * the derivative chain g, gdot, gddot against finite differences.
%
%   WHAT IS CHECKED HERE, on a surrogate:
%     * that the ECBF-CLF-QP of (6.26), carrying all four barrier rows at once,
%       drives a swing foot to a touchdown inside BOTH windows, from a nominal
%       profile that would have missed, under an input bound.
%
%   WHAT IS NOT CHECKED ANYWHERE, and is not claimed:
%     * that DURUS's 23-DoF dynamics can supply the required foot acceleration,
%       or the specific ranges quoted in Section 6.3.3 (l_d in [22, 50] cm etc.).
%       Those are properties of that robot, and this repository cannot speak to
%       them.
%
% ------------------------------------------------------------ what the surrogate is
% The swing foot as a double integrator in its own coordinates:
%
%       r = (l_f, w_f, h_f),     rddot = u,     |u_i| <= spec.a_max
%
% This is the subsystem the Section 6.3 constraints are WRITTEN IN -- (6.17)
% through (6.20) mention nothing but (l_f, w_f, h_f). A real robot differs by
% supplying rddot through LgLf r, a configuration-dependent 3 x nu map, which
% changes how much acceleration is available and in which directions; it does
% not change the barrier construction, the QP structure, or which sets are
% forward invariant. Calling the input bound an "acceleration limit" rather than
% a torque limit is the honest name for that difference.
%
% The nominal swing profile is deliberately WRONG for the commanded foothold:
% it aims at (spec.l_nom, spec.w_nom) while the stone is at (spec.l_d,
% spec.w_d). That is the situation Section 6.3 is about -- one nominal gait,
% many footholds -- and it is what makes the barrier rows do work rather than
% ride along inactive.
%
% Inputs
%   p    : parameter struct (uses p.cbf, p.stones, p.width, p.eps, p.Q_clf,
%          p.clf_construction, p.clf_slack_penalty, p.control_dt)
%   spec : struct
%          .T        nominal step duration [s]
%          .l_nom .w_nom   where the nominal profile aims [m]
%          .l_d .w_d       where the stone actually is [m]
%          .h_apex   nominal mid-step foot clearance [m]
%          .w0       initial (previous) step width [m]
%          .a_max    input bound [m/s^2]
%          .cases    'length' | 'width' | 'both'   (Cases 1, 2, 3)
%
% Output
%   out : struct .t .r (3 x n) .u .l_s .w_s .h_end .in_l .in_w .feasible
%         .log (h, margins, active flags) .spec .geom
%
% See also CH6_BAR_WIDTH, CH6_BAR_STONES, CH6_CBF_ROW, CH6_TEST_BARRIER.

spec = fill_defaults(spec);

%% ---------------------------------------------------------------- the plant
nu = 3;
dt = p.control_dt;
if dt <= 0, dt = 1e-3; end

% CLF on eta = [y; ydot], y = r - r_nom(t). Same RES-CLF construction as
% Chapter 3, at ny = 3 instead of 4 -- ch3_res_clf reads its dimension off p.ny.
p3 = p;
p3.ny    = nu;
p3.Kp    = eye(nu) * 100;
p3.Kd    = eye(nu) *  20;
p3.Q_clf = eye(2*nu);
clf = ch3_res_clf(p3);

%% ------------------------------------------------------ the foothold windows
stone = ch6_resolve_stone(p.stones, ...
                          spec.l_d - spec.stone_sz/2, ...
                          spec.l_d + spec.stone_sz/2);

wp = p.width;
wp.w_min = spec.w_d - spec.stone_sz/2;
wp.w_max = spec.w_d + spec.stone_sz/2;
wp.w0    = spec.w0;
wp.l_f0  = 0;

use_len = any(strcmpi(spec.cases, {'length', 'both'}));
use_wid = any(strcmpi(spec.cases, {'width',  'both'}));

%% ------------------------------------------------------------------ rollout
% START ON THE NOMINAL, POSITION AND VELOCITY. The arch leaves the ground at
% dh/dt = h_apex*pi/T, so a foot starting at REST is already 0.7 m/s below its
% own reference and spends the first tens of milliseconds catching up -- during
% which it dips below z = 0 and the touchdown guard fires at l_f = 0.04 m.
% Measured; it looked like a controller failure and was an initialisation one.
%
% It also matters for the barrier, not just the guard: circle O2 passes through
% the origin of the swing-foot frame, so g_ST2(0) = 0 EXACTLY here (unlike on
% the robot, where the swing foot starts behind the stance foot at -L_prev).
% Admissibility then reduces to gdot(0) >= 0, i.e. the foot must be leaving the
% circle -- which is true at the nominal lift-off velocity and false at rest.
[r_n0, v_n0] = nominal(0, spec);
x  = [r_n0; v_n0];                  % [l w h, dl dw dh]
t  = 0;
T_cap = 2 * spec.T;

n_max = ceil(T_cap/dt) + 2;
Lg = struct('t', nan(1,n_max), 'r', nan(3,n_max), 'u', nan(3,n_max), ...
            'h', nan(4,n_max), 'margin', nan(4,n_max), ...
            'active', false(4,n_max), 'feasible', true(1,n_max), 'n', 0);

fired = false;

while t < T_cap - eps(T_cap)

    [B, nb] = barriers(x, stone, wp, use_len, use_wid);
    [u, ok, marg, act] = solve_qp(x, t, spec, clf, p3, B, nb);

    Lg.n = Lg.n + 1;  i = Lg.n;
    Lg.t(i)      = t;
    Lg.r(:,i)    = x(1:3);
    Lg.u(:,i)    = u;
    Lg.feasible(i) = ok;
    for j = 1:nb
        Lg.h(j,i)      = B(j).g;
        Lg.margin(j,i) = marg(j);
        Lg.active(j,i) = act(j);
    end

    % Touchdown: h_f crosses zero downward, after the foot has actually lifted.
    s = ode45(@(~, xx) [xx(4:6); u], [t, min(t+dt, T_cap)], x, ...
              odeset('RelTol', 1e-9, 'AbsTol', 1e-10, ...
                     'Events', @(tt, xx) guard_event(tt, xx, spec)));

    t = s.x(end);
    x = s.y(:,end);

    if ~isempty(s.ie), fired = true; break; end
    if ~all(isfinite(x)), break; end
end

%% ------------------------------------------------------------------ report
n = Lg.n;
out = struct();
out.t = Lg.t(1:n);
out.r = Lg.r(:,1:n);
out.u = Lg.u(:,1:n);
out.log = struct('h', Lg.h(:,1:n), 'margin', Lg.margin(:,1:n), ...
                 'active', Lg.active(:,1:n), 'feasible', Lg.feasible(1:n));

out.l_s   = x(1);
out.w_s   = x(2);
out.h_end = x(3);
out.T     = t;
out.ok    = fired;
out.feasible = all(Lg.feasible(1:n));

% TOLERANCE. A working barrier lands the foot ON the edge of the window -- that
% is what "outside O2" means when the stone is short -- so the acceptance test
% is decided at the boundary and a 1e-9 tolerance is asking floating point to
% settle it. 1 micron is far below anything physical and far above the
% integrator's own noise.
TOL = 1e-6;
lv = ch6_stone_level(stone, t);
out.in_l = use_len && (out.l_s >= lv.l_min - TOL) && (out.l_s <= lv.l_max + TOL);
out.in_w = use_wid && (out.w_s >= wp.w_min  - TOL) && (out.w_s <= wp.w_max  + TOL);
out.window_l = [lv.l_min, lv.l_max];
out.window_w = [wp.w_min,  wp.w_max];

k_dummy = struct('r', [0;spec.w0], 'rdot', [0;0], 'Lf2r', [0;0], ...
                 'LgLfr', zeros(2,3));
[~, out.geom] = ch6_bar_width(k_dummy, wp);
out.spec  = spec;
out.stone = stone;

end

% ===========================================================================
function [B, nb] = barriers(x, stone, wp, use_len, use_wid)
%BARRIERS  The 2D slices the Chapter-6 barrier functions live in.
%
% Step length is a constraint in (l_f, h_f); step width is a constraint in
% (l_f, w_f). Each gets its own kinematics struct picking the right two of the
% three coordinates, and then the SAME ch6_bar_stones / ch6_bar_width used on
% the robot apply unchanged -- which is the point of ch6_kin returning a plain
% (r, rdot, Lf2r, LgLfr) contract rather than something RABBIT-shaped.
%
% For a double integrator Lf2r = 0 and LgLfr is a selection matrix.
B  = [];
nb = 0;

if use_len
    k = struct('r', [x(1); x(3)], 'rdot', [x(4); x(6)], 'Lf2r', [0;0], ...
               'LgLfr', [1 0 0; 0 0 1]);
    B  = [B, ch6_bar_stones(k, stone, 0)];
    nb = nb + 2;
end

if use_wid
    k = struct('r', [x(1); x(2)], 'rdot', [x(4); x(5)], 'Lf2r', [0;0], ...
               'LgLfr', [1 0 0; 0 1 0]);
    B  = [B, ch6_bar_width(k, wp)];
    nb = nb + 2;
end
end

% ===========================================================================
function [u, ok, marg, act] = solve_qp(x, t, spec, clf, p3, B, nb)
%SOLVE_QP  (6.26) for the surrogate: CLF row relaxed, barrier rows hard, box.

nu = 3;
[rn, vn, an] = nominal(t, spec);

y    = x(1:3) - rn;
ydot = x(4:6) - vn;
eta  = [y; ydot];

[V, LfV, LgV] = ch3_clf_eval(eta, clf, p3.eps); %#ok<ASGLU>
psi = LfV + (clf.c3 / p3.eps) * V;

% ydd = u - an, so mu = u - an and the CLF row is linear in u.
H  = 2 * blkdiag(eye(nu), p3.clf_slack_penalty);
fq = [-2*an; 0];

Aineq = [LgV, -1];
bineq = -psi + LgV * an;

marg = nan(1, nb);
act  = false(1, nb);
rows = cell(1, nb);

for i = 1:nb
    r = ch6_cbf_row(B(i), p3);
    rows{i} = r;
    if r.live
        s = max(1, norm([r.A, r.b], inf));
        Aineq = [Aineq; r.A/s, 0];    %#ok<AGROW>
        bineq = [bineq; r.b/s];       %#ok<AGROW>
    end
end

lb = [-spec.a_max*ones(nu,1); 0];
ub = [ spec.a_max*ones(nu,1); inf];

opts = optimoptions('quadprog','Display','off','Algorithm','interior-point-convex');
[z, ~, ef] = quadprog(H, fq, Aineq, bineq, [], [], lb, ub, [an; 0], opts);

if ef <= 0 || isempty(z)
    u  = max(min(an, spec.a_max), -spec.a_max);
    ok = false;
else
    u  = z(1:nu);
    ok = true;
end

for i = 1:nb
    marg(i) = rows{i}.b - rows{i}.A * u;
    act(i)  = rows{i}.live && marg(i) <= 1e-7 * max(1, abs(rows{i}.b));
end
end

% ===========================================================================
function [r, v, a] = nominal(t, spec)
%NOMINAL  The one nominal swing profile, aimed at the WRONG place on purpose.
%
% Cubic in l and w (zero end velocities), a sine arch in h so the foot leaves
% and returns to the ground at t = 0 and t = T. The barrier's job is to bend
% this into the commanded window; if the nominal already hit it, the run would
% prove nothing.
T = spec.T;
s = min(max(t/T, 0), 1);

sp  = 3*s^2 - 2*s^3;            % cubic ramp, sdot = 0 at both ends
dsp = (6*s - 6*s^2)/T;
ddsp= (6 - 12*s)/T^2;

l  = spec.l_nom * sp;      dl = spec.l_nom * dsp;   ddl = spec.l_nom * ddsp;
w  = spec.w0 + (spec.w_nom - spec.w0)*sp;
dw = (spec.w_nom - spec.w0)*dsp;
ddw= (spec.w_nom - spec.w0)*ddsp;

h   =  spec.h_apex * sin(pi*s);
dh  =  spec.h_apex * (pi/T) * cos(pi*s);
ddh = -spec.h_apex * (pi/T)^2 * sin(pi*s);

r = [l; w; h];
v = [dl; dw; dh];
a = [ddl; ddw; ddh];
end

% ===========================================================================
function [value, isterminal, direction] = guard_event(t, x, spec)
% Hold positive until the foot has had time to lift, exactly as ch6_step does
% for the robot: h_f = 0 at t = 0 by construction and a guard checked there
% would fire immediately.
if t < 0.2*spec.T
    value = 1;
else
    value = x(3);
end
isterminal = 1;
direction  = -1;
end

% ===========================================================================
function s = fill_defaults(s)
d = struct('T', 0.4, 'l_nom', 0.366, 'w_nom', 0.233, ...
           'l_d', 0.30, 'w_d', 0.28, 'h_apex', 0.10, 'w0', 0.233, ...
           'a_max', 40, 'cases', 'both', 'stone_sz', 0.05);
f = fieldnames(d);
for i = 1:numel(f)
    if ~isfield(s, f{i}) || isempty(s.(f{i})), s.(f{i}) = d.(f{i}); end
end
end
