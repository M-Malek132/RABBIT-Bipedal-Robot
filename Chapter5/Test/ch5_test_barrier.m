function pass = ch5_test_barrier()
%CH5_TEST_BARRIER  The Lie stack eta_b of (5.10), against an independent route.
%
% THE ONE TEST THIS FILE EXISTS FOR is check 3. The pendulum's barrier
% derivatives come from committed symbolic code in Model/Generated, and
% symbolic code that is wrong is wrong SILENTLY -- it produces smooth,
% plausible, finite numbers and the ECBF row then enforces a condition on
% something that is not the fourth derivative of anything. Comparing it against
% the source it was generated from would prove only that the generator ran.
%
% So it is checked against a genuinely different computation: integrate the
% DRIFT flow xdot = f(x) to tight tolerance, sample h along it, and
% differentiate with a 9-point stencil. Nothing symbolic, nothing generated,
% and a different discretization of the same mathematical object.
%
% Tolerances there are RELATIVE and loose by design. An order-8 stencil at
% h = 2e-3 recovers a fourth derivative to roughly 1e-7 relative, and the
% quantities involved reach 1e5, so demanding 1e-12 absolute would be
% demanding the finite differences be exact. What the check catches is a wrong
% derivative, which is wrong by O(1) relative, not by 1e-7.
%
% Checks:
%   1. eta_b's leading entry is h, and h carries the right sign convention
%   2. springmass: eta_b against the matrix powers, computed independently
%   3. pendulum: eta_b and L_f^rb h against the drift-flow finite differences
%   4. L_g L_f^(rb-1) h against a direct finite difference in u
%   5. the reciprocal barrier's (5.6) identity  Bdot <= gamma/B  <=>  hdot >= -gamma h^3
%   6. L_g h vanishes exactly for rb > 1, and does not for rb == 1
%
% See also CH5_BARRIER, CH5_GEN_PENDULUM, CH5_LIE_ROWS.

fprintf('\n=== ch5_test_barrier ===\n');
pass = true;

%% -------------------------------------------------- 1. sign convention
p = ch5_params();
xt = ch5_x0(p); xt(3) = 3.4;                 % push x3 past x3max = 3.15
b  = ch5_barrier(xt, p);
ok = (b.h < 0) && abs(b.h - (3.15 - 3.4)) < 1e-14 && (b.eta_b(1) == b.h);
fprintf('  [%s] %-34s h = %+.4f at x3 = 3.4 (limit 3.15)\n', tf(ok), ...
        'h < 0 outside, eta_b(1) == h', b.h);
pass = pass && ok;

q  = ch5_params('system','pendulum');
zt = ch5_x0(q); zt(1) = 0; zt(2) = 0;        % arm straight DOWN, py = -2
bq = ch5_barrier(zt, q);
ok = (bq.h < 0) && abs(bq.q + 2) < 1e-12;
fprintf('  [%s] %-34s py = %+.4f, h = %+.4f (limit %.1f)\n', tf(ok), ...
        'lower bound flips sign the same way', bq.q, bq.h, q.constraint.value);
pass = pass && ok;

%% ------------------------------- 2. springmass eta_b vs independent powers
sm = ch5_springmass(p);
c  = -[0 0 1 0 0 0];
rng(11);
e = 0; eg = 0; ef = 0;
for j = 1:6
    xs = randn(6,1);
    bs = ch5_barrier(xs, p);
    ref = zeros(6,1);
    for i = 0:5
        ref(i+1) = c * (sm.A^i) * xs;        % matrix power, built here
    end
    ref(1) = ref(1) + p.constraint.value;
    e  = max(e,  norm(bs.eta_b - ref, inf) / max(1, norm(ref, inf)));
    ef = max(ef, abs(bs.Lfrb - c*(sm.A^6)*xs));
    eg = max(eg, abs(bs.LgLfrb1 - c*(sm.A^5)*sm.B));
end
pass = rep('springmass eta_b == c A^i x', e,  1e-12, pass);
pass = rep('springmass L_f^6 h',          ef, 1e-11, pass);
pass = rep('springmass L_g L_f^5 h',      eg, 1e-14, pass);

%% ------------------- 3. pendulum eta_b vs finite differences on the drift flow
% JUDGED AGAINST THE STENCIL'S OWN NOISE FLOOR, not a hand-picked tolerance.
%
% Each derivative is a sum of sampled values divided by h^k, so any error in
% those samples -- here the ODE solver's, ~1e-12 -- is amplified by
% sum|c_k|/h^k. At h = 2e-3 that is 1e-9 for the first derivative and 9e-4 for
% the third. A fixed tolerance of 1e-5 on the third derivative is therefore not
% a strict test, it is an IMPOSSIBLE one: it asks the finite differences to be
% three orders of magnitude better than their own arithmetic allows, and it
% fails on correct code. (It did.)
%
% So the floor is computed and reported next to the error. A discrepancy at the
% floor means the two routes agree as well as this comparison can resolve; a
% wrong derivative is wrong by O(1) RELATIVE and clears the floor by many
% orders, which is what this check is actually able to detect.
rng(12);
[floors, sc] = fd_floors(2e-3, 1e-12);
worst = zeros(1,4);
mags  = zeros(1,4);
for j = 1:3
    xs = [randn(2,1); randn(2,1); 0.5*randn(2,1); 0.5*randn(2,1)];
    bs = ch5_barrier(xs, q);

    fd  = fd_drift(xs, q);
    gen = [bs.eta_b(2:4); bs.Lfrb].';

    worst = max(worst, abs(gen - fd));
    mags  = max(mags,  abs(fd));
end

% Relative tolerances, loosening with order because the resolution does. A
% WRONG derivative is wrong by O(1) relative, so even the loosest row below
% still catches one with four orders of margin. Check 3b beneath this is the
% sharp check; this one establishes that the theta derivatives the sharp check
% builds on are themselves right.
names   = {'hdot','hddot','h^(3)','L_f^4 h'};
rel_tol = [1e-8, 1e-7, 1e-4, 1e-4];
for k = 1:4
    tol = max(rel_tol(k) * mags(k), 5 * floors(k));
    ok  = worst(k) <= tol;
    fprintf('  [%s] pendulum %-8s vs FD        err = %.2e  tol = %.2e ', ...
            tf(ok), names{k}, worst(k), tol);
    fprintf('(rel %.1e, floor %.1e)\n', worst(k)/max(mags(k),eps), floors(k));
    pass = pass && ok;
end
fprintf('        (stencil sum|c|/h^k = %s)\n', ...
        strjoin(arrayfun(@(v) sprintf('%.1e', v), sc, 'Uni', 0), ' '));

%% -------- 3b. the SHARP check: py's Lie stack against the chain rule, exactly
% Finite differences cannot resolve a fourth derivative better than a few parts
% in 1e5. This can, to machine precision, because it never differentiates
% numerically at all.
%
% py depends on theta only through two angles, phi = theta1 and psi =
% theta1+theta2:
%
%       py = -l1 cos(phi) - l2 cos(psi)
%
% so its time derivatives follow from phi's and psi's by the scalar Faa di
% Bruno expansion of -l cos(.), written out below. Feeding it the theta
% derivatives from ch5_pend_theta_lie must reproduce ch5_pend_py_lie exactly.
%
% This checks the COMPOSITION -- the part with tensor contractions in it and by
% far the likeliest place for a symbolic derivation to be subtly wrong -- and
% it checks the two generated files against each other, so a fault in either
% one alone shows up. It does not re-check the theta derivatives themselves;
% that is what check 3 above is for.
rng(15);
e_chain = zeros(1,5);
for j = 1:6
    xs = [randn(2,1); randn(2,1); 0.6*randn(2,1); 0.6*randn(2,1)];

    [t2, t3, t4f, t4g] = ch5_pend_theta_lie(xs, q.plant.pv);
    [py, py1, py2, py3, py4f, py4g] = ch5_pend_py_lie(xs, q.plant.pv);

    l1 = q.plant.l(1); l2 = q.plant.l(2);
    t1d = xs(5:6);

    % phi = theta1 and psi = theta1 + theta2, with their four derivatives
    ph = [xs(1),        t1d(1), t2(1), t3(1), t4f(1)];
    ps = [xs(1)+xs(2),  sum(t1d), sum(t2), sum(t3), sum(t4f)];

    r1 = mcos(l1, ph) + mcos(l2, ps);           % [py, py1, py2, py3, py4]

    % the input part: only the phi^(4) / psi^(4) terms carry tau
    in_ref = l1*sin(ph(1)) * t4g(1,:) + l2*sin(ps(1)) * (t4g(1,:) + t4g(2,:));

    got = [py, py1, py2, py3, py4f];
    e_chain(1:5) = max(e_chain(1:5), abs(got - r1) ./ max(1, abs(r1)));
    e_chain(5)   = max(e_chain(5), ...
                       norm(py4g(:).' - in_ref, inf) / max(1, norm(in_ref, inf)));
end
pass = rep('py stack == chain rule (exact)', max(e_chain), 1e-12, pass);

%% ---------------------------- 4. the input coefficient, by difference in u
% h^(4) is affine in tau, so differencing two DIFFERENT constant inputs isolates
% L_g L_f^3 h exactly -- no derivative approximation in the u direction at all,
% only in the time derivative. And the two h^(4) values share the same noise
% floor, so differencing them cancels most of it: this is the sharpest of the
% pendulum checks despite being the highest order.
rng(13);
e_in = 0; mag_in = 0;
for j = 1:3
    xs = [randn(2,1); randn(2,1); 0.4*randn(2,1); 0.4*randn(2,1)];
    bs = ch5_barrier(xs, q);
    h4a = fd_h4_with_u(xs, q, zeros(2,1));
    for col = 1:2
        du = zeros(2,1); du(col) = 1;
        h4b = fd_h4_with_u(xs, q, du);
        e_in = max(e_in, abs((h4b - h4a) - bs.LgLfrb1(col)) / ...
                          max(1, abs(bs.LgLfrb1(col))));
        mag_in = max(mag_in, abs(bs.LgLfrb1(col)));
    end
end
fprintf('  [%s] %-34s err = %.3e (tol %.0e, |L_gL_f^3h| ~ %.1e)\n', ...
        tf(e_in <= 1e-4), 'pendulum L_gL_f^3h vs FD (rel)', e_in, 1e-4, mag_in);
pass = pass && (e_in <= 1e-4);

%% ------------------------------------- 5. the (5.6) reciprocal identity
% Bdot <= gamma/B with B = 1/h says exactly hdot >= -gamma h^3. Checked on the
% relative-degree-1 constraint, where both sides are nonvacuous.
pv1 = ch5_params('constraint', struct('type','v1_max','value',1.2), ...
                 'ecbf.poles', 2.0);
gam = pv1.cbf.gamma;
rng(14);
e_rec = 0;
for j = 1:6
    xs = 0.4*randn(6,1);
    bs = ch5_barrier(xs, pv1);
    if ~bs.rec_defined, continue; end
    us = randn();
    hdot  = bs.Lfh + bs.Lgh*us;
    Bdot  = bs.LfB_rec + bs.LgB_rec*us;
    % Bdot - gamma/B  and  -(hdot + gamma h^3)/h^2  must agree
    e_rec = max(e_rec, abs((Bdot - gam*bs.h) + (hdot + gam*bs.h^3)/bs.h^2) / ...
                        max(1, abs(Bdot)));
end
pass = rep('B=1/h: Bdot<=g/B == hdot>=-g h^3', e_rec, 1e-11, pass);

%% -------------------------------------- 6. L_g h and the relative degree
cases = {p, 'x3_max (rb 6)', false; q, 'py_min (rb 4)', false; ...
         pv1, 'v1_max (rb 1)', true};
for i = 1:size(cases,1)
    pc = cases{i,1};
    bc = ch5_barrier(ch5_x0(pc), pc);
    ok = (bc.rec_controllable == cases{i,3});
    fprintf('  [%s] %-34s L_g h = %s, controllable = %d\n', tf(ok), ...
            cases{i,2}, mat2str(bc.Lgh, 3), bc.rec_controllable);
    pass = pass && ok;
end

fprintf('--- ch5_test_barrier: %s ---\n', tf(pass));

end

% ---------------------------------------------------------------------------
function d = mcos(l, a)
%MCOS  The four time derivatives of  -l*cos(a(t))  given a and its derivatives.
%
%   a = [a, adot, addot, a^(3), a^(4)]
%
% Faa di Bruno for g(a(t)) with g = -l cos:
%
%   g'    =  l sin a       d/dt   = g' adot
%   g''   =  l cos a       d2/dt2 = g'' adot^2 + g' addot
%   g'''  = -l sin a       d3/dt3 = g''' adot^3 + 3 g'' adot addot + g' a3
%   g'''' = -l cos a       d4/dt4 = g'''' adot^4 + 6 g''' adot^2 addot
%                                   + 3 g'' addot^2 + 4 g'' adot a3 + g' a4
%
% The 6, 3, 4 in the last line are the Bell-polynomial coefficients and are the
% part worth writing down carefully -- they are what a hand derivation of a
% relative-degree-4 barrier gets wrong.
s = sin(a(1)); c = cos(a(1));
g1 =  l*s;  g2 =  l*c;  g3 = -l*s;  g4 = -l*c;

ad = a(2); add = a(3); a3 = a(4); a4 = a(5);

d = [ -l*c, ...
      g1*ad, ...
      g2*ad^2 + g1*add, ...
      g3*ad^3 + 3*g2*ad*add + g1*a3, ...
      g4*ad^4 + 6*g3*ad^2*add + 3*g2*add^2 + 4*g2*ad*a3 + g1*a4 ];
end

% ---------------------------------------------------------------------------
function C = fd_stencils()
%FD_STENCILS  Central 9-point coefficients for the first four derivatives.
%
% Rows 1..4 give f^(k)(0) ~ (C(k,:) * v) / h^k on the grid (-4:4)*h. Declared
% once so that fd_drift, fd_h4_with_u and fd_floors cannot disagree about which
% stencil is in use -- if they did, the floor would be computed for one method
% and the error for another, and the test would be meaningless in a way that
% still printed numbers.
%
% Sanity: the odd rows are antisymmetric and the even rows symmetric, and every
% row sums to zero. fd_floors asserts both.
C = [ 1/280, -4/105,    1/5,   -4/5,      0,    4/5,   -1/5,  4/105, -1/280 ; ...
     -1/560,  8/315,   -1/5,    8/5, -205/72,   8/5,   -1/5,  8/315, -1/560 ; ...
     -7/240,   3/10, -169/120, 61/30,     0, -61/30, 169/120, -3/10,  7/240 ; ...
      7/240,   -2/5,  169/60, -122/15,  91/8, -122/15, 169/60, -2/5,  7/240 ];
end

function [floors, scaled] = fd_floors(hs, eps_traj)
%FD_FLOORS  The smallest discrepancy this comparison can resolve, per order.
%
% A derivative formed as (sum_k c_k v_k)/h^K amplifies any error in v by
% sum|c_k|/h^K. Here v comes from a tightly integrated trajectory, so eps_traj
% is the ODE solver's own error -- not machine epsilon, which would understate
% it by three orders.
%
% This is why check 3 reports a floor at all: at h = 2e-3 the third derivative
% carries a floor near 1e-3, and no correct implementation can agree with the
% finite differences more closely than that.
C = fd_stencils();

assert(all(abs(sum(C, 2)) < 1e-12), 'fd_stencils rows must sum to zero');
assert(norm(C(1,:) + fliplr(C(1,:)), inf) < 1e-12, 'row 1 must be antisymmetric');
assert(norm(C(2,:) - fliplr(C(2,:)), inf) < 1e-12, 'row 2 must be symmetric');

scaled = sum(abs(C), 2).' ./ hs.^(1:4);
floors = eps_traj * scaled;
end

function d = fd_drift(x0, p)
%FD_DRIFT  h', h'', h''' and h'''' along xdot = f(x), by a 9-point stencil.
%
% The independent route: no symbolic code, no generated code, just the model's
% own f integrated tightly and h sampled on a uniform grid around t = 0.
d = fd_stack(x0, p, @(~, z) drift(z, p));
end

function d4 = fd_h4_with_u(x0, p, u)
%FD_H4_WITH_U  h'''' along the flow with a CONSTANT input u held.
d  = fd_stack(x0, p, @(tt, zz) ch5_ode_rhs(tt, zz, u, p));
d4 = d(4);
end

function d = fd_stack(x0, p, odef)
%FD_STACK  Sample h along a given flow and apply all four stencils.
hs   = 2e-3;
ts   = (-4:4) * hs;
opts = odeset('RelTol', 1e-13, 'AbsTol', 1e-14);

v = zeros(1, numel(ts));
for i = 1:numel(ts)
    if ts(i) == 0
        z = x0;
    else
        [~, Z] = ode113(odef, [0 ts(i)], x0, opts);
        z = Z(end,:).';
    end
    bb = ch5_barrier(z, p);
    v(i) = bb.h;
end

C = fd_stencils();
d = (C * v.').' ./ hs.^(1:4);
end

function dz = drift(z, p)
[f, ~] = ch5_control_affine(z, p);
dz = f;
end

function ok = rep(name, err, tol, ok_in)
ok = ok_in && (err <= tol);
fprintf('  [%s] %-34s err = %.3e (tol %.0e)\n', tf(err <= tol), name, err, tol);
end

function s = tf(b)
if b, s = 'PASS'; else, s = 'FAIL'; end
end
