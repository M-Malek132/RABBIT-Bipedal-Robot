function pass = ch6_test_barrier()
%CH6_TEST_BARRIER  The geometry of Chapter 6, and the claims it is supposed to make.
%
% Two kinds of check, and the first kind is the important one.
%
% GEOMETRIC IMPLICATIONS. Each barrier exists to make a statement -- "inside O1
% implies l_s <= l_max at h_f = 0", "inside O3 implies w_f <= w_max". Those are
% statements about SETS and can be checked exactly by sampling, with no robot,
% no dynamics and no tolerance. If the implication is wrong the barrier can be
% perfectly enforced and still fail to deliver what the chapter says it does,
% and no dynamic test would notice.
%
% DERIVATIVE CHAIN. g, gdot and gddot against finite differences, including the
% explicit time dependence of Section 6.2 -- which is where the easy mistakes
% are, since a static-stone test cannot see a missing dg/dt term at all.
%
% Checks:
%   1. (6.8): inside O1 at h_f = 0  <=>  l_f <= l_max
%   2. (6.8): outside O2 at h_f = 0, l_f > 0  <=>  l_f >= l_min
%   3. O2's arch height is the clearance ch6_R2_from_clearance was asked for
%   4. (6.19): O3 is tangent to w = w_max and passes through the initial point
%   5. inside O3 => w <= w_max ; inside O4 => w >= w_min ; and the thesis's
%      literal "outside O4" form does NOT imply w >= w_min
%   6. gdot and gddot vs finite differences on a moving stone (Section 6.2)
%   7. Section 6.1.1 == the rb = 2 ECBF: ch6_cbf_row's exponential row equals
%      gddot >= -Kb eta_b with Kb from ch5_pole_gain
%   8. the reciprocal and exponential rows agree at h_CBF = 1 and differ the
%      way h^3 vs h does elsewhere
%
% See also CH6_BAR_STONES, CH6_BAR_WIDTH, CH6_CBF_ROW, CH6_R2_FROM_CLEARANCE.

fprintf('\n=== ch6_test_barrier ===\n');
pass = true;
rng(6);

p = ch6_params();

%% ------------------------------------------------ 1 & 2. (6.8) at touchdown
l_min = 0.25; l_max = 0.55;
stone = ch6_resolve_stone(p.stones, l_min, l_max);

% The comparison needs a dead band. g is a difference of two square roots, so
% at l exactly on the boundary it is zero up to rounding and its SIGN is a coin
% flip -- which would show up as a "mismatch" that is really the test asking
% floating point to decide an equality. Points within tol of the boundary are
% skipped and counted, so a dead band that quietly swallowed everything would
% be visible.
tol  = 1e-9;
bad1 = 0; bad2 = 0; skipped = 0;
for l = linspace(-0.1, 0.9, 401)
    k = kin2(l, 0);                       % ON the ground: h_f = 0
    B = ch6_bar_stones(k, stone, 0);

    if abs(l - l_max) > tol
        if (B(1).g >= 0) ~= (l <= l_max), bad1 = bad1 + 1; end
    else
        skipped = skipped + 1;
    end

    if l > 1e-9
        if abs(l - l_min) > tol
            if (B(2).g >= 0) ~= (l >= l_min), bad2 = bad2 + 1; end
        else
            skipped = skipped + 1;
        end
    end
end
ok = bad1 == 0;
fprintf('  [%s] %-42s %d/401 mismatches\n', tf(ok), ...
        'inside O1 at h=0  <=>  l_s <= l_max', bad1);
pass = pass && ok;

ok = bad2 == 0;
fprintf('  [%s] %-42s %d mismatches (%d on the boundary, skipped)\n', tf(ok), ...
        'outside O2 at h=0, l>0  <=>  l_s >= l_min', bad2, skipped);
pass = pass && ok;

%% -------------------------------------------------- 3. R2 sets the clearance
% The lowest h_f for which the foot is still outside O2, swept over l_f, must
% have its maximum at l_f = l_min/2 and equal the requested clearance.
c_req = 0.05;
[R2, c_used] = ch6_R2_from_clearance(l_min, c_req);
st2 = stone; st2.R2 = R2;

h_needed = zeros(1, 201);
ls = linspace(0, l_min, 201);
for i = 1:numel(ls)
    % solve g_ST2 = 0 for h >= 0 -- closed form: the circle's upper arc
    a   = l_min/2;
    rho = sqrt(R2^2 + a^2);
    d2  = rho^2 - (ls(i) - a)^2;
    h_needed(i) = max(0, sqrt(max(d2,0)) - R2);
end
err = abs(max(h_needed) - c_used);
% cross-check that the barrier agrees with that closed form
kk = kin2(l_min/2, c_used);
Bk = ch6_bar_stones(kk, st2, 0);
ok = err < 1e-12 && abs(Bk(2).g) < 1e-12;
fprintf('  [%s] %-42s arch %.4f m, asked %.4f m, g = %.1e\n', tf(ok), ...
        'R2 delivers the requested clearance', max(h_needed), c_used, Bk(2).g);
pass = pass && ok;

%% ------------------------------------------------ 4. (6.19) tangency + point
wp = p.width;
wp.l_f0 = 0;
[~, geom] = ch6_bar_width(kin2(0, wp.w0), wp);

tang_hi = abs((geom.c3(2) + geom.R3) - wp.w_max);      % centre + R3 == w_max
tang_lo = abs((geom.c4(2) - geom.R4) - wp.w_min);      % centre - R4 == w_min
on3 = abs(norm([0; wp.w0] - geom.c3) - geom.R3);
on4 = abs(norm([0; wp.w0] - geom.c4) - geom.R4);
ok  = max([tang_hi tang_lo on3 on4]) < 1e-12;
fprintf('  [%s] %-42s tangency %.1e, through x0 %.1e\n', tf(ok), ...
        '(6.19): O3/O4 tangent and through x0', max(tang_hi,tang_lo), max(on3,on4));
pass = pass && ok;

%% --------------------------- 5. containment implies the bound; outside does not
bad3 = 0; bad4 = 0;
for l = linspace(-0.5, 1.2, 120)
    for w = linspace(-0.1, 0.6, 120)
        B = ch6_bar_width(kin2(l, w), wp);
        if B(1).g >= 0 && w > wp.w_max + 1e-12, bad3 = bad3 + 1; end
        if B(2).g >= 0 && w < wp.w_min - 1e-12, bad4 = bad4 + 1; end
    end
end
ok = bad3 == 0 && bad4 == 0;
fprintf('  [%s] %-42s O3 %d, O4 %d counterexamples\n', tf(ok), ...
        'inside O3/O4 => w in [w_min, w_max]', bad3, bad4);
pass = pass && ok;

% Now the thesis's literal reading. This is expected to FAIL as a bound, and
% the test asserts that it fails -- so the discrepancy documented in
% ch6_bar_width is reproducible rather than an opinion.
wp_lit = wp;  wp_lit.o4_sense = 'outside';
found = 0;
for l = linspace(-0.5, 1.2, 120)
    for w = linspace(-0.1, wp.w_min - 1e-3, 60)
        B = ch6_bar_width(kin2(l, w), wp_lit);
        if B(2).g >= 0, found = found + 1; end
    end
end
ok = found > 0;
fprintf('  [%s] %-42s %d points admitted with w < w_min\n', tf(ok), ...
        'literal "O4F >= R4" does NOT bound w', found);
pass = pass && ok;

%% ------------------------------ 6. derivative chain on a MOVING stone (6.2)
sp = p.stones;
sp.motion  = 'sinusoidal';
sp.amp     = 0.06;
sp.freq    = 1.5;
mov = ch6_resolve_stone(sp, l_min, l_max);

% STENCILS AND STEP SIZES. A second difference divides by h^2, so cancellation
% sets its noise floor at eps/h^2 -- at h = 1e-6 that is 1e-4 RELATIVE, and a
% tolerance tighter than that would be testing the stencil rather than the
% barrier. The fix is a higher-order stencil at a larger step, not a looser
% tolerance: the 5-point second difference is O(h^4), so at h = 1e-3 truncation
% and roundoff are both ~1e-10 and the check can be tight enough to catch a
% missing term rather than merely a wrong one.
hh1 = 1e-5;      % 3-point first difference,  optimal near eps^(1/3)
hh2 = 1e-3;      % 5-point second difference, O(h^4)
e_gd = 0; e_gdd = 0;
for trial = 1:20
    r0 = [0.05 + 0.5*rand; 0.02 + 0.15*rand];
    v0 = 1.5*randn(2,1);
    a0 = 4*randn(2,1);                        % the "u" of this toy flow
    t0 = 0.3*rand;

    k0 = kin2v(r0, v0, a0);
    B0 = ch6_bar_stones(k0, mov, t0);

    for j = 1:2
        gd_fd  = (traj(r0, v0, a0, mov, t0,  hh1, j) - ...
                  traj(r0, v0, a0, mov, t0, -hh1, j)) / (2*hh1);

        gdd_fd = ( -    traj(r0, v0, a0, mov, t0, -2*hh2, j) ...
                   + 16*traj(r0, v0, a0, mov, t0,   -hh2, j) ...
                   - 30*traj(r0, v0, a0, mov, t0,      0, j) ...
                   + 16*traj(r0, v0, a0, mov, t0,    hh2, j) ...
                   -    traj(r0, v0, a0, mov, t0,  2*hh2, j) ) / (12*hh2^2);

        e_gd  = max(e_gd,  abs(gd_fd  - B0(j).gdot) / max(1, abs(B0(j).gdot)));
        gdd   = B0(j).Lf2g;      % this toy flow has rddot = a0 folded into Lf2r
        e_gdd = max(e_gdd, abs(gdd_fd - gdd) / max(1, abs(gdd)));
    end
end
% 1e-7 on the second derivative is the 5-point stencil's OWN truncation on a
% 1.5 Hz sinusoidal stone (measured 1.7e-08), not slack for the barrier: a
% missing term in Lf2g -- a dropped cross partial, a factor of 2 -- is wrong by
% O(1) relative, six orders above this.
ok = e_gd < 1e-8 && e_gdd < 1e-7;
fprintf('  [%s] %-42s gdot %.1e, gddot %.1e (rel)\n', tf(ok), ...
        'moving-stone derivatives vs FD', e_gd, e_gdd);
pass = pass && ok;

%% --------------------------------- 7. Section 6.1.1 IS the rb = 2 ECBF
q = ch6_params('cbf.form', 'exponential', 'cbf.gamma_b', 7, 'cbf.gamma', 3);
Kb_ref = ch5_pole_gain([7 3], 2);
err = 0;
for trial = 1:10
    k = kin2v([0.2+0.3*rand; 0.05+0.1*rand], randn(2,1), randn(2,1));
    B = ch6_bar_stones(k, stone, 0);
    for j = 1:2
        row = ch6_cbf_row(B(j), q);
        % row is  -LgLfg u <= b, i.e. gddot >= -(b - Lf2g) ... expand it:
        % the ECBF statement is gddot >= -Kb eta_b, so the rhs must satisfy
        %     b - Lf2g == Kb * eta_b
        err = max(err, abs((row.b - B(j).Lf2g) - Kb_ref * B(j).eta_b));
        err = max(err, norm(row.Kb - Kb_ref, inf));
    end
end
ok = err < 1e-10;
fprintf('  [%s] %-42s max residual %.1e\n', tf(ok), ...
        '(6.1)+linear decrease == rb=2 ECBF', err);
pass = pass && ok;

%% ------------------------------- 8. reciprocal vs exponential class-K
qr = ch6_params('cbf.form', 'reciprocal', 'cbf.gamma_b', 7, 'cbf.gamma', 3);
k  = kin2v([0.3; 0.08], [1.0; 0.4], [0; 0]);
B  = ch6_bar_stones(k, stone, 0);
r_e = ch6_cbf_row(B(1), q);
r_r = ch6_cbf_row(B(1), qr);
h   = r_e.h_cbf;
% b_exp - b_rec = gamma (h - h^3)
ok = abs((r_e.b - r_r.b) - 3*(h - h^3)) < 1e-9 * max(1, abs(h)^3);
fprintf('  [%s] %-42s h_CBF = %.4f, gap = %.4f\n', tf(ok), ...
        'cubic vs linear class-K, exactly', h, r_e.b - r_r.b);
pass = pass && ok;

fprintf('  --> %s\n', upper(tf(pass)));

end

% ---------------------------------------------------------------------------
function k = kin2(a, b)
%KIN2  A static point: position only, used for the set-membership checks.
k = struct('r', [a; b], 'rdot', [0;0], 'Lf2r', [0;0], 'LgLfr', zeros(2,4));
end

function k = kin2v(r, v, a)
%KIN2V  A moving point with prescribed acceleration folded into Lf2r, so the
% barrier's gddot is the full second derivative of the toy flow.
k = struct('r', r(:), 'rdot', v(:), 'Lf2r', a(:), 'LgLfr', zeros(2,4));
end

function g = gval(r, stone, t, j)
B = ch6_bar_stones(kin2(r(1), r(2)), stone, t);
g = B(j).g;
end

function g = traj(r0, v0, a0, stone, t0, dt, j)
%TRAJ  g along the toy flow r(t) = r0 + dt v0 + dt^2/2 a0, at time t0 + dt.
%
% One place that knows the trajectory, so the first and second differences
% cannot disagree about what curve they are differentiating.
g = gval(r0 + dt*v0 + 0.5*dt^2*a0, stone, t0 + dt, j);
end

function s = tf(b)
if b, s = 'pass'; else, s = 'FAIL'; end
end
