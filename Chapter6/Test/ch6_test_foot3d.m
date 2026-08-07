function pass = ch6_test_foot3d()
%CH6_TEST_FOOT3D  Section 6.3's three cases on the swing-foot surrogate.
%
% ch6_foot3d's header sets out what this can and cannot establish: DURUS is not
% modelled here, so what is tested is the CONTROLLER carrying the Section 6.3
% barriers, not DURUS's ability to supply the accelerations it asks for.
%
% Checks, mirroring Section 6.3.3's three cases:
%   1. Case 1 (step length only): the foot lands inside [l_min, l_max] from a
%      nominal profile aimed elsewhere
%   2. Case 2 (step width only): lands inside [w_min, w_max], and the step
%      LENGTH is left alone -- Fig. 6.21's "the small vertical bars show the
%      step length doesn't change"
%   3. Case 3 (both): four barrier rows carried simultaneously, both windows hit
%   4. the barriers hold along the whole run, not only at touchdown
%   5. with the rows removed the nominal profile MISSES -- otherwise cases 1-3
%      prove nothing
%   6. the RECOVERY case: a commanded width window that does NOT contain the
%      previous step width, where (6.19) has no circle through the start point
%      and the barrier begins negative
%
% ON THE WIDTH COMMANDS. Cases 2 and 3 move the width by 1.7 cm, not by the
% 5-10 cm Section 6.3.3 quotes, and that is forced by the geometry rather than
% chosen: (6.19) needs w_min < w0 < w_max, so a +-2.5 cm window can only be
% centred within 2.5 cm of the previous width. Check 6 is what covers the rest
% of Section 6.3.3's commanded range, and it is a weaker claim -- recovery, not
% invariance. See ch6_bar_width.
%
% See also CH6_FOOT3D, CH6_BAR_WIDTH, CH6_BAR_STONES.

fprintf('\n=== ch6_test_foot3d ===\n');
pass = true;

p = ch6_params();
p.cbf.gamma_b = 30;
p.cbf.gamma   = 30;
p.control_dt  = 1e-3;

spec = struct('T', 0.45, 'l_nom', 0.366, 'w_nom', 0.233, 'w0', 0.233, ...
              'h_apex', 0.10, 'a_max', 60, 'stone_sz', 0.05);

W_CMD = 0.25;          % window [0.225, 0.275]; contains w0 = 0.233

%% -------------------------------------------------- 5. the nominal misses
% ONLY THE LENGTH CLAIM CAN BE MADE HERE, and the reason is the geometry rather
% than the tuning. (6.19) needs w_min < w0 < w_max, and the window is w_d +-
% 2.5 cm, so a containment-legal width command is within 2.5 cm of the previous
% width -- which means the window ALWAYS contains w0, and a nominal profile that
% holds the width constant is always already inside it. Case 2's containment
% version is therefore a maintenance problem, never a reaching one; the reaching
% version is check 6, and it is the weaker recovery claim.
s0 = spec;  s0.l_d = 0.28;  s0.w_d = W_CMD;  s0.cases = 'none';
r0 = ch6_foot3d(p, s0);
miss_l  = abs(r0.l_s - s0.l_d) > s0.stone_sz/2;
w_inside = abs(r0.w_s - s0.w_d) <= s0.stone_sz/2;
ok = miss_l && w_inside;
fprintf('  [%s] %-42s l_s %.4f (want %.3f, missed), w_s %.4f (already inside)\n', ...
        tf(ok), 'unconstrained nominal misses on length', ...
        r0.l_s, s0.l_d, r0.w_s);
pass = pass && ok;

%% ---------------------------------------------------------- 1. Case 1
s1 = spec;  s1.l_d = 0.28;  s1.w_d = 0.233;  s1.cases = 'length';
r1 = ch6_foot3d(p, s1);
ok = r1.ok && r1.in_l;
fprintf('  [%s] %-42s l_s = %.4f in [%.3f %.3f]\n', tf(ok), ...
        'Case 1: step length only', r1.l_s, r1.window_l(1), r1.window_l(2));
pass = pass && ok;

%% ---------------------------------------------------------- 2. Case 2
s2 = spec;  s2.l_d = 0.366;  s2.w_d = W_CMD;  s2.cases = 'width';
r2 = ch6_foot3d(p, s2);
dl = abs(r2.l_s - r0.l_s);
ok = r2.ok && r2.in_w && dl < 0.05;
fprintf('  [%s] %-42s w_s = %.4f in [%.3f %.3f], l_s moved %.4f m\n', tf(ok), ...
        'Case 2: step width only, length free', r2.w_s, ...
        r2.window_w(1), r2.window_w(2), dl);
pass = pass && ok;

%% ---------------------------------------------------------- 3. Case 3
s3 = spec;  s3.l_d = 0.30;  s3.w_d = W_CMD;  s3.cases = 'both';
r3 = ch6_foot3d(p, s3);
ok = r3.ok && r3.in_l && r3.in_w;
fprintf('  [%s] %-42s l_s = %.4f, w_s = %.4f, 4 rows\n', tf(ok), ...
        'Case 3: length and width together', r3.l_s, r3.w_s);
pass = pass && ok;

%% ------------------------------- 4. the barriers hold along the whole run
% Per barrier, not a single minimum over all four: which one dips is the whole
% diagnostic. Both width circles pass through the start point BY CONSTRUCTION,
% so they begin at g = 0 exactly with zero lateral velocity -- marginally
% admissible in the sense of Corollary 5.2 -- and the tolerance has to admit
% that without admitting a real violation.
worst = min(r3.log.h, [], 2).';
lbl = {'l<=lmax', 'l>=lmin', 'w<=wmax', 'w>=wmin'};
ok = all(worst > -1e-9);
fprintf('  [%s] %-42s min per barrier %s\n', tf(ok), ...
        'g_i >= 0 for the whole swing', ...
        strjoin(arrayfun(@(i) sprintf('%s %+.1e', lbl{i}, worst(i)), 1:4, ...
                         'UniformOutput', false), ', '));
pass = pass && ok;

%% ---------------------- 6. the recovery case, where (6.19) has no circle
% w_d = 0.30 puts the window at [0.275, 0.325], which does NOT contain the
% previous width 0.233. There is no circle tangent to w = w_max through the
% start point, so g starts NEGATIVE: the barrier is being asked to recover, not
% to maintain. Two things must hold, and they are different claims:
%   - the construction reports what it is doing (geom.recovering)
%   - the ECBF still drives g up and the foot still arrives, which is what an
%     exponential CBF can do and a reciprocal one cannot (B = 1/h is not even
%     defined at h <= 0)
kd = struct('r', [0; spec.w0], 'rdot', [0;0], 'Lf2r', [0;0], 'LgLfr', zeros(2,3));
wp = p.width;  wp.w_min = 0.275;  wp.w_max = 0.325;  wp.w0 = spec.w0;
[Bw, gw] = ch6_bar_width(kd, wp);
% B(2) is the LOWER bound (O4). w0 = 0.233 is below w_min = 0.275, so it is O4
% whose circle cannot contain the start point; O3 (the upper bound) is fine and
% would show g = 0, which is why checking the wrong one passes for free.
ok = gw.recovering && Bw(2).g < -1e-6 && abs(Bw(1).g) < 1e-9;
fprintf('  [%s] %-42s g_lower = %+.4f, g_upper = %+.1e, flagged %d\n', tf(ok), ...
        'w0 outside the window is flagged, g < 0', Bw(2).g, Bw(1).g, gw.recovering);
pass = pass && ok;

s6 = spec;  s6.l_d = 0.366;  s6.w_d = 0.30;  s6.cases = 'width';
r6 = ch6_foot3d(p, s6);
recovered = r6.log.h(2, end) > r6.log.h(2, 1);
ok = r6.ok && recovered;
fprintf('  [%s] %-42s g: %+.4f -> %+.4f, w_s = %.4f (want %.3f)\n', tf(ok), ...
        'the ECBF recovers from g < 0', r6.log.h(2,1), r6.log.h(2,end), ...
        r6.w_s, s6.w_d);
pass = pass && ok;

fprintf('  --> %s\n', upper(tf(pass)));

end

function s = tf(b)
if b, s = 'pass'; else, s = 'FAIL'; end
end
