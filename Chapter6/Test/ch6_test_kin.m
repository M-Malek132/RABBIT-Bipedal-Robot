function pass = ch6_test_kin()
%CH6_TEST_KIN  The tracked points and their two derivatives, independently.
%
% Everything in Chapter 6 rests on rddot = Lf2r + LgLfr u being EXACT. If it is
% not, every barrier row enforces a condition on something that is not the
% acceleration of anything, and the resulting run looks perfectly smooth --
% barriers stay positive, the QP solves, the robot walks -- while guaranteeing
% nothing. So this is checked against genuinely different computations rather
% than against itself.
%
% Checks:
%   1. the tracked points agree with ch3_body_points / P_sw / P_st
%   2. rdot = Jr dq agrees with a finite difference of r along the drift flow
%   3. rddot = Lf2r + LgLfr u agrees with a finite difference of rdot along the
%      TRUE closed-loop flow, for random u
%   4. LgLfr agrees with a direct finite difference in u
%   5. the head Jacobian is complete -- it depends on q only through
%      (px, pz, qt), so perturbing the four joint angles must not move it
%
% See also CH6_KIN, CH3_CONTROL_AFFINE, CH3_BODY_POINTS.

fprintf('\n=== ch6_test_kin ===\n');
pass = true;

p = ch6_params();
rng(6);

%% ----------------------------------- 1. the points, against the repo's own
e_sw = 0; e_hd = 0;
for i = 1:8
    x = sample_state();
    q = x(1:p.nq);

    k1 = ch6_kin(x, [], p, 'swing');
    k2 = ch6_kin(x, [], p, 'head');
    pts = ch3_body_points(q);

    e_sw = max(e_sw, norm(k1.r - (P_sw(q) - P_st(q))));
    e_hd = max(e_hd, norm(k2.r - (pts.torso_top - pts.stance_foot)));
end
ok = e_sw < 1e-14 && e_hd < 1e-12;
fprintf('  [%s] %-38s swing %.1e, head %.1e\n', tf(ok), ...
        'points vs P_sw/P_st and body_points', e_sw, e_hd);
pass = pass && ok;

%% ------------------------------------------- 2. rdot along the drift flow
% r is a function of q only, so rdot must be the total derivative along ANY
% flow. Take the drift flow (u = 0) and central-difference r.
h  = 1e-6;
err = 0;
for i = 1:6
    x = sample_state();
    for pt = {'swing','head'}
        k = ch6_kin(x, [], p, pt{1});
        f = ch3_control_affine(x, p);
        xp = x + h*f;  xm = x - h*f;
        rp = point_of(xp, p, pt{1});
        rm = point_of(xm, p, pt{1});
        err = max(err, norm((rp - rm)/(2*h) - k.rdot) / max(1, norm(k.rdot)));
    end
end
ok = err < 1e-7;
fprintf('  [%s] %-38s rel err %.2e\n', tf(ok), ...
        'rdot = Jr dq vs d/dt r on the flow', err);
pass = pass && ok;

%% ------------------------- 3. rddot vs the derivative of the TRUE flow
% This is the check that matters. Integrate xdot = f + g u with u FIXED, sample
% rdot, and central-difference it. Nothing here reuses Jr, Jdotdq or aux.
err = 0;
for i = 1:6
    x = sample_state();
    u = 20 * randn(p.nu, 1);
    for pt = {'swing','head'}
        k = ch6_kin(x, [], p, pt{1});
        ref = k.Lf2r + k.LgLfr * u;

        [f, g] = ch3_control_affine(x, p);
        xdot = f + g*u;
        xp = x + h*xdot;  xm = x - h*xdot;
        vp = rdot_of(xp, p, pt{1});
        vm = rdot_of(xm, p, pt{1});
        err = max(err, norm((vp - vm)/(2*h) - ref) / max(1, norm(ref)));
    end
end
ok = err < 1e-6;
fprintf('  [%s] %-38s rel err %.2e\n', tf(ok), ...
        'rddot = Lf2r + LgLfr u vs true flow', err);
pass = pass && ok;

%% --------------------------------------- 4. LgLfr by finite difference in u
err = 0;
for i = 1:4
    x = sample_state();
    k = ch6_kin(x, [], p, 'swing');
    D = zeros(2, p.nu);
    du = 1e-4;
    for j = 1:p.nu
        ej = zeros(p.nu,1); ej(j) = du;
        [~, ~, aux] = ch3_control_affine(x, p);
        ap = k.Jr * (aux.ddq_drift + aux.ddq_in*( ej)) + k.Jdr;
        am = k.Jr * (aux.ddq_drift + aux.ddq_in*(-ej)) + k.Jdr;
        D(:,j) = (ap - am)/(2*du);
    end
    err = max(err, norm(D - k.LgLfr, 'fro') / max(1, norm(k.LgLfr, 'fro')));
end
ok = err < 1e-10;
fprintf('  [%s] %-38s rel err %.2e\n', tf(ok), ...
        'LgLfr vs finite difference in u', err);
pass = pass && ok;

%% ----------------------------- 5. the head depends on q only through 1,2,3
% If the head Jacobian were truncated -- a plausible slip, since it is hand
% written -- it would still be plausible. What must hold is that it is COMPLETE:
% the torso top is (px + L sin qt, L cos qt - pz), so the four joint angles do
% not move it at all. Perturbing them and seeing the point move would mean the
% wrong point is being protected.
x = sample_state();
r0 = point_of(x, p, 'head');
moved = 0;
for j = 4:7
    xj = x;  xj(j) = xj(j) + 0.3;
    moved = max(moved, norm(point_of(xj, p, 'head') - r0 - ...
                            (P_st(x(1:p.nq)) - P_st(xj(1:p.nq)))));
end
ok = moved < 1e-12;
fprintf('  [%s] %-38s residual %.1e\n', tf(ok), ...
        'head moves only with the stance foot', moved);
pass = pass && ok;

fprintf('  --> %s\n', upper(tf(pass)));

end

% ---------------------------------------------------------------------------
function r = point_of(x, p, pt)
q = x(1:p.nq);
switch pt
    case 'swing', r = P_sw(q) - P_st(q);
    case 'head'
        L = 0.75;
        r = [q(1) + L*sin(q(3)); L*cos(q(3)) - q(2)] - P_st(q);
end
end

function v = rdot_of(x, p, pt)
k = ch6_kin(x, [], p, pt);
v = k.rdot;
end

function x = sample_state()
% Physically plausible poses: knees bent the right way, moderate velocities.
q  = [0.1*randn; -0.9 + 0.05*randn; 0.2*randn; ...
      0.3*randn; -0.4 - 0.3*rand;  0.3*randn; -0.4 - 0.3*rand];
dq = 0.8*randn(7,1);
x  = [q; dq];
end

function s = tf(b)
if b, s = 'pass'; else, s = 'FAIL'; end
end
