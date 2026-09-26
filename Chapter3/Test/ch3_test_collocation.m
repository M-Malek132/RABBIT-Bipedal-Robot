function ch3_test_collocation()
%CH3_TEST_COLLOCATION  Stage-3 verification: transcription, seed quality, timing.
%
% Checks the transcription BEFORE spending minutes in fmincon:
%   1. pack/unpack round-trips exactly.
%   2. The seed rollout lands on the zero dynamics surface and stays there.
%   3. Seed defects are small -- the rollout nearly satisfies the dynamics
%      already, which is the whole point of seeding by simulation.
%   4. Constraint/variable counts leave a positive number of degrees of
%      freedom (an over-determined transcription cannot be solved).
%   5. Hermite-Simpson defects converge at high order in h -- ENFORCED on the
%      Bezier basis every stored gait is solved with, and recorded as a
%      characterization (an expected failure) on the cubic B-spline default,
%      whose interior knots cap the order the maximum defect can show.
%   6. With p.free_theta the phase endpoints ride in z: pack/unpack round-trip
%      them, a fixed-theta gait augments to the same point, and the effective
%      parameters read them back.
%   7. One cost + constraint evaluation is timed, so the cost of a full
%      gradient is known before launching a solve.

fprintf('\n=== ch3_test_collocation ===\n');
p = ch3_params();
pass = true;

%% 1. pack / unpack
rng(7);
N  = p.N_nodes;
Xr = randn(p.nx, N); Tr = 0.7; ar = randn(p.ny, p.n_ctrl);
zr = ch3_col_pack(Xr, Tr, ar, p);
[X2, T2, a2] = ch3_col_unpack(zr, p);
err = max([norm(X2 - Xr, inf), abs(T2 - Tr), norm(a2 - ar, inf)]);
pass = report('pack/unpack round trip', err, 0, pass);

%% 2-3. seed quality
t0 = tic;
[z0, si] = ch3_col_seed(p);
t_seed = toc(t0);
fprintf('        seed: T = %.4f s, guard fired = %d, theta_end = %.4f (target %.4f), %.2f s\n', ...
        si.T, si.guard_fired, si.theta_reached, p.theta_plus, t_seed);

E = ch3_col_eval(z0, p);

% eta along the seed rollout: does the zero dynamics surface stay invariant?
eta_max = 0;
for k = 1:E.N
    [yk, ydk] = ch3_outputs(E.X(:,k), E.alpha, p);
    eta_max = max(eta_max, norm([yk; ydk], inf));
end
pass = report('seed stays on Z (max|eta|)', eta_max, 1e-4, pass);

def_max = max(abs(E.defect(:)));
fprintf('        seed max|defect| = %.3e   (h = %.4f s)\n', def_max, E.h);
pass = pass && def_max < 1e-2;

fprintf('        seed gait: L_step = %.4f m, T = %.4f s, speed = %.4f m/s (v_des %.2f)\n', ...
        E.L_step, E.T, E.L_step/E.T, p.v_des);
fprintf('        seed peak |u| = %.1f Nm, Fz in [%.1f %.1f] N, |impulse| = %.2f Ns\n', ...
        max([abs(E.u(:)); abs(E.um(:))]), min([E.lam(2,:) E.lamm(2,:)]), ...
        max([E.lam(2,:) E.lamm(2,:)]), norm(E.impulse));

%% 4. degrees of freedom
[c, ceq] = ch3_col_constraints(z0, p);
n_var = numel(z0);
n_eq  = numel(ceq);
dof   = n_var - n_eq;
fprintf('        %d variables, %d equalities, %d inequalities -> %d DOF\n', ...
        n_var, n_eq, numel(c), dof);
ok = dof > 0;
fprintf('  [%s] %-30s\n', tf(ok), 'transcription not over-determined');
pass = pass && ok;

% which equalities are actually violated at the seed?
i = 0;
blocks = {'node-1 (on Z, phase, contact)', 13; ...
          'Hermite-Simpson defects', p.nx*(E.N-1); ...
          'node-N (on guard S)', 2; ...
          'periodicity through Delta', 13; ...
          'NEC1 walking rate', 1};
fprintf('        seed equality residuals by block:\n');
for b = 1:size(blocks,1)
    n = blocks{b,2};
    fprintf('           %-32s max = %.3e\n', blocks{b,1}, max(abs(ceq(i+1:i+n))));
    i = i + n;
end

%% 5. Hermite-Simpson order of accuracy
% Refining the node count must drop the defect like a high power of h. This is
% the sharpest single check that the transcription formula is right: a
% mis-placed h/8 or h/6 still gives small defects but degrades the observed
% order.
%
% ENFORCED ON THE BEZIER BASIS. That is the basis every stored gait is solved
% with (posture_195 included), so it is the one whose accuracy the chapter's
% numbers rest on. The Bezier curve is a polynomial in s, smooth everywhere,
% and its max defect falls at order 3.5-4.5 over this ladder.
[ord_bez, ~] = defect_orders(ch3_params('basis', 'bezier'), 'Bezier, degree 5');
ok = ~isempty(ord_bez) && median(ord_bez) > 3.0;
fprintf('  [%s] %-30s median order = %.2f (expect > 3)\n', tf(ok), ...
        'HS order, Bezier basis', median(ord_bez));
pass = pass && ok;

% A CHARACTERIZATION ON THE CUBIC B-SPLINE DEFAULT, EXPECTED TO FAIL. Its
% interior knots make yd only C^2 there, so the MAX defect stalls at the knot
% phases and the observed order collapses (the median defect still falls at
% order ~5, see the Chapter-3 report). That is a property of the basis, not of
% the transcription, so it is reported as xfail and does not fail the suite.
% If it ever PASSES, the explanation above has stopped being true -- that is
% reported as XPASS and fails the suite, so the note cannot go stale silently.
[ord_bsp, ~] = defect_orders(ch3_params('basis', 'bspline', 'bsp_deg', 3), ...
                             'cubic B-spline (default)');
ok_bsp = ~isempty(ord_bsp) && median(ord_bsp) > 3.0;
if ~ok_bsp
    fprintf('  [%s] %-30s median order = %.2f (expect > 3) -- expected\n', ...
            'xfail', 'HS order, cubic B-spline', median(ord_bsp));
    fprintf('        interior knots cap what the MAX defect can show; the\n');
    fprintf('        enforced check above is on the basis the gaits use.\n');
else
    fprintf('  [%s] %-30s median order = %.2f -- the expected failure now\n', ...
            'XPASS', 'HS order, cubic B-spline', median(ord_bsp));
    fprintf('        passes: update this check and the report that cites it.\n');
    pass = false;
end

%% 6. free phase endpoints (p.free_theta)
pf = ch3_params('free_theta', true);
thr = [-0.17; 0.33];
zf  = ch3_col_pack(Xr, Tr, ar, pf, thr);
[X3, T3, a3, th3] = ch3_col_unpack(zf, pf);
err = max([norm(X3 - Xr, inf), abs(T3 - Tr), norm(a3 - ar, inf), norm(th3 - thr, inf)]);
pass = report('free-theta pack/unpack', err, 0, pass);

za  = ch3_col_theta_augment(zr, pf);          % a fixed-theta vector, augmented
[~, ~, ~, th4] = ch3_col_unpack(za, pf);
err = norm(th4 - [pf.theta_minus; pf.theta_plus], inf) + ...
      norm(za(1:numel(zr)) - zr, inf);
pass = report('fixed-theta z augments exactly', err, 0, pass);

pe  = ch3_col_effective_params(zf, pf);
err = norm([pe.theta_minus; pe.theta_plus] - thr, inf);
pass = report('effective params read theta', err, 0, pass);

%% 7. evaluation timing
t0 = tic; for k = 1:3, ch3_col_cost(z0 + 1e-9*k, p); end; t_eval = toc(t0)/3;
n_grad = 2 * n_var;
fprintf('        one cost+constraint eval = %.3f s -> one central-difference\n', t_eval);
fprintf('        gradient (%d vars) ~ %.0f s\n', n_var, n_grad * t_eval);

fprintf('--- ch3_test_collocation: %s ---\n\n', tf(pass));
end

function ok = report(name, err, tol, ok_in)
ok = ok_in && (err <= tol);
fprintf('  [%s] %-30s err = %.3e (tol %.0e)\n', tf(err <= tol), name, err, tol);
end

function [orders, dmax] = defect_orders(p0, label)
%DEFECT_ORDERS  Max Hermite-Simpson defect of the seed under mesh refinement.
fprintf('        defect convergence, %s:\n', label);
Ns = [9 13 21 33];
orders = []; dmax = nan(size(Ns));
prev_d = NaN; prev_h = NaN;
for i = 1:numel(Ns)
    pn = p0;
    pn.N_nodes = Ns(i);
    zn = ch3_col_seed(pn);
    En = ch3_col_eval(zn, pn);
    d  = max(abs(En.defect(:)));
    dmax(i) = d;
    if ~isnan(prev_d)
        ord = log(prev_d/d) / log(prev_h/En.h);
        orders(end+1) = ord; %#ok<AGROW>
        fprintf('           N = %2d  h = %.4f  max|defect| = %.3e   order = %.2f\n', ...
                Ns(i), En.h, d, ord);
    else
        fprintf('           N = %2d  h = %.4f  max|defect| = %.3e\n', Ns(i), En.h, d);
    end
    prev_d = d; prev_h = En.h;
end
end

function s = tf(b)
if b, s = 'PASS'; else, s = 'FAIL'; end
end
