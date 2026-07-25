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
%   5. Hermite-Simpson defects converge at the expected 4th order in h.
%   6. One cost + constraint evaluation is timed, so the cost of a full
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
% Refining the node count must drop the defect like h^4. This is the sharpest
% single check that the transcription formula is right: a mis-placed h/8 or
% h/6 still gives small defects but degrades the observed order.
fprintf('        defect convergence (expect ~4th order):\n');
prev_d = NaN; prev_h = NaN; orders = [];
for Nn = [9 13 21 33]
    pn = ch3_params('N_nodes', Nn);
    zn = ch3_col_seed(pn);
    En = ch3_col_eval(zn, pn);
    d  = max(abs(En.defect(:)));
    if ~isnan(prev_d)
        ord = log(prev_d/d) / log(prev_h/En.h);
        orders(end+1) = ord; %#ok<AGROW>
        fprintf('           N = %2d  h = %.4f  max|defect| = %.3e   order = %.2f\n', ...
                Nn, En.h, d, ord);
    else
        fprintf('           N = %2d  h = %.4f  max|defect| = %.3e\n', Nn, En.h, d);
    end
    prev_d = d; prev_h = En.h;
end
ok = ~isempty(orders) && median(orders) > 3.0;
fprintf('  [%s] %-30s median order = %.2f (expect > 3)\n', tf(ok), ...
        'Hermite-Simpson is 4th order', median(orders));
pass = pass && ok;

%% 6. evaluation timing
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

function s = tf(b)
if b, s = 'PASS'; else, s = 'FAIL'; end
end
