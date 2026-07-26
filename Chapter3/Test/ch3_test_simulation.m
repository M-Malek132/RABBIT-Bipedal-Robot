function ch3_test_simulation()
%CH3_TEST_SIMULATION  Verify the sampled-data (zero-order-hold) integrator.
%
%   ch3_test_simulation
%
% p.control_dt turns ch3_step from "let ode45 call the controller whenever it
% likes" into "solve for u once per control period and hold it". That is the
% mode the constrained CLF-QP has to run in -- its u(x) kinks whenever the QP's
% active set changes, and an adaptive explicit solver responds to a kink by
% shrinking h without bound (measured: 51908 RHS evaluations to cover 0.5% of a
% step). So this path carries real weight and needs its own check.
%
% The independent reference is the CONTINUOUS rollout of the SAME controller,
% using a controller for which both modes are tractable (iolin_pd). A correct
% zero-order hold is a first-order approximation to continuous feedback, so the
% test is not "the two agree" -- they must not, exactly -- but rather:
%
%   1. sampling error DECREASES with the sample period, at first order:
%      halving control_dt must roughly halve the deviation. This is what
%      distinguishes a correct hold from one that is merely close by luck, and
%      it would catch an off-by-one in the period loop or a stale-state bug
%      that a single-tolerance comparison would pass.
%   2. the guard still fires and the step still completes;
%   3. the constrained QP, which cannot be run continuously at all, completes a
%      step under sampling and respects its torque box exactly.
%
% See also CH3_STEP, CH3_PARAMS, CH3_COMPARE_CONTROLLERS.

fprintf('=== ch3_test_simulation ===\n');
ok = true;

p = ch3_params();
[alpha, x0] = ch3_seed(p);

%% ---- 1. ZOH converges to the continuous rollout, at first order ---------
pc = p; pc.controller = 'iolin_pd'; pc.control_dt = 0;
ref = ch3_step(x0, alpha, pc);

ok = report('continuous rollout reaches guard', double(~ref.ok), 0, ok);

dts  = [2e-3 1e-3 5e-4];
devs = zeros(size(dts));
for i = 1:numel(dts)
    pc.control_dt = dts(i);
    s = ch3_step(x0, alpha, pc);
    devs(i) = norm(s.x_end - ref.x_end);
    if ~s.ok, ok = report(sprintf('ZOH %g Hz reaches guard', 1/dts(i)), 1, 0, ok); end
end

% First-order hold => halving dt halves the error. Allow a generous band: the
% ratio is contaminated by the guard firing at a slightly different instant, so
% require only that it is clearly decreasing and not accidentally superlinear.
ratios = devs(1:end-1) ./ max(devs(2:end), realmin);
fprintf('        dev vs continuous: %s\n', sprintf('%.3e ', devs));
fprintf('        halving ratios   : %s (expect ~2)\n', sprintf('%.2f ', ratios));

ok = report('ZOH error shrinks with dt', double(any(diff(devs) >= 0)), 0, ok);
ok = report('ZOH is first order', max(abs(ratios - 2)), 1.0, ok);

%% ---- 2. sampling must not change the step qualitatively -----------------
pc.control_dt = 1e-3;
s = ch3_step(x0, alpha, pc);
ok = report('ZOH step duration matches',  abs(s.T - ref.T),           5e-3, ok);
ok = report('ZOH step length matches',    abs(s.L_step - ref.L_step), 5e-3, ok);

%% ---- 3. the constrained QP runs, and honours its box --------------------
% Box deliberately BELOW what the gait needs, so it binds and the slack has to
% work. This is the configuration that stalls outright without sampling.
u_ff_peak = max(abs(ch3_control(x0, alpha, pc)));
pq = p;
pq.controller = 'clfqp_con';
pq.control_dt = 1e-3;
pq.limits.u_max = max(0.6 * u_ff_peak, 5);
pq.limits.enable.torque = true;

t0 = tic;
sq = ch3_step(x0, alpha, pq);
el = toc(t0);

ok = report('clfqp_con completes a step', double(~sq.ok), 0, ok);
F  = ch3_forces(sq.t, sq.x, alpha, pq);
% Allow a hair of numerical slop on the box: quadprog satisfies bounds to its
% own feasibility tolerance, not to machine precision.
ok = report('clfqp_con respects torque box', ...
            max(F.torque_max - pq.limits.u_max, 0), 1e-6, ok);
fprintf('        box %.2f Nm, peak |u| %.2f Nm, max delta %.3e, %.1f s wall\n', ...
        pq.limits.u_max, F.torque_max, F.delta_max, el);

%% ------------------------------------------------------------------------
if ok
    fprintf('--- ch3_test_simulation: PASS ---\n');
else
    error('ch3_test_simulation:FAIL', 'ch3_test_simulation had failures');
end
end

% ---------------------------------------------------------------------------
function ok = report(name, err, tol, ok_in)
ok = ok_in && (err <= tol);
fprintf('  [%s] %-30s err = %.3e (tol %.0e)\n', tf(err <= tol), name, err, tol);
end

function s = tf(b)
if b, s = 'PASS'; else, s = 'FAIL'; end
end
