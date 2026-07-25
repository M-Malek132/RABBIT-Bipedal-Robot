function ch3_test_vc()
%CH3_TEST_VC  Stage-2 and stage-4 verification: basis, outputs, Lie derivatives.
%
% Checks:
%   1. Bezier endpoint interpolation: yd(0)=alpha_0, yd(1)=alpha_M.
%   2. Bezier dyd/ds and d2yd/ds2 against central differences.
%   3. Bezier == clamped B-spline when degree+1 control points are used
%      (they are the same curve; this cross-validates the two backends).
%   4. ydot returned by ch3_outputs equals d/dt y along the true flow.
%   5. ydd == Lf^2 y + LgLf y u along the true flow, for random u
%      (this is the relative-degree-2 claim, verified numerically).
%   6. u = u_ff makes ydd vanish, i.e. the feedforward really does cancel.
%   7. The decoupling matrix is well conditioned at the test poses.

fprintf('\n=== ch3_test_vc ===\n');
p = ch3_params();
pass = true;

rng(1);
alpha = 0.3*randn(p.ny, p.n_ctrl);

%% 1. endpoint interpolation
y0 = ch3_yd(alpha, 0, p);
y1 = ch3_yd(alpha, 1, p);
pass = report('bezier yd(0) = alpha_0', norm(y0 - alpha(:,1),   inf), 1e-14, pass);
pass = report('bezier yd(1) = alpha_M', norm(y1 - alpha(:,end), inf), 1e-14, pass);

%% 2. analytic vs finite-differenced s-derivatives
h = 1e-6; e1 = 0; e2 = 0;
for s = 0.1:0.1:0.9
    [~, dyd, d2yd] = ch3_bezier(alpha, s);
    [ym, dym] = ch3_bezier(alpha, s-h);
    [yp, dyp] = ch3_bezier(alpha, s+h);
    e1 = max(e1, norm(dyd  - (yp - ym)/(2*h),   inf));
    e2 = max(e2, norm(d2yd - (dyp - dym)/(2*h), inf));
end
pass = report('bezier dyd/ds vs FD',    e1, 1e-6, pass);
pass = report('bezier d2yd/ds2 vs FD',  e2, 1e-5, pass);

%% 3. Bezier == clamped B-spline of matching degree
% A clamped B-spline of degree M with exactly M+1 control points IS the
% degree-M Bezier curve, so the two backends must describe the same curve.
% We assert that on the VALUES, which is the property that actually defines
% the curves.
pb = p;                     pb.basis = 'bezier';
ps = ch3_params('basis','bspline','bsp_deg',p.bez_deg);
eb = 0; ebd = 0;
for s = 0:0.05:1
    ya = ch3_yd(alpha, s, pb);
    yc = ch3_yd(alpha, s, ps);
    eb = max(eb, norm(ya - yc, inf));

    % Bezier's analytic slope must match a finite difference of the B-SPLINE
    % value -- i.e. of the same curve computed by the other backend. This is
    % the cross-backend derivative check, and it deliberately does NOT call
    % BSpline_derivative (see the note below).
    h  = 1e-6;
    sm = min(max(s-h,0),1); sp = min(max(s+h,0),1);
    [~, da] = ch3_bezier(alpha, s);
    dnum = (ch3_yd(alpha, sp, ps) - ch3_yd(alpha, sm, ps)) / (sp - sm);
    if sp > sm && s > 0.02 && s < 0.98
        ebd = max(ebd, norm(da - dnum, inf));
    end
end
pass = report('bezier == bspline (value)', eb,  1e-10, pass);
pass = report('bezier slope vs bspline FD', ebd, 1e-5, pass);

%% 3b. BSpline_derivative self-consistency, at the degree the repo uses
% KNOWN LIMITATION of the inherited Trajectory_Optimization/BSpline_derivative.m:
% its analytic recurrence agrees with a finite difference of its own curve at
% degree 3 (the degree the existing pipeline runs, p.bsp_deg = 3) but NOT at
% degree 5, where the clamped knot vector collapses to the Bezier one and the
% interior recurrence mis-indexes -- the error grows to ~0.4 near s = 1.
% That is why 'bezier' is the default basis here: its derivatives are analytic
% and degree-independent. We assert only the degree-3 behaviour, which is what
% the inherited code is actually relied upon for.
a1 = alpha(1,:); n = numel(a1)-1;
e3 = 0;
for s = 0.05:0.05:0.95
    d_ana = a1 * BSpline_derivative(n, 3, s).';
    h = 1e-7;
    d_fd  = (a1*BSpline(n,3,s+h).' - a1*BSpline(n,3,s-h).')/(2*h);
    e3 = max(e3, abs(d_ana - d_fd));
end
pass = report('BSpline_derivative @ deg 3', e3, 1e-6, pass);

%% test states: pick poses at mid-phase so the phase clamp is never active
x0 = [ +0.0000; -0.9310; +0.2999; -0.7934; +0.6869; -0.6318; +0.9810; ...
       +0.3952; -0.0419; +0.1847; +0.1847; +0.1045; +0.0000; +0.0000];
xa = x0;      xa(3) = xa(3) + 0.225;      % theta -> midway, s ~ 0.5
xb = x0;      xb(3) = xb(3) + 0.115;      % s ~ 0.25
xb(8:14) = xb(8:14) + 0.2*randn(7,1);
states = [xa, xb];

%% 4-7. Lie derivatives along the true flow
e_ydot = 0; e_ydd = 0; e_ff = 0; worst_rc = Inf;
hh = 1e-6;
for k = 1:size(states,2)
    x = states(:,k);
    [Lf2y, LgLfy, u_ff, info] = ch3_io_lin(x, alpha, p);
    worst_rc = min(worst_rc, info.rcond);

    for trial = 1:3
        u = 30*randn(p.nu,1);
        [f, g] = ch3_control_affine(x, p);
        xdot = f + g*u;

        % ydot must equal the directional derivative of y along xdot
        yp = ch3_outputs(x + hh*xdot, alpha, p);
        ym = ch3_outputs(x - hh*xdot, alpha, p);
        e_ydot = max(e_ydot, norm(info.ydot - (yp - ym)/(2*hh), inf));

        % ydd must equal the directional derivative of ydot along xdot
        [~, ydp] = ch3_outputs(x + hh*xdot, alpha, p);
        [~, ydm] = ch3_outputs(x - hh*xdot, alpha, p);
        ydd_num = (ydp - ydm)/(2*hh);
        ydd_ana = Lf2y + LgLfy*u;
        e_ydd = max(e_ydd, norm(ydd_ana - ydd_num, inf));
    end

    % feedforward: u = u_ff should give ydd = 0
    e_ff = max(e_ff, norm(Lf2y + LgLfy*u_ff, inf));
end
pass = report('ydot vs flow derivative',   e_ydot, 1e-5, pass);
pass = report('ydd = Lf2y + LgLfy*u',      e_ydd,  1e-4, pass);
pass = report('u_ff cancels: ydd = 0',     e_ff,   1e-9, pass);
fprintf('  [%s] %-30s rcond = %.3e (min 1e-6)\n', tf(worst_rc > 1e-6), ...
        'decoupling well conditioned', worst_rc);
pass = pass && worst_rc > 1e-6;

fprintf('--- ch3_test_vc: %s ---\n\n', tf(pass));
end

function ok = report(name, err, tol, ok_in)
ok = ok_in && (err <= tol);
fprintf('  [%s] %-30s err = %.3e (tol %.0e)\n', tf(err <= tol), name, err, tol);
end

function s = tf(b)
if b, s = 'PASS'; else, s = 'FAIL'; end
end
