function pass = ch5_test_qp()
%CH5_TEST_QP  The quadratic programs: equivalences, and the motivating failure.
%
% Checks:
%   1. Remark 5.4: the VIOL program (5.15) == the direct program (5.7)
%   2. Remark 5.4 again, for the ECBF pair (5.31) with and without mu_b
%   3. an inactive barrier row reduces (5.31) to the plain CLF-QP closed form
%   4. Definition 5.2 / Remark 5.3: at rb = 1 the ECBF condition IS the ZCBF
%      condition hdot >= -kb h, and it enforces the same set
%   5. THE MOTIVATION: at rb > 1 the Section 5.1 barrier row contains no
%      control, on both plants
%   6. the barrier row is never slacked, and the CLF row always is
%   7. the input box of (5.8)/(5.31) is respected when switched on
%   8. the QPs return finite, feasible answers over a scatter of states
%
% CHECK 5 IS THE POINT OF THE CHAPTER. Everything Section 5.2 builds is
% justified by the claim that nominal CBFs "are difficult to be extended to
% high relative-degree systems". That claim is checkable and it is checked
% here, in the strongest available form: the coefficient of u in the Section
% 5.1 barrier row is EXACTLY zero on both validation plants, so the row cannot
% be satisfied OR violated by any choice of control.
%
% See also CH5_CTRL_CBF_CLF_QP, CH5_CTRL_ECBF_CLF_QP, CH5_CTRL_CLF_QP.

fprintf('\n=== ch5_test_qp ===\n');
pass = true;

rng(31);
p  = ch5_params();
q  = ch5_params('system','pendulum');
XS = ch5_x0(p) + [0.4*randn(3,6); 0.4*randn(3,6)];
ZS = ch5_x0(q) + [0.4*randn(4,6); 0.3*randn(4,6)];

%% ------------------------------------------- 1. Remark 5.4 for Section 5.1
% Only meaningful where the reciprocal barrier has a control in it, i.e.
% relative degree 1. This is why 'v1_max' exists.
%
% The level is 2.0 rather than 0.8 so that every sampled state has h > 0. That
% is a real restriction, not test hygiene: B = 1/h is a barrier only on
% Int(C) = {h > 0}, so Remark 5.4 has nothing to compare outside it. The
% behaviour outside is checked separately below.
pv1 = ch5_params('constraint', struct('type','v1_max','value',2.0), ...
                 'ecbf.poles', 2.0);
e_r54 = 0;
for k = 1:size(XS,2)
    pa = pv1; pa.controller = 'cbfclfqp';
    pb = pv1; pb.controller = 'cbfclfqp_viol';
    ua = ch5_control(XS(:,k), pa);
    ub = ch5_control(XS(:,k), pb);
    e_r54 = max(e_r54, norm(ua-ub, inf) / max(1, norm(ua, inf)));
end
pass = rep('Remark 5.4: (5.15) == (5.7)', e_r54, 1e-7, pass);

% ... and outside the set the reciprocal barrier is UNDEFINED, in both forms.
% An ECBF at the same state is not: its condition is a polynomial in the
% derivatives of h and survives h < 0 perfectly well, which is the difference
% that matters if the state is ever perturbed out of C.
x_out = ch5_x0(pv1); x_out(4) = 3.0;            % xdot1 = 3 > 2, so h < 0
b_out = ch5_barrier(x_out, pv1);

% Both Section 5.1 forms must report the barrier as undefined, return a finite
% control, and stop constraining anything -- i.e. reduce to the CLF-QP.
n_undef = 0;
mu_clf  = ch5_ctrl_clf_qp(ch5_io_lin(x_out, pv1), pv1);
for c = {'cbfclfqp','cbfclfqp_viol'}
    pc = pv1; pc.controller = c{1};
    [uu, in_] = ch5_control(x_out, pc);
    % 1e-4, not 1e-12: ch5_ctrl_clf_qp takes the hard-constrained closed form
    % while these solve the slacked QP, and those differ by 1/slack_penalty by
    % construction. The claim being tested is "the barrier stopped
    % participating", not "the two code paths are bit-identical".
    same_as_clf = norm(in_.mu - mu_clf, inf) <= 1e-4 * max(1, norm(mu_clf, inf));
    n_undef = n_undef + (~in_.qp.barrier_defined && all(isfinite(uu)) && same_as_clf);
end

% The ECBF at the same state is defined, feasible, and still pushing back.
pe = pv1; pe.controller = 'ecbfclfqp';
[ue, ie] = ch5_control(x_out, pe);
ok = (n_undef == 2) && ie.qp.feasible && all(isfinite(ue)) && (b_out.h < 0);
fprintf('  [%s] %-34s h = %+.2f: 1/h undefined (%d/2), ECBF still acts (%d)\n', ...
        tf(ok), 'outside C: 5.1 undefined, 5.2 not', b_out.h, n_undef, ...
        ie.qp.feasible);
pass = pass && ok;

%% -------------------------------------- 2. Remark 5.4 for the ECBF pair
% SAMPLED ALONG AN ACTUAL ROLLOUT, not from a Gaussian cloud around x0, and
% the difference matters. Perturbing all eight pendulum states by 0.4 puts the
% joint velocities far outside anything the arm reaches, where L_f^4 h runs to
% 1e6 and the barrier row demands ||mu|| ~ 2e7. Both formulations still
% describe the same program there, but a QP whose solution is 1e7 with a 1e6
% slack penalty is at the edge of what double precision resolves, and the
% comparison stops measuring the formulations and starts measuring the solver.
%
% Remark 5.4 is an algebraic identity; the useful question is whether it
% survives in the states the controller actually visits.
for pc = {p, q}
    pp = pc{1};
    XX = rollout_states(pp, 6);
    e = 0;
    for k = 1:size(XX,2)
        pa = pp; pa.controller = 'ecbfclfqp';
        pb = pp; pb.controller = 'ecbfclfqp_viol';
        ua = ch5_control(XX(:,k), pa);
        ub = ch5_control(XX(:,k), pb);
        e = max(e, norm(ua-ub, inf) / max(1, norm(ua, inf)));
    end
    pass = rep(sprintf('ECBF direct == VIOL (%s)', pp.system), e, 1e-3, pass);
end

%% --------------------------- 3. an inactive barrier row leaves the CLF alone
% Push the constraint far away so the ECBF row cannot bind; (5.31) must then
% return the Chapter-3 min-norm closed form exactly.
p_far = ch5_params('constraint', 1e4, 'controller', 'ecbfclfqp');
e_far = 0; n_act = 0;
for k = 1:size(XS,2)
    io = ch5_io_lin(XS(:,k), p_far);
    b  = ch5_barrier(XS(:,k), p_far);
    e  = ch5_ecbf_gain(p_far, b.rb);
    mu_e = ch5_ctrl_ecbf_clf_qp(io, b, e, p_far, false);
    mu_c = ch5_ctrl_clf_qp(io, p_far);
    e_far = max(e_far, norm(mu_e - mu_c, inf) / max(1, norm(mu_c, inf)));
    [~,~,qp] = ch5_ctrl_ecbf_clf_qp(io, b, e, p_far, false);
    n_act = n_act + qp.cbf_active;
end
pass = rep('slack barrier -> CLF-QP exactly', e_far, 1e-5, pass);
ok = (n_act == 0);
fprintf('  [%s] %-34s active on %d/%d states\n', tf(ok), ...
        'and the row never binds there', n_act, size(XS,2));
pass = pass && ok;

%% ------------------------- 4. Definition 5.2 / Remark 5.3: rb = 1 is a ZCBF
% With rb = 1 the ECBF row reads h^(1) >= -kb h, which is the zeroing-CBF
% condition verbatim. Check the row the QP builds against that expression, and
% then check that running it enforces the set.
p1 = ch5_params('constraint', struct('type','v1_max','value',0.8), ...
                'ecbf.poles', 2.0, 'controller', 'ecbfclfqp');
kb = p1.ecbf.poles;
e_z = 0;
for k = 1:size(XS,2)
    io = ch5_io_lin(XS(:,k), p1);
    b  = ch5_barrier(XS(:,k), p1);
    e  = ch5_ecbf_gain(p1, b.rb);
    [~, u, qp] = ch5_ctrl_ecbf_clf_qp(io, b, e, p1, false);
    hdot = b.Lfh + b.Lgh * u;
    % the ZCBF residual and the reported y_rb must be the same number
    e_z = max(e_z, abs((hdot + kb*b.h) - qp.y_rb) / max(1, abs(qp.y_rb)));
end
pass = rep('rb=1 ECBF row == ZCBF row', e_z, 1e-9, pass);

p1.T = 20; p1.integrator = 'rk4';
s1 = ch5_simulate(p1);
p1b = p1; p1b.controller = 'clfqp';
s1b = ch5_simulate(p1b);
ok = (s1.h_min >= 0) && (s1b.h_min < s1.h_min);
fprintf('  [%s] %-34s ECBF h_min = %+.4f, baseline %+.4f\n', tf(ok), ...
        'rb=1 ECBF enforces the set', s1.h_min, s1b.h_min);
pass = pass && ok;

%% ============ 5. THE MOTIVATION: Section 5.1 has no row at high rel. degree
fprintf('  --- Section 5.1 at high relative degree ---\n');
for pc = {p, q}
    pp = pc{1};
    XX = ch5_x0(pp) + 0.3*randn(pp.sys.nx, 6);
    b0 = ch5_barrier(ch5_x0(pp), pp);

    max_Lgh = 0; n_unc = 0;
    for k = 1:size(XX,2)
        io = ch5_io_lin(XX(:,k), pp);
        b  = ch5_barrier(XX(:,k), pp);
        max_Lgh = max(max_Lgh, norm(b.Lgh, inf));
        [~,~,qp] = ch5_ctrl_cbf_clf_qp(io, b, pp, false);
        n_unc = n_unc + ~qp.barrier_controllable;
    end

    ok = (max_Lgh == 0) && (n_unc == size(XX,2)) && (b0.rb > 1);
    fprintf('  [%s] %-34s rb = %d, max|L_g h| = %g, uncontrollable %d/%d\n', ...
            tf(ok), sprintf('%s: no control in the row', pp.system), ...
            b0.rb, max_Lgh, n_unc, size(XX,2));
    pass = pass && ok;
end

% ... and the consequence: the reciprocal CBF lets the constraint go
pcb = ch5_params('controller','cbfclfqp','constraint',3.00, ...
                 'integrator','rk4','T',30);
scb = ch5_simulate(pcb);
pec = pcb; pec.controller = 'ecbfclfqp';
sec = ch5_simulate(pec);
ok = (scb.h_min < -1e-3) && (sec.h_min >= 0);
fprintf('  [%s] %-34s Sec 5.1 h_min = %+.4f, Sec 5.2 h_min = %+.4f\n', ...
        tf(ok), 'and so it violates, while 5.2 does not', ...
        scb.h_min, sec.h_min);
pass = pass && ok;

%% ------------------------------------ 6. what is slacked and what is not
% The barrier row must be satisfied exactly; the CLF row may be relaxed. Force
% a conflict with a tight input box and confirm which one gives.
pbox = ch5_params('controller','ecbfclfqp','constraint',3.00, 'u_max', 0.05);
n_dslack = 0; n_bviol = 0; n_ok = 0;
for k = 1:size(XS,2)
    io = ch5_io_lin(XS(:,k), pbox);
    b  = ch5_barrier(XS(:,k), pbox);
    e  = ch5_ecbf_gain(pbox, b.rb);
    [~, u, qp] = ch5_ctrl_ecbf_clf_qp(io, b, e, pbox, false);
    if ~qp.feasible, continue; end
    n_ok = n_ok + 1;
    n_dslack = n_dslack + (qp.delta > 1e-9);
    n_bviol  = n_bviol  + (qp.y_rb < -1e-6 * max(1, abs(qp.Kb_eta)));
    if max(abs(u)) > pbox.u_max + 1e-7, n_bviol = n_bviol + 100; end
end
ok = (n_bviol == 0);
fprintf('  [%s] %-34s CLF slacked on %d/%d, barrier violated %d\n', tf(ok), ...
        'barrier hard, CLF soft', n_dslack, n_ok, n_bviol);
pass = pass && ok;

%% ------------------------------------------------- 7. the input box holds
e_box = 0;
for pc = {ch5_params('controller','ecbfclfqp','u_max',0.5), ...
          ch5_params('system','pendulum','controller','ecbfclfqp','u_max',15)}
    pp = pc{1};
    XX = ch5_x0(pp) + 0.4*randn(pp.sys.nx, 6);
    for k = 1:size(XX,2)
        [u, in] = ch5_control(XX(:,k), pp);
        if in.qp.feasible
            e_box = max(e_box, max(abs(u)) - pp.u_max);
        end
    end
end
pass = rep('|u| <= u_max when boxed', max(e_box, 0), 1e-7, pass);

%% ---------------------------------------- 8. everything stays finite
n_bad = 0; n_all = 0;
for pc = {p, q}
    pp = pc{1};
    XX = ch5_x0(pp) + 0.5*randn(pp.sys.nx, 8);
    for c = {'clfqp','cbfclfqp','ecbfclfqp','ecbfclfqp_viol'}
        pp.controller = c{1};
        for k = 1:size(XX,2)
            u = ch5_control(XX(:,k), pp);
            n_all = n_all + 1;
            n_bad = n_bad + ~all(isfinite(u));
        end
    end
end
ok = (n_bad == 0);
fprintf('  [%s] %-34s %d/%d finite\n', tf(ok), ...
        'all controllers return finite u', n_all - n_bad, n_all);
pass = pass && ok;

fprintf('--- ch5_test_qp: %s ---\n', tf(pass));

end

% ---------------------------------------------------------------------------
function XX = rollout_states(p, n)
%ROLLOUT_STATES  n states spread along a short closed-loop run.
%
% Cheap (rk4, a few seconds of simulated time) and it produces exactly the
% states the chapter's claims are about: on the boundary, off it, barrier
% active and inactive.
pr = p;
pr.controller = 'ecbfclfqp';
pr.integrator = 'rk4';
pr.T = min(p.T, 4);
pr.ecbf.admissibility = 'off';

s  = ch5_simulate(pr);
ok = find(all(isfinite(s.x), 1));
idx = ok(round(linspace(1, numel(ok), n)));
XX  = s.x(:, idx);
end

% ---------------------------------------------------------------------------
function ok = rep(name, err, tol, ok_in)
ok = ok_in && (err <= tol);
fprintf('  [%s] %-34s err = %.3e (tol %.0e)\n', tf(err <= tol), name, err, tol);
end

function s = tf(b)
if b, s = 'PASS'; else, s = 'FAIL'; end
end
