function pass = ch5_test_ecbf()
%CH5_TEST_ECBF  Section 5.2's theory, checked rather than quoted.
%
% Every claim here is one the chapter states as a theorem, corollary or remark,
% and each is checked in the form that would actually catch it being wrong:
%
%   1. pole placement: eig(Fb - Gb Kb) == -poles  (not merely "Hurwitz")
%   2. Remark 5.6: the ECBF row IS y_rb(x) >= 0 for the family (5.28)
%   3. (5.30) built recursively == (5.28) built from the polynomial
%   4. Theorem 5.2's total-negative condition is enforced, not just documented
%   5. Corollary 5.2's pole rule agrees with y_i(x0) >= 0 directly
%   6. Proposition 5.1: forward invariance propagates DOWN the chain of sets
%   7. THE GUARANTEE: h stays >= 0 along closed-loop runs on both plants
%   8. the sampled-data excursion is O(dt) -- i.e. discretization, not the row
%   9. the exponential envelope of Definition 5.1 actually lower-bounds h
%
% Check 1 is not pedantry. Kb = fliplr(poly(-p)) and Kb = poly(-p) both give a
% Hurwitz Ab; only one gives the requested poles, and the wrong one produces a
% controller that works and enforces the wrong exponential rate.
%
% See also CH5_ECBF_GAIN, CH5_ECBF_ADMISSIBLE, CH5_CTRL_ECBF_CLF_QP.

fprintf('\n=== ch5_test_ecbf ===\n');
pass = true;

%% ---------------------------------------------------- 1. pole placement
e_pole = 0;
for rb = 1:6
    poles = 0.7 + 0.4*(1:rb);
    e = ch5_ecbf_gain(struct('ecbf', struct('poles', poles)), rb);
    got = sort(real(e.eig_Ab));
    e_pole = max(e_pole, norm(got - sort(-poles), inf) / max(poles));
    % and Cb picks out h from eta_b, eq (5.11)
    e_pole = max(e_pole, norm(e.Cb - [1 zeros(1,rb-1)], inf));
end
pass = rep('eig(Fb - Gb Kb) == -poles', e_pole, 1e-9, pass);

%% -------------------------------------- 2. Remark 5.6, on both plants
% The row the QP enforces is  h^(rb) + Kb eta_b >= 0. Remark 5.6 says that
% expression IS y_rb of (5.28). Built here from the polynomial coefficients,
% which is a different route from the pole-placement gain the QP uses.
rng(21);
for pc = {ch5_params(), ch5_params('system','pendulum')}
    p  = pc{1};
    b0 = ch5_barrier(ch5_x0(p), p);
    e  = ch5_ecbf_gain(p, b0.rb);
    a  = poly(-e.poles);                 % [1, a1, ..., a_rb]

    err = 0;
    for j = 1:5
        xs = ch5_x0(p) + 0.3*randn(p.sys.nx, 1);
        b  = ch5_barrier(xs, p);
        hd = [b.eta_b(:); b.Lfrb];       % h^(0..rb), drift value on top

        % y_rb = sum_j a(j+1) h^(rb-j)
        y_poly = 0;
        for j = 0:b.rb
            y_poly = y_poly + a(j+1) * hd(b.rb - j + 1);
        end

        % what the controller computes: h^(rb) + Kb eta_b, at u = 0
        y_ctrl = b.Lfrb + e.Kb * b.eta_b;

        err = max(err, abs(y_poly - y_ctrl) / max(1, abs(y_poly)));
    end
    pass = rep(sprintf('Remark 5.6 on %s', p.system), err, 1e-11, pass);
end

%% --------------------------------- 3. recursion (5.30) == polynomial (5.28)
% y_i = ydot_(i-1) + p_i y_(i-1), applied to a random derivative stack.
rng(22);
rb = 5;
poles = [0.8 1.1 1.4 1.9 2.6];
err = 0;
for j = 1:8
    hd = randn(rb+2, 1);                 % h^(0..rb+1), arbitrary

    % recursive form: y_i and its derivative, carried together
    yc = [1; zeros(rb, 1)];              % coefficients of y_0 = h
    for i = 1:rb
        % ydot_(i-1): shift the coefficient vector up one derivative
        ycd = [0; yc(1:end-1)];
        yc  = ycd + poles(i)*yc;
    end
    y_rec = sum(yc .* hd(1:rb+1));

    % polynomial form
    a = poly(-poles);
    y_pol = 0;
    for i = 0:rb
        y_pol = y_pol + a(i+1) * hd(rb - i + 1);
    end
    err = max(err, abs(y_rec - y_pol) / max(1, abs(y_pol)));
end
pass = rep('(5.30) recursion == (5.28)', err, 1e-11, pass);

%% ------------------------------- 4. Theorem 5.2's total-negative condition
% A complex pole pair gives a Hurwitz Ab with the right characteristic
% polynomial and breaks Proposition 5.1's induction. It must be REFUSED.
pbad = struct('ecbf', struct('poles', [1+1i, 1-1i, 2, 3]));
threw = false; msg = '';
try
    ch5_ecbf_gain(pbad, 4);
catch err
    threw = true; msg = err.identifier;
end
ok = threw;
fprintf('  [%s] %-34s %s\n', tf(ok), 'complex poles are refused', msg);
pass = pass && ok;

% ... and a genuinely real, positive set is accepted
ok = true;
try
    eg = ch5_ecbf_gain(struct('ecbf', struct('poles', [5 5.5 6 10])), 4);
    ok = eg.total_negative && eg.hurwitz;
catch
    ok = false;
end
fprintf('  [%s] %-34s\n', tf(ok), 'real positive poles accepted');
pass = pass && ok;

%% ------------------------------ 5. Corollary 5.2 agrees with y_i(x0) >= 0
% The corollary's pole rule p_i >= -ydot_(i-1)/y_(i-1) is derived from
% y_i(x0) >= 0. Compute both independently and require them to say the same
% thing on states drawn to make BOTH answers occur.
rng(23);
p = ch5_params();
n_agree = 0; n_tot = 0; disagree = 0;
for j = 1:40
    xs = ch5_x0(p) + [1.2*randn(3,1); 1.2*randn(3,1)];
    adm = ch5_ecbf_admissible(xs, p);

    % the exactly-checkable rungs
    by_sets  = all(adm.y(1:adm.rb) >= 0);
    by_poles = all(adm.slack(1:adm.rb-1) >= -1e-9) && adm.y(1) >= 0;

    n_tot = n_tot + 1;
    if by_sets == by_poles, n_agree = n_agree + 1; else, disagree = disagree + 1; end
end
ok = (disagree == 0) && (n_agree == n_tot);
fprintf('  [%s] %-34s %d/%d states agree\n', tf(ok), ...
        'Cor 5.2 == y_i(x0) >= 0', n_agree, n_tot);
pass = pass && ok;

% and the default x0 is admissible on both plants, as the studies assume
for pc = {ch5_params(), ch5_params('system','pendulum')}
    a = ch5_ecbf_admissible(ch5_x0(pc{1}), pc{1});
    fprintf('  [%s] %-34s %s\n', tf(a.ok), ...
            sprintf('%s x0 admissible', pc{1}.system), a.msg);
    pass = pass && a.ok;
end

%% --------------------------- 5b. Remark 5.7: the relative-degree-2 special case
% Remark 5.7 relates ECBFs to the "modified CBF" of [82]: for a position
% constraint g(x) >= 0 of relative degree 2, enforcing the STANDARD relative-
% degree-1 CBF condition on
%
%       h(x) = (d/dt + gamma_b) o g(x)
%
% yields g(x) >= 0, so g is a relative-degree-2 ECBF. Checked as an algebraic
% identity on the derivative stack: the rb = 2 ECBF row with poles
% [gamma_b, p2] must equal the rb = 1 zeroing condition  hdot + p2 h >= 0
% applied to h = gdot + gamma_b g.
%
% Worth checking because it is the bridge between this chapter and the prior
% literature: if the two disagreed, one of the two constructions would be
% enforcing a different set than advertised.
rng(25);
err = 0;
for j = 1:8
    gb = 0.5 + 2*rand();          % gamma_b
    p2 = 0.5 + 2*rand();
    gd = randn(3,1);              % [g, gdot, gddot]

    % route A: the rb = 2 ECBF row, straight from ch5_pole_gain
    Kb2 = ch5_pole_gain([gb p2], 2);
    yA  = gd(3) + Kb2 * gd(1:2);

    % route B: h = gdot + gb g, treated as a relative-degree-1 barrier
    h    = gd(2) + gb*gd(1);
    hdot = gd(3) + gb*gd(2);
    yB   = hdot + p2*h;

    err = max(err, abs(yA - yB) / max(1, abs(yB)));
end
pass = rep('Remark 5.7: rb=2 == modified CBF', err, 1e-12, pass);

%% ------------------------------------------- 6. Proposition 5.1, numerically
% If C_i is forward invariant then so is C_(i-1). Simulate the LINEAR chain
% etadot_b = Fb eta_b + Gb mu_b directly with mu_b = -Kb eta_b (the equality
% case, which makes C_rb invariant by construction) and check that every lower
% set stays entered once entered.
rb = 4;
poles = [1.5 2 2.5 3];
eg = ch5_ecbf_gain(struct('ecbf', struct('poles', poles)), rb);
rng(24);
worst_drop = 0;
for j = 1:20
    % start strictly inside every C_i by construction: pick y_i > 0 and invert
    yv = 0.5 + rand(rb+1, 1);
    eb = y_to_eta(yv, poles);
    for i = 1:4000
        eb = eb + 1e-3 * (eg.Ab * eb);
        yv2 = eta_to_y(eb, poles);
        worst_drop = min(worst_drop, min(yv2(1:rb)));
    end
end
ok = worst_drop >= -1e-9;
fprintf('  [%s] %-34s min y_i over all runs = %+.3e\n', tf(ok), ...
        'Prop 5.1: sets stay entered', worst_drop);
pass = pass && ok;

%% ------------------------------------ 7-9. the guarantee, closed loop
fprintf('  --- closed loop ---\n');

% springmass: cheap enough to run in the suite
ps = ch5_params('controller','ecbfclfqp','constraint',3.00, ...
                'integrator','rk4','T',30);
ss = ch5_simulate(ps);
ok = ss.h_min >= 0;
fprintf('  [%s] %-34s h_min = %+.3e\n', tf(ok), ...
        'springmass h >= 0 (rb = 6)', ss.h_min);
pass = pass && ok;

yr = ss.y_rb(~isnan(ss.y_rb));
sc = max(1, abs(ss.Kb_eta(~isnan(ss.y_rb))));
ok = min(yr ./ sc) >= -1e-6;
fprintf('  [%s] %-34s min y_rb/scale = %+.3e\n', tf(ok), ...
        'springmass y_rb >= 0', min(yr ./ sc));
pass = pass && ok;

% the exponential envelope of Definition 5.1: h(t) >= Cb exp(Ab t) eta_b(0)
eg = ch5_ecbf_gain(ps, 6);
env = zeros(1, ss.n);
for i = 1:ss.n
    env(i) = eg.Cb * expm(eg.Ab * ss.t(i)) * ss.eta_b(:,1);
end
viol = min(ss.h(1:ss.n) - env);
ok = viol >= -1e-6 * max(abs(env));
fprintf('  [%s] %-34s min(h - envelope) = %+.3e\n', tf(ok), ...
        'Def 5.1 envelope lower-bounds h', viol);
pass = pass && ok;

%% ---------------- WHAT IS REPRODUCIBLE ACROSS INTEGRATORS, AND WHAT IS NOT
% Worth stating plainly, because it decides what any single plot from this
% chapter is evidence of.
%
% The springmass closed loop is reproducible to machine precision: rk4 and
% ode45 agree to ~1e-14 over the whole run.
%
% THE PENDULUM CLOSED LOOP IS NOT. Once the barrier row is active the
% controller holds the state on the constraint boundary with a zero-order
% hold, which is a discrete-time sliding mode: the state leaves the boundary
% between samples and is pushed back at the next one. Which side of an
% active-set switch a given sample lands on is decided far below the
% integration error, so two integrators that agree to 1e-10 at t = 0 separate
% steadily -- measured, ~0.5 rad after 3 s.
%
% So an individual pendulum TRAJECTORY is not a reproducible object, and no
% claim in this chapter should rest on one. What IS reproducible is the thing
% the chapter actually asserts: the safe set is never left, by either
% integrator, by more than the sampled-data O(dt) margin. That is what gets
% asserted here, and the trajectory divergence is printed next to it rather
% than hidden.
fprintf('  --- integrator dependence ---\n');
for pc = {ch5_params('controller','ecbfclfqp','constraint',3.00,'T',6), ...
          ch5_params('system','pendulum','controller','ecbfclfqp', ...
                     'constraint',-1.0,'T',3)}
    pa = pc{1}; pa.integrator = 'rk4';
    pb = pc{1}; pb.integrator = 'ode45';
    sa = ch5_simulate(pa);
    sb = ch5_simulate(pb);

    div = max(abs(sa.x(:) - sb.x(:))) / max(1, max(abs(sb.x(:))));
    tol = 10 * pa.control_dt;
    ok  = (sa.h_min >= -tol) && (sb.h_min >= -tol);

    fprintf('  [%s] %-34s h_min rk4 %+.2e / ode45 %+.2e  (traj diverges %.1e)\n', ...
            tf(ok), sprintf('%s: SAFETY is integrator-free', pa.system), ...
            sa.h_min, sb.h_min, div);
    pass = pass && ok;
end

% pendulum: short window, coarse-but-honest, and the O(dt) scaling
fprintf('  --- sampled-data excursion is O(dt) ---\n');
hmins = zeros(1,3);
dts   = [2e-3 1e-3 5e-4];
for i = 1:3
    qq = ch5_params('system','pendulum','controller','ecbfclfqp', ...
                    'constraint',-1.0,'integrator','rk4','T',7, ...
                    'control_dt',dts(i));
    sq = ch5_simulate(qq);
    hmins(i) = sq.h_min;
end
fprintf('        dt      = %s\n', mat2str(dts));
fprintf('        h_min   = %s\n', mat2str(round(hmins, 8)));

% Either strictly safe, or shrinking as dt does. The claim under test is that
% any excursion is the DISCRETIZATION rather than the barrier row failing.
neg = min(hmins, 0);
ok  = all(hmins >= 0) || (abs(neg(3)) <= abs(neg(1)) + 1e-12);
fprintf('  [%s] %-34s worst excursion %.2e -> %.2e\n', tf(ok), ...
        'excursion shrinks with dt', neg(1), neg(3));
pass = pass && ok;

fprintf('--- ch5_test_ecbf: %s ---\n', tf(pass));

end

% ---------------------------------------------------------------------------
function eb = y_to_eta(yv, poles)
%Y_TO_ETA  Build eta_b = [h; hdot; ...] from the desired y_0..y_(rb-1).
%
% y_i = sum_j c_i(j+1) h^(i-j) is lower triangular in the derivatives, so this
% is a forward substitution -- which is also the constructive proof that a
% point strictly inside every C_i exists.
rb = numel(poles);
eb = zeros(rb, 1);
eb(1) = yv(1);
for i = 1:rb-1
    c = poly(-poles(1:i));
    acc = 0;
    for j = 1:i
        acc = acc + c(j+1) * eb(i - j + 1);
    end
    eb(i+1) = yv(i+1) - acc;
end
end

function yv = eta_to_y(eb, poles)
rb = numel(poles);
yv = zeros(rb, 1);
for i = 0:rb-1
    c = poly(-poles(1:i));
    acc = 0;
    for j = 0:i
        acc = acc + c(j+1) * eb(i - j + 1);
    end
    yv(i+1) = acc;
end
end

function ok = rep(name, err, tol, ok_in)
ok = ok_in && (err <= tol);
fprintf('  [%s] %-34s err = %.3e (tol %.0e)\n', tf(err <= tol), name, err, tol);
end

function s = tf(b)
if b, s = 'PASS'; else, s = 'FAIL'; end
end
