function pass = ch6_test_sim()
%CH6_TEST_SIM  The stepping-stone simulation, against things it cannot fake.
%
% Checks:
%   1. with no barriers and the Chapter-3 controller, one step reproduces the
%      reference gait -- the harness does not perturb what it is measuring
%   2. l_s is the swing-foot position at impact, recomputed from x_end
%   3. the barrier actually MOVES the placement: on a window that excludes the
%      nominal step, the CBF run lands somewhere the baseline does not
%   4. a moving stone is scored against the window AT IMPACT, not at t = 0
%   5. the sampled-data hold is consistent -- halving control_dt changes the
%      trajectory by O(dt), not O(1)
%   6. a run that cannot complete reports why, rather than returning a short
%      result that looks like a success
%
% See also CH6_STEP, CH6_SIMULATE, CH6_TERRAIN.

fprintf('\n=== ch6_test_sim ===\n');
pass = true;

[x0, alpha] = reference_gait();
p = ch6_params();
p.limits.enable.torque   = true;  p.limits.u_max  = 300;
p.limits.enable.friction = true;  p.limits.mu_s   = 0.6;
p.limits.enable.grf      = true;  p.limits.Fz_min = 50;

%% ------------------------------------- 1. the baseline reproduces the gait
q = p;  q.controller = 'iolin_pd';  q.cbf.problem = 'none';
s0 = ch6_step(x0, alpha, q);
ok = s0.ok && abs(s0.T - 0.3009) < 5e-3 && abs(s0.l_s - 0.3533) < 5e-3;
fprintf('  [%s] %-42s T = %.4f s, l_s = %.4f m\n', tf(ok), ...
        'no-CBF baseline == reference gait', s0.T, s0.l_s);
pass = pass && ok;

%% ---------------------------------------- 2. l_s is what it says it is
q2 = P_sw(s0.x_end(1:p.nq)) - P_st(s0.x_end(1:p.nq));
ok = abs(q2(1) - s0.l_s) < 1e-12 && abs(q2(2) - s0.h_f_end) < 1e-12;
fprintf('  [%s] %-42s residual %.1e\n', tf(ok), ...
        'l_s recomputed from x_end', abs(q2(1) - s0.l_s));
pass = pass && ok;

%% --------------------------------- 3. the barrier moves the foot placement
% A window that the nominal step misses. The baseline must miss it (it has no
% mechanism not to); the CBF run must land somewhere different. This is the
% weakest form of "the barrier does something" that cannot be satisfied by
% accident.
win = [0.25 0.30];
qb = p;  qb.controller = 'iolin_pd';  qb.cbf.problem = 'none';
qc = p;  qc.controller = 'cbf_clf_qp';
qb.stone = ch6_resolve_stone(p.stones, win(1), win(2));
qc.stone = qb.stone;

sb = ch6_step(x0, alpha, qb);
sc = ch6_step(x0, alpha, qc);
moved = abs(sc.l_s - sb.l_s);
ok = ~sb.in_window && moved > 0.01 && (sc.l_s < sb.l_s);
fprintf('  [%s] %-42s baseline %.4f -> CBF %.4f (window [%.2f %.2f])\n', ...
        tf(ok), 'barrier shortens the step toward the window', ...
        sb.l_s, sc.l_s, win(1), win(2));
pass = pass && ok;

%% --------------------------- 4. a moving stone is scored at IMPACT
qm = p;
qm.stones.motion  = 'linear';
qm.stones.v_stone = 0.20;                    % 20 cm/s drift
qm.stone = ch6_resolve_stone(qm.stones, 0.30, 0.40);
sm = ch6_step(x0, alpha, qm);
expect = 0.30 + 0.20*sm.T;
ok = abs(sm.stone_end.l_min - expect) < 1e-9 && ...
     abs(sm.stone_end.l_min - 0.30) > 1e-3;
fprintf('  [%s] %-42s l_min: %.4f at t=0 -> %.4f at T=%.3f\n', tf(ok), ...
        'moving window read at the impact time', 0.30, sm.stone_end.l_min, sm.T);
pass = pass && ok;

%% ------------------------------------- 5. the zero-order hold is consistent
% Halving the control period must change the trajectory by O(dt). An O(1)
% change means the hold is not converging to anything and every number in the
% chapter is an artifact of the sample rate. ch3_test_simulation makes the same
% check for Chapter 3.
%
% COMPARE AT A FIXED TIME, NOT AT THE GUARD. x_end is sampled at the impact,
% and the impact TIME itself moves with dt, so comparing x_end measures the
% state at three different instants of a trajectory travelling at metres per
% second -- an O(1) difference that says nothing about the integrator. The
% comparison here is at t* = 0.15 s, mid-swing and well before any guard.
t_star = 0.15;
d1 = state_at(qc, x0, alpha, 2e-3, t_star);
d2 = state_at(qc, x0, alpha, 1e-3, t_star);
d3 = state_at(qc, x0, alpha, 2.5e-4, t_star);
e1 = norm(d1 - d3);  e2 = norm(d2 - d3);
ratio = e1 / max(e2, 1e-14);
ok = e2 < e1 && ratio > 1.3;
fprintf('  [%s] %-42s |dt=2ms| %.2e, |dt=1ms| %.2e, ratio %.2f\n', tf(ok), ...
        'halving control_dt shrinks the error', e1, e2, ratio);
pass = pass && ok;

%% ------------------------------------- 6. a failing run says why
qf = p;
qf.limits.u_max = 0.5;                       % no authority at all
qf.limits.enable.torque = true;
terr = ch6_terrain(4, [0.30 0.32], 0.05, 1);
sf = ch6_simulate(x0, alpha, qf, terr);
ok = sf.failed && ~isempty(sf.reason) && ~strcmp(sf.reason, 'completed');
fprintf('  [%s] %-42s "%s"\n', tf(ok), ...
        'an impossible run reports a reason', trunc(sf.reason, 44));
pass = pass && ok;

fprintf('  --> %s\n', upper(tf(pass)));

end

% ---------------------------------------------------------------------------
function xt = state_at(p, x0, alpha, dt, t_star)
%STATE_AT  The state at a fixed time, interpolated off the solver's own grid.
q = p;  q.control_dt = dt;
s = ch6_step(x0, alpha, q);
if s.t(end) < t_star
    error('ch6_test_sim:tooShort', ...
          'The step ended at %.4f s, before t* = %.3f s.', s.t(end), t_star);
end
% The output grid is dense (every control period plus the solver's own points),
% so linear interpolation between neighbours is well inside the sampling error
% being measured.
xt = interp1(s.t.', s.x.', t_star, 'linear').';
end

function [x0, alpha] = reference_gait()
persistent X A
if isempty(X)
    root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
    S = load(fullfile(root, 'Results', 'ch3_reference_gait.mat'));
    [Xn, ~, A] = ch3_col_unpack(S.z_opt, S.p);
    X = Xn(:,1);
end
x0 = X;  alpha = A;
end

function s = trunc(s, n)
if numel(s) > n, s = [s(1:n-3) '...']; end
end

function s = tf(b)
if b, s = 'pass'; else, s = 'FAIL'; end
end
