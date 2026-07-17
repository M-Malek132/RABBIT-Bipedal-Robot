function report = inspect_solution(z, p)
%INSPECT_SOLUTION Break down cost/status/constraint violations for a
% decision vector z (e.g. z_opt returned by fmincon), so you can see
% which term is actually driving infeasibility instead of guessing.

    [coeffs, x_start] = unpack_z(z, p);
    [x_end, total_torque_sq, max_penetration, status] = simulate_hzd_gait(coeffs, x_start, p);
    [c, ceq] = hzd_constraints(z, p);

    nq = p.nq;
    idx_periodic = 2:14;   % matches hzd_constraints.m: excludes px (translation gauge)

    report.status          = status;
    report.total_torque_sq = total_torque_sq;
    report.max_penetration = max_penetration;
    report.ground_c        = c;
    report.periodicity     = ceq(1:13);
    report.foot_height     = ceq(14);
    report.foot_velocity   = ceq(15:16);
    report.theta_end       = ceq(17);

    fprintf('\n=== Solution Inspection ===\n');
    if status < 0
        fprintf('STATUS: FAILED (impact event never fired, or NaNs) -- everything below is meaningless.\n');
    else
        fprintf('STATUS: ok (impact event fired)\n');
    end
    fprintf('total torque^2 cost:        %.4f\n', total_torque_sq);
    fprintf('ground inequality c (<=0 wanted): %.6f  %s\n', c, ternary(c <= p.ground_tol, '[OK]', '[VIOLATED]'));
    fprintf('  (max_penetration = %.6f; positive means something went below ground)\n', max_penetration);

    fprintf('\nperiodicity ceq (want ~0), state idx 2:14 = [pz qt q1 q2 q3 q4  dpx dpz dqt dq1 dq2 dq3 dq4]:\n');
    labels = {'pz','qt','q1','q2','q3','q4','dpx','dpz','dqt','dq1','dq2','dq3','dq4'};
    for i = 1:numel(labels)
        fprintf('  %-4s : %10.4f\n', labels{i}, report.periodicity(i));
    end
    fprintf('  max |periodicity| = %.4f\n', max(abs(report.periodicity)));

    fprintf('\nfoot_height ceq (want ~0):   %10.4f\n', report.foot_height);
    fprintf('foot_velocity ceq (want ~0): %10.4f  %10.4f\n', report.foot_velocity(1), report.foot_velocity(2));
    fprintf('theta_end ceq (want ~0):     %10.4f\n', report.theta_end);

    fprintf('\nlargest single violation: ');
    all_viol = [max(0, c); abs(report.periodicity(:)); abs(report.foot_height); ...
                abs(report.foot_velocity(:)); abs(report.theta_end)];
    all_names = [{'ground_c'}, labels, {'foot_height','foot_vel_1','foot_vel_2','theta_end'}];
    [worst_val, worst_idx] = max(all_viol);
    fprintf('%s = %.4f\n\n', all_names{worst_idx}, worst_val);
end

function out = ternary(cond, a, b)
    if cond, out = a; else, out = b; end
end
