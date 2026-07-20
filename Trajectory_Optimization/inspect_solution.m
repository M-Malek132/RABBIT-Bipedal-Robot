function report = inspect_solution(z, p)
%INSPECT_SOLUTION Break down cost/status/constraint violations for a
% decision vector z (e.g. z_opt returned by fmincon), so you can see
% which term is actually driving infeasibility instead of guessing.

    [coeffs, x_start] = unpack_z(z, p);
    [x_end, total_torque_sq, max_penetration, status, swing_clearance, T_step] = simulate_hzd_gait(coeffs, x_start, p);
    [c, ceq] = hzd_constraints(z, p);

    nq = p.nq;

    % Step length and the normalized cost fmincon actually minimizes
    q_end  = x_end(1:nq);
    sw     = P_sw(q_end);
    st     = P_st(q_end);
    L_step = sw(1) - st(1);
    v_avg  = L_step / max(T_step, eps);   % achieved average walking speed

    report.status          = status;
    report.total_torque_sq = total_torque_sq;
    report.step_length     = L_step;
    report.step_duration   = T_step;
    report.walking_speed   = v_avg;
    report.cost_normalized = total_torque_sq / max(L_step, eps);
    report.max_penetration = max_penetration;
    report.swing_clearance = swing_clearance;
    report.ground_c        = c(1);
    report.clearance_c     = c(2);
    report.step_length_c   = c(3);
    report.periodicity     = ceq(1:13);
    report.foot_height     = ceq(14);
    report.foot_velocity   = ceq(15:16);
    report.theta_end       = ceq(17);
    report.velocity_ceq    = ceq(18);

    fprintf('\n=== Solution Inspection ===\n');
    if status < 0
        fprintf('STATUS: FAILED (impact event never fired, or NaNs) -- everything below is meaningless.\n');
    else
        fprintf('STATUS: ok (impact event fired)\n');
    end
    fprintf('theta: start %+.4f -> end %+.4f   (theta_plus target %+.4f)\n', ...
            theta_of_q(x_start(1:nq)), theta_of_q(q_end), p.theta_plus);
    fprintf('step length:                %.4f m\n', L_step);
    fprintf('step duration:              %.4f s\n', T_step);
    fprintf('walking speed L/T:          %.4f m/s   (v_des %.4f, NEC1 resid %+.4f)\n', ...
            v_avg, p.v_des, ceq(18));
    fprintf('raw   torque^2:             %.4f\n', total_torque_sq);
    fprintf('COST (torque^2 / L_step):   %.4f\n', report.cost_normalized);
    fprintf('ground inequality  (<=0 wanted): %+.6f  %s\n', c(1), ternary(c(1) <= 0, '[OK]', '[VIOLATED]'));
    fprintf('  (max_penetration = %+.6f; positive means something went below ground)\n', max_penetration);
    fprintf('swing clearance    (<=0 wanted): %+.6f  %s\n', c(2), ternary(c(2) <= 0, '[OK]', '[VIOLATED]'));
    fprintf('  (min swing-foot height mid-step = %+.6f m, need >= %.3f)\n', swing_clearance, p.swing_clearance_min);
    fprintf('step length        (<=0 wanted): %+.6f  %s\n', c(3), ternary(c(3) <= 0, '[OK]', '[VIOLATED]'));
    fprintf('  (L_step = %+.4f m, need >= %.3f)\n', L_step, p.step_length_min);

    fprintf('\nperiodicity ceq (want ~0), state idx 2:14 = [pz qt q1 q2 q3 q4  dpx dpz dqt dq1 dq2 dq3 dq4]:\n');
    labels = {'pz','qt','q1','q2','q3','q4','dpx','dpz','dqt','dq1','dq2','dq3','dq4'};
    for i = 1:numel(labels)
        fprintf('  %-4s : %10.4f\n', labels{i}, report.periodicity(i));
    end
    fprintf('  max |periodicity| = %.4f\n', max(abs(report.periodicity)));

    fprintf('\nfoot_height ceq (want ~0):   %10.4f\n', report.foot_height);
    fprintf('foot_velocity ceq (want ~0): %10.4f  %10.4f\n', report.foot_velocity(1), report.foot_velocity(2));
    fprintf('theta_end ceq (want ~0):     %10.4f\n', report.theta_end);
    fprintf('velocity  ceq (want ~0):     %10.4f\n', report.velocity_ceq);

    fprintf('\nlargest single violation: ');
    all_viol = [max(0, c(:)); abs(report.periodicity(:)); abs(report.foot_height); ...
                abs(report.foot_velocity(:)); abs(report.theta_end); abs(report.velocity_ceq)];
    all_names = [{'ground_c','swing_clearance','step_length_c'}, labels, ...
                 {'foot_height','foot_vel_1','foot_vel_2','theta_end','velocity'}];
    [worst_val, worst_idx] = max(all_viol);
    fprintf('%s = %.4f\n\n', all_names{worst_idx}, worst_val);
end

function out = ternary(cond, a, b)
    if cond, out = a; else, out = b; end
end
