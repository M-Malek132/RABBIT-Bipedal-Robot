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

    % STANCE-FOOT DRIFT: the stance foot is held only by an acceleration-level
    % holonomic constraint, so ode45 lets it creep. Re-simulate the full
    % trajectory and measure how far P_st wanders from its t=0 position.
    if status > 0
        [~, x_traj] = simulate_hzd_gait_full(coeffs, x_start, p);
        foot0 = P_st(x_traj(1, 1:nq)');
        drift = 0;
        for kk = 1:size(x_traj, 1)
            fk    = P_st(x_traj(kk, 1:nq)');
            drift = max(drift, norm(fk - foot0));
        end
        report.stance_foot_drift = drift;
    else
        report.stance_foot_drift = NaN;
    end
    report.cost_normalized = total_torque_sq / max(L_step, eps);
    report.max_penetration = max_penetration;
    report.swing_clearance = swing_clearance;
    report.ground_c        = c(1);
    report.clearance_c     = c(2);
    report.step_length_c   = c(3);
    report.stability_c     = c(4);
    report.torque_c        = c(5);
    report.friction_c      = c(6);
    report.grf_c           = c(7);
    report.impulse_c       = c(8);
    report.periodicity     = ceq(1:13);
    report.foot_height     = ceq(14);
    report.foot_velocity   = ceq(15:16);
    report.theta_end       = ceq(17);
    report.velocity_ceq    = ceq(18);

    % ORBITAL STABILITY: spectral radius of the step-to-step Poincare map.
    % rho < 1 => the gait is attracting (walkable); rho > 1 => it falls after a
    % few steps. Computed here as a DIAGNOSTIC regardless of p.enforce_stability
    % (this is ~26 step sims; acceptable for a one-off inspection).
    if status > 0
        report.rho = poincare_stability(coeffs, x_start, p);
    else
        report.rho = NaN;
    end

    % PHYSICAL QUANTITIES for the Table 3.1 constraints (torque, GRF, friction,
    % impact impulse). Reported only -- not yet enforced. Use these to pick
    % RABBIT-appropriate limits (the thesis values are ATRIAS-specific).
    if status > 0
        report.forces = gait_forces(coeffs, x_start, p);
    else
        report.forces = struct('peak_torque_max',NaN,'min_vertical_grf',NaN, ...
            'mean_vertical_grf',NaN,'max_friction_ratio',NaN,'impact_impulse',NaN);
    end

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
    fprintf('stance-foot drift in step:  %.4f m   (should be ~0; large => constraint drift)\n', ...
            report.stance_foot_drift);
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
    fprintf('STABILITY rho (Poincare):   %.4f   %s\n', report.rho, ...
            ternary(report.rho < 1, '[STABLE, <1]', '[UNSTABLE, >=1 -- falls over]'));
    if isfield(p,'enforce_stability') && p.enforce_stability
        fprintf('  stability inequality (<=0 wanted): %+.6f  %s\n', c(4), ...
                ternary(c(4) <= 0, '[OK]', '[VIOLATED]'));
    else
        fprintf('  (stability NOT enforced in the solve; p.enforce_stability = false)\n');
    end

    % --- Physical quantities (Table 3.1); each enforced independently -------
    Fq = report.forces;
    ef_t = isfield(p,'enforce_torque')   && p.enforce_torque;
    ef_f = isfield(p,'enforce_friction') && p.enforce_friction;
    ef_g = isfield(p,'enforce_grf')      && p.enforce_grf;
    ef_i = isfield(p,'enforce_impulse')  && p.enforce_impulse;
    fprintf('\n--- physical quantities (Table 3.1;  * = enforced this solve) ---\n');
    fprintf('%s peak joint torque |u|:  %8.3f Nm   c=%+.3f %s\n', ...
            star(ef_t), Fq.peak_torque_max, c(5), verdict(ef_t, c(5)));
    fprintf('%s max friction |Fx/Fz|:   %8.3f      c=%+.3f %s\n', ...
            star(ef_f), Fq.max_friction_ratio, c(6), verdict(ef_f, c(6)));
    fprintf('%s min vertical GRF F_z:   %8.1f N    c=%+.3f %s  (mean %.1f N, weight ~294 N)\n', ...
            star(ef_g), Fq.min_vertical_grf, c(7), verdict(ef_g, c(7)), Fq.mean_vertical_grf);
    fprintf('%s impact impulse |F_e|:   %8.3f Ns   c=%+.3f %s\n', ...
            star(ef_i), Fq.impact_impulse, c(8), verdict(ef_i, c(8)));

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
    all_names = [{'ground_c','swing_clearance','step_length_c','stability_c', ...
                  'torque_c','friction_c','grf_c','impulse_c'}, labels, ...
                 {'foot_height','foot_vel_1','foot_vel_2','theta_end','velocity'}];
    [worst_val, worst_idx] = max(all_viol);
    fprintf('%s = %.4f\n\n', all_names{worst_idx}, worst_val);
end

function out = ternary(cond, a, b)
    if cond, out = a; else, out = b; end
end

function s = star(enforced)
    if enforced, s = '*'; else, s = ' '; end
end

function s = verdict(enforced, cval)
% Verdict tag for a physical inequality: only meaningful when enforced
% (when disabled the inequality is a placeholder -1).
    if ~enforced
        s = '';
    elseif cval <= 0
        s = '[OK]';
    else
        s = '[VIOLATED]';
    end
end
