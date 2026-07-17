%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [c, ceq] = hzd_constraints(z, p)
    [coeffs, x_start] = unpack_z(z, p);
    nq = p.nq;

    [x_end, ~, max_penetration, status] = simulate_hzd_gait(coeffs, x_start, p);

    if status < 0
        % Failed sim: return large-but-finite violations, never NaN,
        % so fmincon treats this as "badly infeasible" instead of crashing.
        ceq = 1e3 * ones(2*nq+3, 1);   % matches size computed below
        c   = 1e3;
        return;
    end

    x_next = rabbit_reset_map(rabbit_impact_map(x_end));
    idx_periodic = 2:14;
    ceq_periodicity = x_next(idx_periodic) - x_start(idx_periodic);

    q_start  = x_start(1:nq);
    dq_start = x_start(nq+1:end);
    foot0 = P_st(q_start);
    ceq_foot_height = foot0(2);
    Jst = J_st(q_start);
    ceq_foot_vel = Jst * dq_start;

    q_end = x_end(1:nq);
    ceq_theta_end = theta_of_q(q_end) - p.theta_plus;

    ceq = [ceq_periodicity; ceq_foot_height; ceq_foot_vel; ceq_theta_end];  % 13+1+2+1 = 17 = 2*nq+3

    c = max_penetration - p.ground_tol;
end
