%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [c, ceq] = hzd_constraints(z, p)
    [coeffs, x_start] = unpack_z(z, p);
    nq = p.nq;

    [x_end, ~, max_penetration, status, swing_clearance, T_step] = simulate_hzd_gait(coeffs, x_start, p);

    if status < 0
        % Failed sim: return a finite "badly infeasible" marker, never NaN,
        % so fmincon steers away instead of crashing.
        %
        % Keep this the same ORDER OF MAGNITUDE as real constraint values.
        % A huge sentinel (this used to be 1e3) poisons fmincon's constraint
        % scaling when the very first evaluation fails: everything afterwards
        % gets divided by that magnitude, and genuine violations of ~1 then
        % look like ~1e-3 and are falsely reported as "constraints satisfied".
        ceq = 10 * ones(2*nq+4, 1);   % matches size computed below
        c   = 10 * ones(3,1);         % matches the 3 inequalities below
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

    % Step length = horizontal stance-foot -> swing-foot distance at impact
    % (Westervelt's p^h_2(q0^-)).
    sw_end = P_sw(q_end);
    st_end = P_st(q_end);
    L_step = sw_end(1) - st_end(1);

    % NEC1 -- average walking rate: step length / step duration = desired
    % speed (Westervelt eq. 6.50 / thesis v_des). This PINS the walking
    % speed; without it, speed was an uncontrolled byproduct of the dynamics.
    % T_step > 0 always (a step took finite time), so the ratio is well posed.
    ceq_velocity = L_step / T_step - p.v_des;

    ceq = [ceq_periodicity; ceq_foot_height; ceq_foot_vel; ceq_theta_end; ceq_velocity];  % 13+1+2+1+1 = 18 = 2*nq+4

    % Inequalities (fmincon wants c <= 0):
    %  1) nothing penetrates the ground at any point in the step
    %  2) NIC3 / mid-step swing-foot clearance: the swing foot must stay at
    %     least p.swing_clearance_min above the ground through the middle of
    %     the step, so the impact surface S is reached ONLY at the end.
    %  3) DESIRED STEP LENGTH (Westervelt NEC1 / thesis "desired step length"):
    %     the robot must actually travel. Without this the optimizer collapses
    %     to a degenerate "step in place" (L_step ~ 0.015 m, swing foot
    %     scuffing) because a tiny step costs little torque.
    c = [ max_penetration - p.ground_tol
          p.swing_clearance_min - swing_clearance
          p.step_length_min - L_step ];
end
