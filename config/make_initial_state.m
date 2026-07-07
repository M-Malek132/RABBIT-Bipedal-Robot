function x0 = make_initial_state()
    nq = 7;
    x0 = zeros(2*nq, 1);

    %-----------------------
    % Configuration (q)
    %-----------------------
    q0 = [0; -1; 0.1; -0.3; 0.6; -1.0; 0.6]; % [px, pz, qt, q1, q2, q3, q4]

    % Enforce stance foot on ground (z = 0)
    tmp = P_st(q0);        % expected: scalar stance foot height
    q0(2) = q0(2) + tmp(2);

    x0(1:nq) = q0;

    %-----------------------
    % Velocity (dq)
    %-----------------------
    dq0 = zeros(nq, 1);
    dq0(3) = 0.3;

    % Enforce stance foot zero velocity: J_st(q)*dq = 0
    J = J_st(q0);
    dq0 = (eye(nq) - pinv(J)*J) * dq0;

    x0(nq+1:end) = dq0;
end
