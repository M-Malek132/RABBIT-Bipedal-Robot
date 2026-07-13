function x0 = make_random_initial_state(p)
    nq = p.nq;
    x0 = zeros(2*nq, 1);

    %-----------------------
    % Configuration (q) — randomized within plausible joint ranges
    %-----------------------
    px = 0;                          % arbitrary — periodicity excludes this anyway
    pz = -1 + 0.1*randn();
    qt = 0.1 + 0.1*randn();
    q1 = -0.3 + 0.2*randn();
    q2 =  0.6 + 0.2*randn();
    q3 = -1.0 + 0.2*randn();
    q4 =  0.6 + 0.2*randn();

    q0 = [px; pz; qt; q1; q2; q3; q4];

    % Enforce stance foot on ground (z = 0) — same projection as make_initial_state
    tmp = P_st(q0);
    q0(2) = q0(2) + tmp(2);          % NOTE: sign — see caution below
    x0(1:nq) = q0;

    %-----------------------
    % Velocity (dq) — randomized, then projected to satisfy J_st*dq = 0
    %-----------------------
    dq0 = 0.3*randn(nq,1);

    J = J_st(q0);
    dq0 = (eye(nq) - pinv(J)*J) * dq0;
    x0(nq+1:end) = dq0;
end