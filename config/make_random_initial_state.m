%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% RANDOM INITIAL STATE — rejection-sampled: ground-valid AND far enough
% from theta_plus to avoid the phase-window denominator blowing up
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function x0 = make_random_initial_state(p)
    nq = p.nq;
    ok = false;

    for attempt = 1:50
        px = 0;
        pz = -1 + 0.1*randn();
        qt = 0.1 + 0.1*randn();
        q1 = -0.3 + 0.2*randn();
        q2 =  0.6 + 0.2*randn();
        q3 = -1.0 + 0.2*randn();
        q4 =  0.6 + 0.2*randn();
        q0 = [px; pz; qt; q1; q2; q3; q4];

        tmp = P_st(q0);
        q0(2) = q0(2) + tmp(2);   % project stance foot onto ground

        th = theta_of_q(q0);
        r  = check_ground_validity(q0, []);

        if r.is_valid && abs(p.theta_plus - th) > p.min_theta_gap
            ok = true;
            break;
        end
    end

    if ~ok
        warning('make_random_initial_state: no valid sample found in 50 attempts, using last draw anyway.');
    end

    x0 = zeros(2*nq,1);
    x0(1:nq) = q0;

    dq0 = 0.3*randn(nq,1);
    J = J_st(q0);
    dq0 = (eye(nq) - pinv(J)*J) * dq0;
    x0(nq+1:end) = dq0;
end
