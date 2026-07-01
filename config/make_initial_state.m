function [x0] = make_initial_state()
    nq = 7;
    x0 = zeros(2 * nq, 1);
    
    % Configuration
    q0 = [0; 0; 0.1; -0.3; 0.6; -1.0; 0.6]; % [px, pz, qt, q1, q2, q3, q4]
    x0(1:nq) = q0;
    
    % Velocity
    dq0 = zeros(nq, 1);
    dq0(3) = 0.3; 
    
    % Constraint enforcement
    J = J_st(q0); 
    dq0_corrected = (eye(nq) - pinv(J) * J) * dq0;
    x0(nq + 1:end) = dq0_corrected;
end
