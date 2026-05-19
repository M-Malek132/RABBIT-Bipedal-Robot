function u = rabbit_clf_cbf_controller(x, params)
    % 1. Get dynamics
    [D, C, G, B] = rabbit_dynamics(x, params);
    
    % Use feedback linearization matrices
    % xdot = f(x) + g(x)u 
    % (Calculate LfV, LgV, Lfh, Lgh based on your outputs and barrier functions)
    
    % --- Example Placeholders ---
    LfV = 0; LgV = zeros(1, 4); V = 0; % From your Lyapunov function
    Lfh = 0; Lgh = zeros(1, 4); h = 1; % From your Barrier function (e.g., foot height > 0)
    
    % 2. Setup QP Cost Function: min 0.5 * z^T H z + f^T z
    % z = [u1; u2; u3; u4; delta]
    weight_u = 1; 
    weight_delta = 1000; % High penalty for violating stability to ensure safety
    
    H = diag([weight_u, weight_u, weight_u, weight_u, weight_delta]);
    f_qp = zeros(5, 1);
    
    % 3. Setup Inequality Constraints: A_ineq * z <= b_ineq
    
    % CLF constraint: LgV * u - delta <= -LfV - c*V
    c_clf = 10;
    A_clf = [LgV, -1];
    b_clf = -LfV - c_clf * V;
    
    % CBF constraint: -Lgh * u <= Lfh + gamma*h
    gamma_cbf = 10;
    A_cbf = [-Lgh, 0];
    b_cbf = Lfh + gamma_cbf * h;
    
    % Combine inequalities
    A_ineq = [A_clf; A_cbf];
    b_ineq = [b_clf; b_cbf];
    
    % 4. Setup Input Saturation Bounds (Torque limits)
    u_max = 150; % Max torque Nm
    lb = [-u_max; -u_max; -u_max; -u_max; -Inf]; % Lower bounds
    ub = [ u_max;  u_max;  u_max;  u_max;  Inf]; % Upper bounds
    
    % 5. Solve QP using quadprog
    options = optimoptions('quadprog', 'Display', 'off');
    [z_opt, ~, exitflag] = quadprog(H, f_qp, A_ineq, b_ineq, [], [], lb, ub, [], options);
    
    if exitflag == 1
        u = z_opt(1:4); % Extract control inputs
    else
        warning('QP failed to find a solution!');
        u = zeros(4,1); % Fallback
    end
end
