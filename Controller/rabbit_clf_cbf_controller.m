function [u, delta] = rabbit_clf_cbf_controller(t, x, params, target_stone_idx)
    % ---------------------------------------------------------------------
    % 1. Robot Dynamics: M(q)*ddq + H(q,dq) = B*u
    % ---------------------------------------------------------------------
    p = packParameters(params);
    
    q = x(1:7,1);
    dq = x(8:14,1);
    
    D = D_matrix(q,p);
    C = C_vector(q,dq,p);
    G = G_vector(q,p);
    H = C + G;
    B = input_matrix();
    D_inv = inv(D);

    % ---------------------------------------------------------------------
    % 2. CLF Formulation (Trajectory Tracking)
    % ---------------------------------------------------------------------
    % Get actual outputs (e.g., actuated joints or virtual constraints)
    [y, J_y, Jdot_y_dq] = get_ActualOutputs(q, dq); 
    [y_d, dy_d, ddy_d] = get_DesiredOutputs(t); 

    e = y - y_d;             
    de = J_y * dq - dy_d;      
    eta = [e; de];           

    % Error dynamics drift and control terms
    Lfe_drift = Jdot_y_dq - J_y * D_inv * H - ddy_d;
    Lge_ctrl  = J_y * D_inv * B;

    F_drift = [de; Lfe_drift];
    G_ctrl  = [zeros(length(e), size(B,2)); Lge_ctrl];

    % CLF Matrix P 
    num_outputs = length(e);
    Kp_mat = 100 * eye(num_outputs); 
    Kd_mat = 20 * eye(num_outputs);
    P = [Kp_mat + 0.5*Kd_mat^2, 0.5*Kd_mat; 
         0.5*Kd_mat,            0.5*eye(num_outputs)]; 

    % Lie derivatives of V(eta) = eta' * P * eta
    V = eta' * P * eta;
    LfV = 2 * eta' * P * F_drift;
    LgV = 2 * eta' * P * G_ctrl;
            
    c_clf = 10; % CLF decay rate
        
    % ---------------------------------------------------------------------
    % 3. HOCBF Formulation (Stepping on Stones)
    % ---------------------------------------------------------------------
    % Get Target Stone Parameters directly from the main params struct
    stone = params.stones(target_stone_idx, :); % Fixed: Changed p.stones{...} to params.stones(...)
    start_x = stone(1);
    end_x = stone(2);
    stone_width = end_x - start_x;
    center_x = start_x + stone_width / 2;

    % Forward Kinematics for Swing Foot
    [~,p_swing,~,~,~,~] = rabbit_kinematics(q,p); % Fixed typo: param -> p
    J_sw = J_swing(q,p); % Fixed signature: removed dq
    Jdot_sw_dq = Jdotdq_swing(q,dq,p); % Returns 2x1 drift vector [ax_drift; az_drift]

    p_swing_x = p_swing(1);
    p_swing_z = p_swing(2);
    v_swing_x = J_sw(1,:) * dq; % Linear horizontal velocity of swing foot

    % Gaussian Safe Boundary
    A = 0.15; 
    sigma = stone_width / 4;
    d_horizontal = p_swing_x - center_x;
    
    exp_term = exp(-(d_horizontal^2) / (2 * sigma^2));
    z_boundary = A * (1 - exp_term);
    
    % Boundary gradients w.r.t x
    dz_bound_dx = A * (d_horizontal / sigma^2) * exp_term;
    d2z_bound_dx2 = A * exp_term * ( (1/sigma^2) - (d_horizontal^2 / sigma^4) );

    % CBF Function h(q)
    h = p_swing_z - z_boundary;

    % First Derivative of h (psi_1)
    Jh = J_sw(2,:) - dz_bound_dx * J_sw(1,:); % Jacobian of h (1x7 row vector)
    h_dot = Jh * dq;
    
    gamma_1 = 10;
    psi_1 = h_dot + gamma_1 * h;

    % Second Derivative preparation (Direct scalar isolation of drift terms)
    a_drift_x = Jdot_sw_dq(1);
    a_drift_z = Jdot_sw_dq(2);
    
    % Combined scalar drift acceleration for h_ddot
    h_drift = a_drift_z - (d2z_bound_dx2 * (v_swing_x^2) + dz_bound_dx * a_drift_x);
    
    % Final Lie derivatives for the HOCBF acceleration constraint
    Lf_psi1 = h_drift - Jh * D_inv * H + gamma_1 * h_dot;
    Lg_psi1 = Jh * D_inv * B;

    gamma_2 = 10;

    % ---------------------------------------------------------------------
    % 4. Quadratic Program (QP) Setup
    % ---------------------------------------------------------------------
    weight_u = 1;
    weight_delta = 1000;
    H_qp = diag([weight_u, weight_u, weight_u, weight_u, weight_delta]);
    f_qp = zeros(5, 1);

    % Inequality Constraints: A_ineq * z <= b_ineq
    % 1. CLF: LgV * u - delta <= -LfV - c_clf * V
    A_clf = [LgV, -1];
    b_clf = -LfV - c_clf * V;

    % 2. CBF: -Lg_psi1 * u <= Lf_psi1 + gamma_2 * psi1
    A_cbf = [-Lg_psi1, 0];
    b_cbf = Lf_psi1 + gamma_2 * psi_1;

    A_ineq = [A_clf; A_cbf];
    b_ineq = [b_clf; b_cbf];

    % Actuator Limits (Bounds)
    u_max = 150;
    lb = [-u_max; -u_max; -u_max; -u_max; 0]; % delta >= 0
    ub = [ u_max;  u_max;  u_max;  u_max; Inf];

    % ---------------------------------------------------------------------
    % 5. Solve QP
    % ---------------------------------------------------------------------
    options = optimoptions('quadprog', 'Display', 'off');
    [z_opt, ~, exitflag] = quadprog(H_qp, f_qp, A_ineq, b_ineq, [], [], lb, ub, [], options);

    if exitflag == 1
        u = z_opt(1:4);
        delta = z_opt(5);
    else
        warning('QP Failed to find a solution. Outputting zero torques.');
        u = zeros(4,1);
        delta = 0;
    end
end


function [y, J_y, Jdot_y_dq] = get_ActualOutputs(q, dq)
    % GET_ACTUALOUTPUTS Isolates the actuated degrees of freedom.
    %
    % Outputs:
    %   y          - 4x1 vector of actual controlled outputs
    %   J_y        - 4x7 Jacobian matrix of the outputs w.r.t q
    %   Jdot_y_dq  - 4x1 vector representing d/dt(J_y)*dq

    % Selection matrix H_0 to isolate the 4 actuated joints (last 4 states)
    H_0 = [zeros(4, 3), eye(4)];
    
    % Actual output values
    y = H_0 * q;
    
    % The output Jacobian is simply the constant selection matrix
    J_y = H_0;
    
    % Because J_y is constant, its time derivative is zero
    Jdot_y_dq = zeros(4, 1);
end


function [y_d, dy_d, ddy_d] = get_DesiredOutputs(t)
    % GET_DESIREDOUTPUTS Generates a smooth joint reference trajectory.
    %
    % Outputs:
    %   y_d   - 4x1 vector of desired joint angles [rad]
    %   dy_d  - 4x1 vector of desired joint velocities [rad/s]
    %   ddy_d - 4x1 vector of desired joint accelerations [rad/s^2]

    % Nominal standing/walking joint targets (Stance Hip, Stance Knee, Swing Hip, Swing Knee)
    y0 = [0.2;  0.4; -0.2;  0.4]; 
    
    % Oscillatory tracking amplitudes and tracking frequency
    Amp = [0.1;  0.15;  0.1;  0.15]; 
    omega = 2 * pi * 1.5; % 1.5 Hz walking frequency cyclic cadence
    
    % Compute analytic trajectories, velocities, and accelerations
    y_d   = y0 + Amp .* sin(omega * t);
    dy_d  = Amp .* omega .* cos(omega * t);
    ddy_d = -Amp .* (omega^2) .* sin(omega * t);
end