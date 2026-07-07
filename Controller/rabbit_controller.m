function u = rabbit_controller(x)
    % 1. Extract States
    nq = 7;
    q  = x(1:nq);
    dq = x(nq+1:end);

    % 2. Get Robot Dynamics Matrices
    % You need these from your model (e.g., via Symbolic Toolbox or Euler-Lagrange)
    M_matrix = M(q);       % 7x7 Mass Matrix
    C = V([q; dq]); % 7x1 Coriolis/Gravity vector
    
    % 3. Define Phase Variable (s)
    theta = q(4) + 0.5*q(6);
    theta0 = -0.3; thetaf = 0.3;
    s = (theta - theta0)/(thetaf - theta0);
    s = min(max(s,0),1);

    % 4. Virtual Constraints
    % hd should be the desired joint positions
    hd = desired_gait(s);
    dhd = desired_gait_velocity(s) * (dq(4) + 0.5*dq(6)); % d(hd)/dt
    
    % Tracking Error
    y  = q(4:7) - hd;
    dy = dq(4:7) - dhd;

    % 5. Computed Torque Control
    % We want ddq_desired = -Kp*y - Kd*dy
    % So, u = M * (-Kp*y - Kd*dy) + C
    Kp = diag([400, 400, 300, 300]);
    Kd = diag([40, 40, 30, 30]);

    % We assume only joints 4,5,6,7 are actuated
    % Map the control input to the actuated joints
    B = [zeros(3,4); eye(4)]; 
    
    % Compute required acceleration for tracking
    ddq_ref = -Kp*y - Kd*dy;
    
    % The control law
    u = M_matrix(4:7, :) * [zeros(3,1); ddq_ref] + C(4:7);
end
